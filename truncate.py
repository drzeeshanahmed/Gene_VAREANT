import argparse
import os
import sys
import multiprocessing as mp
import itertools


def strip_text(text, start, end):
    # type: (str, str, str) -> str
    start = len(start) if text.startswith(start) else 0
    end = -len(end) if text.endswith(end) else len(text)
    return text[start:end]


class MetaLine:
    def __init__(self, line):
        # type: (MetaLine, str) -> None
        tokens = strip_text(line, "##", "\n").split("=", 1)
        assert (
            len(tokens) == 2
        ), "Invalid VCF File. Meta-Information line must be of form ##Key=Value."
        self.key = tokens[0]
        self.value = tokens[1]

    def __str__(self):
        return "##{}={}\n".format(self.key, self.value)

    def get_tokens(self):
        # type: (MetaLine) -> dict[str, str]
        d = dict()
        if self.value.startswith("<") and self.value.endswith(">"):
            in_quotes = False
            last = 0
            text = self.value[1:-2] + ","
            for i in range(len(text)):
                if text[i] == "," and not in_quotes:
                    token = text[last:i].strip().split("=", 1)
                    last = i + 1
                    d[token[0]] = token[1].strip('"') if len(token) == 2 else None
                elif text[i] == '"':
                    in_quotes = not in_quotes
        return d


class Header:
    def __init__(self, line):
        # type: (Header, str) -> None
        tokens = strip_text(line, "#", "\n").split("\t")
        assert (
            len(tokens) >= 8
        ), "Invalid VCF File. Header line must have 8 fixed columns."
        assert (
            tokens[0] == "CHROM"
        ), "Invalid VCF File. First header column must be CHROM."
        assert tokens[1] == "POS", "Invalid VCF File. First header column must be POS."
        assert tokens[2] == "ID", "Invalid VCF File. First header column must be ID."
        assert tokens[3] == "REF", "Invalid VCF File. First header column must be REF."
        assert tokens[4] == "ALT", "Invalid VCF File. First header column must be ALT."
        assert (
            tokens[5] == "QUAL"
        ), "Invalid VCF File. First header column must be QUAL."
        assert (
            tokens[6] == "FILTER"
        ), "Invalid VCF File. First header column must be FILTER."
        assert (
            tokens[7] == "INFO"
        ), "Invalid VCF File. First header column must be INFO."
        if len(tokens) > 8:
            assert (
                tokens[8] == "FORMAT"
            ), "Invalid VCF File. Optional 9th column must be FORMAT."

        self.hasSamples = False
        if len(tokens) > 8:
            self.hasSamples = True
            self.format = tokens[8]
            self.samples = tokens[9:]
            self.samples_mask = [True for _ in self.samples]

        self.tokens = tokens

    def __str__(self):
        text = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        if self.hasSamples:
            text += "\tFORMAT"
            if len(self.samples) > 0:
                text += "\t" + "\t".join(
                    [s for s, m in zip(self.samples, self.samples_mask) if m]
                )
        return text + "\n"


class Variant:
    def __init__(self, line, header):
        # type: (Variant, str, Header) -> None
        tokens = strip_text(line, "", "\n").split("\t")
        assert len(tokens) == len(header.tokens), (
            "Invalid VCF File. Missing columns for variant:" + line
        )

        self.chrom = tokens[0]
        self.pos = int(tokens[1])
        self.id = tokens[2]
        self.ref = tokens[3]
        self.alt = tokens[4]
        self.qual = float(tokens[5])
        self.filter = tokens[6]
        self.info = dict()
        for pair in tokens[7].split(";"):
            p = pair.split("=", 1)
            self.info[p[0]] = p[1] if len(p) == 2 else True

        self.info_csq = self.info.get("CSQ", "")
        self.info_ann = self.info.get("ANN", "")

        self.hasSamples = False
        if len(tokens) > 8:
            self.hasSamples = True
            self.format = tokens[8]
            self.samples = tokens[9:]

        self._header = header

    def __str__(self):
        tokens = [
            self.chrom,
            str(self.pos),
            self.id,
            self.ref,
            self.alt,
            str(self.qual),
            self.filter,
            ";".join(
                [
                    (k + "=" + v) if isinstance(v, str) else k
                    for k, v in self.info.items()
                ]
            ),
        ]
        if self.format is not None:
            tokens.append(self.format)
        if self.samples is not None:
            tokens.extend(
                [s for s, m in zip(self.samples, self._header.samples_mask) if m]
            )

        return "\t".join(tokens) + "\n"


class VCFFilterer:
    def __init__(self, config):
        # type: (VCFFilterer, dict[str, str]) -> None
        c = config
        self.retain_info_entries = (
            set(c.get("Retain-Info-Entries").split("|"))
            if "Retain-Info-Entries" in c
            else None
        )
        self.retain_variant_by_gene_symbol = (
            set(c.get("Retain-Variant-By-Gene-Symbol").split("|"))
            if "Retain-Variant-By-Gene-Symbol" in c
            else None
        )
        self.retain_variant_by_ensembl_id = (
            set(c.get("Retain-Variant-By-Ensembl-ID").split("|"))
            if "Retain-Variant-By-Ensembl-ID" in c
            else None
        )
        self.retain_variant_by_rsid = (
            set(c.get("Retain-Variant-By-rsID").split("|"))
            if "Retain-Variant-By-rsID" in c
            else None
        )
        self.retain_if_passed_filters = (
            set(c.get("Retain-If-Passed-Filters").split("|"))
            if "Retain-If-Passed-Filters" in c
            else None
        )
        self.reformat_genotype = (
            c.get("Reformat-Genotype").split(":") if "Reformat-Genotype" in c else None
        )
        self.minimum_qual = (
            float(c.get("Minimum-Quality-Threshold"))
            if "Minimum-Quality-Threshold" in c
            else None
        )
        self.retain_samples = (
            set(c.get("Retain-Samples").split("|")) if "Retain-Samples" in c else None
        )

        if self.reformat_genotype is not None:
            assert (
                self.reformat_genotype[0] == "GT"
            ), "Invalid configuration. First Genotype entry must be 'GT'."

    def should_keep_meta(self, meta):
        # type: (VCFFilterer, MetaLine) -> bool
        if (
            meta.key == "FORMAT"
            and self.reformat_genotype is not None
            and meta.get_tokens().get("ID") not in self.reformat_genotype
        ):
            return False
        elif (
            meta.key == "INFO"
            and self.retain_info_entries is not None
            and meta.get_tokens().get("ID") not in self.retain_info_entries
        ):
            return False
        elif meta.key == "FILTER" and self.retain_if_passed_filters is not None:
            return True
            # id = meta.get_tokens().get("ID")
            # if "PASS" in self.retain_if_passed_filters and id != "PASS":
            #     return False
            # elif id in self.retain_if_passed_filters:
            #     return False
        return True

    def transform_header(self, header):
        # type: (VCFFilterer, Header) -> Header
        if self.retain_samples is not None:
            for i, sample in enumerate(header.samples):
                header.samples_mask[i] = sample in self.retain_samples

    def transform_variant(self, variant):
        # type: (VCFFilterer, Variant) -> None
        if self.retain_info_entries:
            variant.info = {
                k: v for k, v in variant.info.items() if k in self.retain_info_entries
            }
        if self.reformat_genotype and variant.hasSamples:
            old_format = variant.format.split(":")
            indexes = [
                old_format.index(f) for f in self.reformat_genotype if f in old_format
            ]
            new_format = [f for f in self.reformat_genotype if f in old_format]

            variant.format = ":".join(new_format)
            for i in range(len(variant.samples)):
                s = variant.samples[i].split(":")
                assert len(s) == len(
                    old_format
                ), "Invalid VCF. Each sample must adhere to FORMAT column."
                variant.samples[i] = ":".join([s[idx] for idx in indexes])

    def should_keep_variant(self, variant):
        # type: (VCFFilterer, Variant) -> bool
        if self.minimum_qual is not None and variant.qual < self.minimum_qual:
            return False

        if self.retain_if_passed_filters is not None:
            if "PASS" in self.retain_if_passed_filters:
                if variant.filter != "PASS":
                    return False
            else:
                for failed in variant.filter.split(";"):
                    if failed in self.retain_if_passed_filters:
                        return False
        # Genes
        if (
            self.retain_variant_by_rsid is not None
            or self.retain_variant_by_gene_symbol is not None
            or self.retain_variant_by_ensembl_id is not None
        ):
            matches_gene_symbol = False
            matches_ensembl_id = False
            matches_rsid = False

            gene_symbols = (
                self.retain_variant_by_gene_symbol
                if self.retain_variant_by_gene_symbol is not None
                else set()
            )
            ensembl_ids = (
                self.retain_variant_by_ensembl_id
                if self.retain_variant_by_ensembl_id is not None
                else set()
            )
            rsids = (
                self.retain_variant_by_rsid
                if self.retain_variant_by_rsid is not None
                else set()
            )

            for i in variant.info_ann.split(","):
                ann = i.split("|")
                if len(ann) == 16:
                    matches_gene_symbol = (
                        matches_gene_symbol
                        or len(gene_symbols.intersection(ann[3].split("&"))) > 0
                    )
                    matches_ensembl_id = (
                        matches_ensembl_id
                        or len(ensembl_ids.intersection(ann[4].split("&"))) > 0
                    )

            # for i in variant.info_csq.split(","):
            #     csq = i.split("|")
            #     if len(csq) == 42:
            #         matches_gene_symbol = matches_gene_symbol or len(gene_symbols.intersection(csq[3].split("&"))) > 0
            #         matches_ensembl_id = matches_ensembl_id or len(ensembl_ids.intersection(csq[4].split("&"))) > 0
            #         matches_rsid = matches_rsid or len(rsids.intersection(csq[17].split("&"))) > 0

            if not matches_gene_symbol and not matches_rsid and not matches_ensembl_id:
                return False

        return True


def job(data):
    line, header, filterer = data
    variant = Variant(line, header)
    filterer.transform_variant(variant)
    if filterer.should_keep_variant(variant):
        return str(variant)
    return None


def truncate_vcf(filterer, input, output, keep_meta):
    line = input.readline()
    assert (
        line == "##fileformat=VCFv4.2\n"
    ), "Invalid VCF File. First line must be ##fileformat=VCFv4.2"
    output.write(line)

    for line in input:
        if not line.startswith("##"):
            break
        meta = MetaLine(line)
        if keep_meta or filterer.should_keep_meta(meta):
            output.write(str(meta))

    output.write("##IntelliGenes-Truncate-VCF='" + " ".join(sys.argv) + "'\n")

    assert line.startswith(
        "#"
    ), "Invalid VCF File. Header line must follow meta-information lines."
    header = Header(line)
    filterer.transform_header(header)
    output.write(str(header))

    pool = mp.Pool()  # uses cpu_count
    for line in pool.imap(
        job,
        zip(input, itertools.repeat(header), itertools.repeat(filterer)),
        chunksize=40000,
    ):  # 40000 seems to work well
        if line is not None:
            output.write(line)
    pool.close()


def main(config_path, input_path, output_path, keep_meta):
    print("Truncating VCF...")
    with open(config_path, "r") as config_file:
        config = dict()
        for line in config_file:
            if line.startswith("#"):
                continue
            
            config_parameter = line.split(": ")
            assert (
                len(config_parameter) == 2
            ), "Config parameter must be of the form 'Key: Value'. See configuration specification."
            config[config_parameter[0]] = strip_text(config_parameter[1], "", "\n")

        filterer = VCFFilterer(config)

    output_file = os.path.join(output_path, "truncated.vcf")
    with open(input_path, "r") as input, open(output_file, "w") as output:
        truncate_vcf(filterer, input, output, keep_meta)
    
    print("Finished truncating VCF")
    return output_file
