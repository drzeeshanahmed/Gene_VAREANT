import os
import sqlite3
import pandas as pd
import re
import warnings
import io
import pysam

from log import print_text

warnings.filterwarnings("ignore")


def get_tokens(text: str) -> dict[str, str]:
    d = dict()
    inQuotes = False
    last = 0
    text += ","  # end with comma just so that we get all tokens including last

    for i in range(len(text)):
        if text[i] == "," and not inQuotes:
            token = text[last:i].strip().split("=", maxsplit=1)
            last = i + 1
            if len(token) == 1:
                d[token[0]] = None
            else:  # len == 2
                d[token[0]] = token[1].strip('"')
        elif text[i] == '"':
            inQuotes = not inQuotes

    return d


class Filter:
    def __init__(self, d: dict[str, str]):
        self.id = d.get("ID")
        self.description = d.get("Description")

    def __str__(self) -> str:
        return f"ID({self.id}); Description({self.description})"


class Format:
    def __init__(self, d: dict[str, str]):
        self.id = d.get("ID")
        self.number = d.get("Number")
        self.type = d.get("Type")
        self.description = d.get("Description")

    def __str__(self) -> str:
        return f"ID({self.id}); Number({self.number}); Type({self.type}); Description({self.description})"


class Info:
    def __init__(self, d: dict[str, str]):
        self.id = d.get("ID")
        self.number = d.get("Number")
        self.type = d.get("Type")
        self.description = d.get("Description")
        self.source = d.get("Source")
        self.version = d.get("Version")

    def __str__(self) -> str:
        return f"ID({self.id}); Number({self.number}); Type({self.type}); Description({self.description}); Source({self.source}); Version({self.version})"


class Contig:
    def __init__(self, d: dict[str, str]):
        self.id = d.get("ID")
        self.length = d.get("length")

    def __str__(self) -> str:
        return f"ID({self.id}); Length({self.length})"


filters: list[Filter] = []
formats: list[Format] = []
infos: list[Info] = []
contigs: list[Contig] = []
cols: list[str] = []


class Variant:
    def __init__(
        self, cols: list[str], tokens: list[str], extraction_value: str | None
    ):
        # If VCF is missing format column, add the format columns
        if len(cols) == 8:
            [chrom, pos, id, ref, alt, qual, filter, info] = tokens
            format = None
            samples = None
        else:
            [chrom, pos, id, ref, alt, qual, filter, info, format, *samples] = tokens

        self.chrom = chrom
        try:
            self.pos = int(pos)
        except Exception:
            self.pos = None

        self.id = id
        self.ref = ref
        self.alt = alt
        try:
            self.qual = float(qual)
        except Exception:
            self.qual = None

        self.filter = None if filter == "PASS" else filter
        self.gene_symbol = None

        self.info = dict()
        for token in info.split(";"):
            split = token.split("=", maxsplit=1)
            if len(split) == 1:
                self.info[split[0]] = True
            elif len(split) == 2:
                self.info[split[0]] = split[1]

        self.format = format.split(":")
        self.samples = [None if s.startswith("./.") else s for s in samples]
        if extraction_value not in self.format:
            # binary matrix
            self.values = [int(s is not None) for s in self.samples]
        else:
            value_index = self.format.index(extraction_value)
            self.values = [
                s.split(":")[value_index] if s is not None else None
                for s in self.samples
            ]

    def __str__(self) -> str:
        return f"""VARIANT:
    CHROM({self.chrom});
    POS({self.pos});
    ID({self.id});
    REF({self.ref});
    ALT({self.alt});
    QUAL({self.qual});
    FILTER({self.filter});
    INFO({self.info});
    FORMAT({self.format});
    Samples({self.samples});
    Values({self.values});"""


variants: list[Variant] = []


def main(config_path, input_path, output_path, log_file):
    print_text(log_file, "Extracting VCF")
    extract_unidentified_variants = True
    extraction_value = None
    print_text(log_file, "Reading from configuration " + config_path)
    with open(config_path, "r") as config_file:
        for line in config_file:
            print_text(log_file, line)
            if line.startswith("Extract-Unidentified-Variants: False"):
                extract_unidentified_variants = False
            if line.startswith("Extract-Genotype-Value: "):
                extraction_value = line[len("Extract-Genotype-Value: ") :].strip()

    db = os.path.join(output_path, "variant_db.sqlite")
    cigt_path = os.path.join(output_path, "cigt_matrix.cigt.csv")

    print_text(log_file, "Saving SQLite to " + db)
    print_text(log_file, "Saving CIGT to " + cigt_path)

    print_text(log_file, "Reading input VCF at " + input_path)
    if not (input_path.endswith(".vcf.gz") or input_path.endswith(".vcf")):
        raise ValueError(
            "Must provide a compressed (.vcf.gz) or uncompressed (.vcf) input file."
        )

    count = 0

    with (
        io.TextIOWrapper(pysam.BGZFile(input_path, "r"), "utf-8")
        if input_path.endswith(".vcf.gz")
        else open(input_path, "r")
    ) as f:
        while line := f.readline():
            line = line.strip()  # remove new line characters

            if line.startswith("##"):
                if (m := re.fullmatch("##FILTER=<(.*)>", line)) is not None:
                    if not line.startswith("##FILTER=<ID=PASS"):
                        filters.append(Filter(get_tokens(m.group(1))))
                elif (m := re.fullmatch("##FORMAT=<(.*)>", line)) is not None:
                    formats.append(Format(get_tokens(m.group(1))))
                elif (m := re.fullmatch("##INFO=<(.*)>", line)) is not None:
                    infos.append(Info(get_tokens(m.group(1))))
                elif (m := re.fullmatch("##contig=<(.*)>", line)) is not None:
                    contigs.append(Contig(get_tokens(m.group(1))))
            elif line.startswith("#"):  # main heading line
                cols = line.strip("#").split("\t")
                for i in range(len(cols)):
                    m = re.match("DNA-BR2-(.+)_filtered_snps_dip1", cols[i])
                    if m is not None:
                        cols[i] = m.group(1)
            else:  # variant entries
                var = Variant(cols, line.split("\t"), extraction_value)
                if var.id != "." or extract_unidentified_variants:
                    variants.append(var)

                count += 1
                if count % 5000 == 0:
                    print_text(log_file, "# variants read: " + str(count))

    print_text(log_file, "# variants read: " + str(count))

    print_text(log_file, "Opening DB connection")

    sq = sqlite3.connect(db)
    cur = sq.cursor()

    info_columns = []
    print_text(log_file, "Preparing INFO headers")
    for info in infos:
        info_columns.append(f"WESI_{info.id} {info.type.upper()}")
        info_columns.append(f"WESI_{info.id}_Description TEXT")

    sample_columns = []
    print_text(log_file, "Preparing FORMAT definitions")
    for f in formats:
        sample_columns.append(f"WESS_{f.id} TEXT")
        sample_columns.append(f"WESS_{f.id}_Description TEXT")

    print_text(log_file, "Creating WES_VARIANT table")
    cur.execute("DROP TABLE IF EXISTS WES_VARIANT;")
    print_text(log_file, "Creating WES_SAMPLE table")
    cur.execute("DROP TABLE IF EXISTS WES_SAMPLE;")
    print_text(log_file, "Creating WES_INFO table")
    cur.execute("DROP TABLE IF EXISTS WES_INFO;")

    cur.execute("""
    CREATE TABLE WES_VARIANT(
        WESV_Id INTEGER PRIMARY KEY,
        WESV_Crom_Number TEXT,
        WESV_Crom_Position INTEGER,
        WESV_Identifier TEXT,
        WESV_Ref_Base TEXT,
        WESV_Alt_Base TEXT,
        WESV_Quality INTEGER,
        WESV_Filter TEXT,
        WESV_Filter_Desc TEXT
    );""")

    cur.execute(f"""
    CREATE TABLE WES_SAMPLE(
        WESS_Id INTEGER PRIMARY KEY,
        WESS_WESV_Id INTEGER,
        {', '.join(sample_columns)},
        FOREIGN KEY(WESS_WESV_Id) REFERENCES WES_VARIANT(WESV_Id)
    );""")

    cur.execute(f"""
    CREATE TABLE WES_INFO(
        WESI_Id INTEGER PRIMARY KEY,
        WESI_WESV_Id INTEGER,
        {', '.join(info_columns)},
        FOREIGN KEY(WESI_WESV_Id) REFERENCES WES_VARIANT(WESV_Id)
    );""")

    filter_desc_map = {f.id: f.description for f in filters}
    format_desc_map = {f.id: f.description for f in formats}

    print_text(
        log_file, "Inserting variant data (chrom, pos, id, ref, alt, qual, filter)"
    )
    print_text(
        log_file,
        f"Inserting sample data (variant_id, {', '.join([f.id for f in formats])})",
    )
    print_text(
        log_file,
        f"Inserting info annotations (variant_id, {', '.join([i.id for i in infos])})",
    )
    print_text(log_file, f"Processing {len(variants)} total variants")
    count = 0
    for i, v in enumerate(variants):
        index = i + 1
        filter_desc = (
            None
            if v.filter is None
            else ";".join([filter_desc_map.get(s) for s in v.filter.split(";")])
        )
        cur.execute(
            "INSERT INTO WES_VARIANT VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?);",
            (index, v.chrom, v.pos, v.id, v.ref, v.alt, v.qual, v.filter, filter_desc),
        )

        info_entries = [None, index]
        for info in infos:
            info_entries.append(v.info.get(info.id))
            info_entries.append(info.description)

        cur.execute(
            f"INSERT INTO WES_INFO VALUES(?, ?{', ?' * len(info_columns)});",
            tuple(info_entries),
        )

        # sample_entries
        sample_entries = []
        for sample in v.samples:
            if sample is None:
                continue

            entries = [None, index]
            mapping = {k: s for k, s in zip(v.format, sample.split(":"))}
            for k in formats:
                entries.append(mapping.get(k.id))
                entries.append(k.description)

            sample_entries.append(tuple(entries))

        cur.executemany(
            f"INSERT INTO WES_SAMPLE VALUES(?, ?{', ?' * len(sample_columns)})",
            sample_entries,
        )

        count += 1
        if count % 5000 == 0:
            print_text(log_file, "# variants finished: " + str(count))

    print_text(log_file, "# variants finished: " + str(count))
    print_text(log_file, "Closing DB connection")
    sq.commit()
    sq.close()

    print_text(log_file, "Creating DataFrame")
    # CIGT
    df = pd.DataFrame([v.values for v in variants])
    df.columns = cols[9:]
    df.index = [v.id for v in variants]

    print_text(log_file, "Saving CIGT DataFrame")
    df.transpose().to_csv(cigt_path, index=True, index_label="SID")

    print_text(log_file, "Finished Extracting VCF")
