import subprocess
import os
import io
import pysam

from log import print_text


def main(input_path, output_path, snpsift, dbsnp, dbnsfp, clinvar, log_file):
    print_text(log_file, "Annotating VCF")

    if dbsnp is None and dbnsfp is None and clinvar is None:
        print_text(
            log_file, "No annotation databases were specified. Skipping annotation..."
        )
        return input_path

    if input_path.endswith(".vcf.gz"):
        basename = "annotated.vcf.gz"
    elif input_path.endswith(".vcf"):
        basename = "annotated.vcf"
    else:
        raise ValueError(
            "Must provide a compressed (.vcf.gz) or uncompressed (.vcf) input file."
        )

    if dbsnp is not None:
        print_text(log_file, "Annotating with dbSNP at " + dbsnp)
        basename = "dbsnp_" + basename
        output = os.path.join(output_path, basename)
        print_text(log_file, "Saving to " + output)
        print_text(
            log_file,
            " ".join(["java", "-jar", snpsift, "annotate", dbsnp, input_path]),
        )

        with (
            io.TextIOWrapper(pysam.BGZFile(output, "w"), "utf-8")
            if output.endswith(".vcf.gz")
            else open(output, "w")
        ) as f:
            process = subprocess.Popen(
                ["java", "-jar", snpsift, "annotate", dbsnp, input_path],
                stdout=subprocess.PIPE,
                text=True,
            )
            for line in process.stdout:
                f.write(line)

            process.stdout.close()
            process.wait()

        print_text(log_file, "Finished dbSNP")
        input_path = output

    if dbnsfp is not None:
        print_text(log_file, "Annotating with dbNSFP at " + dbnsfp)
        basename = "dbnsfp_" + basename
        output = os.path.join(output_path, basename)
        print_text(log_file, "Saving to " + output)
        print_text(
            log_file,
            " ".join(["java", "-jar", snpsift, "dbnsfp", "-db", dbnsfp, input_path]),
        )

        with (
            io.TextIOWrapper(pysam.BGZFile(output, "w"), "utf-8")
            if output.endswith(".vcf.gz")
            else open(output, "w")
        ) as f:
            process = subprocess.Popen(
                ["java", "-jar", snpsift, "dbnsfp", "-db", dbnsfp, input_path],
                stdout=subprocess.PIPE,
                text=True,
            )
            for line in process.stdout:
                f.write(line)

            process.stdout.close()
            process.wait()

        print_text(log_file, "Finished dbNSFP")
        input_path = output

    if clinvar is not None:
        print_text(log_file, "Annotating with ClinVar at " + clinvar)
        basename = "clinvar_" + basename
        output = os.path.join(output_path, basename)
        print_text(log_file, "Saving to " + output)
        print_text(
            log_file,
            " ".join(["java", "-jar", snpsift, "annotate", clinvar, input_path]),
        )

        with (
            io.TextIOWrapper(pysam.BGZFile(output, "w"), "utf-8")
            if output.endswith(".vcf.gz")
            else open(output, "w")
        ) as f:
            process = subprocess.Popen(
                ["java", "-jar", snpsift, "annotate", clinvar, input_path],
                stdout=subprocess.PIPE,
                text=True,
            )
            for line in process.stdout:
                f.write(line)

            process.stdout.close()
            process.wait()

        print_text(log_file, "Finished ClinVar")
        input_path = output

    print_text(log_file, "Finished Annotating VCF")

    return input_path
