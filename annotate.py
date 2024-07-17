import subprocess
import os


def main(input_path, output_path, snpsift, dbsnp, dbnsfp, clinvar):
    print("Annotating VCF...")
    basename = "annotated.vcf"
    if dbsnp is not None:
        print("Annotating with dbsnp...")
        basename = "dbsnp_" + basename
        output = os.path.join(output_path, basename)
        with open(output, "w") as f:
            subprocess.run(
                ["java", "-jar", snpsift, "annotate", dbsnp, input_path],
                stdout=f,
            )
        print("Saved to " + output)
        input_path = output

    if dbnsfp is not None:
        print("Annotating with dbnsfp...")
        basename = "dbnsfp_" + basename
        output = os.path.join(output_path, basename)
        with open(output, "w") as f:
            subprocess.run(
                ["java", "-jar", snpsift, "dbnsfp", "-db", dbnsfp, input_path], stdout=f
            )
        print("Saved to " + output)
        input_path = output

    if clinvar is not None:
        print("Annotating with clinvar...")
        basename = "clinvar_" + basename
        output = os.path.join(output_path, basename)
        with open(output, "w") as f:
            subprocess.run(
                ["java", "-jar", snpsift, "annotate", clinvar, input_path], stdout=f
            )
        print("Saved to " + output)
        input_path = output

    print("Finished Annotating VCF")
    return input_path