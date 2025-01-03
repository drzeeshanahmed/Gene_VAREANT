import traceback
import sys
import argparse
import os
from datetime import datetime
from log import print_text
from truncate import main as truncate_pipeline
from annotate import main as annotate_pipeline
from extract import main as extract_pipeline

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(dest="command")
    truncate = subparsers.add_parser("truncate")
    annotate = subparsers.add_parser("annotate")
    extract = subparsers.add_parser("extract")
    all = subparsers.add_parser("all")

    # Truncate
    truncate.add_argument(
        "--config",
        required=True,
        help="Configuration rulesets for how to truncate the input VCF.",
    )
    truncate.add_argument(
        "--input", required=True, help="Input VCF. Must be compliant with VCF v4.2"
    )
    truncate.add_argument(
        "--output", required=True, help="Output directory to store filtered VCFs."
    )
    truncate.add_argument(
        "--preserve-meta",
        help="Output directory to store filtered VCFs.",
        action="store_true",
    )

    # Annotate
    annotate.add_argument(
        "--input",
        required=True,
        help="Input VCF that gets annotated with output files.",
    )
    annotate.add_argument(
        "--output", required=True, help="Output directory to store annotated VCFs."
    )
    annotate.add_argument(
        "--snpsift", required=True, help="Path to SnpSift.jar executable"
    )
    annotate.add_argument(
        "--dbsnp",
        help="Path to dbsnp vcf.gz file. A corresponding tabix-indexed file must exist in the same folder.",
    )
    annotate.add_argument(
        "--dbnsfp",
        help="Path to dbnsfp vcf.gz file. A corresponding tabix-indexed file must exist in the same folder.",
    )
    annotate.add_argument(
        "--clinvar",
        help="Path to clinvar vcf.gz file. A corresponding tabix-indexed file must exist in the same folder.",
    )

    # Extract
    extract.add_argument(
        "--config",
        required=True,
        help="Configuration rulesets for how to truncate the input VCF.",
    )
    extract.add_argument(
        "--input", required=True, help="Input VCF that gets extracted into DB and CIGT."
    )
    extract.add_argument(
        "--output",
        required=True,
        help="Output directory to store SQLite3 DB and CIGT formatted file.",
    )

    # Combined
    all.add_argument(
        "--input",
        required=True,
        help="Input VCF File",
    )
    all.add_argument(
        "--output",
        required=True,
        help="Output directory to store results.",
    )
    all.add_argument(
        "--config",
        required=True,
        help="Configuration rulesets for how to truncate the input VCF.",
    )
    all.add_argument(
        "--preserve-meta",
        help="Output directory to store filtered VCFs.",
        action="store_true",
    )
    all.add_argument("--snpsift", required=True, help="Path to SnpSift.jar executable")
    all.add_argument(
        "--dbsnp",
        help="Path to dbsnp vcf.gz file. A corresponding tabix-indexed file must exist in the same folder.",
    )
    all.add_argument(
        "--dbnsfp",
        help="Path to dbnsfp vcf.gz file. A corresponding tabix-indexed file must exist in the same folder.",
    )
    all.add_argument(
        "--clinvar",
        help="Path to clinvar vcf.gz file. A corresponding tabix-indexed file must exist in the same folder.",
    )

    args = parser.parse_args()

    if os.path.exists(args.output):
        args.output = os.path.join(
            args.output, datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        )

    log_path = os.path.join(args.output, "log.txt")
    os.mkdir(args.output)

    log_file = open(log_path, "a")

    print_text(log_file, "Starting VAREANT with " + " ".join(sys.argv))

    try:
        if args.command == "truncate":
            truncate_pipeline(
                args.config,
                args.input,
                args.output,
                args.preserve_meta,
                log_file,
            )
        elif args.command == "annotate":
            annotate_pipeline(
                args.input,
                args.output,
                args.snpsift,
                args.dbsnp,
                args.dbnsfp,
                args.clinvar,
                log_file,
            )
        elif args.command == "extract":
            extract_pipeline(
                args.config,
                args.input,
                args.output,
                log_file,
            )
        elif args.command == "all":
            input = truncate_pipeline(
                args.config,
                args.input,
                args.output,
                args.preserve_meta,
                log_file,
            )
            input = annotate_pipeline(
                input,
                args.output,
                args.snpsift,
                args.dbsnp,
                args.dbnsfp,
                args.clinvar,
                log_file,
            )
            extract_pipeline(
                args.config,
                input,
                args.output,
                log_file,
            )
    except Exception as e:
        print_text(log_file, "Execution failed with Exception:")
        print_text(log_file, "".join(traceback.format_exception(e)))

    print_text(log_file, "Finished VAREANT")
    log_file.flush()
    log_file.close()
