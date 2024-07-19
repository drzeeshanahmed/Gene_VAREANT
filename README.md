# Gene_VAREANT

## Setup

This pipeline requires git, Python, and Java. Please verify that they are installed:

```bash
git --version
python --version # >= 2.7
java --version # >= 11
```

First clone this repository using git:

```bash
git clone https://github.com/drzeeshanahmed/Gene_VAREANT
```

OPTIONAL: Set up a virtual environment to install all python dependencies. Note that
the `venv` may not exist for earlier version of python (before 3.6). If this is the case, you
can skip this step of creating a virtual environment.

```bash
python -m venv venv
source venv/bin/activate
```

Install python dependencies:

```bash
pip install -r requirements.txt
```

## Third-party dependencies

The annotate phase of this pipeline requires that `SnpSift.jar` is available on system,
as well as any of the supported annotation databases. To install these, please refer
to their documentations. As the annotation databases can be quite large (tens of GB), please
install only those which you require.

1. SnpSift can be downloaded here: https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
1. ClinVar and its index file can be downloaded at https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
1. dbSNP and its index file can be downloaded from https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/
1. dbnsfp and its index file can be downloaded following the instructions at https://pcingola.github.io/SnpEff/snpsift/dbnsfp/

## Pipeline

The main.py file contains all the utilities to execute the analysis

### Truncation

```bash
python main.py truncate --input {INPUT_VCF} --output {OUTPUT_DIR} --config {CONFIG_FILE}
```

Replace `INPUT_VCF` with the path to the VCF file to filter, `OUTPUT_DIR` with the directory to store the
resulting `truncated.vcf` file, and `CONFIG_FILE` with the configuration file as specified in the user guide.

### Annotation

```bash
python main.py annotate --input {INPUT_VCF} --output {OUTPUT_DIR} --snpsift {SNPSIFT_JAR} --dbsnp {DBSNP_PATH} --dbnsfp {DBNSFP_PATH} --clinvar {CLINVAR_PATH}
```

Replace `INPUT_VCF` with the path to the VCF file to annotate, `OUTPUT_DIR` with the directory to store the
resulting annotated VCF file, `SNPSIFT_JAR` with the path the to SnpSift.jar executable,
`DBSNP_PATH` to the vcf.gz file path, `DBNSFP_PATH` to the txt.gz file path, `CLIVNAR_PATH` to the vcf.gz file path.

Note that for each of `DBSNP_PATH`, `DBNSFP_PATH`, and `CLINVAR_PATH`, those flags are optional and ignored if omitted.
If included, a corresponding tabix-indexed file with the same basename (`{*}.tbi`) must exist in the same directory. Please refer
to the SnpSift documentation for this.

### Extraction

```bash
python main.py extract --input {INPUT_VCF} --output {OUTPUT_DIR}
```

Replace `INPUT_VCF` with the path to the VCF file to extract, `OUTPUT_DIR` with the directory to store the
resulting `variant_db.sqlite` and `cigt_matrix.cigt.csv` files.

### All

Lastly, as a convenience script to execute all three stages together and save results to the same directory, you can execute the following command:

```bash
python main.py all --input {INPUT_VCF} --output {OUTPUT_DIR} --config {CONFIG_FILE} --snpsift {SNPSIFT_JAR} --dbsnp {DBSNP_PATH} --dbnsfp {DBNSFP_PATH} --clinvar {CLINVAR_PATH}
```

Refer to previous commands for a description of each flag.
