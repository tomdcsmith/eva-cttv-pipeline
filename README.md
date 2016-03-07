## README ##

Building
-------
1. "git clone git@github.com:EBIvariation/eva-cttv-pipeline.git"
2. "cd eva-cttv-pipeline"
3. and then:
	4. to install: "python3 setup.py install"
	5. to install to develop: "python3 setup.py develop"
	6. to build a source distribution: "python3 setup.py sdist"

Setting up virtual environment
-------
For a Python virtual environment to work with the pipeline:

1. "cd eva-cttv-pipeline"
2. "virtualenv -p python3.5 venv"
3. "source venv/bin/activate" ("venv/bin/deactivate" to deactivate virtualenv)
4. pip install -r /path/to/requirements.txt

Usage
-------

python3 clinvar_to_evidence_strings.py --out \<OUTPUT_FILE\> -e \<EFO_MAPPING_FILE\> -g \<SNP_2_GENE_MAPPING_FILE\> -v \<VARIANT_SUMMARY_FILE\> [--clinSig \<CLINICAL_SIGNIFICANCE_LIST\>] [--ignore \<TERM_URL_IGNORE_FILE\>]

OUTPUT_FILE: path to output directory

EFO_MAPPING_FILE: file with mappings of Clinvar trait names to URLS

SNP_2_GENE_MAPPING_FILE: file from CTTV with mappings from rs IDs to Ensembl gene IDs, and to functional consequences

VARIANT_SUMMARY_FILE: variant_summary file from Clinvar

CLINICAL_SIGNIFICANCE_LIST: comma separated (no spaces) list of clinical significances allowed to generate evidence strings (defaults to "pathogenic,likely pathogenic")

TERM_URL_IGNORE_FILE: path to file containing list of invalid EFO URLs which will be used to filter out evidence strings with matching URLs; file should contain list of URLs, one on each line
