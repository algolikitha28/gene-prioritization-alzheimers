# src/config.py

# Disease we are studying
DISEASE_NAME = "Alzheimer disease"

# Genome-wide significance threshold
GWAS_P_VALUE_THRESHOLD = 5e-8

# GWAS Catalog API
GWAS_BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2"

# Number of records requested per API page
PAGE_SIZE = 100

# Output files
RAW_GWAS_FILE = "data/raw/gwas_ad_raw.csv"
PROCESSED_GWAS_FILE = "data/processed/gwas_gene_scores.csv"