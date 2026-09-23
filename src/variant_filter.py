import pandas as pd
import numpy as np

# --------------------------------------------------
# File paths
# --------------------------------------------------

input_file = "data/raw/gwas_ad_raw.csv"
output_file = "data/processed/variants_filtered.csv"

# --------------------------------------------------
# Load GWAS data
# --------------------------------------------------

df = pd.read_csv(input_file)

print("Original shape:", df.shape)

# --------------------------------------------------
# Basic validation
# --------------------------------------------------

required_columns = [
    "gene",
    "variant",
    "p_value",
    "trait",
    "study",
    "pubmed_id",
    "author"
]

missing_columns = [
    col for col in required_columns
    if col not in df.columns
]

if missing_columns:
    raise ValueError(
        f"Missing required columns: {missing_columns}"
    )

# --------------------------------------------------
# Remove rows with missing essential information
# --------------------------------------------------

df = df.dropna(
    subset=["gene", "variant", "p_value"]
).copy()

# --------------------------------------------------
# Keep only valid positive p-values
# --------------------------------------------------

df = df[df["p_value"] > 0].copy()

# --------------------------------------------------
# Classify statistical significance
# --------------------------------------------------

def classify_pvalue(p):
    if p < 5e-8:
        return "genome_wide_significant"
    elif p < 1e-5:
        return "suggestive"
    else:
        return "weak"


df["significance_class"] = df["p_value"].apply(
    classify_pvalue
)

# --------------------------------------------------
# Calculate -log10(p)
# --------------------------------------------------

df["minus_log10_p"] = -np.log10(df["p_value"])

# --------------------------------------------------
# Save processed variant table
# --------------------------------------------------

df.to_csv(output_file, index=False)

# --------------------------------------------------
# Summary
# --------------------------------------------------

print("\nProcessed shape:", df.shape)

print("\nSignificance classes:")
print(
    df["significance_class"]
    .value_counts()
)

print("\nUnique variants:",
      df["variant"].nunique())

print("Unique genes:",
      df["gene"].nunique())

print("\nTop 10 strongest records:")
print(
    df.sort_values(
        "p_value"
    )[
        [
            "gene",
            "variant",
            "p_value",
            "minus_log10_p",
            "significance_class"
        ]
    ].head(10).to_string(index=False)
)

print(
    f"\nSaved to: {output_file}"
)