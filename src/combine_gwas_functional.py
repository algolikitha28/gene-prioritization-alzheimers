import pandas as pd


GWAS_FILE = "data/processed/variants_filtered.csv"
FUNCTIONAL_FILE = "data/processed/variant_functional_summary.csv"
OUTPUT_FILE = "data/processed/gwas_functional_annotations.csv"


print("Reading GWAS data...")

gwas = pd.read_csv(GWAS_FILE)

print(f"GWAS rows: {len(gwas)}")
print(f"Unique GWAS variants: {gwas['variant'].nunique()}")
print(f"Unique GWAS genes: {gwas['gene'].nunique()}")


print()
print("Reading functional annotation summary...")

functional = pd.read_csv(FUNCTIONAL_FILE)

print(f"Functional rows: {len(functional)}")
print(f"Unique functional variants: {functional['variant'].nunique()}")


# Keep only the functional columns needed for the merge
functional_selected = functional[
    [
        "variant",
        "most_severe_consequence",
        "functional_impact",
        "annotation_count"
    ]
]


print()
print("Merging GWAS and functional annotations...")


merged = gwas.merge(
    functional_selected,
    on="variant",
    how="left"
)


merged.to_csv(
    OUTPUT_FILE,
    index=False
)


print()
print("Combined table created.")
print(f"Rows: {len(merged)}")
print(f"Columns: {merged.columns.tolist()}")


print()
print("Missing functional annotations:")

print(
    merged["most_severe_consequence"]
    .isna()
    .sum()
)


print()
print("Functional impact distribution:")

print(
    merged["functional_impact"]
    .value_counts(dropna=False)
    .to_string()
)


print()
print("First 20 rows:")

print(
    merged[
        [
            "gene",
            "variant",
            "p_value",
            "most_severe_consequence",
            "functional_impact",
            "annotation_count"
        ]
    ]
    .head(20)
    .to_string(index=False)
)