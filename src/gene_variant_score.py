import pandas as pd

input_file = "data/processed/variants_filtered.csv"
output_file = "data/processed/gene_variant_scores.csv"

# Load filtered variant-level data
df = pd.read_csv(input_file)

print("Input shape:", df.shape)

# Aggregate variant-level evidence to gene level
gene_df = df.groupby("gene").agg(
    variant_count=("variant", "nunique"),

    significant_variant_count=(
        "significance_class",
        lambda x: (x == "genome_wide_significant").sum()
    ),

    best_p_value=("p_value", "min"),

    max_minus_log10_p=("minus_log10_p", "max")
).reset_index()


# Min-Max normalization function
def min_max_normalize(series):
    min_value = series.min()
    max_value = series.max()

    if max_value == min_value:
        return pd.Series(0.0, index=series.index)

    return (series - min_value) / (max_value - min_value)


# Normalize statistical strength
gene_df["normalized_p_score"] = min_max_normalize(
    gene_df["max_minus_log10_p"]
)

# Normalize significant variant support
gene_df["normalized_variant_score"] = min_max_normalize(
    gene_df["significant_variant_count"]
)


# Combine the two evidence components
gene_df["gwas_variant_score"] = (
    0.7 * gene_df["normalized_p_score"]
    +
    0.3 * gene_df["normalized_variant_score"]
)


# Rank genes
gene_df = gene_df.sort_values(
    "gwas_variant_score",
    ascending=False
)


# Save results
gene_df.to_csv(output_file, index=False)


print("\nGene-level shape:", gene_df.shape)

print("\nTop 20 genes:")
print(
    gene_df.head(20).to_string(index=False)
)

print(
    f"\nSaved to: {output_file}"
)