import pandas as pd


INPUT_FILE = "data/processed/gwas_functional_annotations.csv"
OUTPUT_FILE = "data/processed/gene_functional_summary.csv"


print("Reading GWAS + functional annotation table...")

df = pd.read_csv(INPUT_FILE)

print(f"Rows: {len(df)}")
print(f"Unique genes: {df['gene'].nunique()}")
print(f"Unique variants: {df['variant'].nunique()}")


# --------------------------------------------------
# Count unique variants per gene
# --------------------------------------------------

variant_counts = (
    df.groupby("gene")["variant"]
    .nunique()
    .rename("variant_count")
)


# --------------------------------------------------
# Count functional impact categories per gene
# --------------------------------------------------

impact_counts = pd.crosstab(
    df["gene"],
    df["functional_impact"]
)


# Make sure all four columns exist
for impact in ["HIGH", "MODERATE", "LOW", "MODIFIER"]:

    if impact not in impact_counts.columns:
        impact_counts[impact] = 0


impact_counts = impact_counts[
    ["HIGH", "MODERATE", "LOW", "MODIFIER"]
]


# --------------------------------------------------
# Combine
# --------------------------------------------------

summary = variant_counts.to_frame()

summary = summary.join(
    impact_counts,
    how="left"
)


# --------------------------------------------------
# Find the most severe consequence observed
# --------------------------------------------------

severity_order = {
    "stop_gained": 13,
    "frameshift_variant": 12,
    "splice_donor_variant": 11,
    "splice_polypyrimidine_tract_variant": 10,
    "inframe_deletion": 9,
    "inframe_insertion": 8,
    "missense_variant": 7,
    "synonymous_variant": 6,
    "3_prime_UTR_variant": 5,
    "5_prime_UTR_variant": 4,
    "non_coding_transcript_exon_variant": 3,
    "intron_variant": 2,
    "upstream_gene_variant": 1,
    "downstream_gene_variant": 1,
    "intergenic_variant": 0
}


def get_most_severe_consequence(group):

    consequences = (
        group["most_severe_consequence"]
        .dropna()
        .astype(str)
        .tolist()
    )

    if not consequences:
        return None

    return max(
        consequences,
        key=lambda x: severity_order.get(x, -1)
    )


most_severe = (
    df.groupby("gene")
    .apply(get_most_severe_consequence)
    .rename("most_severe_functional_consequence")
)


summary = summary.join(
    most_severe,
    how="left"
)


summary = summary.reset_index()


summary.to_csv(
    OUTPUT_FILE,
    index=False
)


print()
print("Gene-level functional summary created.")

print(f"Rows: {len(summary)}")

print(
    f"Columns: {summary.columns.tolist()}"
)

print()
print("First 20 genes:")

print(
    summary.head(20)
    .to_string(index=False)
)