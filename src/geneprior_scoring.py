
import pandas as pd
from pathlib import Path


# --------------------------------------------------
# Paths
# --------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent

GWAS_FILE = BASE_DIR / "data" / "processed" / "gene_variant_scores.csv"
FUNCTIONAL_FILE = BASE_DIR / "data" / "processed" / "gene_functional_summary.csv"

OUTPUT_FILE = BASE_DIR / "data" / "processed" / "geneprior_ranked.csv"


# --------------------------------------------------
# Read input files
# --------------------------------------------------

print("Reading GenePrior evidence tables...")

gwas = pd.read_csv(GWAS_FILE)
functional = pd.read_csv(FUNCTIONAL_FILE)

print("GWAS table shape:", gwas.shape)
print("Functional table shape:", functional.shape)


# --------------------------------------------------
# Calculate functional score
# --------------------------------------------------

print("\nCalculating functional scores...")

impact_weights = {
    "HIGH": 1.0,
    "MODERATE": 0.7,
    "LOW": 0.3,
    "MODIFIER": 0.1
}

functional["functional_score_raw"] = (
    functional["HIGH"] * impact_weights["HIGH"]
    + functional["MODERATE"] * impact_weights["MODERATE"]
    + functional["LOW"] * impact_weights["LOW"]
    + functional["MODIFIER"] * impact_weights["MODIFIER"]
)

functional["functional_score"] = (
    functional["functional_score_raw"]
    / functional["variant_count"]
)

functional["functional_score"] = (
    functional["functional_score"]
    .clip(0, 1)
)


# --------------------------------------------------
# Select GWAS evidence
# --------------------------------------------------

gwas_evidence = gwas[
    [
        "gene",
        "gwas_variant_score",
        "best_p_value",
        "variant_count",
        "significant_variant_count"
    ]
].copy()

gwas_evidence = gwas_evidence.rename(
    columns={
        "gwas_variant_score": "gwas_score"
    }
)


# --------------------------------------------------
# Select functional evidence
# --------------------------------------------------

functional_evidence = functional[
    [
        "gene",
        "functional_score",
        "most_severe_functional_consequence",
        "HIGH",
        "MODERATE",
        "LOW",
        "MODIFIER"
    ]
].copy()


# --------------------------------------------------
# Merge evidence
# --------------------------------------------------

print("\nMerging GWAS and functional evidence...")

evidence = pd.merge(
    gwas_evidence,
    functional_evidence,
    on="gene",
    how="outer"
)


# --------------------------------------------------
# Check missing values
# --------------------------------------------------

print("\nMissing values after merge:")

print(evidence.isnull().sum())


# --------------------------------------------------
# Calculate GenePrior v1 score
# --------------------------------------------------

print("\nCalculating GenePrior v1 score...")

GWAS_WEIGHT = 0.70
FUNCTIONAL_WEIGHT = 0.30

evidence["geneprior_score"] = (
    GWAS_WEIGHT * evidence["gwas_score"]
    + FUNCTIONAL_WEIGHT * evidence["functional_score"]
)

# Convert score to 0–100 scale
evidence["geneprior_score_100"] = (
    evidence["geneprior_score"] * 100
)


# --------------------------------------------------
# Rank genes
# --------------------------------------------------

print("\nRanking genes...")

evidence = evidence.sort_values(
    by="geneprior_score",
    ascending=False
).reset_index(drop=True)

evidence["geneprior_rank"] = (
    evidence.index + 1
)


# --------------------------------------------------
# Reorder columns
# --------------------------------------------------

column_order = [
    "gene",
    "geneprior_rank",
    "geneprior_score",
    "geneprior_score_100",
    "gwas_score",
    "functional_score",
    "best_p_value",
    "variant_count",
    "significant_variant_count",
    "most_severe_functional_consequence",
    "HIGH",
    "MODERATE",
    "LOW",
    "MODIFIER"
]

evidence = evidence[column_order]


# --------------------------------------------------
# Save final ranked table
# --------------------------------------------------

evidence.to_csv(
    OUTPUT_FILE,
    index=False
)


# --------------------------------------------------
# Summary
# --------------------------------------------------

print("\nGenePrior v1 ranking created.")

print("Shape:", evidence.shape)

print("\nColumns:")
print(evidence.columns.tolist())

print("\nTop 20 GenePrior genes:")

print(
    evidence.head(20).to_string(index=False)
)

print("\nSaved to:")
print(OUTPUT_FILE)

