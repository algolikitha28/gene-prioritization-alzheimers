import pandas as pd

INPUT_FILE = "data/processed/geneprior_ranked.csv"

print("Reading GenePrior ranking...")
df = pd.read_csv(INPUT_FILE)

print("\n========== BASIC VALIDATION ==========")

print(f"Total rows: {len(df)}")
print(f"Unique genes: {df['gene'].nunique()}")

print("\nMissing values:")
print(df.isnull().sum())

print("\nDuplicate genes:")
print(df['gene'].duplicated().sum())

print("\nRank range:")
print(f"Minimum rank: {df['geneprior_rank'].min()}")
print(f"Maximum rank: {df['geneprior_rank'].max()}")

print("\nScore range:")
print(f"Minimum GenePrior score: {df['geneprior_score'].min()}")
print(f"Maximum GenePrior score: {df['geneprior_score'].max()}")

print("\nTop 10 genes:")
print(
    df[
        [
            "geneprior_rank",
            "gene",
            "geneprior_score",
            "gwas_score",
            "functional_score",
            "gwas_contribution",
            "functional_contribution"
        ]
    ]
    .head(10)
    .to_string(index=False)
)

print("\n========== MATHEMATICAL VALIDATION ==========")

# Check 1: GWAS contribution
expected_gwas = 0.70 * df["gwas_score"]

gwas_check = (
    (df["gwas_contribution"] - expected_gwas).abs() < 1e-10
).all()

print(f"GWAS contribution formula correct: {gwas_check}")


# Check 2: Functional contribution
expected_functional = 0.30 * df["functional_score"]

functional_check = (
    (df["functional_contribution"] - expected_functional).abs() < 1e-10
).all()

print(f"Functional contribution formula correct: {functional_check}")


# Check 3: Final GenePrior score
expected_geneprior = (
    df["gwas_contribution"]
    + df["functional_contribution"]
)

geneprior_check = (
    (df["geneprior_score"] - expected_geneprior).abs() < 1e-10
).all()

print(f"GenePrior formula correct: {geneprior_check}")


# Check 4: Score × 100
expected_score_100 = df["geneprior_score"] * 100

score_100_check = (
    (df["geneprior_score_100"] - expected_score_100).abs() < 1e-10
).all()

print(f"Score × 100 correct: {score_100_check}")


# Check 5: Descending order
ranking_check = df["geneprior_score"].is_monotonic_decreasing

print(f"Ranking sorted correctly: {ranking_check}")


# Check 6: Score range
score_range_check = (
    (df["geneprior_score"] >= 0)
    & (df["geneprior_score"] <= 1)
).all()

print(f"All GenePrior scores between 0 and 1: {score_range_check}")


# Check 7: Tied scores
ties = df[
    df["geneprior_score"].duplicated(keep=False)
].sort_values("geneprior_score", ascending=False)

print("\nTied GenePrior scores:")

if len(ties) == 0:
    print("No ties found.")
else:
    print(
        ties[
            ["geneprior_rank", "gene", "geneprior_score"]
        ].to_string(index=False)
    )