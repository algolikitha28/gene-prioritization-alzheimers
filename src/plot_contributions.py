import pandas as pd
import matplotlib.pyplot as plt


# --------------------------------------------------
# 1. Paths
# --------------------------------------------------

INPUT_FILE = "data/processed/geneprior_ranked.csv"
OUTPUT_FILE = "data/processed/geneprior_contributions_top20.png"


# --------------------------------------------------
# 2. Read ranked GenePrior table
# --------------------------------------------------

print("Reading GenePrior ranking...")

df = pd.read_csv(INPUT_FILE)

print(f"Total genes: {len(df)}")


# --------------------------------------------------
# 3. Select top 20
# --------------------------------------------------

top20 = df.head(20).copy()

# Reverse order so rank 1 appears at the top
top20 = top20.iloc[::-1]


# --------------------------------------------------
# 4. Create figure
# --------------------------------------------------

plt.figure(figsize=(12, 8))


# --------------------------------------------------
# 5. GWAS contribution
# --------------------------------------------------

plt.barh(
    top20["gene"],
    top20["gwas_contribution"],
    label="GWAS contribution"
)


# --------------------------------------------------
# 6. Functional contribution
# --------------------------------------------------

plt.barh(
    top20["gene"],
    top20["functional_contribution"],
    left=top20["gwas_contribution"],
    label="Functional contribution"
)


# --------------------------------------------------
# 7. Labels
# --------------------------------------------------

plt.xlabel("Contribution to GenePrior Score")

plt.ylabel("Gene")

plt.title(
    "GenePrior Score Contributions: Top 20 Genes"
)

plt.legend()


# --------------------------------------------------
# 8. Save
# --------------------------------------------------

plt.tight_layout()

plt.savefig(
    OUTPUT_FILE,
    dpi=300
)

print(f"\nFigure saved to: {OUTPUT_FILE}")

plt.show()