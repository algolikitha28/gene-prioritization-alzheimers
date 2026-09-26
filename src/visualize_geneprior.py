import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------
# 1. File paths
# --------------------------------------------------

INPUT_FILE = "data/processed/geneprior_ranked.csv"
OUTPUT_FILE = "data/processed/geneprior_top20.png"


# --------------------------------------------------
# 2. Read GenePrior ranking
# --------------------------------------------------

print("Reading GenePrior ranking...")

df = pd.read_csv(INPUT_FILE)

print(f"Total genes: {len(df)}")


# --------------------------------------------------
# 3. Select top 20 genes
# --------------------------------------------------

top20 = df.head(20).copy()

print("\nTop 20 genes:")
print(
    top20[
        [
            "geneprior_rank",
            "gene",
            "geneprior_score_100"
        ]
    ].to_string(index=False)
)


# --------------------------------------------------
# 4. Create plot
# --------------------------------------------------

plt.figure(figsize=(12, 7))

plt.barh(
    top20["gene"][::-1],
    top20["geneprior_score_100"][::-1]
)

plt.xlabel("GenePrior Score (0–100)")
plt.ylabel("Gene")
plt.title("Top 20 Genes Ranked by GenePrior Score")

plt.tight_layout()


# --------------------------------------------------
# 5. Save figure
# --------------------------------------------------

plt.savefig(OUTPUT_FILE, dpi=300)

print(f"\nFigure saved to: {OUTPUT_FILE}")

plt.show()