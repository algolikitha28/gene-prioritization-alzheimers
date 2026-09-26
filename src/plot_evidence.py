import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------
# 1. File paths
# --------------------------------------------------

INPUT_FILE = "data/processed/geneprior_ranked.csv"
OUTPUT_FILE = "data/processed/geneprior_gwas_vs_functional.png"


# --------------------------------------------------
# 2. Read data
# --------------------------------------------------

print("Reading GenePrior ranking...")

df = pd.read_csv(INPUT_FILE)

print(f"Total genes: {len(df)}")


# --------------------------------------------------
# 3. Create scatter plot
# --------------------------------------------------

plt.figure(figsize=(10, 7))

plt.scatter(
    df["gwas_score"],
    df["functional_score"],
    alpha=0.7
)


# --------------------------------------------------
# 4. Important genes to label
# --------------------------------------------------

genes_to_label = [
    "BIN1",
    "TREM2",
    "IL34",
    "ARHGAP45",
    "APH1B",
    "PLCG2",
    "MME"
]


# --------------------------------------------------
# 5. Add labels
# --------------------------------------------------

for _, row in df[df["gene"].isin(genes_to_label)].iterrows():

    plt.annotate(
        row["gene"],
        (row["gwas_score"], row["functional_score"]),
        xytext=(6, 6),
        textcoords="offset points"
    )


# --------------------------------------------------
# 6. Labels and title
# --------------------------------------------------

plt.xlabel("GWAS Score")
plt.ylabel("Functional Score")

plt.title(
    "GWAS Evidence vs Functional Evidence"
)


# --------------------------------------------------
# 7. Save figure
# --------------------------------------------------

plt.tight_layout()

plt.savefig(
    OUTPUT_FILE,
    dpi=300
)

print(f"\nFigure saved to: {OUTPUT_FILE}")

plt.show()