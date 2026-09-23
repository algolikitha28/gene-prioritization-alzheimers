import pandas as pd


INPUT_FILE = "data/processed/variant_annotations.csv"
OUTPUT_FILE = "data/processed/variant_functional_summary.csv"


# Impact severity order
IMPACT_ORDER = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1
}


print("Reading VEP annotation table...")

df = pd.read_csv(INPUT_FILE)

print(f"Annotation rows: {len(df)}")
print(f"Unique variants: {df['variant'].nunique()}")


summary_rows = []


for variant, group in df.groupby("variant"):

    # VEP's most severe consequence for this variant
    consequences = (
        group["most_severe_consequence"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    if len(consequences) > 0:
        most_severe_consequence = consequences[0]
    else:
        most_severe_consequence = None


    # Find rows associated with the most severe consequence
    matching_rows = group[
        group["consequence"].fillna("").str.contains(
            str(most_severe_consequence),
            regex=False
        )
    ]


    # Determine impact
    impacts = (
        matching_rows["impact"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    if not impacts:
        impacts = (
            group["impact"]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )


    if impacts:
        functional_impact = max(
            impacts,
            key=lambda x: IMPACT_ORDER.get(x, 0)
        )
    else:
        functional_impact = None


    summary_rows.append({
        "variant": variant,
        "most_severe_consequence": most_severe_consequence,
        "functional_impact": functional_impact,
        "annotation_count": len(group)
    })


summary_df = pd.DataFrame(summary_rows)


summary_df = summary_df.sort_values(
    by="variant"
).reset_index(drop=True)


summary_df.to_csv(
    OUTPUT_FILE,
    index=False
)


print()
print("Functional summary created.")
print(f"Rows: {len(summary_df)}")
print(f"Columns: {summary_df.columns.tolist()}")
print()
print("Impact distribution:")
print(
    summary_df["functional_impact"]
    .value_counts(dropna=False)
    .to_string()
)
print()
print("Most severe consequence distribution:")
print(
    summary_df["most_severe_consequence"]
    .value_counts(dropna=False)
    .to_string()
)
print()
print("First 20 variants:")
print(
    summary_df.head(20)
    .to_string(index=False)
)