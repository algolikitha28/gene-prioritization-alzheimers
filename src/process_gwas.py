# src/process_gwas.py

import pandas as pd
import numpy as np

from config import (
    GWAS_P_VALUE_THRESHOLD,
    RAW_GWAS_FILE,
    PROCESSED_GWAS_FILE
)


def load_data():

    print("Loading raw GWAS data...")

    df = pd.read_csv(RAW_GWAS_FILE)

    print(f"Rows loaded: {len(df)}")

    return df


def clean_data(df):

    print("\n--- Cleaning data ---")

    # Remove rows without gene
    before = len(df)

    df = df.dropna(
        subset=["gene"]
    )

    print(
        f"Removed {before - len(df)} rows without gene."
    )

    # Convert p-value to numeric
    df["p_value"] = pd.to_numeric(
        df["p_value"],
        errors="coerce"
    )

    # Remove rows where p-value is missing
    before = len(df)

    df = df.dropna(
        subset=["p_value"]
    )

    print(
        f"Removed {before - len(df)} rows without p-value."
    )

    # Keep valid p-values
    before = len(df)

    df = df[
        (df["p_value"] > 0) &
        (df["p_value"] <= 1)
    ]

    print(
        f"Removed {before - len(df)} invalid p-values."
    )

    # Clean gene names
    df["gene"] = (
        df["gene"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    # Remove duplicate records
    before = len(df)

    df = df.drop_duplicates()

    print(
        f"Removed {before - len(df)} duplicate rows."
    )

    return df


def filter_significant(df):

    print("\n--- Filtering significant associations ---")

    significant = df[
        df["p_value"] <= GWAS_P_VALUE_THRESHOLD
    ].copy()

    print(
        f"Significant associations: {len(significant)}"
    )

    print(
        f"Threshold: {GWAS_P_VALUE_THRESHOLD}"
    )

    return significant


def calculate_gene_summary(df):

    print("\n--- Creating gene-level summary ---")

    gene_summary = (
        df
        .groupby("gene")
        .agg(
            best_p_value=(
                "p_value",
                "min"
            ),

            significant_variants=(
                "variant",
                "nunique"
            )
        )
        .reset_index()
    )

    return gene_summary


def calculate_gwas_strength(df):

    print("\n--- Calculating GWAS strength ---")

    df["gwas_strength"] = (
        -np.log10(
            df["best_p_value"]
        )
    )

    return df


def normalize_score(df):

    print("\n--- Normalizing GWAS score ---")

    minimum = df["gwas_strength"].min()

    maximum = df["gwas_strength"].max()

    if maximum == minimum:

        df["gwas_score"] = 1.0

    else:

        df["gwas_score"] = (
            (
                df["gwas_strength"]
                - minimum
            )
            /
            (
                maximum
                - minimum
            )
        )

    return df


def main():

    print("=" * 60)
    print("GenePrior - GWAS Processing")
    print("=" * 60)

    # Step 1
    df = load_data()

    print("\nRaw data:")
    print(df.head())

    # Step 2
    df = clean_data(df)

    # Step 3
    significant = filter_significant(df)

    if significant.empty:

        print(
            "\nNo genome-wide significant "
            "associations found."
        )

        return

    # Step 4
    gene_summary = calculate_gene_summary(
        significant
    )

    # Step 5
    gene_summary = calculate_gwas_strength(
        gene_summary
    )

    # Step 6
    gene_summary = normalize_score(
        gene_summary
    )

    # Sort strongest genes first
    gene_summary = gene_summary.sort_values(
        by="gwas_score",
        ascending=False
    )

    # Step 7
    gene_summary.to_csv(
        PROCESSED_GWAS_FILE,
        index=False
    )

    print("\n" + "=" * 60)
    print("TOP GENES")
    print("=" * 60)

    print(
        gene_summary.head(20).to_string(
            index=False
        )
    )

    print("\nSaved:")
    print(PROCESSED_GWAS_FILE)


if __name__ == "__main__":
    main()