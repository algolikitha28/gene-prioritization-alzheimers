import requests
import pandas as pd

from config import (
    RAW_GWAS_FILE,
    PAGE_SIZE
)


BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2"

# Validated Alzheimer's disease GWAS
# Bellenguez et al. 2022
ALZHEIMER_STUDY = "GCST90027158"


def get_json(url, params=None):

    response = requests.get(
        url,
        params=params,
        timeout=60
    )

    print(f"\nGET {response.url}")
    print("Status:", response.status_code)

    response.raise_for_status()

    return response.json()


def get_study():

    print("\n" + "=" * 60)
    print("STEP 1 - CHECK ALZHEIMER GWAS STUDY")
    print("=" * 60)

    url = f"{BASE_URL}/studies/{ALZHEIMER_STUDY}"

    data = get_json(url)

    print("\nStudy information:")

    print(
        "Accession:",
        data.get("accession_id")
    )

    print(
        "First author:",
        data.get("first_author")
    )

    print(
        "PubMed ID:",
        data.get("pubmed_id")
    )

    print(
        "Reported trait:",
        data.get("reported_trait")
    )

    return data


def get_associations():

    print("\n" + "=" * 60)
    print("STEP 2 - DOWNLOAD ALZHEIMER ASSOCIATIONS")
    print("=" * 60)

    url = f"{BASE_URL}/associations"

    params = {
        "accession_id": ALZHEIMER_STUDY,
        "size": PAGE_SIZE,
        "page": 0
    }

    data = get_json(
        url,
        params
    )

    embedded = data.get(
        "_embedded",
        {}
    )

    associations = embedded.get(
        "associations",
        []
    )

    print(
        "\nAssociations received:",
        len(associations)
    )

    return associations


def extract_records(associations):

    records = []

    for association in associations:

        p_value = association.get(
            "p_value"
        )

        mapped_genes = association.get(
            "mapped_genes",
            []
        )

        snp_alleles = association.get(
            "snp_allele",
            []
        )

        # -------------------------
        # Extract variant
        # -------------------------

        variant_id = None

        if (
            isinstance(snp_alleles, list)
            and snp_alleles
        ):

            first_snp = snp_alleles[0]

            if isinstance(
                first_snp,
                dict
            ):

                variant_id = first_snp.get(
                    "rs_id"
                )

        # -------------------------
        # Extract genes
        # -------------------------

        if not isinstance(
            mapped_genes,
            list
        ):

            mapped_genes = []

        # -------------------------
        # Extract trait
        # -------------------------

        reported_trait = association.get(
            "reported_trait",
            []
        )

        if isinstance(
            reported_trait,
            list
        ):

            trait = (
                reported_trait[0]
                if reported_trait
                else ""
            )

        else:

            trait = str(
                reported_trait
            )

        # -------------------------
        # Create records
        # -------------------------

        for gene in mapped_genes:

            if not gene:
                continue

            records.append(
                {
                    "gene": str(
                        gene
                    ).strip().upper(),

                    "variant": variant_id,

                    "p_value": p_value,

                    "trait": trait,

                    "study": association.get(
                        "accession_id"
                    ),

                    "pubmed_id": association.get(
                        "pubmed_id"
                    ),

                    "author": association.get(
                        "first_author"
                    )
                }
            )

    return records


def main():

    print("=" * 60)
    print("GenePrior - GWAS Data Download")
    print("=" * 60)

    print(
        f"Alzheimer GWAS: {ALZHEIMER_STUDY}"
    )

    # ============================================
    # STEP 1
    # Verify study
    # ============================================

    study = get_study()

    # ============================================
    # STEP 2
    # Get associations
    # ============================================

    associations = get_associations()

    if not associations:

        print(
            "\nNo associations returned."
        )

        return

    # ============================================
    # STEP 3
    # Extract genes
    # ============================================

    records = extract_records(
        associations
    )

    print(
        "\nGene-association records:",
        len(records)
    )

    if not records:

        print(
            "\nNo mapped genes found."
        )

        print(
            "\nFirst association:"
        )

        print(
            associations[0]
        )

        return

    # ============================================
    # STEP 4
    # DataFrame
    # ============================================

    df = pd.DataFrame(
        records
    )

    print("\nColumns:")

    print(
        df.columns.tolist()
    )

    print("\nFirst 10 records:")

    print(
        df.head(10).to_string(
            index=False
        )
    )

    # ============================================
    # STEP 5
    # Validation
    # ============================================

    print("\n" + "=" * 60)
    print("VALIDATION")
    print("=" * 60)

    print(
        "\nStudy IDs:"
    )

    print(
        df["study"].value_counts()
    )

    print(
        "\nTraits:"
    )

    print(
        df["trait"].value_counts().head(20)
    )

    # ============================================
    # STEP 6
    # Save
    # ============================================

    df.to_csv(
        RAW_GWAS_FILE,
        index=False
    )

    print(
        f"\nSaved to: {RAW_GWAS_FILE}"
    )

    print("\n" + "=" * 60)
    print("GWAS DOWNLOAD COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()