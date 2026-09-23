import time
import json
import requests
import pandas as pd


VEP_URL = "https://rest.ensembl.org/vep/human/id"

INPUT_FILE = "data/processed/unique_variants.csv"
OUTPUT_FILE = "data/processed/variant_annotations.csv"
RAW_OUTPUT_FILE = "data/processed/vep_raw.json"


def annotate_variants(rsids, max_retries=5):

    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }

    payload = {
        "ids": rsids
    }

    for attempt in range(1, max_retries + 1):

        try:
            response = requests.post(
                VEP_URL,
                headers=headers,
                json=payload,
                timeout=120
            )

            print(f"Attempt {attempt}: Status {response.status_code}")

            if response.status_code == 200:
                return response.json()

            elif response.status_code in [429, 500, 502, 503, 504]:

                wait_time = attempt * 10

                print(
                    f"Ensembl temporarily unavailable. "
                    f"Waiting {wait_time} seconds..."
                )

                time.sleep(wait_time)

            else:

                print(response.text)
                response.raise_for_status()

        except requests.RequestException as e:

            print(f"Request error: {e}")

            if attempt < max_retries:

                wait_time = attempt * 10

                print(
                    f"Retrying in {wait_time} seconds..."
                )

                time.sleep(wait_time)

            else:
                raise

    raise RuntimeError(
        "VEP annotation failed after multiple attempts."
    )


def parse_vEP_results(data):

    annotations = []

    for record in data:

        variant = record.get("id")

        most_severe = record.get(
            "most_severe_consequence"
        )

        # -----------------------------
        # Transcript-level consequences
        # -----------------------------

        transcript_consequences = record.get(
            "transcript_consequences",
            []
        )

        for consequence in transcript_consequences:

            annotations.append({
                "variant": variant,
                "vep_gene": consequence.get("gene_symbol"),
                "gene_id": consequence.get("gene_id"),
                "transcript_id": consequence.get("transcript_id"),
                "consequence": ",".join(
                    consequence.get(
                        "consequence_terms",
                        []
                    )
                ),
                "impact": consequence.get("impact"),
                "distance": consequence.get("distance"),
                "most_severe_consequence": most_severe
            })

        # -----------------------------
        # Intergenic consequences
        # -----------------------------

        intergenic_consequences = record.get(
            "intergenic_consequences",
            []
        )

        for consequence in intergenic_consequences:

            annotations.append({
                "variant": variant,
                "vep_gene": None,
                "gene_id": None,
                "transcript_id": None,
                "consequence": ",".join(
                    consequence.get(
                        "consequence_terms",
                        []
                    )
                ),
                "impact": consequence.get("impact"),
                "distance": None,
                "most_severe_consequence": most_severe
            })

    return annotations


if __name__ == "__main__":

    print("Reading unique variants...")

    df = pd.read_csv(INPUT_FILE)

    rsids = (
        df["variant"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    print(f"Unique variants found: {len(rsids)}")

    print()
    print("Sending variants to Ensembl VEP...")

    data = annotate_variants(rsids)

    print()
    print(f"VEP returned {len(data)} variant records.")

    # Save raw VEP response
    with open(
        RAW_OUTPUT_FILE,
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            data,
            f,
            indent=2
        )

    print("Raw VEP response saved to:")
    print(RAW_OUTPUT_FILE)

    # Parse VEP results
    annotations = parse_vEP_results(data)

    annotation_df = pd.DataFrame(annotations)

    # Remove exact duplicate annotation rows
    annotation_df = annotation_df.drop_duplicates()

    annotation_df.to_csv(
        OUTPUT_FILE,
        index=False
    )

    print()
    print("Duplicate rows removed.")
    print(f"Clean annotation rows: {len(annotation_df)}")

    print()
    print("Annotation table created.")
    print(f"Rows: {len(annotation_df)}")
    print(
        f"Columns: {annotation_df.columns.tolist()}"
    )

    print()
    print("Variants successfully represented in annotation table:")

    print(
        annotation_df["variant"]
        .nunique()
    )

    print()
    print(
        annotation_df.head(10)
        .to_string(index=False)
    )