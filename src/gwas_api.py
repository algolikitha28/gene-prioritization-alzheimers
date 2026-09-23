# src/gwas_api.py

import requests


def get_json(url, params=None):
    """
    Send a GET request to the GWAS Catalog API
    and return the JSON response.
    """

    response = requests.get(
        url,
        params=params,
        timeout=60
    )

    response.raise_for_status()

    return response.json()


def get_studies(disease_name, page_size=100):
    """
    Retrieve GWAS studies related to a disease.
    """

    url = "https://www.ebi.ac.uk/gwas/rest/api/v2/studies"

    params = {
        "disease_trait": disease_name,
        "size": page_size
    }

    return get_json(url, params)


def get_associations(
    disease_name=None,
    mapped_gene=None,
    page_size=100
):
    """
    Retrieve GWAS associations.

    We can filter by disease and/or gene.
    """

    url = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"

    params = {
        "size": page_size
    }

    if disease_name:
        params["disease_trait"] = disease_name

    if mapped_gene:
        params["mapped_gene"] = mapped_gene

    return get_json(url, params)