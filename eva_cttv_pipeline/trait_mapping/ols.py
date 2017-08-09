from functools import lru_cache
import requests
import urllib

from eva_cttv_pipeline.trait_mapping.utils import request_retry_helper


def ols_query_helper(url: str) -> str:
    """
    Given a url for OLS, make a get request and return the label for the term, from the response
    from OLS.

    :param url: OLS url to which to make a get request to query for a term.
    :return: The ontology label of the term specified in the url.
    """
    try:
        json_response = requests.get(url).json()
        for term in json_response["_embedded"]["terms"]:
            if term["is_defining_ontology"]:
                return term["label"]
    except:
        return None


@lru_cache(maxsize=16384)
def get_ontology_label_from_ols(ontology_uri: str) -> str:
    """
    Using provided ontology uri, build an OLS url with which to make a request for the uri to find
    the term label for this uri.

    :param ontology_uri: A uri for a term in an ontology.
    :return: Term label for the ontology uri provided in the parameters.
    """
    url = build_ols_query(ontology_uri)
    label = request_retry_helper(ols_query_helper, 4, url)
    return label


def build_ols_query(ontology_uri: str) -> str:
    """Build a url to query OLS for a given ontology uri."""
    url = "http://www.ebi.ac.uk/ols/api/terms?iri={}".format(ontology_uri)
    return url


def double_encode_uri(uri: str) -> str:
    """Double encode a given uri."""
    return urllib.parse.quote(urllib.parse.quote(uri, safe=""), safe="")


def ols_efo_query(uri: str) -> requests.Response:
    """
    Query EFO using OLS for a given ontology uri, returning the response from the request.

    :param uri: Ontology uri to use in querying EFO using OLS
    :return: Response from OLS
    """
    double_encoded_uri = double_encode_uri(uri)
    return requests.get(
        "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/{}".format(double_encoded_uri))


@lru_cache(maxsize=16384)
def is_current_and_in_efo(uri: str) -> bool:
    """
    Checks whether given ontology uri is a valid and non-obsolete term in EFO.

    :param uri: Ontology uri to use in querying EFO using OLS
    :return: Boolean value, true if ontology uri is valid and non-obsolete term in EFO
    """
    response = ols_efo_query(uri)
    if response.status_code != 200:
        return False
    response_json = response.json()
    return not response_json["is_obsolete"]


@lru_cache(maxsize=16384)
def is_in_efo(uri: str) -> bool:
    """
    Checks whether given ontology uri is a valid term in EFO.

    :param uri: Ontology uri to use in querying EFO using OLS
    :return: Boolean value, true if ontology uri is valid and non-obsolete term in EFO
    """
    response = ols_efo_query(uri)
    return response.status_code == 200
