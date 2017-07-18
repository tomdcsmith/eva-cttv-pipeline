import json
import requests
import urllib




##
## Zooma functions
##

def build_zooma_query(trait_name, filters, zooma_host):
    url = "{}/spot/zooma/v2/api/services/annotate?propertyValue={}".format(zooma_host, trait_name)
    url_filters = [
                    "required:[{}]".format(filters["required"]),
                    "ontologies:[{}]".format(filters["ontologies"]),
                    "preferred:[{}]".format(filters["preferred"])
                  ]
    url += "&filter={}".format(",".join(url_filters))
    return url






##
## OLS functions
##

def double_encode_uri(uri):
    return urllib.parse.quote(urllib.parse.quote(uri, safe=""), safe="")


def ols_efo_query(uri):
    double_encoded_uri = double_encode_uri(uri)
    return requests.get("http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/{}".format(double_encoded_uri))


def is_current_and_in_efo(uri):
    response = ols_efo_query(uri)
    if response.status_code == 400:
        return False
    response_json = response.json()
    return not response_json["is_obsolete"]


def is_in_efo(uri):
    response = ols_efo_query(uri)
    return response.status_code == 200
