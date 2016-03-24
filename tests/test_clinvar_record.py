import json
import unittest

from eva_cttv_pipeline import clinvar_record


class ClinvarRecord(unittest.TestCase):
    def setUp(self):
        self.clinvar_record = get_test_record()




def get_test_record():
    test_record = """{
  "chromosome": "16",
  "start": 2138765,
  "end": 2138765,
  "reference": "G",
  "alternate": "GG",
  "clinvarSet": {
    "recordStatus": "current",
    "title": "NM_000548.3(TSC2):c.*154dup AND Tuberous sclerosis syndrome",
    "referenceClinVarAssertion": {
      "clinVarAccession": {
        "acc": "RCV000055062",
        "version": 1,
        "type": "RCV",
        "dateUpdated": 1412982000000
      },
      "recordStatus": "current",
      "clinicalSignificance": {
        "reviewStatus": "NOT_CLASSIFIED_BY_SUBMITTER",
        "description": "not provided"
      },
      "assertion": {
        "type": "VARIATION_TO_DISEASE"
      },
      "observedIn": [
        {
          "sample": {
            "origin": "germline",
            "species": {
              "value": "human",
              "taxonomyId": 9606
            },
            "affectedStatus": "yes"
          },
          "method": [
            {
              "methodType": "LITERATURE_ONLY"
            }
          ],
          "observedData": [
            {
              "attribute": {
                "integerValue": 1,
                "type": "VariantAlleles"
              },
              "id": 3619513
            }
          ]
        }
      ],
      "measureSet": {
        "measure": [
          {
            "name": [
              {
                "elementValue": {
                  "value": "NM_000548.3(TSC2):c.*154dup",
                  "type": "Preferred"
                }
              }
            ],
            "attributeSet": [
              {
                "attribute": {
                  "value": "NM_000548.3:c.*154dup",
                  "type": "HGVS, coding, RefSeq",
                  "change": "c.*154dup"
                }
              },
              {
                "attribute": {
                  "value": "NM_001009944.2:c.*963dupC",
                  "type": "HGVS, coding, RefSeq",
                  "change": "c.*963dupC"
                }
              },
              {
                "attribute": {
                  "value": "NG_005895.1:g.44459dupG",
                  "type": "HGVS, genomic, RefSeqGene",
                  "change": "g.44459dupG"
                }
              },
              {
                "attribute": {
                  "value": "NC_000016.10:g.2088764dupG",
                  "integerValue": 38,
                  "type": "HGVS, genomic, top level",
                  "change": "g.2088764dupG"
                }
              },
              {
                "attribute": {
                  "value": "NC_000016.9:g.2138765dupG",
                  "integerValue": 37,
                  "type": "HGVS, genomic, top level, previous",
                  "change": "g.2138765dupG"
                }
              },
              {
                "attribute": {
                  "value": "p.(=)",
                  "type": "HGVS, non-validated"
                }
              },
              {
                "attribute": {
                  "value": "Exon 41",
                  "type": "Location"
                }
              },
              {
                "attribute": {
                  "value": "3 prime UTR variant",
                  "type": "MolecularConsequence"
                },
                "xref": [
                  {
                    "db": "Sequence Ontology",
                    "id": "SO:0001624",
                    "status": "CURRENT"
                  },
                  {
                    "db": "RefSeq",
                    "id": "NM_001009944.2:c.*962_*963insC",
                    "status": "CURRENT"
                  }
                ]
              },
              {
                "attribute": {
                  "value": "500B downstream variant",
                  "type": "MolecularConsequence"
                },
                "xref": [
                  {
                    "db": "Sequence Ontology",
                    "id": "SO:0001634",
                    "status": "CURRENT"
                  },
                  {
                    "db": "RefSeq",
                    "id": "NM_000548.3:c.*154_*155insG",
                    "status": "CURRENT"
                  }
                ]
              }
            ],
            "cytogeneticLocation": [
              "16p13.3"
            ],
            "sequenceLocation": [
              {
                "assembly": "GRCh37",
                "chr": "16",
                "accession": "NC_000016.9",
                "start": 2138765,
                "stop": 2138765,
                "variantLength": 2,
                "referenceAllele": "G",
                "alternateAllele": "GG"
              },
              {
                "assembly": "GRCh38",
                "chr": "16",
                "accession": "NC_000016.10",
                "start": 2088764,
                "stop": 2088764,
                "variantLength": 2,
                "referenceAllele": "G",
                "alternateAllele": "GG"
              }
            ],
            "measureRelationship": [
              {
                "name": [
                  {
                    "elementValue": {
                      "value": "polycystic kidney disease 1 (autosomal dominant)",
                      "type": "Preferred"
                    }
                  }
                ],
                "symbol": [
                  {
                    "elementValue": {
                      "value": "PKD1",
                      "type": "Preferred"
                    }
                  }
                ],
                "attributeSet": [
                  {
                    "attribute": {
                      "value": "Sufficient evidence for dosage pathogenicity",
                      "dateValue": 1329868800000,
                      "type": "Haploinsufficiency"
                    },
                    "citation": [
                      {
                        "url": "http://www.ncbi.nlm.nih.gov/projects/dbvar/ISCA/isca_gene.cgi?sym=PKD1"
                      }
                    ]
                  },
                  {
                    "attribute": {
                      "value": "No evidence available",
                      "dateValue": 1329868800000,
                      "type": "Triplosensitivity"
                    },
                    "citation": [
                      {
                        "url": "http://www.ncbi.nlm.nih.gov/projects/dbvar/ISCA/isca_gene.cgi?sym=PKD1"
                      }
                    ]
                  }
                ],
                "sequenceLocation": [
                  {
                    "assembly": "GRCh37",
                    "chr": "16",
                    "accession": "NC_000016.9",
                    "start": 2138710,
                    "stop": 2185898,
                    "strand": "-"
                  },
                  {
                    "assembly": "GRCh38",
                    "chr": "16",
                    "accession": "NC_000016.10",
                    "start": 2088707,
                    "stop": 2135897,
                    "strand": "-"
                  }
                ],
                "type": "variant in gene",
                "xref": [
                  {
                    "db": "Gene",
                    "id": "5310",
                    "status": "CURRENT"
                  },
                  {
                    "db": "OMIM",
                    "id": "601313",
                    "type": "MIM",
                    "status": "CURRENT"
                  }
                ]
              },
              {
                "name": [
                  {
                    "elementValue": {
                      "value": "tuberous sclerosis 2",
                      "type": "Preferred"
                    }
                  }
                ],
                "symbol": [
                  {
                    "elementValue": {
                      "value": "TSC2",
                      "type": "Preferred"
                    }
                  }
                ],
                "attributeSet": [
                  {
                    "attribute": {
                      "value": "Sufficient evidence for dosage pathogenicity",
                      "dateValue": 1329868800000,
                      "type": "Haploinsufficiency"
                    },
                    "citation": [
                      {
                        "url": "http://www.ncbi.nlm.nih.gov/projects/dbvar/ISCA/isca_gene.cgi?sym=TSC2"
                      }
                    ]
                  },
                  {
                    "attribute": {
                      "value": "No evidence available",
                      "dateValue": 1329868800000,
                      "type": "Triplosensitivity"
                    },
                    "citation": [
                      {
                        "url": "http://www.ncbi.nlm.nih.gov/projects/dbvar/ISCA/isca_gene.cgi?sym=TSC2"
                      }
                    ]
                  }
                ],
                "sequenceLocation": [
                  {
                    "assembly": "GRCh37",
                    "chr": "16",
                    "accession": "NC_000016.9",
                    "start": 2097989,
                    "stop": 2138712,
                    "strand": "+"
                  },
                  {
                    "assembly": "GRCh38",
                    "chr": "16",
                    "accession": "NC_000016.10",
                    "start": 2047801,
                    "stop": 2088711,
                    "strand": "+"
                  }
                ],
                "comment": [
                  {
                    "value": "This gene is cited in the ACMG recommendations of 2013 (PubMed 23788249) for reporting incidental findings in exons.",
                    "dataSource": "NCBI curation"
                  }
                ],
                "type": "variant in gene",
                "xref": [
                  {
                    "db": "Gene",
                    "id": "7249",
                    "status": "CURRENT"
                  },
                  {
                    "db": "OMIM",
                    "id": "191092",
                    "type": "MIM",
                    "status": "CURRENT"
                  }
                ]
              }
            ],
            "type": "Duplication",
            "id": 75791,
            "xref": [
              {
                "db": "Tuberous sclerosis database (TSC2)",
                "id": "TSC2_02318",
                "status": "CURRENT"
              },
              {
                "db": "dbSNP",
                "id": "397514891",
                "type": "rs",
                "status": "CURRENT"
              }
            ]
          }
        ],
        "name": [
          {
            "elementValue": {
              "value": "NM_000548.3(TSC2):c.*154dup",
              "type": "preferred name"
            }
          }
        ],
        "type": "Variant",
        "id": 64862
      },
      "traitSet": {
        "trait": [
          {
            "name": [
              {
                "elementValue": {
                  "value": "Tuberous sclerosis syndrome",
                  "type": "Preferred"
                },
                "xref": [
                  {
                    "db": "SNOMED CT",
                    "id": "7199000",
                    "status": "CURRENT"
                  }
                ]
              }
            ],
            "symbol": [
              {
                "elementValue": {
                  "value": "TSC",
                  "type": "Preferred"
                },
                "xref": [
                  {
                    "db": "OMIM",
                    "id": "191100",
                    "type": "MIM",
                    "status": "CURRENT"
                  }
                ]
              },
              {
                "elementValue": {
                  "value": "TS",
                  "type": "Alternate"
                },
                "xref": [
                  {
                    "db": "OMIM",
                    "id": "191100",
                    "type": "MIM",
                    "status": "CURRENT"
                  }
                ]
              }
            ],
            "attributeSet": [
              {
                "attribute": {
                  "value": "Tuberous sclerosis complex (TSC) involves abnormalities of the skin (hypomelanotic macules, facial angiofibromas, shagreen patches, fibrous facial plaques, ungual fibromas); brain (cortical tubers, subependymal nodules [SENs] and subependymal giant cell astrocytomas [SEGAs], seizures, intellectual disability/developmental delay); kidney (angiomyolipomas, cysts, renal cell carcinomas); heart (rhabdomyomas, arrhythmias); and lungs (lymphangioleiomyomatosis [LAM]). CNS tumors are the leading cause of morbidity and mortality; renal disease is the second leading cause of early death.",
                  "type": "public definition"
                },
                "xref": [
                  {
                    "db": "GeneReviews",
                    "id": "NBK1220",
                    "status": "CURRENT"
                  }
                ]
              },
              {
                "attribute": {
                  "value": "Neoplasm",
                  "type": "keyword"
                }
              },
              {
                "attribute": {
                  "value": "Hereditary cancer syndrome",
                  "type": "keyword"
                }
              }
            ],
            "citation": [
              {
                "url": "https://www.orpha.net/data/patho/Pro/en/Emergency_TuberousSclerosis.pdf",
                "citationText": "Orphanet, Tuberous sclerosis, 2007",
                "type": "practice guideline",
                "abbrev": "Orphanet, 2007"
              },
              {
                "id": {
                  "value": "20301399",
                  "source": "PubMed"
                },
                "type": "review",
                "abbrev": "GeneReviews"
              }
            ],
            "type": "Disease",
            "id": 15993,
            "xref": [
              {
                "db": "MedGen",
                "id": "C0041341",
                "status": "CURRENT"
              }
            ]
          }
        ],
        "type": "Disease",
        "id": 8139
      },
      "dateCreated": 1379286000000,
      "dateLastUpdated": 1412982000000,
      "id": 144533
    },
    "clinVarAssertion": [
      {
        "clinVarSubmissionID": {
          "submitter": "Tuberous sclerosis database (TSC2)",
          "title": "NM_000548.3:c.*154dup AND TSC",
          "localKey": "TSC2_02318_TSC",
          "submitterDate": 1376002800000
        },
        "clinVarAccession": {
          "acc": "SCV000083280",
          "version": 1,
          "type": "SCV",
          "orgID": 500074,
          "dateUpdated": 1403478000000
        },
        "recordStatus": "current",
        "clinicalSignificance": {
          "reviewStatus": "CLASSIFIED_BY_SINGLE_SUBMITTER",
          "description": [
            "not provided"
          ]
        },
        "assertion": {
          "type": "variation to disease"
        },
        "externalID": {
          "db": "Tuberous sclerosis database (TSC2)",
          "id": "TSC2_02318",
          "status": "CURRENT"
        },
        "observedIn": [
          {
            "sample": {
              "origin": "germline",
              "species": {
                "value": "human"
              },
              "affectedStatus": "yes"
            },
            "method": [
              {
                "methodType": "LITERATURE_ONLY"
              }
            ],
            "observedData": [
              {
                "attribute": {
                  "value": "1",
                  "type": "VariantAlleles"
                }
              }
            ]
          }
        ],
        "measureSet": {
          "measure": [
            {
              "attributeSet": [
                {
                  "attribute": {
                    "value": "NM_000548.3:c.*154dup",
                    "type": "HGVS"
                  }
                },
                {
                  "attribute": {
                    "value": "p.(=)",
                    "type": "HGVS"
                  }
                },
                {
                  "attribute": {
                    "value": "Exon 41",
                    "type": "Location"
                  }
                }
              ],
              "measureRelationship": [
                {
                  "symbol": [
                    {
                      "elementValue": {
                        "value": "TSC2",
                        "type": "Preferred"
                      }
                    }
                  ],
                  "type": "variant in gene"
                }
              ],
              "type": "Variation"
            }
          ],
          "type": "Variant"
        },
        "traitSet": {
          "trait": [
            {
              "name": [
                {
                  "elementValue": {
                    "value": "TSC",
                    "type": "Preferred"
                  }
                }
              ],
              "type": "Disease"
            }
          ],
          "type": "Disease"
        },
        "id": 143786
      }
    ],
    "id": 3756609
  },
  "annot": {
    "chromosome": "16",
    "start": 2138765,
    "end": 2138765,
    "referenceAllele": "G",
    "alternativeAllele": "GG",
    "consequenceTypes": [
      {
        "geneName": "TSC2",
        "ensemblGeneId": "ENSG00000103197",
        "ensemblTranscriptId": "ENST00000439673",
        "strand": "+",
        "biotype": "protein_coding",
        "aaChange": "-",
        "codon": "-",
        "soTerms": [
          {
            "soName": "downstream_gene_variant",
            "soAccession": "SO:0001632"
          }
        ]
      },
      {
        "geneName": "MIR1225",
        "ensemblGeneId": "ENSG00000221656",
        "ensemblTranscriptId": "ENST00000408729",
        "strand": "-",
        "biotype": "miRNA",
        "aaChange": "-",
        "codon": "-",
        "soTerms": [
          {
            "soName": "downstream_gene_variant",
            "soAccession": "SO:0001632"
          }
        ]
      },
      {
        "geneName": "PKD1",
        "ensemblGeneId": "ENSG00000008710",
        "ensemblTranscriptId": "ENST00000262304",
        "strand": "-",
        "biotype": "protein_coding",
        "cDnaPosition": 14084,
        "aaChange": "-",
        "codon": "-",
        "soTerms": [
          {
            "soName": "3_prime_UTR_variant",
            "soAccession": "SO:0001624"
          },
          {
            "soName": "feature_elongation",
            "soAccession": "SO:0001907"
          }
        ]
      },
      {
        "geneName": "TSC2",
        "ensemblGeneId": "ENSG00000103197",
        "ensemblTranscriptId": "ENST00000219476",
        "strand": "+",
        "biotype": "protein_coding",
        "aaChange": "-",
        "codon": "-",
        "soTerms": [
          {
            "soName": "downstream_gene_variant",
            "soAccession": "SO:0001632"
          }
        ]
      },
      {
        "geneName": "TSC2",
        "ensemblGeneId": "ENSG00000103197",
        "ensemblTranscriptId": "ENST00000353929",
        "strand": "+",
        "biotype": "protein_coding",
        "aaChange": "-",
        "codon": "-",
        "soTerms": [
          {
            "soName": "downstream_gene_variant",
            "soAccession": "SO:0001632"
          }
        ]
      },
      {
        "geneName": "TSC2",
        "ensemblGeneId": "ENSG00000103197",
        "ensemblTranscriptId": "ENST00000350773",
        "strand": "+",
        "biotype": "protein_coding",
        "aaChange": "-",
        "codon": "-",
        "soTerms": [
          {
            "soName": "downstream_gene_variant",
            "soAccession": "SO:0001632"
          }
        ]
      },
      {
        "geneName": "TSC2",
        "ensemblGeneId": "ENSG00000103197",
        "ensemblTranscriptId": "ENST00000568454",
        "strand": "+",
        "biotype": "protein_coding",
        "aaChange": "-",
        "codon": "-",
        "soTerms": [
          {
            "soName": "downstream_gene_variant",
            "soAccession": "SO:0001632"
          }
        ]
      },
      {
        "geneName": "TSC2",
        "ensemblGeneId": "ENSG00000103197",
        "ensemblTranscriptId": "ENST00000382538",
        "strand": "+",
        "biotype": "protein_coding",
        "aaChange": "-",
        "codon": "-",
        "soTerms": [
          {
            "soName": "downstream_gene_variant",
            "soAccession": "SO:0001632"
          }
        ]
      },
      {
        "geneName": "TSC2",
        "ensemblGeneId": "ENSG00000103197",
        "ensemblTranscriptId": "ENST00000401874",
        "strand": "+",
        "biotype": "protein_coding",
        "aaChange": "-",
        "codon": "-",
        "soTerms": [
          {
            "soName": "downstream_gene_variant",
            "soAccession": "SO:0001632"
          }
        ]
      },
      {
        "geneName": "PKD1",
        "ensemblGeneId": "ENSG00000008710",
        "ensemblTranscriptId": "ENST00000423118",
        "strand": "-",
        "biotype": "protein_coding",
        "cDnaPosition": 14081,
        "aaChange": "-",
        "codon": "-",
        "soTerms": [
          {
            "soName": "3_prime_UTR_variant",
            "soAccession": "SO:0001624"
          },
          {
            "soName": "feature_elongation",
            "soAccession": "SO:0001907"
          }
        ]
      },
      {
        "geneName": "RP11-304L19.1",
        "ensemblGeneId": "ENSG00000259933",
        "ensemblTranscriptId": "ENST00000570072",
        "strand": "+",
        "biotype": "sense_overlapping",
        "aaChange": "-",
        "codon": "-",
        "soTerms": [
          {
            "soName": "upstream_gene_variant",
            "soAccession": "SO:0001631"
          }
        ]
      },
      {
        "soTerms": [
          {
            "soName": "regulatory_region_variant",
            "soAccession": "SO:0001566"
          }
        ]
      }
    ]
  }
}"""
    record_string = json.load(test_record)
    return record_string