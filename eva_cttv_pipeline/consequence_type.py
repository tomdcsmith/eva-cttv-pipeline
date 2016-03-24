import xlrd

__author__ = 'Javier Lopez: javild@gmail.com'


def _process_consequence_type_file_xls(snp_2_gene_file):

    consequence_type_dict = {}
    one_rs_multiple_genes = set()

    ct_mapping_read_book = xlrd.open_workbook(snp_2_gene_file, formatting_info=True)
    ct_mapping_read_sheet = ct_mapping_read_book.sheet_by_index(0)
    for i in range(1, ct_mapping_read_sheet.nrows):
        if ct_mapping_read_sheet.cell_value(rowx=i, colx=2) != 'Not found':

            rs_id = ct_mapping_read_sheet.cell_value(rowx=i, colx=0)
            ensembl_gene_id = ct_mapping_read_sheet.cell_value(rowx=i, colx=2)
            so_term = ct_mapping_read_sheet.cell_value(rowx=i, colx=1)

            if rs_id in consequence_type_dict:
                if ensembl_gene_id != consequence_type_dict[rs_id].getEnsemblGeneId():
                    print('WARNING (clinvar_record.py): different genes and annotations found for a given gene.')
                    print(' Variant id: ' + rs_id + ', ENSG: ' + so_term + ', ENSG: ' + consequence_type_dict[rs_id].getEnsemblGeneId())
                    print('Skipping')
                    one_rs_multiple_genes.add(rs_id)
                else:
                    consequence_type_dict[rs_id].add_so_term(so_term)
            else:
                consequence_type_dict[rs_id] = ConsequenceType(ensembl_gene_id, [so_term])

    return consequence_type_dict, one_rs_multiple_genes


def _process_gene(consequence_type_dict, rs_id, ensembl_gene_id, so_term):
    if rs_id in consequence_type_dict:
        consequence_type_dict[rs_id].add_ensembl_gene_id(ensembl_gene_id)
        consequence_type_dict[rs_id].add_so_term(so_term)
    else:
        consequence_type_dict[rs_id] = ConsequenceType([ensembl_gene_id], [so_term])


def _process_consequence_type_file_tsv(snp_2_gene_file):

    consequence_type_dict = {}
    one_rs_multiple_genes = set()

    with open(snp_2_gene_file, "rt") as f:
        for line in f:
            line = line.rstrip()
            line_list = line.split("\t")

            rs_id = line_list[0]
            ensembl_gene_id = line_list[2]
            if not ensembl_gene_id or rs_id == "rs":
                continue
            so_term = line_list[4]

            ensembl_gene_ids = ensembl_gene_id.split(",")
            for ensembl_gene_id in ensembl_gene_ids:
                _process_gene(consequence_type_dict, rs_id, ensembl_gene_id, so_term)

    return consequence_type_dict, one_rs_multiple_genes


def process_consequence_type_file(snp_2_gene_file):

    print('Loading mapping rs->ENSG/SOterms')

    if snp_2_gene_file.endswith(".xls"):
        consequence_type_dict, one_rs_multiple_genes = _process_consequence_type_file_xls(snp_2_gene_file)
    else:
        consequence_type_dict, one_rs_multiple_genes = _process_consequence_type_file_tsv(snp_2_gene_file)

    print(str(len(consequence_type_dict)) + ' rs->ENSG/SOterms mappings loaded')
    print(str(len(one_rs_multiple_genes)) + ' rsIds with multiple gene associations')
    print('Done.')

    return consequence_type_dict


class SoTerm(object):

    so_accession_name_dict = {'transcript_ablation': 1893,
                              'splice_donor_variant': 1575,
                              'splice_acceptor_variant': 1574,
                              'stop_gained': 1587,
                              'frameshift_variant': 1589,
                              'stop_lost': 1578,
                              'initiator_codon_variant': 1582,
                              'inframe_insertion': 1821,
                              'inframe_deletion': 1822,
                              'missense_variant': 1583,
                              'transcript_amplification': 1889,
                              'splice_region_variant': 1630,
                              'incomplete_terminal_codon_variant': 1626,
                              'synonymous_variant': 1819,
                              'stop_retained_variant': 1567,
                              'coding_sequence_variant': 1580,
                              'miRNA': 276, 'miRNA_target_site': 934,
                              'mature_miRNA_variant': 1620,
                              '5_prime_UTR_variant': 1623,
                              '3_prime_UTR_variant': 1624,
                              'exon_variant': 1791,
                              'non_coding_transcript_exon_variant': 1792,
                              'non_coding_transcript_variant': 1619,
                              'intron_variant': 1627,
                              'NMD_transcript_variant': 1621,
                              'TFBS_ablation': 1895,
                              'TFBS_amplification': 1892,
                              'TF_binding_site_variant': 1782,
                              'regulatory_region_variant': 1566,
                              'regulatory_region_ablation': 1894,
                              'regulatory_region_amplification': 1891,
                              'feature_elongation': 1907,
                              'feature_truncation': 1906,
                              'intergenic_variant': 1628,
                              'lincRNA': 1463,
                              'downstream_gene_variant': 1632,
                              '2KB_downstream_gene_variant': 1632,
                              'upstream_gene_variant': 1631,
                              '2KB_upstream_gene_variant': 1631,
                              'SNV': 1483,
                              'SNP': 694,
                              'RNA_polymerase_promoter': 1203,
                              'CpG_island': 307,
                              'DNAseI_hypersensitive_site': 685,
                              'polypeptide_variation_site': 336,
                              'start_lost': 2012,
                              'protein_altering_variant': 1818}

    ranked_so_names_list = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained',
                         'frameshift_variant', 'stop_lost', 'initiator_codon_variant', 'transcript_amplification',
                         'inframe_insertion', 'inframe_deletion', 'missense_variant', 'splice_region_variant',
                         'incomplete_terminal_codon_variant', 'stop_retained_variant', 'synonymous_variant',
                         'coding_sequence_variant', 'mature_miRNA_variant', '5_prime_UTR_variant',
                            '3_prime_UTR_variant', 'non_coding_transcript_exon_variant', 'intron_variant',
                            'NMD_transcript_variant', 'non_coding_transcript_variant', 'upstream_gene_variant',
                            'downstream_gene_variant', 'TFBS_ablation', 'TFBS_amplification',
                            'TF_binding_site_variant', 'regulatory_region_ablation',
                            'regulatory_region_amplification', 'regulatory_region_variant', 'feature_elongation',
                            'feature_truncation', 'intergenic_variant']

    def __init__(self, so_name):
        self._so_name = so_name
        if so_name in SoTerm.so_accession_name_dict:
            self._so_accession = SoTerm.so_accession_name_dict[so_name]
        else:
            self._so_accession = None

    def get_name(self):
        return self._so_name

    def get_accession(self):
        if self._so_accession is not None:
            accession_number_str = str(self._so_accession)
            return 'SO:' + accession_number_str.rjust(7, '0')
        else:
            return None

    def get_rank(self):
        # If So name not in Ensembl's ranked list, return the least severe rank
        if self.get_name() not in SoTerm.ranked_so_names_list:
            return len(SoTerm.ranked_so_names_list)
        else:
            return SoTerm.ranked_so_names_list.index(self.get_name())

    @staticmethod
    def get_ranked_so_names():
        return SoTerm.ranked_so_names_list

    def __eq__(self, other):
        return self.get_accession() == other.get_accession()

    def __hash__(self):
        return hash(self._so_accession)


class ConsequenceType(object):

    def __init__(self, ensembl_gene_ids=None, so_names=None):
        if ensembl_gene_ids:
            self._ensembl_gene_ids = set(ensembl_gene_ids)
        else:
            self._ensembl_gene_ids = set()
        self._ensemblTranscriptId = None

        if so_names is not None:
            self._so_terms = set([SoTerm(so_name) for so_name in so_names])
        else:
            self._so_terms = None

    # def getEnsemblGeneId(self):
    #     return self._ensemblGeneId

    def get_ensembl_gene_ids(self):
        return self._ensembl_gene_ids

    def add_so_term(self, soName):
        if self._so_terms is None:
            self._so_terms = {SoTerm(soName)}
        else:
            self._so_terms.add(SoTerm(soName))

    def get_so_accessions(self):
        return [soTerm.get_accession() for soTerm in self._so_terms]

    def getMostSevereSo(self):
        return min(list(self._so_terms), key=lambda x: x.get_rank())

    def get_so_terms(self):
        return self._so_terms

    def add_ensembl_gene_id(self, new_ensembl_gene_id):
        self._ensembl_gene_ids.add(new_ensembl_gene_id)
