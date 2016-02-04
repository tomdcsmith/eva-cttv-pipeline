__author__ = 'Javier Lopez: javild@gmail.com'


class ConsequenceType:

    class SoTerm:

        soAccessionToNameDict = {}
        soAccessionToNameDict['transcript_ablation'] = 1893
        soAccessionToNameDict['splice_donor_variant'] = 1575
        soAccessionToNameDict['splice_acceptor_variant'] = 1574
        soAccessionToNameDict['stop_gained'] = 1587
        soAccessionToNameDict['frameshift_variant'] = 1589
        soAccessionToNameDict['stop_lost'] = 1578
        soAccessionToNameDict['initiator_codon_variant'] = 1582
        soAccessionToNameDict['inframe_insertion'] = 1821
        soAccessionToNameDict['inframe_deletion'] = 1822
        soAccessionToNameDict['missense_variant'] = 1583
        soAccessionToNameDict['transcript_amplification'] = 1889
        soAccessionToNameDict['splice_region_variant'] = 1630
        soAccessionToNameDict['incomplete_terminal_codon_variant'] = 1626
        soAccessionToNameDict['synonymous_variant'] = 1819
        soAccessionToNameDict['stop_retained_variant'] = 1567
        soAccessionToNameDict['coding_sequence_variant'] = 1580
        soAccessionToNameDict['miRNA'] = 276
        soAccessionToNameDict['miRNA_target_site'] = 934
        soAccessionToNameDict['mature_miRNA_variant'] = 1620
        soAccessionToNameDict['5_prime_UTR_variant'] = 1623
        soAccessionToNameDict['3_prime_UTR_variant'] = 1624
        soAccessionToNameDict['exon_variant'] = 1791
        soAccessionToNameDict['non_coding_transcript_exon_variant'] = 1792
        soAccessionToNameDict['non_coding_transcript_variant'] = 1619
        soAccessionToNameDict['intron_variant'] = 1627
        soAccessionToNameDict['NMD_transcript_variant'] = 1621
        soAccessionToNameDict['TFBS_ablation'] = 1895
        soAccessionToNameDict['TFBS_amplification'] = 1892
        soAccessionToNameDict['TF_binding_site_variant'] = 1782
        soAccessionToNameDict['regulatory_region_variant'] = 1566
        soAccessionToNameDict['regulatory_region_ablation'] = 1894
        soAccessionToNameDict['regulatory_region_amplification'] = 1891
        soAccessionToNameDict['feature_elongation'] = 1907
        soAccessionToNameDict['feature_truncation'] = 1906
        soAccessionToNameDict['intergenic_variant'] = 1628
        soAccessionToNameDict['lincRNA'] = 1463
        soAccessionToNameDict['downstream_gene_variant'] = 1632
        soAccessionToNameDict['2KB_downstream_gene_variant'] = 1632
        soAccessionToNameDict['upstream_gene_variant'] = 1631
        soAccessionToNameDict['2KB_upstream_gene_variant'] = 1631
        soAccessionToNameDict['SNV'] = 1483
        soAccessionToNameDict['SNP'] = 694
        soAccessionToNameDict['RNA_polymerase_promoter'] = 1203
        soAccessionToNameDict['CpG_island'] = 307
        soAccessionToNameDict['DNAseI_hypersensitive_site'] = 685
        soAccessionToNameDict['polypeptide_variation_site'] = 336

        rankedSoNamesList = ['transcript_ablation','splice_acceptor_variant','splice_donor_variant','stop_gained',
                                  'frameshift_variant','stop_lost','initiator_codon_variant','transcript_amplification',
                                  'inframe_insertion','inframe_deletion','missense_variant','splice_region_variant',
                                  'incomplete_terminal_codon_variant','stop_retained_variant','synonymous_variant',
                                  'coding_sequence_variant','mature_miRNA_variant','5_prime_UTR_variant',
                                  '3_prime_UTR_variant','non_coding_transcript_exon_variant','intron_variant',
                                  'NMD_transcript_variant','non_coding_transcript_variant','upstream_gene_variant',
                                  'downstream_gene_variant','TFBS_ablation','TFBS_amplification',
                                  'TF_binding_site_variant','regulatory_region_ablation',
                                  'regulatory_region_amplification','regulatory_region_variant','feature_elongation',
                                  'feature_truncation','intergenic_variant']

        def __init__(self, soName):
            self._soName = soName
            if(soName in ConsequenceType.SoTerm.soAccessionToNameDict):
                self._soAccession = ConsequenceType.SoTerm.soAccessionToNameDict[soName]
            else:
                self._soAccession = None

        def getName(self):
            return self._soName

        def getAccession(self):
            if(self._soAccession != None):
                accessionNumberStr = str(self._soAccession)
                return 'SO:'+accessionNumberStr.rjust(7,'0')
            else:
                return None

        def getRank(self):
            # If So name not in Ensembl's ranked list, return the least severe rank
            if(self.getName() not in ConsequenceType.SoTerm.rankedSoNamesList):
                return len(ConsequenceType.SoTerm.rankedSoNamesList)
            else:
                return ConsequenceType.SoTerm.rankedSoNamesList.index(self.getName())

        def getRankedSoNames(self):
            return ConsequenceType.SoTerm.rankedSoNamesList

        def __eq__(self, other):
            return self.getAccession()==other.getAccession()

        def __hash__(self):
            return hash((self._soAccession))

    def __init__(self, ensemblGeneId=None, soNames=None):
        self._ensemblGeneId = ensemblGeneId
        self._ensemblTranscriptId = None

        if(soNames != None):
            self._soTerms= set([ConsequenceType.SoTerm(soName) for soName in soNames])
        else:
            self._soTerms = None

    def getEnsemblGeneId(self):
        return self._ensemblGeneId

    def addSoTerm(self, soName):
        if(self._soTerms == None):
            self._soTerms = set([ConsequenceType.SoTerm(soName)])
        else:
            self._soTerms.add(ConsequenceType.SoTerm(soName))

    def getSoAccessions(self):
        return [soTerm.getAccession() for soTerm in self._soTerms]

    def getMostSevereSo(self):
        return min(list(self._soTerms), key=lambda x:x.getRank())
