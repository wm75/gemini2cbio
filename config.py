#!/usr/bin/env python
# -*- coding: utf-8 -*- 


inputfile = 'beispiel2.tabular'
outputfile = 'output_panda.maf'
metadatafile = 'meta.txt'


mutationStatus = 'somatisch'
tumorSampleBarcode = 'abx'
transcriptId = 0
ncbiBuild = 'GRCh37'
center = 'Sequenz Zentrum'
strand = '+'

#Meta data
studyName = 'study_1'
geneticAlterationType = 'COPY_NUMBER_ALTERATION'
datatype = 'SEG'
referenceGenomeId = 'hg19'
description = 'Somatic CNA data (copy number ratio from tumor samples minus ratio from matched normals) from TCGA.'
dataFile = 'brca_tcga_data_cna_hg19.seg'
