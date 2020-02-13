#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import pandas as pd

import config


#read input and save as df
df = pd.read_csv(config.inputfile, sep='\t')


#replace empty cells in columns impact and impact_so with Targeted Region 
df['impact'] = df.impact.fillna('Targeted_Region')
df['impact_so'] = df.impact_so.fillna('Targeted_Region')


#replace cells in column impact
df['impact'] = df.impact.replace({
    #VEP terms (uses SO by default)
	'splice_acceptor_variant' : 'Splice_Site',
	'splice_donor_variant' : 'Splice_Site',
	'stop_gained' : 'Nonsense_Mutation',
	'stop_lost' : 'Nonstop_Mutation',
	'frameshift_variant' : 'Frame_Shift_',
	'initiator_codon_variant' : 'Translation_Start_Site',
	'missense_variant' : 'Missense_Mutation',
	'inframe_insertion' : 'In_Frame_Ins',
	'inframe_deletion' : 'In_Frame_Del',
	'splice_region_variant' : 'Splice_Region',
	'mature_miRNA_variant' : 'RNA',
	'regulatory_region_variant' : 'IGR',
	'TF_binding_site_variant' : 'IGR',
	'regulatory_region_ablation' : '',
	'regulatory_region_amplification' : '',
	'TFBS_ablation' : '',
	'TFBS_amplification' : '',
	'stop_retained_variant' : 'Silent',
	'synonymous_variant' : 'Silent',
	'5_prime_UTR_variant' : "5'UTR",
	'3_prime_UTR_variant' : "3'UTR",
	'intron_variant' : 'Intron',
	'coding_sequence_variant' : 'Missense_Mutation',
	'upstream_gene_variant' : "5'Flank",
	'downstream_gene_variant' : "3'Flank",
	'intergenic_variant' : 'RNA',
	'nc_transcript_variant' : 'RNA',
	'NMD_transcript_variant' : 'Silent',
	'incomplete_terminal_codon_variant' : 'Silent',
	'non_coding_exon_variant' : 'RNA',
	'transcript_ablation' : 'Splice_Site',
	'transcript_amplification' : 'Intron',
	'feature_elongation' : '',
	'feature_truncation' : ''
})


#replace cells in column impact_so
df['impact_so'] = df.impact_so.replace({
    #snpEff terms
    'SPLICE_SITE_ACCEPTOR' : 'Splice_Site',
	'SPLICE_SITE_DONOR' : 'Splice_Site', 
	'STOP_GAINED' : 'Nonsense_Mutation',
	'STOP_LOST' : 'Nonstop_Mutation', 
	'FRAME_SHIFT' : 'Frame_Shift_', 
	'START_LOST' : '', 
	'EXON_DELETED' : '',
	'NON_SYNONYMOUS_START' : '',
	'CHROMOSOME_LARGE_DELETION' : '',
	'RARE_AMINO_ACID' : '',
	'NON_SYNONYMOUS_CODING' : 'Missense_Mutation',
	'CODON_INSERTION' : 'In_Frame_Ins',
	'CODON_DELETION' : 'In_Frame_Del',
	'CODON_CHANGE' : '',
	'CODON_CHANGE_PLUS_CODON_DELETION' : '',
	'CODON_CHANGE_PLUS_CODON_INSERTION' : '',
	'UTR_5_DELETED' : '',
	'UTR_3_DELETED' : '',
	'SPLICE_SITE_REGION' : 'Splice_Region',
	'SYNONYMOUS_STOP' : 'Silent',
	'SYNONYMOUS_CODING' : 'Silent',
	'UTR_5_PRIME' : "5'UTR",
	'UTR_3_PRIME' : "3'UTR",
	'INTRON' : 'Intron',
	'CDS' : 'Missense_Mutation',
	'UPSTREAM' : "5'Flank",
	'DOWNSTREAM' : "3'Flank",
	'INTERGENIC' : 'RNA',
	'INTERGENIC_CONSERVED' : '',
	'INTRAGENIC' : 'Intron',
	'GENE' : '',
	'TRANSCRIPT' : '',
	'EXON' : 'RNA',
	'START_GAINED' : '',
	'SYNONYMOUS_START' : '',
	'INTRON_CONSERVED' : ''
})


#add missing columns
if 'Mutation_Status' not in df:
    df['Mutation_Status'] = config.mutationStatus

if 'Tumor_Sample_Barcode' not in df:
    df['Tumor_Sample_Barcode'] = config.tumorSampleBarcode

if 'NCBI_Build' not in df:
    df['NCBI_Build'] = config.ncbiBuild

if 'Center' not in df:
    df['Center'] = config.center


#merge impact and impact_so to column impact (replace impact value with impact_so value in case that impact value is empty)
df.loc[df.impact == 'Targeted_Region', 'impact'] = df.impact_so


#add ins or del from sub_type to impact in case of Frame_Shift
df.loc[df.impact == 'Frame_Shift_', 'impact'] = df.impact.astype(str) + df.sub_type.astype(str).apply(lambda x: x.capitalize())


#select columns and order
df1 = df[['gene', 'entrez_id', 'Center', 'NCBI_Build', 'start', 'end', 'strand', 'impact', 'type', 'ref', 'alt', 'Tumor_Sample_Barcode', 'Mutation_Status', 'aa_change', 'chrom']]


#rename columns (gemini --> cBioPortal)
df2 = df1.rename(columns={
    'chrom' : 'Chromosome',
    'start' : 'Start_Position',
	'end' : 'End_Position',
	'ref' : 'Reference_Allele',
	'alt' : 'Tumor_Seq_Allele1',
	'type' : 'Variant_Type',
	'gene' : 'Hugo_Symbol',
	'aa_change' : 'HGVSp_Short',
	'impact' : 'Variant_Classification',
	'transcript' : 'Transcript_ID',
	'entrez_id' : 'Entrez_Gene_Id',
	'strand' : 'Strand'
})


#generate output
df2.to_csv(config.outputfile, sep='\t', encoding='utf-8', index=False)
