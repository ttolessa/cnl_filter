
if __name__ == '__main__':
    
    import shutil
    from pybedtools import BedTool
    import pandas as pd
    import os
    from Bio import SeqIO
    import argparse
    import subprocess
    import numpy as np
    
    description = """This program parses out the gene models on the identified NLR loci.\n
It filters them for NB-ARC domain containing proteins only.\n
Works with augustus.gff3 or braker.gtf.\n
It needs 'gene' in the feature column. It needs 'transcript_id' and 'gene_id' in the featur column."""
    
    
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('anno_fn', type=str, nargs='+',
    help='Please provide a gtf or gff3 annotation file.')
    
    parser.add_argument('pfam_fn', type=str, nargs='+',
    help='Please provide interpro tsv output file for the corresponding proteins.\nInterpro should be run with Pfam, Coils, and Gene3D.')
    
    parser.add_argument('fa_fn', type=str, nargs='+',
    help='Please provide NLR loci filename used for gene prediction.')
    
    parser.add_argument('pro_fn', type=str, nargs='+',
    help='Please provide protein filename of the annotated gene models.')
    
    parser.add_argument('genome_fn', type=str, nargs='+',
    help='Please provide protein filename of original genome used for NLR prediction.')


    augustus_gff_fn = os.path.abspath(parser.parse_args().anno_fn[0])
    pfam_fn = os.path.abspath(parser.parse_args().pfam_fn[0])
    NLR_loci_fa_fn = os.path.abspath(parser.parse_args().fa_fn[0])
    protein_fn = os.path.abspath(parser.parse_args().pro_fn[0])
    genome_fn = os.path.abspath(parser.parse_args().genome_fn[0])
    
    #print(augustus_gff_fn)
    #print(pfam_fn)
    
    prefix_list = os.path.basename(protein_fn).split('.')

    #if len(prefix_list) > 2:  
        #prefix = '.'.join(prefix_list[:-1])
    #else:
    prefix = os.path.basename(protein_fn).split('.')[0]
    suffix = os.path.basename(augustus_gff_fn).split('.')[-1]
    
    #input sorted move to setting up the tmp directory
    TMP_DIR = os.path.join(os.getcwd(), 'tmp')
    if not os.path.exists(TMP_DIR):
        os.makedirs(TMP_DIR)
    
    ##get initial chrom names from NLR file
    NLR_chrom_list = []
    with open(NLR_loci_fa_fn, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('>'):
                line = line.strip('\n')
                NLR_chrom = line[1:]
                NLR_chrom_list.append(NLR_chrom)
                
    ##split up braker.gtf or augustus.gff3 into tmp dir files one for each NLR_locus
    
    
    tmp_beds = []
    for NLR_chrom in NLR_chrom_list:
        tmp_fn = os.path.join(TMP_DIR, F"{NLR_chrom}.tmp.{suffix}")
        tmp_beds.append(tmp_fn)
        with open (tmp_fn, 'w') as tmp_out_fh:
            with open(augustus_gff_fn, 'r') as  augfh:
                for line in augfh.readlines():
                    line = line.strip('\n')
                    if line.startswith(NLR_chrom) and 'gene' == line.split('\t')[2]:
                        print(line, file = tmp_out_fh)
                        
    #check for bedfiles with overlap
    bed_with_overlap = []
    bedtool_list = []
    for tmp_gff in tmp_beds:
        try:
            bedtool_list.append(BedTool(tmp_gff).sort().merge(c=1, o='count').to_dataframe())
        except:
            print('This file tmp_gff file raises an error', tmp_gff)
           # bedtool_list.append(float('NaN'))
    
    #check if there are any overlapping gene models
    for tmp_df in bedtool_list:
        try:
            if sum(tmp_df.name > 1) > 0:
                bed_with_overlap.append(tmp_df.chrom.unique())
                print(F"Overlapping gene models in {tmp_df.chrom.unique()}")

        except:
            print(F'Having and issue with {tmp_df.chrom.unique()}')
            
    if len(bed_with_overlap):
        print('Script may want updating, as some gene models overlap. Might also stick with PFAM filtering for now.')
    
    
    #read in interpro results
    #this should be done with try and except really.

    interpro_df = pd.read_csv(pfam_fn, sep='\t', header = None, names=[x for x in range(0,13)])
    interpro_df = interpro_df.iloc[:,0:11].copy()
    pfam_df = interpro_df[interpro_df[3] ==  'Pfam']
    Gene3d_df = interpro_df[interpro_df[3] ==  'Gene3D']
    Coils_df = interpro_df[interpro_df[3] ==  'Coils']

    #domains to filter
    PFAM_LRR_domains = ['PF08263', 'PF12799', 'PF13306', 'PF13855', 'PF13516', 'PF00560', 'PF07725', 'PF07723', 'PF18837', 'PF01463', 'PF01462', 'PF18831', 'PF18805']
    Rx_N_domain = ['PF18052']
    TIR_domain = ['PF01582']
    NB_ARC_domain = ['PF00931' 'PF17862' 'PF00004']
    
    #individual gene containing domain lists
    NB_ARCs = np.concatenate((pfam_df[pfam_df[4].isin(NB_ARC_domain)][0].unique(),\
			       Gene3d_df[Gene3d_df[4]=='G3DSA:1.10.10.10'][0].unique(),\
			       Gene3d_df[Gene3d_df[4]=='G3DSA:3.40.50.300'][0].unique(),\
			       Gene3d_df[Gene3d_df[4]=='G3DSA:1.10.8.430'][0].unique()), axis=None)
    LRRs = np.concatenate((pfam_df[pfam_df[4].isin(PFAM_LRR_domains)][0].unique(),\
			       Gene3d_df[Gene3d_df[4]=='G3DSA:3.80.10.10'][0].unique(),\
			       Gene3d_df[Gene3d_df[4]=='G3DSA:1.25.10.10'][0].unique()), axis=None)
    TIRs = np.concatenate((pfam_df[pfam_df[4].isin(TIR_domain)][0].unique(),\
			       Gene3d_df[Gene3d_df[4]=='G3DSA:3.40.50.10140'][0].unique()), axis=None)
    RPWs = pfam_df[pfam_df[4] == 'PF05659'][0].unique()
    Coils = np.concatenate((pfam_df[pfam_df[4].isin(Rx_N_domain)][0].unique(),\
			       Coils_df[0].unique()), axis=None)
    JACs = pfam_df[pfam_df[4] == 'PF01419'][0].unique()
     
    #combinations
    NBLRRs = np.intersect1d(NB_ARCs, LRRs)	# Centrally conserved NBS with C-terminal lucin-rich repeat(LRR) domains
    NBS = list((((set(NB_ARCs) - set(LRRs)) - set(TIRs)) - set(Coils)) - set(RPWs))	# Only the centrally conserved NBS without N and C treminal domains
    NL = list((((set(NB_ARCs) & set(LRRs)) - set(TIRs)) - set(Coils)) - set(RPWs))  # Conserved NBS with lucin-rich repeat(LRR) but no N-terminal TIR, Coil and RPW8 domains
    NB_ARC = list(set(NB_ARCs))		# All the NB-ARC doamins
    TNLs = np.intersect1d(NBLRRs, TIRs) # Centrally conserved NBS with lucin-rich repeat(LRR) and N-terminal TIR only
    TNJs = list((set(NB_ARCs) & set(JACs)) & set(TIRs)) # Centrally conserved NBS with C-terminal Jacalin and N-terminal TIR only
    TN = list(((set(TIRs) - set(LRRs)) & set(NB_ARCs)) - set(Coils)) # Centrally conserved NBS with N-terminal TIR domain only
    NJ = list((set(NB_ARCs) & set(JACs)) - set(TIRs)) # Centrally conserved NBS with C-terminal Jacalin domain only
    TL = list((((set(TIRs) & set(LRRs)) - set(NB_ARCs)) -set(Coils)) -set(RPWs)) # TIR domain with lucin-rich repeat but without centrally conserved NBS domain
    TIR = list((((set(TIRs) - set(LRRs)) - set(NB_ARCs)) -set(Coils)) -set(RPWs)) # TIRonly domain 
    CNLs = np.intersect1d(NBLRRs, Coils) # Centrally conserved NBS with lucin-rich repeat(LRR) and N-terminal Coiled-coil domain
    CN = list((set(Coils) - set(LRRs)) & set(NB_ARCs)) # Centrally conserved NBS with N-terminal Coiled-coil domain 
    CL = list(((set(Coils) & set(LRRs)) - set(NB_ARCs)) -set(TIRs)) # Coiled-coil domain with lucin-rich repeat but without centrally conserved NBS domain
    COIL = list((((set(Coils) - set(LRRs)) - set(NB_ARCs)) -set(TIRs)) -set(RPWs)) # Coiled-coil only domain
    CTNL = list((set(NBLRRs) & set(TIRs)) & set(Coils)) # Centrally conserved NBS domain with C-terminal LRR and N-terminal TIR and Coiled-coil
    CT = list((set(TIRs) & set(Coils)) - set(NBLRRs)) # Domains containing only TIR and Coiled-Coil 
    RNLs = np.intersect1d(NBLRRs, RPWs) # Centrally conserved NBS with lucin-rich repeat(LRR) and N-terminal RPW8 domain
    RN = list((set(RPWs) - set(LRRs)) & set(NB_ARCs)) # Centrally conserved NBS with N-terminal RPW8 but no C-terminal LRR domain
    
#generate the rename dicts and put them in a dict for easy usage downstream
    NBLRR_rename = dict(zip(NBLRRs, [F'{prefix}_{x}' for x in NBLRRs]))
    NBS_rename = dict(zip(NBS, [F'{prefix}_{x}' for x in NBS]))
    NL_rename = dict(zip(NL, [F'{prefix}_{x}' for x in NL]))
    NB_ARC_rename = dict(zip(NB_ARC, [F'{prefix}_{x}' for x in NB_ARC]))
    TNL_rename = dict(zip(TNLs, [F'{prefix}_{x}' for x in TNLs]))
    TNJ_rename = dict(zip(TNJs, [F'{prefix}_{x}' for x in TNJs]))
    TN_rename = dict(zip(TN, [F'{prefix}_{x}' for x in TN]))
    NJ_rename = dict(zip(NJ, [F'{prefix}_{x}' for x in NJ]))
    TL_rename = dict(zip(TL, [F'{prefix}_{x}' for x in TL]))
    TIR_rename = dict(zip(TIR, [F'{prefix}_{x}' for x in TIR]))
    CNL_rename = dict(zip(CNLs, [F'{prefix}_{x}' for x in CNLs]))
    CN_rename = dict(zip(CN, [F'{prefix}_{x}' for x in CN]))
    CL_rename = dict(zip(CL, [F'{prefix}_{x}' for x in CL]))
    COIL_rename = dict(zip(COIL, [F'{prefix}_{x}' for x in COIL]))
    CTNL_rename = dict(zip(CTNL, [F'{prefix}_{x}' for x in CTNL]))
    CT_rename = dict(zip(CT, [F'{prefix}_{x}' for x in CT]))
    RNL_rename = dict(zip(RNLs, [F'{prefix}_{x}' for x in RNLs]))
    RN_rename = dict(zip(RN, [F'{prefix}_{x}' for x in RN]))    
    rename_dicts = {}
    rename_dicts['NBLRR'] = NBLRR_rename
    rename_dicts['NBSonly'] = NBS_rename
    rename_dicts['NLonly'] = NL_rename
    rename_dicts['all_NB_ARC'] = NB_ARC_rename
    rename_dicts['TNL'] = TNL_rename
    rename_dicts['TNJ'] = TNJ_rename
    rename_dicts['TNonly'] = TN_rename
    rename_dicts['NJonly'] = NJ_rename
    rename_dicts['TLonly'] = TN_rename
    rename_dicts['TIRonly'] = TIR_rename
    rename_dicts['CNL'] = CNL_rename
    rename_dicts['CNonly'] = CN_rename
    rename_dicts['CLonly'] = CN_rename
    rename_dicts['COILonly'] = COIL_rename
    rename_dicts['CTNL'] = CTNL_rename
    rename_dicts['CTonly'] = CT_rename
    rename_dicts['RNL'] = RNL_rename
    rename_dicts['RNonly'] = RN_rename
    
#report back
    for key, values in rename_dicts.items():
        print(F"This is the number of gene models with {key} domains containing genes: {len(values)}")


#pull out gene models and save them
    for key, rename_dict in rename_dicts.items():
        protein_filter_fn = protein_fn.replace('.aa', F'.{key}_filtered.aa')
        NLR_protein_seqs = []
        for seq in SeqIO.parse(protein_fn, 'fasta'):
            if seq.id in rename_dict.keys():
                seq.id = rename_dict[seq.id]
                NLR_protein_seqs.append(seq)
        SeqIO.write(NLR_protein_seqs, protein_filter_fn, 'fasta')
        #report back
        print(F"Pulled out filtered gene models for {key} and saved corresponding proteins \
    as {protein_filter_fn}\n")
    
    #now filter the augusts gff file and rename all the gene and transcript ids with file the prefix
    #also adjust the coordinates to fit the original genome again. This was a bit tricky based on the 1 off errors.
    for key, rename_dict in rename_dicts.items():
        augustus_gff_filter_fn = os.path.join(os.getcwd(), F'{prefix}.{key}.{suffix}')
        with open(augustus_gff_filter_fn, 'w') as out_fh:
            with open(augustus_gff_fn) as aug_fh:
                for line in aug_fh.readlines():
                    if line.startswith('#'):
                        continue
                    else:
                        try:
                            line = line.strip('\n')
                            line_list = line.split('\t')

                        except:
                            print(F"Funny line {line}")
                    #print(len(line_list), line)


                    if any([ x for x in rename_dict.keys() if F"transcript_id \"file_1_file_1_{x}\"" in  line]): #catches old transcript id lines
                        try: 
                            for old, new in rename_dict.items():
                                if line_list[8] == F"transcript_id \"file_1_file_1_{old}\"; gene_id \"file_1_file_1_{old.split('.')[0]}\";":
                                    line_list[8] = F"transcript_id \"{new}\"; gene_id \"{new.split('.')[0]}\";"
                                    #line = F'{line_list[0]}\t{line_list[1]}\t{line_list[2]}\t{line_list[3]}\t{line_list[4]}\t{line_list[5]}\t{line_list[6]}\t{line_list[7]}\t{line_list[8]}'
                                    line = F'{line_list[0].split(":")[0]}\t{line_list[1]}\t{line_list[2]}\t{int(line_list[0].split(":")[1].split("-")[0])+int(line_list[3]) }\t{int(line_list[0].split(":")[1].split("-")[0])+int(line_list[4])}\t{line_list[5]}\t{line_list[6]}\t{line_list[7]}\t{line_list[8]}'
                                    print(line, file = out_fh)

                        except:
                            print('Issue_1')
                        #print(line, file = out_fh)
                    elif line_list[2] == 'gene': #catches old gene id lines
                        #print('gene')
                        try: 
                            for old, new in rename_dict.items():
                                if line_list[8] ==  old.split('.')[0]:
                                    line_list[8] = new.split('.')[0]
                                    #line = F'{line_list[0]}\t{line_list[1]}\t{line_list[2]}\t{line_list[3]}\t{line_list[4]}\t{line_list[5]}\t{line_list[6]}\t{line_list[7]}\t{line_list[8]}'
                                    line = F'{line_list[0].split(":")[0]}\t{line_list[1]}\t{line_list[2]}\t{int(line_list[0].split(":")[1].split("-")[0])+int(line_list[3]) }\t{int(line_list[0].split(":")[1].split("-")[0])+int(line_list[4])}\t{line_list[5]}\t{line_list[6]}\t{line_list[7]}\t{line_list[8]}'
                                    print(line, file = out_fh)
                        except:
                            print('Issue')

                    elif line_list[2] == 'transcript':
                        for old, new in rename_dict.items():
                            if line_list[8] ==  old:
                                line_list[8] = new
                                #line = F'{line_list[0]}\t{line_list[1]}\t{line_list[2]}\t{line_list[3]}\t{line_list[4]}\t{line_list[5]}\t{line_list[6]}\t{line_list[7]}\t{line_list[8]}'
                                line = F'{line_list[0].split(":")[0]}\t{line_list[1]}\t{line_list[2]}\t{int(line_list[0].split(":")[1].split("-")[0])+int(line_list[3]) }\t{int(line_list[0].split(":")[1].split("-")[0])+int(line_list[4])}\t{line_list[5]}\t{line_list[6]}\t{line_list[7]}\t{line_list[8]}'
                                print(line, file = out_fh)

        command = F"/home/tamene/anaconda3/bin/getAnnoFasta.pl {augustus_gff_filter_fn} --seqfile={genome_fn}"
        try:
            print(command)
            output = subprocess.check_output(
                command,
                stderr=subprocess.STDOUT,
                shell=True)
        except:
            print(F'{command}.\n This did not work please check.')


        #translate coding sequence
        NLR_protein_fn = augustus_gff_filter_fn.replace(suffix, 'aa')
        NLR_proteins = []
        for seq in SeqIO.parse(augustus_gff_filter_fn.replace(suffix, 'codingseq'), 'fasta'):
            seq.seq = seq.seq.translate()
            NLR_proteins.append(seq)
        SeqIO.write(NLR_proteins, NLR_protein_fn, 'fasta')

        print(F"Generated new translation of proteins in {NLR_protein_fn}")

    shutil.rmtree(TMP_DIR)
    
    print('All done selecting your NLR gene models.')
