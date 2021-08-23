from Bio import SeqIO
import os
import sys
import pandas as pd


def find_mut_from_aligned_fasta(aligned_fasta):
    # get 1st record as target
    record_iterator = list(SeqIO.parse(aligned_fasta, 'fasta'))
    for record in record_iterator:
        if record.name == 'IGHG1':
            ref = record.seq
        elif record.name in ['sp|P01857|IGHG1_HUMAN', 'sp|P01859|IGHG2_HUMAN', 'sp|P01860|IGHG3_HUMAN', 'sp|P01859|IGHG2_HUMAN', 'sp|P01861|IGHG4_HUMAN', 'sp|P01877|IGHA2_HUMAN', 'sp|P01876|IGHA1_HUMAN', 'sp|P01871|IGHM_HUMAN', 'sp|P01880|IGHD_HUMAN', 'sp|P01854|IGHE_HUMAN']:
            type_query = record.seq
        else:
            query = record.seq
    ignored = [' ', '-']
    num = 0  # res numbering starts from 1, and is based on the ref not query
    mutations = []
    length = min(len(ref), len(type_query), len(query))
    for i in range(length):  # iterate on query
        if ref[i] not in ignored:  # not gap in query
            num += 1
            if query[i] in ignored and type_query[i] not in ignored:
                mutations.append(type_query[i] + str(num+117)+'del')
            elif query[i] not in ignored and type_query[i] not in ignored and query[i] != type_query[i]:  # not match
                mutation = type_query[i] + str(num+117) + query[i]
                # print(mutation)
                mutations.append(mutation)
            elif query[i] not in ignored and type_query[i] in ignored:
                mutations.append(query[i] + str(num+117)+'insert')
    return mutations


def get_Fc(query_file, db):  # return Fc domain
    inputfa = query_file
    result = 'temp.out'
    #cmd_1=""" makeblastdb -in IGHG1.fst   -dbtype prot -parse_seqids -out IGHG1"""
    cmd_2 = """blastp -query %s -db %s -evalue 1000 -num_threads 30  -max_target_seqs 30 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen  slen" -out %s """ % (
        inputfa, db, result)
    cmd_3 = """ sed -i "1i\#Query\tSbjct\tIdent\tAlignment_length\tMsmatch\tGapopen\tQstart\tQend\tSstart\tSend\tEvalue\tBitscore\tQlen\tSlen" %s """ % result
    # os.system(cmd_1)
    os.system(cmd_2)
    os.system(cmd_3)
    df = pd.read_csv(result, sep='\t', engine='python')
    qstart = int(df.iloc[0, 6])
    qend = int(df.iloc[0, 7])

    query = open(inputfa, 'r')
    subject = open('%s.fst' % db, 'r')
    fout = open('align.fst', 'w')
    fout.write(""">IGHG1
ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSS
GLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGG
PSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYN
STYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDE
LTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRW
QQGNVFSCSVMHEALHNHYTQKSLSLSPGK\n""")
    for record in SeqIO.parse(subject, 'fasta'):
        fout.write('>%s\n%s\n' % (record.description, record.seq))

    for record in SeqIO.parse(query, 'fasta'):
        seq = record.seq
        align_seq = seq[qstart-1:qend]
        fout.write('>%s\n%s\n' % (record.description, align_seq))
    query.close()
    subject.close()
    fout.close()
    return qstart, qend, str(align_seq)


# get CH1 Hinge CH2 CH3
def Fc_numbering(query_file, db='./db/IGHG1'):
    qstart, qend, Fc_seq = get_Fc(query_file, db='./db/IGHG1')
    cmd_4 = """muscle -out Align.fst -quiet -in align.fst"""
    os.system(cmd_4)
    aligned_fasta = open('Align.fst', 'r')
    record_iterator = list(SeqIO.parse(aligned_fasta, 'fasta'))
    ignored = [' ', '-']
    num = 0
    for record in record_iterator:
        if record.name == 'IGHG1':
            ref = record.seq
        elif record.name in ['sp|P01857|IGHG1_HUMAN', 'sp|P01859|IGHG2_HUMAN', 'sp|P01860|IGHG3_HUMAN', 'sp|P01859|IGHG2_HUMAN', 'sp|P01861|IGHG4_HUMAN', 'sp|P01877|IGHA2_HUMAN', 'sp|P01876|IGHA1_HUMAN', 'sp|P01871|IGHM_HUMAN', 'sp|P01880|IGHD_HUMAN', 'sp|P01854|IGHE_HUMAN']:
            type_query = record.seq
        else:
            query = record.seq
    length = min(len(ref), len(type_query), len(query))
    for i in range(length):
        if ref[i] not in ignored:
            if num == 97:
                end_CH1 = i
            elif num == 98:
                start_Hinge = i
            elif num == 112:
                end_Hinge = i
            elif num == 113:
                start_CH2 = i
            elif num == 222:
                end_CH2 = i
            elif num == 223:
                start_CH3 = i
            elif num == 329:
                end_CH3 = i
            num += 1
    ch1 = str(query[0:end_CH1+1]).replace('-', '')
    hinge = str(query[start_Hinge:end_Hinge+1]).replace('-', '')
    ch2 = str(query[start_CH2:end_CH2+1]).replace('-', '')
    ch3 = str(query[start_CH3:end_CH3+1]).replace('-', '')
    fc_list = [ch1, hinge, ch2, ch3]
    return fc_list


ha = {'P01857': 'IGHG1',
      'P01859': 'IGHG2',
      'P01860': 'IGHG3',
      'P01861': 'IGHG4',
      'P01877': 'IGHA2',
      'P01876': 'IGHA1',
      'P01871': 'IGHM',
      'P01880': 'IGHD',
      'P01854': 'IGHE'}

inputfa = sys.argv[1]
mutation_set = []
for record in SeqIO.parse(inputfa, 'fasta'):
    temfout = open('temp.fst', 'w')
    db = './db/IG'
    result = 'temp_1.out'
    temfout.write('>%s\n%s\n' % (record.description, record.seq))
    temfout.close()
    cmd_2 = """blastp -query temp.fst -db %s -evalue 1000 -num_threads 30  -max_target_seqs 30 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen  slen" -out %s """ % (
        db, result)
    cmd_3 = """ sed -i "1i\#Query\tSbjct\tIdent\tAlignment_length\tMsmatch\tGapopen\tQstart\tQend\tSstart\tSend\tEvalue\tBitscore\tQlen\tSlen" %s """ % result
    os.system(cmd_2)
    os.system(cmd_3)
    df = pd.read_csv(result, sep='\t', engine='python')
    name = df.iloc[0, 0]
    db = 'db/%s' % ha[df.iloc[0, 1]]
    get_Fc('temp.fst', db=db)
    # get mutation
    cmd_4 = """muscle -out Align.fst -quiet -in align.fst"""
    os.system(cmd_4)
    #from bpkit.mutate import (convert_mutations_to_df, find_mut_from_aligned_fasta, group_mut_df,write_mut_df_to_txt)
    mutations = find_mut_from_aligned_fasta('./Align.fst')
    mutation_set.append([name, df.iloc[0, 2], df.iloc[0, 1], mutations])

df1 = pd.DataFrame(mutation_set)
df1.columns = ['Name', 'Identity', 'uniprot_id', 'mutations']
for key in ha:
    df1.loc[df1[df1['uniprot_id'] == key].index, 'Subtype'] = ha[key]
df1.to_excel('Heavy_chain_mutation.xlsx', index=False)
os.system('rm ./temp_1.out')
os.system('rm ./temp.out')
os.system('rm ./align.fst')
os.system('rm ./Align.fst')
os.system('rm ./temp.fst')
