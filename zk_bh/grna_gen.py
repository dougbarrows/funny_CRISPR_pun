#! /usr/bin/env python3


import argparse, re

parser = argparse.ArgumentParser(description='Imports fasta, extracts possible gRNAs')
parser.add_argument('-f','--fasta',required=True,help='Fasta file to extract gRNAs from')
parser.add_argument('-c','--cas_list',required=True,help='Tab delimitted list of Cas protein parameters')
args = parser.parse_args()


def cas2dict(cas_prot_tsv):

    with open(cas_prot_tsv, 'r') as cas_file:
        cas_dict = {}
        for line in cas_file:
            if line[0] != '#':
                if line[-1] == '\n':
                    line = line.rstrip()
                cas_param = line.split('\t')
                cas_prot = cas_param[0]
                pam = cas_param[1]
                pam_orient = cas_param[2]
                grna_len = cas_param[3]
                cas_dict[cas_prot] = [pam, pam_orient, grna_len]

    return cas_dict


def fasta2dict(fasta_file):
    fasta_dict = {}

    with open(fasta_file, 'r') as fasta:
        str_fasta = ''
        for line in fasta:
            # rstrips if it's not a sequence line and isn't a line with only \n
            if line[0] != '>' and line[0] != '\n':
                line = line.rstrip()
            # concatenates a new line character if it is the first sequence line
            if line[0] == '>' and str_fasta != '':
                str_fasta = str_fasta + '\n'
            # concatenates the prep string with the prepped line from the fasta file
            str_fasta = str_fasta + line

        # extracts 0) seq ID 1) description (or the \n if there is non) 2) sequence
        extracted = re.findall(r'(>\S+)([^\n]\s\w+.+\n|\n)(\S+[^>])', str_fasta)

        # adds a new dictionary for each gene
        for index in range(len(extracted)):
            gene = extracted[index][0]
            gene = gene[1:]

            # remove description if it is not actually present
            descrip = extracted[index][1]
            descrip = descrip[1:]
            seq = extracted[index][2]
            if seq[-1] is '\n' or '\r':
                seq = seq.rstrip()

            # prepares dictionaries for each gene with description, seq, rvcmpl_seq, and codons
            fasta_dict[gene] = {}
            if descrip == '\n':
                fasta_dict[gene]['description'] = None
            else:
                fasta_dict[gene]['description'] = descrip
            fasta_dict[gene]['sequence'] = seq.upper()

    return fasta_dict


def grna_finder(fasta_dict,cas_dict):

    # dict for translating ambiguous nucleotides
    ambiguous_nt = {
            'N': 'ATCGYRBDHV',
            'R': 'AG',
            'Y': 'CT',
            'B': 'CTGY',
            'D': 'ATGR',
            'H': 'ATCY',
            'V': 'ACGR'
            }

    # dict for translating sequences to their reverse complement, including ambiguous nucleotides
    reverse_nt = {
            'N': 'ATCGYRBDHV',
            'R': 'TC',
            'Y': 'GA',
            'B': 'GACVR',
            'D': 'TACY',
            'H': 'GATR',
            'V': 'CGTY',
            'G': 'C',
            'A': 'T',
            'C': 'G',
            'T': 'A'
            }


    grna_dict_list = []

    for cas in cas_dict:

        pam = cas_dict[cas][0]
        pam_orientation = cas_dict[cas][1]
        grna_length = int(cas_dict[cas][2])

        # populating PAM strings for concatenation
        pam_re = ''
        rev_pam_re = ''

        # for each character in the inputted Pam
        for char in pam:
            # if it is a standard nucleotide, add it to 'pam_re'
            if char in ['A','T','C','G']:
                pam_re += char
                rev_pam_re = reverse_nt[char] + rev_pam_re
            else:
                pam_re += '[' + ambiguous_nt[char] + ']'
                rev_pam_re = '[' + reverse_nt[char] + ']' + rev_pam_re

        grna_counter = {}

        for seq_name in fasta_dict:
            sequence = fasta_dict[seq_name]['sequence']
            if pam_orientation == "3'":
                grna_matches = re.finditer(r'(?=(\w{' + str(grna_length) + '}' + pam_re + '))', sequence)
                rev_grna_matches = re.finditer(r'(?=('+ rev_pam_re +'\w{' + str(grna_length) + '}))', sequence)
                match_info = [(match.start(1), match.end(1), match.group(1)[:grna_length], match.group(1)[grna_length:]) for match in grna_matches]
                rev_match_info = [(rev_match.start(1), rev_match.end(1), rev_match.group(1)[len(pam):], rev_match.group(1)[:len(pam)]) for rev_match in rev_grna_matches] 
            elif pam_orientation == "5'":
                grna_matches = re.finditer(r'(?=(' + pam_re + '\w{' + str(grna_length) + '}))', sequence)
                rev_grna_matches = re.finditer(r'(?=(\w{' + str(grna_length) + '}'+rev_pam_re+'))', sequence)
                match_info = [(match.start(1), match.end(1), match.group(1)[len(pam):], match.group(1)[:len(pam)]) for match in grna_matches]
                rev_match_info = [(rev_match.start(1), rev_match.end(1), rev_match.group(1)[:grna_length], rev_match.group(1)[grna_length:]) for rev_match in rev_grna_matches]
            index = -1
            for grna in match_info:
                index += 1
                if grna_counter.get(grna[2]):
                    grna_counter[grna[2]].append([index,'+'])
                else:
                    grna_counter[grna[2]] = [[index,'+']]
                grna_dict_list.append({})
                grna_dict = grna_dict_list[-1]
                grna_dict['chrom'] = seq_name
                grna_dict['name'] = cas
                grna_dict['chrom_start'] = grna[0]
                grna_dict['chrom_stop'] = grna[1]
                grna_dict['grna_seq'] = grna[2]
                grna_dict['strand'] = '+'
                grna_dict['pam'] = grna[3]

            for grna in rev_match_info:
                index += 1
                if grna_counter.get(grna[2]):
                    grna_counter[grna[2]].append([index,'-'])
                else:
                    grna_counter[grna[2]] = [[index,'-']]
                grna_dict_list.append({})
                grna_dict = grna_dict_list[-1]
                grna_dict['chrom'] = seq_name
                grna_dict['name'] = cas
                grna_dict['chrom_start'] = grna[0]
                grna_dict['chrom_stop'] = grna[1]
                grna_dict['grna_seq'] = grna[2]
                grna_dict['strand'] = '-'
                grna_dict['pam'] = grna[3]

            deletion_index = []
            for grna_seq in grna_counter:
                seq_count = grna_counter[grna_seq]
                if len(seq_count) > 1:
                    for index in seq_count:
                        deletion_index.append(index[0])
            deletion_index.sort(reverse=True)

            for index in deletion_index:
                del grna_dict_list[index]

    return grna_dict_list


cas_dict = cas2dict(args.cas_list)
fasta_dict = fasta2dict(args.fasta)
hits = grna_finder(fasta_dict,cas_dict)

for hit in hits:
    print(hit)
