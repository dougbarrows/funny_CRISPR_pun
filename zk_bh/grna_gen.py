#! /usr/bin/env python3


import argparse, re

parser = argparse.ArgumentParser(description='Imports fasta, extracts possible gRNAs')
parser.add_argument('-f','--fasta',required=True,help='Fasta file to extract gRNAs from')
parser.add_argument('-c','--cas_list',required=True,help='Tab delimitted list of Cas protein parameters')
parser.add_argument('--gc_upper',default=.8,help='GC content upper bound (in decimal) of gRNA sequences to not consider, default = .8')
parser.add_argument('--gc_lower',default=.1,help='GC content lower bound (in decimal) of gRNA sequences to not consider, defulat = .1')
args = parser.parse_args()


def cas2dict(cas_prot_tsv):

    # import a cas protein tsv to get cas specific information
    with open(cas_prot_tsv, 'r') as cas_file:
        cas_dict = {}
        for line in cas_file:
            # prep the lines to work with by removing irrelevant white space
            if line[0] != '#':
                if line[-1] == '\n':
                    line = line.rstrip()
                # populate the dictionary of cas protein information
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


def grna_finder(fasta_dict,cas_dict,gc_upper_bound,gc_lower_bound):

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

    # creating an empty list to populate with dictionaries of gRNAs
    grna_dict_list = []

    # for each cas protein in the input tsv, gather the info in each column
    for cas in cas_dict:
        pam = cas_dict[cas][0]
        pam_orientation = cas_dict[cas][1]
        grna_length = int(cas_dict[cas][2])

        # populating PAM strings for concatenation
        pam_re = ''
        rev_pam_re = ''

        # for each character in the inputted PAMs
        for char in pam:
            # if it is a standard nucleotide, add it to 'pam_re'
            if char in ['A','T','C','G']:
                pam_re += char
                rev_pam_re = reverse_nt[char] + rev_pam_re
            # if it isn't a standard nucleotide, refer to the ambiguous dictionary to obtain the \
            # nucleotides the nonstandard nucleotide is composed
            else:
                pam_re += '[' + ambiguous_nt[char] + ']'
                rev_pam_re = '[' + reverse_nt[char] + ']' + rev_pam_re

        # create a dictionary for counting duplicate gRNAs
        grna_counter = {}

        # for each sequence in the fasta dictionary
        for seq_name in fasta_dict:
            sequence = fasta_dict[seq_name]['sequence']
            # if the PAM sequence is on the 3' end of the gRNA
            if pam_orientation == "3'":
                # match a regular expression containing word characters equivalent to the gRNA len \
                # and obtain the PAM nucleotides from the regular expression ready converted sequence
                grna_matches = re.finditer(r'(?=(\w{' + str(grna_length) + '}' + pam_re + '))', sequence)
                # match a regular expression containing the prepared reverse PAM regular expression \
                # and obtain word characters afterward up to the gRNA len
                rev_grna_matches = re.finditer(r'(?=('+ rev_pam_re +'\w{' + str(grna_length) + '}))', sequence)
                # populate the match information from the iterable regular expression into a tuple \
                # index 0 is the start index of the gRNA PAM combo, 1 is the end index, 2 is the gRNA \
                # 3 is the PAM sequence - remember this is with respect to the orientation of the PAM \
                # relative to the gRNA
                match_info = [(match.start(1), match.end(1), match.group(1)[:grna_length], match.group(1)[grna_length:]) for match in grna_matches]
                # same for the reverse sequence
                rev_match_info = [(rev_match.start(1), rev_match.end(1), rev_match.group(1)[len(pam):], rev_match.group(1)[:len(pam)]) for rev_match in rev_grna_matches] 
            # if the PAM sequence is on the 5' end of the gRNA
            elif pam_orientation == "5'":
                # this is essentially the reciprocal of the above 'if' command
                grna_matches = re.finditer(r'(?=(' + pam_re + '\w{' + str(grna_length) + '}))', sequence)
                rev_grna_matches = re.finditer(r'(?=(\w{' + str(grna_length) + '}'+rev_pam_re+'))', sequence)
                match_info = [(match.start(1), match.end(1), match.group(1)[len(pam):], match.group(1)[:len(pam)]) for match in grna_matches]
                rev_match_info = [(rev_match.start(1), rev_match.end(1), rev_match.group(1)[:grna_length], rev_match.group(1)[grna_length:]) for rev_match in rev_grna_matches]

            print('\n'+str(len(match_info)+len(rev_match_info)),'potential',cas,'gRNAs that match its PAM requirements.')
            # create a new index to populate index values for deletion from the compiled dict
            # we will be deleting any gRNA that is a complete match cuz who wants a 100% off \
            # target ( if it's you, then whatever dude, go figure your life out )
            index = -1
            gc_removal_count = 0
            at_removal_count = 0
            for grna in match_info:
                g_count = grna[2].count('G')
                c_count = grna[2].count('C')
                gc_count = g_count + c_count
                gc_content = gc_count/len(grna[2])
                if gc_content > gc_upper_bound:
                    gc_removal_count += 1
                    continue
                if gc_content < gc_lower_bound:
                    at_removal_count += 1
                    continue
                # add 1 to the index for each iteration
                index += 1
                # if the gRNA is a key in the counter dictionary
                # append the index of the current seq to the gRNA counter
                if grna_counter.get(grna[2]):
                    grna_counter[grna[2]].append([index,'+'])
                else:
                # otherwise, add a new key for the gRNA and populate it with the current index
                    grna_counter[grna[2]] = [[index,'+']]
                # append a dictionary to the list of gRNAs and populate it with .bed specific \
                # parameters
                grna_dict_list.append({})
                grna_dict = grna_dict_list[-1]
                grna_dict['chrom'] = seq_name
                grna_dict['chromStart'] = grna[0]
                grna_dict['chromEnd'] = grna[1]
                grna_dict['grna_seq'] = grna[2]
                grna_dict['strand'] = '+'
                grna_dict['pam'] = grna[3]
                grna_dict['name'] = cas

            # now continue with the reverse sequences
            for grna in rev_match_info:
                g_count = grna[2].count('G')
                c_count = grna[2].count('C')
                gc_count = g_count + c_count
                gc_content = gc_count/len(grna[2])
                if gc_content > gc_upper_bound:
                    gc_removal_count += 1
                    continue
                if gc_content < gc_lower_bound:
                    at_removal_count += 1
                    continue
                index += 1
                if grna_counter.get(grna[2]):
                    grna_counter[grna[2]].append([index,'-'])
                else:
                    grna_counter[grna[2]] = [[index,'-']]
                grna_dict_list.append({})
                grna_dict = grna_dict_list[-1]
                grna_dict['chrom'] = seq_name
                grna_dict['chromStart'] = grna[0]
                grna_dict['chromEnd'] = grna[1]
                # time to create the reverse complements of each sequence
                # these will still be reported as the start and stop indices of the sense \
                # sequence in the .bed file, though these are the antisense sequences
                before_com_seq = list(grna[2])
                for nt_index in range(len(before_com_seq)):
                    before_com_seq[nt_index] = reverse_nt[before_com_seq[nt_index]]
                com_grna_seq = ''.join(before_com_seq)
                rev_grna_seq = com_grna_seq[::-1]
                grna_dict['grna_seq'] = rev_grna_seq
                grna_dict['strand'] = '-'
                before_com_pam = list(grna[3])
                # get the antisense for the PAM as well
                for nt_index in range(len(before_com_pam)):
                    before_com_pam[nt_index] = reverse_nt[before_com_pam[nt_index]]
                com_pam_seq = ''.join(before_com_pam)
                rev_pam_seq = com_pam_seq[::-1]
                grna_dict['pam'] = rev_pam_seq
                grna_dict['name'] = cas

            # delete 100% gRNA matches by first populating a list of deletion indices
            print('\nRemoving gRNAs that 100% overlap... \ncuz who needs them, amiright?\n')
            deletion_index = []
            for grna_seq in grna_counter:
                seq_count = grna_counter[grna_seq]
                if len(seq_count) > 1:
                    for index in seq_count:
                        deletion_index.append(index[0])
            # so as to not break the index-list relationship, sort the indices from \
            # highest to lowest
            deletion_index.sort(reverse=True)

            print('Removed',len(deletion_index),cas,'gRNAs that have 100% similarity to another.')
            print('Removed',gc_removal_count,cas,'gRNAs that exceed GC upper boundary:',gc_upper_bound)
            print('Removed',at_removal_count,cas,'gRNAs that exceed GC lower boundary:',gc_lower_bound)

            # delete the largest index first and work through the rest
            for index in deletion_index:
                del grna_dict_list[index]
        print('\n'+str(len(grna_dict_list)),'valid gRNAs.')
    return grna_dict_list


def grna_dict2bed(fasta,hits):

    # prepare the header for the output bed string
    output_string = '#chrom\tchromStart\tchromEnd\tguideSeq\tstrand\tpam\tname'
    # name the output file by concatenating the following to the end of the fasta
    with open(fasta+'_grna_hits.bed', 'w') as output_file:
        for grna in hits:
            output_string += '\n'
            output_string += grna['chrom']
            output_string += '\t'+str(grna['chromStart'])
            output_string += '\t'+str(grna['chromEnd'])
            output_string += '\t'+grna['grna_seq']
            output_string += '\t'+grna['strand']
            output_string += '\t'+grna['pam']
            output_string += '\t'+grna['name']
        output_file.write(output_string)


cas_dict = cas2dict(args.cas_list)
fasta_dict = fasta2dict(args.fasta)
hits = grna_finder(fasta_dict,cas_dict,float(args.gc_upper),float(args.gc_lower))
grna_dict2bed(args.fasta,hits)
