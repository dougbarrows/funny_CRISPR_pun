#!/usr/bin/env python3

import sys
import os
import re

# Function that takes gff3filename as input
# Outputs a bed file for the TSS

def gff3_to_TSSbed(gff3filename):

	
	# Determine the output file name
	outputfilename = os.path.splitext(gff3filename)[0]+".TSS.bed"

	# fields in a gff3 files
	field_str = "seqid source genetype start end score strand phase attributes"
	fields = field_str.split(' ')
	
	# open gff3 file, read in lines with gene information
	with open(gff3filename, 'r') as fo, open(outputfilename, 'w') as outputfile:
		
		for line in fo:
			line = line.strip()
			
			# skip header lines
			if line.startswith('#'):
				continue
			# get info from non-header lines
			else:
				data = dict(zip(fields,line.split('\t')))		
				# only keep lines that are of type "gene"
				# print relevant columns to file

				if data['genetype'] == 'gene':

					attributes = data['attributes'].split(";")
					gene_name = attributes[0]+";"+attributes[1]
					strand = data['strand']
										
					bed_line = "chr{}\t{}\t{}\t{}\t{}\t{}\n"
					outputfile.write(bed_line.format(data['seqid'], data['start'], data['start'], gene_name, '', strand))
					
	

def main():
	
	gff3filename = "Homo_sapiens.GRCh38.98.gff3"
	gff3_to_TSSbed(gff3filename)
	
if __name__ == "__main__":
	main()
