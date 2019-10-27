#!/usr/bin/env python3

import sys
import os
import re
import argparse


def singlestring(line):
	
	if len(re.findall(r"(\w+)",line)) > 1:
		return False
	else:
		return True

###############################################
# FUNCTION: Make a dictionary from a bed file
###############################################
def bedfile_to_dict(bedfile, genetype):
	gene_dict = {}

	if genetype == "ENS":
		value = 1
	elif genetype == "genename":
		value = 2
	else:
		print("ERROR! genetype can only be 'ID' (Ensembl gene ID) or 'genename'.")
		exit(1)
	
	# Open bed file, save all the values to a dictionary
	with open(bedfile, 'r') as fo_bed:
		for line in fo_bed:
			line = line.rstrip()
			
			# Parse through data file
			field_str = "chr start end namefield score strand"
			fields = field_str.split(' ')
			genedata = dict(zip(fields,line.split('\t')))
			
			# Get the gene name type that you want, use as IDs for dictionary
			namefield = genedata.pop('namefield')
			namefield_regex = re.search(r"^ID=gene:(\w+);Name=(\w+)" , namefield)	
			name = namefield_regex.group(value)

			gene_dict[name] = genedata

	# Return the dectionary
	return gene_dict

######################################################################################
# FUNCTION: Given a gene name and a gene dictionary, print out a bed file format line#
######################################################################################
def gene_to_bedline(gene, gene_dict):

	geneinfo = gene_dict[gene]
	bedline_to_print = "{}\t{}\t{}\t{}\t{}\t{}".format(geneinfo['chr'],geneinfo['start'],geneinfo['start'],gene,0,geneinfo['strand'])
	return(bedline_to_print)

#####################
# Main function
#####################

def main():

	bedfile = "Homo_sapiens.GRCh38.98.chromosome.9.TSS.bed"

	genelistfile = "genelist.txt"
	
	genetype = "genename"

	# Make a dictionary of all of the genenames from the TSS bed file
	gene_dict = bedfile_to_dict(bedfile, "genename")

	####################################################################
	# Output the bed file of only genes of interest
	####################################################################
	
	with open(genelistfile, 'r') as fo_genelist:
		for rawline in fo_genelist:
			# The line should contain only one string. If more than one, quit and provide error
			line = rawline.rstrip()
			if not singlestring(line):
				print('\nERROR! Your genelist file is not formatted correctly. Each line should contain one gene name.!!\n')
				exit(1)
			gene = line
			
			
			bedline_to_print = gene_to_bedline(gene, gene_dict)
			print(bedline_to_print)	
		

if __name__ == "__main__":
	main()
