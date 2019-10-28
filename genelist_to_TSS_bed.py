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
	
	# if the gene is in the dictionary, return a bed-type line as a string
	if gene in gene_dict:
		geneinfo = gene_dict[gene]
		bedline_to_print = "{}\t{}\t{}\t{}\t{}\t{}".format(geneinfo['chr'],geneinfo['start'],geneinfo['start'],gene,0,geneinfo['strand'])
		return bedline_to_print
	# if the gene is not in the dictionary, return false
	else:
		return False

############################################################################################
# FUNCTION: Given a TSS bed file and a list of genes, output only TSSs of genes of interest
############################################################################################
def TSSs_of_interest_bed(TSSbedfile, genelistfile):

	# Make a list of genes from all of the files
	genelist = []
	# Read genes in the file if the user provided a file
	if genelistfile.endswith(".txt"):
		with open(genelistfile, 'r') as fo_genelist:
			for rawline in fo_genelist:
				line = rawline.rstrip()
				if not singlestring(line): # Quit if line of file contains more than one string
					print('\nERROR! Your genelist file is not formatted correctly. Each line should contain one gene name.!!\n')
					exit(1)
				genelist.append(line)
	# If it's not a file, then read in the genelistfile
	else:
		genelist.append(genelistfile)

	# Check if the genelist has Ensembl gene IDs to determine what "type" of genelist you have
	if genelist[0].startswith("ENS"):
		genenametype = "ENS"
	else:
		genenametype = "genename"

	# Make a dictionary of all of the genes / TSS sites with the 
	gene_dict = bedfile_to_dict(TSSbedfile , genenametype) # make a dictionary

	# Go through the list of genes, make a new bed file with the TSSs of interest
	outputfilename = os.path.splitext(os.path.basename(TSSbedfile))[0]+".genesofinterest.bed"
	#outputfilename = os.path.splitext(TSSbedfile)[0]+".genesofinterest.bed"
	unmappedlist = [] # list of unmapped genes
	with open(outputfilename, 'w') as fo_output:
		for gene in genelist:
			bedline_to_print = gene_to_bedline(gene, gene_dict)
			if bedline_to_print:
				fo_output.write(bedline_to_print+"\n")
			else:
				unmappedlist.append(gene)
	# Print note to user if there were any unmapped gene names / IDs
	# Print out these genes to a new file so user can check the unmapped names
	unmappedcount = len(unmappedlist)
	if unmappedcount > 0:
		unmappedgenesfile =os.path.splitext(os.path.basename(genelistfile))[0]+".unmappedgenes.txt"
		with open(unmappedgenesfile, 'w') as fo_unmapped:
			for unmappedgene in unmappedlist:
				fo_unmapped.write(gene+"\n")
		print("\nNOTE: Your genelist had {} unmapped genes.\nThese genes are listed in {}.\n".format(unmappedcount, unmappedgenesfile))

#####################
# Main function
#####################

def main():

	gff3_to_TSSbed("Homo_sapiens.GRCh38.98.chromosome.9.bed")

	TSSbedfile = "Homo_sapiens.GRCh38.98.chromosome.9.TSS.bed"
	genelistfile = "folder/genelist.txt"

	# Call function to print out a bed file from a bed file with a TSS of interest
	# and with a list of genes of interest
	TSSs_of_interest_bed(TSSbedfile, "STRA8") # can use a single gene
	TSSs_of_interest_bed(TSSbedfile, genelistfile) # can use a .txt file with a gene list

if __name__ == "__main__":
	main()
