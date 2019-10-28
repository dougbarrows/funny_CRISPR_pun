#!/usr/bin/env python3

import argparse
import os, sys, re, glob, shutil
import subprocess
import gzip
from gff3_to_TSSbed import gff3_to_TSSbed
from genelist_to_TSS_bed import TSSs_of_interest_bed
from overlap import range_overlap

parser = argparse.ArgumentParser(description="""
	This program identifies a list of guide RNAs for CRISPR-based gene attenuation or promotion\n
	for a given list of genes and associated species. The two required inputs are a gene or a \n
	list of genes and the species of interest (see pattern requirements below). The output is a\n
	.bed file that contains a list of gRNAs for each gene. 
	""")
parser.add_argument('geneNames', help='The specific gene(s) of interest entered as a string or the name of the file with the list of genes.')
parser.add_argument('speciesName', help='The name of the target species. Uses underscore for spaces (e.g., homo_sapiens)')
parser.add_argument('-o', '--output', help='Optional: supply output file directory name, otherwise write to program default', dest = 'out')
parser.add_argument('-u', '--upstream', help='Optional: supply distance upstream of the TSS site to check for gRNAs', dest = 'upStream')
parser.add_argument('-d', '--downstream', help='Optional: supply distance downstream of the TSS site to check for gRNAs', dest = 'downStream')
parser.add_argument('-s', '--sortby', help="Optional: if 'de novo' is 'True' sort by is an option for sorting by score using either 'Doench2016_perc', 'Doench2016_score', 'Moreno_Matos_perc', 'Moreno_Matos_score', 'MIT_specificity'", dest = 'sortBy')
parser.add_argument('deNovo', help="Choice of using a de novo table from this program ('True') or using a user-specified table loaded to the working directory ('False').")
args = parser.parse_args()

geneNames = args.geneNames
speciesName = args.speciesName
speciesName = speciesName.lower()
deNovoOpt = args.deNovo

if args.out:
	outputFileDirName = args.out
if args.upStream:
	upStr = args.upStream
if args.downStream:
	downStr = args.downStream

def main():

	currPath = os.getcwd()
	fileDir = currPath + '/' + speciesName	
	if os.path.exists(fileDir):
		contents = os.listdir(fileDir)
		for files in contents:
			files = files.rstrip()
			if files.endswith('.TSS.bed'):
				print('A reference TSS file already exists for your species. Is it OK to use that file (yes/no)?')
				userIn1 = input()
				if userIn1 == 'yes':
					print('The current file to be used as a TSS reference is:', files)
				else:
					print('This program does not yet have the capability for an external user to upload their own reference file. Sorry!')
					sys.exit(0)

	else:
		print('A TSS reference file is not present. Retrieving necessary files from the Ensembl ftp download site.')
		fileList = 'wget ' + '\'' + 'ftp://ftp.ensembl.org/pub/release-98/gff3/' + speciesName + '/CHECKSUMS' + '\'' 		
		fileList_run = subprocess.run(fileList, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if fileList_run.returncode != 0:
			print("FAILED! Unable to locate list file. Check for existence on the Ensembl ftp site.")
			exit(2)
		else:
			with open('CHECKSUMS') as gZList:
				for lines in gZList:
					lines = lines.rstrip()
					if lines.endswith('.98.gff3.gz'):
						temp = lines.split()
						fileGrab = temp[2]
		
		siteCmd = 'wget ' + '\'' + 'ftp://ftp.ensembl.org/pub/release-98/gff3/' + speciesName + '/' + fileGrab + '\''
		siteCmd_run = subprocess.run(siteCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if siteCmd_run.returncode != 0:
			print("FAILED! Either your species is not available or the species did not match. Please check the Ensembl FTP site and try again")
			exit(2)
		else:
			contents=os.listdir(currPath)
			for files in contents:
				if files[-10:] == '98.gff3.gz':
					preGz = files 				# the file name that matches the desired gff3 file					
					with gzip.open(preGz,'rb') as f_in:
						with open(preGz[:-3], 'wb') as f_out:
							shutil.copyfileobj(f_in, f_out) 	# unzipping the gff3 file
			newName = preGz[:-3]		# removing the '.gz' from the name since it is now unzipped
			os.remove(preGz)			# deleting the .gz file
			os.remove('CHECKSUMS')	# deleting the CHECKSUMS file
			gff3_to_TSSbed(newName)	# this is generating a new .bed file of TSSs for the species -- saved so it can be used again
			newerName = newName[:-4]+'TSS.bed' # taking the output from above function and renaming file to append 'TSS.bed'
			newPath = currPath + '/' + speciesName	# creating a new directory using the species name for storage of input files
			os.mkdir(newPath)
			os.rename(newerName,newPath+'/'+newerName)
			os.rename(newName,newPath+'/'+newName)
			print('{} was successfully downloaded.\n{} was successfully generated.\nBoth files are in: {}\n'.format(newName, newerName, newPath)) 
		
		contents = os.listdir(fileDir)
		fastaOK = 1
		while fastaOK == 1:
			for files in contents:
				files = files.rstrip()
				if files.endswith('.fa'):
					print('A fasta file already exists for your species. Is it OK to use that file (yes/no)?')
					userIn2 = input()
					if userIn2 == 'yes':
						fastaOK = 0
						break
					else:
						print('Downloading file from the Ensembl fasta ftp site.')
						fastaOK = 2
						break
				else:
					fastaOK = 2
		if fastaOK == 2:  
			fileList = 'wget ' + '\'' + 'ftp://ftp.ensembl.org/pub/release-98/fasta/' + speciesName + '/dna/CHECKSUMS' + '\''
			fileList_run = subprocess.run(fileList, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
			if fileList_run.returncode != 0:
				print("FAILED! Unable to locate list file. Check for existence on the Ensembl ftp site.")
				exit(2)
			else: 
				with open('CHECKSUMS') as fastaList:
					for lines in fastaList:
						lines = lines.rstrip()
						if lines.endswith('dna_rm.toplevel.fa.gz'):
							temp = lines.split()
							fileGrab = temp[2]
			siteCmd = 'wget ' + '\'' + 'ftp://ftp.ensembl.org/pub/release-98/fasta/' + speciesName + '/dna/' + fileGrab + '\''
			siteCmd_run = subprocess.run(siteCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			if siteCmd_run.returncode != 0:
				print("FAILED! Either your species is not available or the species did not match. Please check the Ensembl FTP site and try again")
				exit(2)
			else:
				contents=os.listdir(currPath)
				for files in content:
					if files.endswith == 'dna_rm.toplevel.fa.gz':
						fastaGz = files
						with gzip.open(fastaGz,'rb') as f_in:
							with open(fastaGz[:-3],'wb') as f_out:
								shutil.copyfileobj(f_in, f_out)
				newName = fastaGz[:-3]
				os.remove(fastaGz)
				os.remove('CHECKSUMS')
				os.rename(newName,newPath+'/'+newName)
				print('{} was successfully downloaded and is located in {}'.format(newName,newPath))	

	
	if args.out:
		outputDir = currPath + '/' + outputFileDirName 
	else:
		outputDir = currPath + '/' + speciesName + '_OutputFiles'
	os.mkdir(outputDir)
	contents=os.listdir(fileDir)
	for files in contents:
		if files.endswith('.TSS.bed'):
			inLocFile = fileDir+'/' + files
			os.chdir(outputDir)
			geneNamesIn = currPath + '/' + geneNames
			TSSs_of_interest_bed(inLocFile, geneNamesIn)	# this is generating a new .bed file of TSSs specific to the desired genes
			
			overlapInput = files[:-3]+'genesofinterest.bed'
			gRNAs = currPath + '/' + speciesName + '/selectGrnas.txt'
			if args.upStream:
				upStrIn = upStr
			else:
				upStrIn = 100
			if args.downStream:
				downStrIn = downStr
			else:
				downStrIn = 50
			if args.sortBy:
				sortByIn = args.sortBy
			else:
				sortByIn = 'Doench2016'
			range_overlap(overlapInput, gRNAs, 'overLapGRnas_OUTPUT', upStrIn, downStrIn, sortByIn, deNovoOpt)
	 
	
	sys.exit(0)


if __name__ == '__main__':
	main()
