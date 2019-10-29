#!/usr/bin/env python3


import os, sys, re, glob, shutil
import subprocess
import gzip
from gff3_to_TSSbed import gff3_to_TSSbed
from genelist_to_TSS_bed import TSSs_of_interest_bed

usage = '''


usage: {}   <gene name>   <species>

# Initial command line checks. The current version allows 2 inputs: gene name
 (this can be a single gene entered at the command line or a list of genes in
 a .txt file), species (this should be entered in lower case and using an
underscore for any spaces -- consistent with the Ensembl file)   

# The main shell for generating a list of transcript start sites for a given CAS, gene, and species. The current version does the following:
# -- Requires 2 inputs (gene name, species)
# -- Checks to see if a TSS for that species exists. If the TSS exists, the existing file is passed to a program that identifies the TSS.
# -- If the TSS does not exist, the Ensembl database is accessed for that species and a new TSS file is generated and then added to the database.
# -- Once the input .bed file is identified containing the TSSs, the file is interesected with the gene list and yields a new .bed file
# -- The next module takes the new .bed file (a list of TSS for each gene) and passes it to a module that identifies the gRNAs for each gene






'''.format(sys.argv[0])

# Initial command line checks. The current version allows 2 inputs: gene name (this can be a single gene entered at the command line or a list of genes in a .txt file), species (this should be entered in lower case and using an underscore for any spaces -- consistent with the Ensembl file)
varIns = 3

if len(sys.argv) < varIns: 				# check number of input files for too few
	print('Input files are too few!')
	sys.stderr.write(usage)
	sys.exit(1)
if len(sys.argv) > varIns:					# check for too many input files, which may be caused by uses of spaces in the species
	print('Too many input files. Or you may have a space in one of your entries. Replace the space with an underscore.')
	sys.stderr.write(usage)
	sys.exit(1)

if re.search(r'\s', sys.argv[1]) and len(sys.argv > varIns):
	print('Your gene name can not contain spaces. If a space is needed, please try an underscore (\'_\').')
	sys.stderr.write(usage)	
	sys.exit(1)
if re.search(r'\s', sys.argv[2]) and len(sys.argv > varIns):
	print('Your species name can not contain spaces. If a space is needed, please use an underscore (e.g., \'homo_sapiens\').')
	sys.stderr.write(usage)	
	sys.exit(1)

def main():

	# capture command line arguments
	geneNames = sys.argv[1]
	speciesName = sys.argv[2]

	currPath = os.getcwd()
	currFile = currPath + '/' + speciesName	
	if os.path.exists(currFile):
		print('A reference file already exists for your species. Is it OK to use that file?')
		userIn1 = input('yes/no: ')
		if userIn1 == 'yes':
			#this is where we will tag for input to Doug/John's function
			print('The current file to be used is in:', currFile)
		else:
			print('This program does not yet have the capability for an external user to upload their own reference file. Sorry!')
			sys.exit(0)
	else:
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
					preGz = files					
					with gzip.open(preGz,'rb') as f_in:
						with open(preGz[:-3], 'wb') as f_out:
							shutil.copyfileobj(f_in, f_out)
			newName = preGz[:-3]
			os.remove(preGz)
			os.remove('CHECKSUMS')
			gff3_to_TSSbed(newName)							# this is generating a new .bed file of TSSs for the species -- saved so it can be used again
			newerName = newName[:-4]+'TSS.bed'
			os.remove(newName)
			newPath = currPath + '/' + speciesName
			os.mkdir(newPath)
			os.rename(newerName,newPath+'/'+newerName)

	fileLocName = currFile + '/'
	contents=os.listdir(fileLocName)
	for files in contents:
		inLocFile = fileLocName + files
		TSSs_of_interest_bed(inLocFile, geneNames)	# this is generating a new .bed file of TSSs specific to the desired genes
		
		dougInput = files[-3:]+'geneofinterest.bed'
	 
	
	sys.exit(0)


if __name__ == '__main__':
	main()
