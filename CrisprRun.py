#!/usr/bin/env python3

# requires 2 inputs (gene name, species)

import os, sys, re, glob, shutil
import subprocess
import gzip

usage = '\n\n\tusage: {} \n\n\n'

if len(sys.argv) < 3:
	print('Input files are too few!')
	sys.stderr.write(usage)
	sys.exit(1)
if len(sys.argv) > 3:
	print('Too many input files. Or you may have a space in one of your entries. Replace the space with an underscore.')
	sys.stderr.write(usage)
	sys.exit(1)

if re.search(r'\s', sys.argv[1]) and len(sys.argv > 3):
	print('Your gene name can not contain spaces. If a space is needed, please try an underscore (\'_\').')
	sys.stderr.write(usage)	
	sys.exit(1)
if re.search(r'\s', sys.argv[2]) and len(sys.argv > 3):
	print('Your species name can not contain spaces. If a space is needed, please use an underscore (e.g., \'homo_sapiens\').')
	sys.stderr.write(usage)	
	sys.exit(1)

def main():

	# capture command line argument
	geneName = sys.argv[1]
	speciesName = sys.argv[2]

	currPath = os.getcwd()
	currFile = currPath + '/' + speciesName + '.bed'	
	if os.path.isfile(currFile):
		print('A reference file already exists for your species. Is it OK to use that file?')
		userIn1 = input('yes/no')
		if userIn1 == 'yes':
			#this is where we will tag for input to Doug/John's function
			print('The current file to be used is:', currFile)
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
			os.remove(preGz)
			print('File ready for next step!')
	
	sys.exit(0)


if __name__ == '__main__':
	main()
