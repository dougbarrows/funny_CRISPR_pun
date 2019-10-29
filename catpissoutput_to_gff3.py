#!/usr/bin/env python3

# Last updated: October 29, 2019
# Author: Mina Kojima


##############################################################################
# FUNCTION: Takes catpiss outputfile as input, outputs a gff3 file for JBrowse
##############################################################################
def convert_catpiss_output_to_gff3(catpiss_output, scoretype, outputfolder='.'):

	import os
	import re

	# Determine the output file name
	gff3_output = outputfolder+"/"+os.path.splitext(os.path.basename(catpiss_output))[0]+".gff3"

	# Read in input file, output gf3 file
	
	with open(catpiss_output, 'r') as fo_in, open(gff3_output, 'w') as fo_out:
		
		linenumber = 0
		for line in fo_in:
			line = line.rstrip()
			linenumber += 1

			# from header line, get the name of the strings
			if linenumber == 1:
				fields_str = line
				fields = fields_str.split('\t')
			# Read in all the other lines and get information
			else:
				data = dict(zip(fields, line.split('\t')))
				# Find all the data for the chromosomes
				col1 = re.search(r"chr([\w+])",data["Chromosome"]).group(1)
				col2 = "catpiss"				
				col3 = "sgRNA"
				col4 = data["Start"]
				col5 = data["End"]
				col6 = data[scoretype]
				col7 = data["strand"]
				col8 = "."
				name_info = data["Gene"].split(" ")
				if len(name_info) == 2:
					ENS_ID = name_info[0]
					genename = name_info[1]
					Name = "Name=sgRNA_"+ENS_ID+"_"+genename+"_score"+data[scoretype]
					Note = "Note="+ENS_ID+"_"+genename
				else:
					genename = name_info[0]
					Name = "Name=sgRNA_"+genename+"_score_"+data[scoretype]
					Note = "Note="+genename
				ID = "ID=sgRNA_"+data["Chromosome"]+"_"+data["Start"]+"_"+data["End"]
				otherscores = "otherscores=Doench2016_perc:"+data["Doench2016_perc"]
				otherscores += ",Doench2016_score:"+data["Doench2016_score"]
				otherscores += ",Moreno_Matos_perc:"+data["Moreno_Matos_perc"]
				otherscores += ",Moreno_Matos_score:"+data["Moreno_Matos_score"]
				otherscores += ",MIT_specificity:"+data["MIT_specificity"]
				distance_to_tss = "distance_to_tss="+data["GuideEnd_to_TSS"]
				col9_list = [ID, Name, Note, otherscores, distance_to_tss]
				col9 = ";".join(col9_list)
				col_list = [col1, col2, col3, col4, col5, col6, col7, col8, col9]
				
				# Figure out the gff3 format line to print out, print to output file
				line_to_print = "\t".join(col_list)
				fo_out.write(line_to_print+"\n")

	# Print note to user about what file was written
	print("\nConverted gff3 file:\n\t{}\n".format(gff3_output))

	# Return output filename to user
	return gff3_output
	


#########
# MAIN
#########

def main():

	import sys

	if len(sys.argv) < 3:
		print("\n\t!!! Program Usage Error !!!\n")
		print("Usage: {} crisproutputfile scoretype".format(sys.argv[0]))
		sys.exit(1)

	infile = sys.argv[1]
	scoretype = sys.argv[2]

	convert_catpiss_output_to_gff3(infile,scoretype)

if __name__ == "__main__":
	main()
