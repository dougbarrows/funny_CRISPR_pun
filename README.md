# PFB 2019 Final Project

## Group members:

- Matthew Waas (TA)

- Doug Barrows
- John Werner
- Hunter Underhill
- Mina Kojima
- Becky Hennessey
- Zach Konkel

Our project is centered around providing users with a set of possible guide RNA locations
for CRISPRa or CRISPRi given a set of gene IDs. 

The approach is split into three parts
1) linking genes and start sites
2) creating a set of possible guides for a given PAM sequence
3) cross-refencing the list of start sites with the list of possible guides



### More information on some functions:

###### FUNCTION: gff3_to_TSSbed

- Usage: ``` gff3_to_TSSbed(gff3filename, outputfolder=".") ```

- INPUT: Takes a gff3 annotation file

- OUTPUT: bed file with the TSS site for all of the genes

  

###### FUNCTION: TSSs_of_interest_bed

- Usage:  ``` TSSs_of_interest_bed(TSSbedfile, genelistfile.txt, outputfolder=".")```

- INPUT: 1) bed file with all TSSs, 2) gene name (or genelist in .txt format)

- OUTPUT: bed file with only the TSSs of interest

- Note about gene names:

  - can take either genenames (such as 'DMRT1') and Ensembl ID names (ENSGxxxxxxxxxx)
  - if a user-provided gene does not match a known genename or gene ID in the annotation file, provides user with a helpful message about the genes without a match

  

###### FUNCTION: convert_catpiss_output_to_gff3

- Usage: ```convert_catpiss_output_to_gff3(catpiss_output, scoretype, outputfolder='.')```

- INPUT: output file from the catpiss program, with all gRNAs

  - also specify the "scoretype" to display in JBrowse

- OUTPUT: gff3 file to use for JBrowse

  

## Visualization of designed sgRNAs in JBrowse

- Go to JBrowse site: [crispr.programmingforbiology.org/jbrowse](crispr.programmingforbiology.org/jbrowse)
- gRNAs colored by score:
  - Blue: score >50
  - Red: score <50




