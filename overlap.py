#!/usr/bin/env python3
# sort_by can be either 'Doench2016_perc', 'Doench2016_score', 'Moreno_Matos_perc', 'Moreno_Matos_score', or 'MIT_specificity'

def range_overlap(user_bed_path, guide_loc_path, output_name, upstream , downstream , sort_by, de_novo):
    print(user_bed_path)
    print(guide_loc_path)
    print(output_name)
    print(upstream)
    print(downstream)
    print(sort_by)
    print(de_novo)    
    #note: if de_novo = True, then sort_by is automatically set to be "mismatch_score"
    
    import pyranges as pr
    import pandas as pd

    #######################
    # read in the bed file from the user (right now this will be a file of TSS's for a specific gene list from Mina's function)
    #######################
    
    
    user_bed = pd.read_csv(user_bed_path, sep = '\t', header = None) # important: this assumes that the bed file has no column names

    # need to have first three columns be called 'Chromosome', 'Start', 'End' for a pyRanges object so here we will change the column names
    column_names_user = user_bed.columns.values
    column_names_user_list = list(column_names_user)
    column_names_user_list_str = [str(i) for i in column_names_user_list]

    column_names_user_list_str[0:6] = ['Chromosome', 'Start', 'End', 'Gene', '5', "Strand"]
    user_bed.columns = column_names_user_list_str
    
    # iterate over the rows and change the start and end positions in the user bed file based on the upstream and downstream arguments
    for index, row in user_bed.iterrows():
        if user_bed.at[index, 'Strand'] == '-':
            if user_bed.at[index, 'Start'] < downstream:
                user_bed.at[index, 'Start'] = 0
                user_bed.at[index, 'End'] += upstream
            else:   
                user_bed.at[index, 'End'] += upstream
                user_bed.at[index, 'Start'] -= downstream
        else:
            if user_bed.at[index, 'Start'] < upstream:
                user_bed.at[index, 'Start'] = 0
                user_bed.at[index, 'End'] += downstream
            else:
                user_bed.at[index, 'Start'] -= upstream
                user_bed.at[index, 'End'] += downstream
    
    user_bed_pyR = pr.PyRanges(user_bed) # convert the panda df to a pyRanges object which is required for the overlap function
    #user_bed_pyR_merge = user_bed_pyR.merge() # if we wanted to collapse overlapping ranges, we could use this, but removed it becuase we lose gene column when we do it
    

    # these are the guides already determined for human genome
    guide_locs = pd.read_csv(guide_loc_path, sep = '\t') # important: this assumes that the guide table has column names, need to change if this isnt true!
    column_names_gloc = guide_locs.columns.values
    column_names_gloc_list = list(column_names_gloc)
    column_names_gloc_list_str = [str(i) for i in column_names_gloc_list]

    column_names_gloc_list_str[0:3] = ['Chromosome', 'Start', 'End']
    guide_locs.columns = column_names_gloc_list_str
    guide_locs_noNaN = guide_locs.fillna(0)
    if de_novo == False:
        guide_locs_noNaN[['Doench2016_perc','Doench2016_score']] = guide_locs_noNaN.fusi.str.split('%', expand=True) 
        guide_locs_noNaN['Doench2016_score'] = guide_locs_noNaN['Doench2016_score'].str.replace(r"[\(\)]","") 
        guide_locs_noNaN[['Moreno_Matos_perc','Moreno_Matos_score']] = guide_locs_noNaN.crisprScan.str.split('%', expand=True)
        guide_locs_noNaN['Moreno_Matos_score'] = guide_locs_noNaN['Moreno_Matos_score'].str.replace(r"[\(\)]","")
        guide_locs_noNaN[['MIT_specificity']] = guide_locs_noNaN[['scoreDesc']]
        guide_locs_noNaN['MIT_specificity'] = guide_locs_noNaN['MIT_specificity'].str.replace(r"[A-Za-z\s\.\-]+","0")

        guide_locs_noNaN[["Doench2016_perc", "Doench2016_score", "Moreno_Matos_perc", "Moreno_Matos_score" ,"MIT_specificity"]] = guide_locs_noNaN[["Doench2016_perc", "Doench2016_score", "Moreno_Matos_perc", "Moreno_Matos_score" ,"MIT_specificity"]].apply(pd.to_numeric)
        
        guide_locs_noNaN_select = guide_locs_noNaN.iloc[ : ,[0,1,2,3,4,5,6,7,11,20,21,22,23,24]]
        
        if sort_by not in ['Doench2016_perc', 'Doench2016_score', 'Moreno_Matos_perc', 'Moreno_Matos_score', 'MIT_specificity']:
            print("The sort_by argument must be 'Doench2016_perc', 'Doench2016_score', 'Moreno_Matos_perc', 'Moreno_Matos_score', 'MIT_specificity', because you failed to comply the default ('Doench2016_perc') will be used")
            sort_by = "Doench2016_perc"
        
        guide_locs_pyR = pr.PyRanges(guide_locs_noNaN_select)
        
    else: 
        sort_by = "mismatch_score"
        
        guide_locs_pyR = pr.PyRanges(guide_locs_noNaN)
        print(type(guide_locs_pyR))
        print(guide_locs_pyR.head())
    print(sort_by)
    
    guide_locs_pyR_overlap = guide_locs_pyR.overlap(user_bed_pyR)
    for k, guide_locs_pyR_overlap_df in guide_locs_pyR_overlap: # convert to pandas with this loop to more easily manipulate the df
        guide_locs_pyR_overlap_df
    print(guide_locs_pyR_overlap_df.head())    
    gene_list = []
    for index_ol, row_ol in guide_locs_pyR_overlap_df.iterrows():
        for index_ub, row_ub in user_bed.iterrows():
            if row_ol[0] == row_ub[0]:
                if row_ol[1] in range(row_ub[1], (row_ub[2] + 1)) or row_ol[2] in range(row_ub[1], (row_ub[2] + 1)):
                    gene_list.append(row_ub[3])
    
    guide_locs_pyR_overlap_df = guide_locs_pyR_overlap_df.assign(Gene = gene_list)
    guide_locs_pyR_overlap_df_sort = guide_locs_pyR_overlap_df.sort_values(["Gene", sort_by], ascending = [True, False])

    guide_locs_pyR_overlap_df_sort.to_csv((output_name + '.txt'), sep = '\t', header = True)
        
        
    return(guide_locs_pyR_overlap_df_sort)

	
