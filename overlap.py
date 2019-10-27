import pyranges as pr
import pandas as pd

user_bed_path = "testTSSfile.bed"
guide_loc_path = "ch9_hg19.txt"
upstream = 100
downstream = 50

range_overlap(user_bed_path = user_bed_path, guide_loc_path = guide_loc_path, output_name = "mina_testFile", upstream = upstream, downstream = downstream)

def range_overlap(user_bed_path, guide_loc_path, output_name, upstream, downstream):
    
    user_bed = pd.read_csv(user_bed_path, sep = '\t', header = None) # important: this assumes that hte bed file has no column names

    # need to have first three columns be called 'Chromosome', 'Start', 'End' for a pyRanges object so here we will change the column names
    column_names_user = user_bed.columns.values
    column_names_user_list = list(column_names_user)
    column_names_user_list_str = [str(i) for i in column_names_user_list]

    column_names_user_list_str[0:6] = ['Chromosome', 'Start', 'End', '4', '5', "Strand"]
    user_bed.columns = column_names_user_list_str
        
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
    
    user_bed_pyR = pr.PyRanges(user_bed)
    user_bed_pyR_merge = user_bed_pyR.merge()
    
    print(user_bed_pyR_merge)

    # these are the guides already determined for human genome
    guide_locs = pd.read_csv(guide_loc_path, sep = '\t') # important: this assumes that the guide table has column names, need to change if this isnt true!
    column_names_gloc = guide_locs.columns.values
    column_names_gloc_list = list(column_names_gloc)
    column_names_gloc_list_str = [str(i) for i in column_names_gloc_list]

    column_names_gloc_list_str[0:3] = ['Chromosome', 'Start', 'End']
    guide_locs.columns = column_names_gloc_list_str
    guide_locs_noNaN = guide_locs.fillna(0)
    guide_locs_pyR = pr.PyRanges(guide_locs_noNaN)


    guide_locs_pyR_overlap = guide_locs_pyR.overlap(user_bed_pyR_merge)

    guide_locs_pyR_overlap.to_csv((output_name + '.txt'), sep = '\t', header = True)
        
        
    return(guide_locs_pyR_overlap)


    # trying to confirm we have the right overlaps (and not intersect)
    guide_locs_start = guide_locs.loc[: ,'Start']
    guide_locs_end = guide_locs.loc[: ,'End']
    
    end_minus_start_guide = guide_locs["End"] - guide_locs["Start"]
    
    guide_locs_pyR_overlap_inp = pd.read_csv("test_overlap.txt", sep = '\t', header = 0)
    end_minus_start_overlap = guide_locs_pyR_overlap_inp["End"] - guide_locs_pyR_overlap_inp["Start"]

    guide_locs_pyR_intersect_inp = pd.read_csv("test_intersect.txt", sep = '\t', header = 0)
    end_minus_start_intersect = guide_locs_pyR_intersect_inp["End"] - guide_locs_pyR_intersect_inp["Start"]
    
    # make own dataframe to test intersect vs overlap
    
    guide = [['chr1' , 1, 24, '+']]
    user_plus = [['chr1', 10, 20, '+'], ['chr1', 20, 30, '+']]
    user_minus = [['chr1', 10, 20, '-'], ['chr1', 20, 30, '-']]
    
    df_guide = pd.DataFrame(guide, columns = ['Chromosome', 'Start', 'End', 'Strand'])
    df_guide.to_csv("guide_locs_testSet.txt", sep = '\t', header = True, index = False)
    df_user_plus = pd.DataFrame(user_plus, columns = ['Chromosome', 'Start', 'End', 'Strand'])
    df_user_plus["empty1"] = 0
    df_user_plus["empt2"] = 0
    df_user_plus = df_user_plus[["Chromosome", "Start", "End", "empty1", "empt2", "Strand"]]
    df_user_plus.to_csv("user_plus_testSet.txt", sep = '\t', header = None, index = False)
    df_user_minus = pd.DataFrame(user_minus, columns = ['Chromosome', 'Start', 'End', 'Strand'])
    df_user_minus["empty1"] = 0
    df_user_minus["empt2"] = 0
    df_user_minus = df_user_minus[["Chromosome", "Start", "End", "empty1", "empt2", "Strand"]]
    df_user_minus.to_csv("user_minus_testSet.txt", sep = '\t', header = None, index = False)
