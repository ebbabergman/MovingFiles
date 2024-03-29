import os
import pandas as pd
import shutil

def make_non_existing_path(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(path)

def make_new_empty_folder(path):
    if  os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)
    print(path)

def use_only_good_images(exclude_images_path,image_number_heading, meta_data_header, df_used, total_flags_column = "total"):
    if(exclude_images_path == ""):
        return df_used
    
    df_bad_images = pd.read_csv(exclude_images_path , delimiter= ";")
    df_bad_images.columns= df_bad_images.columns.str.lower()
    df_do_not_use = pd.merge(df_used,df_bad_images, on = meta_data_header, how = "left" )
    df_do_not_use = df_do_not_use[df_do_not_use[total_flags_column] == 1 ] ## total ==1 means at least one flag has been raised and the image should be excluded
    df_used = df_used[~df_used[image_number_heading].isin(df_do_not_use[image_number_heading])]
    return df_used

def get_included_groups(include_header,included_groups,exclude_groups, df):
    if(len(included_groups) == 0):
        included_groups = df[include_header].unique()
        
    included_groups = [group for group in included_groups if group not in exclude_groups]
    return included_groups