import csv
import os
from pickle import FALSE
import shutil
import numpy as np
import random
import pandas as pd
from pandas.core import indexing
import math

class MakeKFolds:
                   
    def __init__(self,
                labels_path = '/home/jovyan/Data/Specs/Specs_Labels.csv',
                output_dir='/home/jovyan/Inputs/SPECS_Nuclei_Cutoff_K_folds/',
                #include_groups = [], #Empty for everything included,
                include_groups = ["heat shock response signalling agonist", "phosphodiesterase inhibitor", "methyltransferase inhibitor","DILI","HDAC inhibitor","topoisomerase inhibitor", "mTOR inhibitor","NFkB pathway inhibitor","JAK inhibitor","pregnane x receptor agonist"], #Empty for everything included,
                #include_groups = ["negcon","DNA polymerase inhibitor", "mTOR inhibitor", "topoisomerase inhibitor"],
                include_header = "selected_mechanism",
                class_column_header = "selected_mechanism",
                exclude_groups = [],
                exclude_header = "selected_mechanism",
                # Exclude images isn't used. Change? There has to be a simpler way to do this.
                exclude_images_path = "/home/jovyan/Data/Specs/Flaggs/images__outside_nuclei_cut_82_149.csv",
                intact_group_header = 'compound_id',
                intact_control_group_headers = ['plate', 'well'], # NOTE: hard coded for 2 headers to to troubles with dataframe
                meta_data_header = ['plate', 'well', 'site'],
                image_number_heading = "nr",   
                has_controls = True,
                frac_of_controls_to_use = 1,
                k_folds = "3",
                divide_on_header = 'compound_id',
                ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.exclude_images_path = exclude_images_path
        self.included_groups = include_groups
        self.include_header = include_header
        self.exclude_groups = exclude_groups
        self.exclude_header = exclude_header
        self.meta_data_header = meta_data_header
        self.image_number_heading = image_number_heading
        self.class_column_header =  class_column_header 
        self.intact_group_header = intact_group_header
        self.intact_control_group_headers = intact_control_group_headers
        self.has_controls = has_controls
        self.frac_of_controls_to_use = frac_of_controls_to_use
        self.frac_of_controls_to_use = frac_of_controls_to_use
        self.k_folds = int(k_folds)
        self.divide_on_header = divide_on_header


    def main(self):
        print("Started get info.")

        df_base = pd.read_csv(self.labels_path , delimiter= ",")
        df_base.dropna(subset = [self.class_column_header], inplace=True)
        df_base.drop_duplicates(inplace=True)

        self.included_groups = self.get_included_groups(df_base)
        df = df_base[df_base[self.include_header].isin(self.included_groups) & ~df_base[self.exclude_header].isin(self.exclude_groups)]
       
        k_folds = self.get_k_folds(df)

        ##Make some statistics 
        s_statistics = df.groupby(self.class_column_header)[self.divide_on_header].nunique()
        df_statistics = pd.DataFrame(s_statistics.index)
        df_statistics["Total"] = s_statistics.values

        number_of_folds = self.k_folds
        for k_fold in range(0,number_of_folds):
            df_fold = k_folds[k_fold] 
            df_grouped = df_fold.groupby(self.class_column_header)
            statistic_column_header = self.divide_on_header+"s_in_fold_"+str(k_fold)
            df_statistics[statistic_column_header] = df_grouped[self.divide_on_header].nunique().values
   
        print(df_statistics.to_latex())
        ##Write out data 
        if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print("made the output dir")   

        df_statistics.to_csv(self.output_dir + "k_fold_statistics.csv", index = False)
        for k_fold in range(0,number_of_folds):
            df_fold = k_folds[k_fold] 
            df_fold.to_csv(self.output_dir + "k_fold_"+ str(k_fold +1)+".csv", index = False)

        print("Finished. Find output in: " + self.output_dir)

    def get_included_groups(self, df):
        included_groups = self.included_groups
        if(len(included_groups) == 0):
            included_groups = df[self.include_header].unique()
            
        included_groups = [group for group in included_groups if group not in self.exclude_groups]
        return included_groups

    def get_k_folds(self, df):
        number_of_folds = self.k_folds
        k_fold_frac = 1/number_of_folds
        k_folds = [None]*number_of_folds
        group_n = {}
        df_used_wells = pd.DataFrame()
        df_unused = df.copy()

        for group in self.included_groups:
            df_group = df[df[self.include_header].isin([group])]
            if self.has_controls and group == 'negcon':
                unique_controls = len(df_group[self.intact_control_group_headers].groupby(self.intact_control_group_headers))
                group_n[group] = math.floor(unique_controls*k_fold_frac)
            else:
                group_n[group] = math.floor(df_group[self.divide_on_header].nunique()*k_fold_frac)
            if group_n[group] < 1: group_n[group] = 1

        for k_fold in range(0,number_of_folds-1):
            df_fold= pd.DataFrame()
            for group in self.included_groups:
                df_group = df_unused[df_unused[self.include_header].isin([group])]
                if df_group.empty:
                    df_used = df[df[self.include_header].isin([group])]
                    df_unused.append(df_used)
                    df_group = df_unused[df_unused[self.include_header].isin([group])]
                    print("Every unit from group had been used, re-using values for group " + str(group) + ".")
                
                if self.has_controls and group == 'negcon':
                    df_sampled, df_used_wells = self.getControlSampel( df_group, df_used_wells, group_n[group])
                    df_fold = df_fold.append(df_sampled)
                else:
                    unique_entries = df_group[self.intact_group_header].unique()
                    group_choice = np.random.choice(unique_entries, size = group_n[group])
                    df_group_coice = df_group[df_group[self.divide_on_header].isin(group_choice)]
                    df_fold = df_fold.append(df_group_coice)
    #                 pd.merge(df1, df2, on=['A','B'], how='outer', indicator=True).query("_merge != 'both'")
    #    .drop('_merge', axis=1)
    #    .reset_index(drop=True)
            df_unused = pd.concat([df_unused, df_fold, df_fold]).drop_duplicates(keep=False)
            #df_unused =  df_unused[~df_unused.isin(df_fold)].dropna(how = 'all')
            k_folds[k_fold] = df_fold
        k_folds[number_of_folds-1] = df_unused

        # if self.has_controls and df[df[self.class_column_header] == '[dmso]'].empty:
        #     df_group = df[df[self.include_header].isin(['[dmso]'])]
        #     df_sampled, df_used_wells, df = self.getControlSampel( df_group, df_used_wells, df, group_n['[dmso]'])
        #     k_folds[number_of_folds-1] =  k_folds[number_of_folds-1].append(df_sampled)
            
        return k_folds

    def getControlSampel(self, df_group, df_used_wells, n_sample):
        if df_group.empty:
            df_group = df_used_wells

        df_control = pd.DataFrame({'count' : df_group[self.intact_control_group_headers].groupby(self.intact_control_group_headers).size()}).reset_index()
        sampled_well_index = np.random.choice(df_control.index, n_sample)
        sampled_wells = df_control.iloc[sampled_well_index]
        df_sampled= df_group.merge(sampled_wells)
        df_sampled.drop('count', axis = 1, inplace = True)

        df_used_wells = df_used_wells.append(df_sampled)
        return df_sampled, df_used_wells, 

if __name__ == "__main__":
    MakeKFolds().main()