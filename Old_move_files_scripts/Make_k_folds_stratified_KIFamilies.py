import csv
import os
import shutil
import numpy as np
import random
import pandas as pd
from pandas.core import indexing
import math

class MakeKFolds:
                   
    def __init__(self,
                labels_path = "~/Inputs/KinasInhibitors/New_labels/Labels.csv",
                output_dir = '/home/jovyan/Inputs/3_fold_Kinase_3Family_Strict/',
                exclude_images_path = "~/Inputs/Kinase_Flagged_Sites/QC_KinaseInhibitors_OnlyStrictFlags_AllPlates.csv",
                include_groups = ["control","PI3K","EGFR", "PIKK"], #Empty for everything included,
                include_header = 'family',
                exclude_groups = ['P009063','P009083'], #Empty for everything included,
                exclude_header = 'plate',
                class_column_header = 'family',
                intact_group_header = 'compoundname',
                intact_control_group_headers = ['plate', 'well'], # NOTE: hard coded for 2 headers to to troubles with dataframe
                meta_data_header = ['plate', 'well', 'site'],
                image_number_heading = "nr",   
                has_controls = True,
                frac_of_controls_to_use = 1,
                k_folds = "10",
                divide_on_header = 'compoundname',
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
        self.has_controls = has_controls,
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
        df_used_wells =  pd.DataFrame()
        df_unused = df.copy()

        for group in self.included_groups:
            df_group = df[df[self.include_header].isin([group])]
            if self.has_controls and group == 'control':
                group_n[group] = math.floor(df_group[self.intact_group_header].nunique()*k_fold_frac)
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
                
                if self.has_controls and group == 'control':
                    df_sampled, df_used_wells, df = self.getControlSampel( df_group, df_used_wells, df, group_n[group])
                    df_fold = df_fold.append(df_sampled)
                else:
                    unique_entries = df_group[self.intact_group_header].unique()
                    group_choice = np.random.choice(unique_entries, size = group_n[group])
                    df_group_coice = df_group[df_group[self.divide_on_header].isin(group_choice)]
                    df_fold = df_fold.append(df_group_coice)
            df_unused =  df_unused[~df_unused.isin(df_fold)].dropna(how = 'all')
            k_folds[k_fold] = df_fold
        k_folds[number_of_folds-1] = df_unused

        if self.has_controls and df[df[self.class_column_header] == 'control'].empty:
            df_group = df[df[self.include_header].isin(['control'])]
            df_sampled, df_used_wells, df = self.getControlSampel( df_group, df_used_wells, df, group_n['control'])
            k_folds[number_of_folds-1] =  k_folds[number_of_folds-1].append(df_sampled)
            
        return k_folds

    def getControlSampel(self, df_group, df_used_wells, df_used, n_sample):
        if(df_group[self.intact_group_header].count() == 0):
            df_group = df_used_wells
            df_used.append(df_used_wells)
        sampled_well = np.random.choice(df_group[self.intact_group_header].unique(), n_sample)
        df_sampled = df_group[df_group[self.intact_group_header].isin(sampled_well)]
        df_used_wells = df_used_wells.append(df_sampled)
        return df_sampled, df_used_wells, df_used

if __name__ == "__main__":
    MakeKFolds().main()