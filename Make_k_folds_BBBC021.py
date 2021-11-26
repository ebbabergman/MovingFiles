import csv
import os
import shutil
import numpy as np
import random
import pandas as pd
import math

class MakeKFolds:
   
    def __init__(self,
                labels_path = '/home/jovyan/scratch-shared/Ebba/BBBC021_Filtered_Data/Labels.csv',
                output_dir = '/home/jovyan/Inputs/BBBC021_K_folds/',
                include_groups = [], #Empty for everything included,
                include_header = 'moa',
                exclude_groups = ["Cholesterol-lowering","Eg5 inhibitors"], #Empty for everything included,
                exclude_header = 'moa',
                class_column_header = 'moa',
                well_column_header = 'compound',
                k_folds = "3",
                has_controls = False,
                frac_of_controls_to_use = 0.20,
                divide_on_header = 'compound'
                ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.included_groups = include_groups
        self.include_header = include_header
        self.exclude_groups = exclude_groups
        self.exclude_header = exclude_header
        self.class_column_header =  class_column_header 
        self.well_column_header = well_column_header
        self.k_folds = int(k_folds)
        self.has_controls = has_controls
        self.frac_of_controls_to_use = frac_of_controls_to_use
        self.divide_on_header = divide_on_header

    def main(self):
        print("Started get info.")
        entries_list = []

        df = pd.read_csv(self.labels_path , delimiter= ",")
        df.dropna(subset = [self.class_column_header], inplace=True)

        if(len(self.included_groups) == 0):
            self.included_groups = df[self.include_header].unique()
            if self.include_header == self.exclude_header:
                self.included_groups = [group for group in self.include_header if group not in self.exclude_header]

        df_used = df[df[self.include_header].isin(self.included_groups) & ~df[self.exclude_header].isin(self.exclude_groups)]
       
        k_folds = self.get_k_folds(df_used)

        ##Make some statistics 
        df_statistics_base = df_used
        df_statistics_base = df_statistics_base[[self.class_column_header, self.well_column_header]]

        df_statistics =pd.DataFrame(df_statistics_base.groupby(self.class_column_header).count()[[self.well_column_header]].reset_index().values, columns=[self.class_column_header,self.well_column_header])
        df_statistics.rename(columns={self.well_column_header: "total"}, inplace=True)
    
        number_of_folds = self.k_folds
        for k_fold in range(0,number_of_folds):
            df_fold = k_folds[k_fold] 
            df_statistics[str(k_fold)] = df_fold.groupby(self.class_column_header).count().reset_index()[[self.well_column_header]]
       
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


    def get_k_folds(self, df_used):
        number_of_folds = self.k_folds
        k_fold_frac = 1/number_of_folds
        k_folds = [None]*number_of_folds
        group_n = {}
        df_used_wells =  pd.DataFrame()

        for group in self.included_groups:
            df_group = df_used[df_used[self.include_header].isin([group])]
            if self.has_controls and group == 'control':
                group_n[group] = math.floor(df_group[self.well_column_header].nunique()*k_fold_frac)
            else:
                group_n[group] = math.floor(df_group[self.divide_on_header].nunique()*k_fold_frac)
            if group_n[group] < 1: group_n[group] = 1
                ## EBBA TODO  THis groupby header does not seem to be correct, double check and adjust
        for k_fold in range(0,number_of_folds-1):
            df_fold= pd.DataFrame()
            for group in self.included_groups:
                df_group = df_used[df_used[self.include_header].isin([group])]
                if self.has_controls and group == 'control':
                    df_sampled, df_used_wells, df_used = self.getControlSampel( df_group, df_used_wells, df_used, group_n[group])
                    df_fold = df_fold.append(df_sampled)
                    
                else:
                    df_group = df_group.sample(n = group_n[group])
                    df_fold = df_fold.append(df_group)
                df_used = pd.concat([df_used, df_fold, df_fold]).drop_duplicates(keep=False)
            k_folds[k_fold] = df_fold
        k_folds[number_of_folds-1] = df_used

        if self.has_controls and df_used[df_used[self.class_column_header] == 'control'].empty:
            df_group = df_used[df_used[self.include_header].isin(['control'])]
            df_sampled, df_used_wells, df_used = self.getControlSampel( df_group, df_used_wells, df_used, group_n['control'])
            k_folds[number_of_folds-1] =  k_folds[number_of_folds-1].append(df_sampled)
            
        return k_folds

    def getControlSampel(self, df_group, df_used_wells, df_used, n_sample):
        if(df_group[self.well_column_header].count() == 0):
            df_group = df_used_wells
            df_used.append(df_used_wells)
        sampled_well = np.random.choice(df_group[self.well_column_header].unique(), n_sample)
        df_sampled = df_group[df_group[self.well_column_header].isin(sampled_well)]
        df_used_wells = df_used_wells.append(df_sampled)
        return df_sampled, df_used_wells, df_used

if __name__ == "__main__":
    MakeKFolds().main()
