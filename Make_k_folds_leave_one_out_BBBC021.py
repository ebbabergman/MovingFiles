import csv
import os
import shutil
import numpy as np
import random
import pandas as pd
import math

import General

class MakeKFolds:
   

    def __init__(self,
                labels_path = '/home/jovyan/scratch-shared/Ebba/BBBC021_Filtered_Data/Labels.csv',
                exclude_images_path = "", # Empty for none
                output_dir = '/home/jovyan/Inputs/BBBC021_LeaveOneOut_Kfolds',
                include_groups = [], #Empty for everything included,
                include_header = 'moa',
                exclude_groups = ["Cholesterol-lowering","Eg5 inhibitors"], #Empty for everything included,
                exclude_header = 'moa',
                class_column_header = 'moa',
                meta_data_header = ['plate', 'well', 'site'],
                image_number_heading = "image_number",
                intact_group_header = 'compound',
                has_controls = False,
                frac_of_controls_to_use = 0.0,
                intact_control_group_headers = ['plate', 'well'], # NOTE: hard coded for 2 headers to to troubles with dataframe
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

    def main(self):
        print("Started get info.")
      
        df_base = pd.read_csv(self.labels_path , delimiter= ",")
        df_base.dropna(subset = [self.class_column_header], inplace=True)
        df_base.drop_duplicates(inplace=True)

        self.included_groups = self.get_included_groups(df_base)
        df = df_base[df_base[self.include_header].isin(self.included_groups) & ~df_base[self.exclude_header].isin(self.exclude_groups)]
       
        k_folds = self.get_k_folds(df)
        
        if self.has_controls:
            controll_k_folds = self.get_k_folds_control(df)
            k_folds.extend(controll_k_folds)

        print("Made " + str(len(k_folds)) +" k-folds")

        ##Write out data 
        if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print("made the output dir")   

        print("Starting to write to files")
        fold_number = 1
        for df_fold in k_folds:
            df_fold.to_csv(self.output_dir + "k_fold_"+ str(fold_number)+".csv", index = False)
            fold_number = fold_number + 1


        print("Finished. Find output in: " + self.output_dir)

    def get_included_groups(self, df):
        included_groups = self.included_groups
        if(len(included_groups) == 0):
            included_groups = df[self.include_header].unique()
            
        included_groups = [group for group in included_groups if group not in self.exclude_groups]
        return included_groups

    def get_k_folds(self, df):
        k_folds = []
        df_used = df[df[self.include_header].isin(self.included_groups) & ~df[self.exclude_header].isin(self.exclude_groups)]
       
        df_used = General.use_only_good_images(self.exclude_images_path,self.image_number_heading, self.meta_data_header, df_used)
       
        #group by metadataheaders except sites 
        df_used = df_used[df_used[self.include_header].isin(self.included_groups)]
        groups = df_used[self.intact_group_header].unique()

        for group in groups: 
            group_rows = df_used[df_used[self.intact_group_header] == group]
            k_folds.append(group_rows)

        return k_folds

    def get_k_folds_control(self, df):
        k_folds = []
        df_control = df[df['type'] == "control"]
        unique_combos = df_control[df_control.columns & self.intact_control_group_headers].drop_duplicates().to_numpy()
        for combo in unique_combos:
            k_fold = df_control[(df_control[self.intact_control_group_headers[0]] == combo[0]) & (df_control[self.intact_control_group_headers[1]] == combo[1])]
            k_folds.append(k_fold)

        return k_folds


    def getControlSampel(self, df_group, df_used_wells, df_used, n_sample):
        if(df_group[self.intact_control_group_header].count() == 0):
            return False
        sampled_well = np.random.choice(df_group[self.intact_control_group_header].unique(), n_sample)
        df_sampled = df_group[df_group[self.intact_control_group_header].isin(sampled_well)]
        df_used_wells = df_used_wells.append(df_sampled)
        return True, df_sampled, df_used_wells, df_used

if __name__ == "__main__":
    MakeKFolds().main()