import csv
import os
import shutil
import numpy as np
import random
import pandas as pd

class MakeKFolds:
   
    def __init__(self,
                labels_path = '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/dataframe.csv',
                output_dir = '/home/jovyan/Inputs/TEST_Kinase_compound_K_folds_one_by_one/',
                exclude_images_path = "/home/jovyan/Inputs/Kinase_Flagged_Sites/KinaseInhibitor_CP_and_Aut.csv",
                include_groups = ['TK','CMGC','AGC'], #Empty for everything included,
                include_header = 'group',
                exclude_groups = ['P009063','P009083'], #Empty for everything included,
                exclude_header = 'plate',
                class_column_header = 'group',
                intact_group_header = 'compoundname',
                intact_control_group_headers = ['plate', 'well'], # NOTE: hard coded for 2 headers to to troubles with dataframe
                meta_data_header = ['plate', 'well', 'site'],
                image_number_heading = "nr",   
                has_controls = False,
                frac_of_controls_to_use = 0.20
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

    def main(self):
        print("Started get info.")
      
        df = pd.read_csv(self.labels_path , delimiter= ";")
       
        #k_folds = self.get_k_folds(df)

        k_folds = []
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
            df_fold.to_csv(self.output_dir + "k_fold_"+ str(fold_number)+".csv")
            fold_number = fold_number + 1


        print("Finished. Find output in: " + self.output_dir)


    def get_k_folds(self, df):
        k_folds = []

        df_used = df[df[self.include_header].isin(self.included_groups) & ~df[self.exclude_header].isin(self.exclude_groups)]
       
        df_bad_images = pd.read_csv(self.exclude_images_path , delimiter= ";")
        df_bad_images.columns= df_bad_images.columns.str.lower()
        df_do_not_use = pd.merge(df_used,df_bad_images, on = self.meta_data_header, how = "left" ) # rename metadata_well to well in pre prossessing step
        df_do_not_use = df_do_not_use[df_do_not_use["total"] == 1 ]
        df_used = df_used[df_used[self.image_number_heading].isin(df_do_not_use[self.image_number_heading])]
       

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
            k_fold = df_control[(df_control[self.intact_control_group_headers[0]] == combo[0]) & (self.intact_control_group_headers[1] == combo[1])]
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
