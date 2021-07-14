import csv
import os
import shutil
import numpy as np
import random
import pandas as pd

class MakeKFolds:
   
    def __init__(self,
                labels_path = '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/dataframe.csv',
                output_dir = '/home/jovyan/Inputs/Kinase_compound_K_folds_one_by_one/',
                exclude_images_path = '/home/jovyan/Inputs/Kinase_Flagged_Sites/QC_KinaseInhibitors_OnlyFlaggedAut_AllPlates.csv',
                include_groups = ['control', 'TK','CMGC','AGC'], #Empty for everything included,
                include_header = 'group',
                exclude_groups = ['P009063','P009083'], #Empty for everything included,
                exclude_header = 'plate',
                class_column_header = 'group',
                intact_group_header = 'compoundname',
                intact_control_group_header = 'well',
                meta_data_header = ['plate', 'well', 'site'],
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
        self.class_column_header =  class_column_header 
        self.intact_group_header = intact_group_header
        self.intact_control_group_header = intact_control_group_header
        self.has_controls = has_controls,
        self.frac_of_controls_to_use = frac_of_controls_to_use

    def main(self):
        print("Started get info.")
      
        df = pd.read_csv(self.labels_path , delimiter= ";")
        df_used = df[df[self.include_header].isin(self.included_groups) & ~df[self.exclude_header].isin(self.exclude_groups)]
       
        df_bad_images = pd.read_csv(self.exclude_images_path , delimiter= ";")
        df_bad_images.columns= df_bad_images.columns.str.lower()
        df_do_not_use = pd.merge(df_used,df_bad_images, on = self.meta_data_header, how = "left" ) # rename metadata_well to well in pre prossessing step
        df_do_not_use = df_do_not_use[df_do_not_use["total"] == 1 ]
        df_used = df_used[df_used[self.image_number_heading].isin(df_do_not_use[self.image_number_heading])]
       
        k_folds = self.get_k_folds(df_used)

        ##Write out data 
        if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print("made the output dir")   

        fold_number = 1
        for k_fold in k_folds:
            df_fold = k_folds[k_fold] 
            df_fold.to_csv(self.output_dir + "k_fold_"+ str(fold_number)+".csv")
            fold_number = fold_number + 1

        print("Finished. Find output in: " + self.output_dir)


    def get_k_folds(self, df_used):
        k_folds = []
        group_n = {}
        df_used_wells =  pd.DataFrame()

        for group in self.included_groups:
            test = df_used[df_used[self.include_header].isin([group])]
            if self.has_controls and group == 'control':
                test = test.groupby(self.intact_control_group_header)
                group_n[group] = int(int(test[[self.include_header]].count().count())) 
            else:
                group_n[group] = int(test.count()[[self.intact_group_header]]) 
                

        for group in self.included_groups:
            number_of_examples = group_n[group]
            for group_sample_number in number_of_examples:
                df_group = df_used[df_used[self.include_header].isin([group])]
                if self.has_controls and group == 'control':
                    df_sampled, df_used_wells, df_used = self.getControlSampel( df_group, df_used_wells, df_used, 1)
                    df_fold = df_fold.append(df_sampled)
                else:
                    df_group = df_group.sample(n = 1)
                    df_fold = df_fold.append(df_group)
            df_used = pd.concat([df_used, df_fold, df_fold]).drop_duplicates(keep=False)
            k_folds.extend(df_fold)

        if self.has_controls and df_used[df_used[self.class_column_header] == 'control'].empty:
            df_group = df_used[df_used[self.include_header].isin(['control'])].copy()
            controls_left_to_do = True
            while  controls_left_to_do: # for length of df_group
                controls_left_to_do, df_sampled, df_used_wells, df_used = self.getControlSampel(df_group, df_used_wells, df_used, 1)
                k_folds.extend(df_sampled)                        
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
