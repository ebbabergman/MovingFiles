import csv
import os
import shutil
import numpy as np
import random
import pandas as pd

class MakeKFolds:
   
    def __init__(self,
                labels_path = '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/dataframe.csv',
                output_dir = '/home/jovyan/Inputs/Kinase_compound_K_folds/',
                 include_groups = ['control', 'TK','CMGC','AGC'], #Empty for everything included,
                include_header = 'group',
                exclude_groups = ['P009063','P009083'], #Empty for everything included,
                exclude_header = 'plate',
                class_column_header = 'group',
                intact_group_header = 'compound',
                intact_control_group_header = 'well',
                k_folds = "10",
                has_controls = False,
                frac_of_controls_to_use = 0.20
                ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.included_groups = include_groups
        self.include_header = include_header
        self.exclude_groups = exclude_groups
        self.exclude_header = exclude_header
        self.class_column_header =  class_column_header 
        self.intact_group_header = intact_group_header
        self.intact_control_group_header = intact_control_group_header
        self.k_folds = int(k_folds)
        self.has_controls = has_controls,
        self.frac_of_controls_to_use = frac_of_controls_to_use

    def main(self):
        print("Started get info.")
        entries_list = []

        
        df = pd.read_csv(self.labels_path , delimiter= ";")
       
        df_used = df[df[self.include_header].isin(self.included_groups) & ~df[self.exclude_header].isin(self.exclude_groups)]
       
        k_folds = self.get_k_folds(df_used)

        ##Make some statistics 
        df_statistics_base = df_used
        df_statistics_base = df_statistics_base[[self.class_column_header, self.intact_group_header]]

        df_statistics =pd.DataFrame(df_statistics_base.groupby(self.class_column_header).count()[[self.intact_group_header]].reset_index().values, columns=[self.class_column_header,self.intact_group_header])
        df_statistics.rename(columns={self.intact_group_header: "total"}, inplace=True)
    
        number_of_folds = self.k_folds
        for k_fold in range(0,number_of_folds):
            df_fold = k_folds[k_fold] 
            df_statistics[str(k_fold)] = df_fold.groupby(self.class_column_header).count().reset_index()[[self.intact_group_header]]
       
        print(df_statistics)
        ##Write out data 
        if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print("made the output dir")   

        df_statistics.to_csv(self.output_dir + "k_fold_statistics.csv")
        for k_fold in range(0,number_of_folds):
            df_fold = k_folds[k_fold] 
            df_fold.to_csv(self.output_dir + "k_fold_"+ str(k_fold +1)+".csv")

        print("Finished. Find output in: " + self.output_dir)


    def get_k_folds(self, df_used):
        number_of_folds = self.k_folds
        k_fold_frac = 1/number_of_folds
        k_folds = [None]*number_of_folds
        group_n = {}
        df_used_wells =  pd.DataFrame()

        for group in self.included_groups:
            test = df_used[df_used[self.include_header].isin([group])]
            if self.has_controls and group == 'control':
                test = test.groupby(self.intact_group_header)
                group_n[group] = int(int(test[[self.include_header]].count().count())*k_fold_frac) 
                if group_n[group] < 1: group_n[group] = 1
            else:
                group_n[group] = int(test.count()[[self.intact_group_header]]*k_fold_frac) 

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
        if(df_group[self.intact_control_group_header].count() == 0):
            df_group = df_used_wells
            df_used.append(df_used_wells)
        sampled_well = np.random.choice(df_group[self.intact_control_group_header].unique(), n_sample)
        df_sampled = df_group[df_group[self.intact_control_group_header].isin(sampled_well)]
        df_used_wells = df_used_wells.append(df_sampled)
        return df_sampled, df_used_wells, df_used

if __name__ == "__main__":
    MakeKFolds().main()
