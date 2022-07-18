import csv
import os
from pickle import FALSE
import shutil
from xxlimited import Str
import numpy as np
import random
import pandas as pd
from pandas.core import indexing
import math

# Make K-folds, including their train and validation parts, with csv files as outpus


class MakeKFolds:
                   
    def __init__(self,
                labels_path = '/home/jovyan/Data/Specs/Labels.csv',
                output_dir='/home/jovyan/Inputs/SPECS_Nuclei_Cutoff_CP_AUT_K_folds/',
                include_groups = [], #Empty for everything included,
                #include_groups = ["negcon","heat shock response signalling agonist", "phosphodiesterase inhibitor", "methyltransferase inhibitor","DILI","HDAC inhibitor","topoisomerase inhibitor", "mTOR inhibitor","NFkB pathway inhibitor","JAK inhibitor","pregnane x receptor agonist"], #Empty for everything included,
                #include_groups = ["negcon","DNA polymerase inhibitor", "mTOR inhibitor", "topoisomerase inhibitor"],
                include_header = "selected_mechanism",
                class_column_header = "selected_mechanism",
                exclude_groups = ["poscon","empty"],
                exclude_groups_header = "pert_type",
                exclude_images_path = "/home/jovyan/Data/Specs/Flaggs/Cell_Profiler_Flagged_images_outside_nuclei_cut_82_149.csv",
                intact_group_header = 'compound_id',
                unique_sample_headers = ['plate', 'well', 'site'],
                image_number_heading = "nr",   
                k_folds = "3",
                divide_on_header = 'compound_id',
                make_train_valid = True,
                valid_fraction = 0.25 # 1 = 100%,  Percentage of images remaining afte the test set has been excluded
                ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.exclude_images_path = exclude_images_path
        self.included_groups = include_groups
        self.include_header = include_header
        self.exclude_groups = exclude_groups
        self.exclude_groups_header = exclude_groups_header
        self.unique_sample_headers = unique_sample_headers
        self.image_number_heading = image_number_heading
        self.class_column_header =  class_column_header 
        self.intact_group_header = intact_group_header
        self.k_folds = int(k_folds)
        self.divide_on_header = divide_on_header
        self.make_train_valid = make_train_valid
        self.valid_fraction = valid_fraction

    def main(self):
        print("Started make k-folds.")

        df_base = pd.read_csv(self.labels_path , delimiter= ",")
        df_base.dropna(subset = [self.class_column_header], inplace=True)
        df_base.drop_duplicates(inplace=True)

        all_groups_that_could_be_included = self.get_included_groups(df_base)
        df = df_base[df_base[self.include_header].isin(all_groups_that_could_be_included) & ~df_base[self.exclude_groups_header].isin(self.exclude_groups)]
        self.included_groups = self.get_included_groups(df)

        print("Included and exlcuded groups")

        k_folds_test = self.get_k_folds_test(df)

        if self.make_train_valid:
            k_fold_train, k_fold_validation=  self.get_k_folds_tv(df, k_fold_test, self.valid_fraction)

        ##Make some statistics 
        s_statistics = df.groupby(self.class_column_header)[self.divide_on_header].nunique()
        df_statistics = pd.DataFrame(s_statistics.index)
        df_statistics["Total"] = s_statistics.values

        number_of_folds = self.k_folds
        for k_fold in range(0,number_of_folds):
            df_fold = k_folds_test[k_fold] 
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
            df_fold = k_folds_test[k_fold] 
            df_fold.to_csv(self.output_dir + "k_fold_"+ str(k_fold +1)+".csv", index = False)

        print("Finished. Find output in: " + self.output_dir)

    def get_included_groups(self, df):
        included_groups = self.included_groups
        if(len(included_groups) == 0):
            included_groups = df[self.include_header].unique()
            
        included_groups = [group for group in included_groups if group not in self.exclude_groups]
        return included_groups

    def get_k_folds_test(self, df):
        number_of_folds = self.k_folds
        k_fold_frac = 1/number_of_folds
        k_folds = [None]*number_of_folds
        group_n = {}
        df_unused = df.copy()

        for group in self.included_groups:
            df_group = df[df[self.include_header].isin([group])]
            # TODO START HERE
            ## Is there an error here?
            group_n[group] = math.floor(df_group[self.divide_on_header].nunique()*k_fold_frac)
            if group_n[group] < 1: 
                group_n[group] = 1
                print("A group did not have enough unique groupings to have 1 unique entry per fold. Group:" + str(group) )


        for k_fold in range(1,number_of_folds +1):
            df_fold= pd.DataFrame()
            for group in self.included_groups:
                df_group = df_unused[df_unused[self.include_header].isin([group])]
                if df_group.empty:
                    df_used = df[df[self.include_header].isin([group])]
                    df_unused.append(df_used)
                    df_group = df_unused[df_unused[self.include_header].isin([group])]
                    print("Every unit from group had been used, re-using values for group " + str(group) + ".")

                unique_entries = df_group[self.intact_group_header].unique()
                group_choice = np.random.choice(unique_entries, size = group_n[group])
                df_group_coice = df_group[df_group[self.divide_on_header].isin(group_choice)]
                df_fold = df_fold.append(df_group_coice)
            df_unused = pd.concat([df_unused, df_fold, df_fold]).drop_duplicates(keep=False)
            k_folds[k_fold] = df_fold

       # TODO use  df_unused but put them into the other fold iterativel ( I feel like I already did this somewhere..)

        print("Made test sets for k-folds")    
        return k_folds

    def get_k_folds_tv(self, df, k_fold_test):
        number_of_folds = self.k_folds
        k_fold_train = [None]*number_of_folds
        k_fold_validation = [None]*number_of_folds

        # TODO implement the actual function
        return k_fold_train, k_fold_validation


if __name__ == "__main__":
    MakeKFolds().main()