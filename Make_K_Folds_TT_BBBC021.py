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

# Make K-folds, Test and Train. Validation will be duplicates of Train and not really show anything


class MakeKFolds:
                   
    def __init__(self,
                labels_path = "/home/jovyan/Data/BBBC021/BBBC021_Labels.csv",
                output_dir='/home/jovyan/Inputs/BBBC021_All_Leave_One_Out_Test_Train_Only/',
                #
                include_groups = [], #Empty for everything included,
                include_header = "moa",
                class_column_header = "moa",
                #exclude_groups = [["DMSO","Cholesterol-lowering","Eg5 inhibitors"]],
                exclude_groups = [["DMSO"]],
                exclude_groups_headers = ["moa"],
                exclude_images_path = "",
                intact_group_header = 'compound',
                unique_sample_headers = ["ImageNumber"],
                image_number_heading = "ImageNumber",   
                k_folds = "3",
                divide_on_header = 'compound',
                make_train_valid = True,
                valid_fraction = 0.25, # 1 = 100%,  Percentage of images remaining afte the test set has been excluded
                leave_one_out = True,
                ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.exclude_images_path = exclude_images_path
        self.included_groups = include_groups
        self.include_header = include_header
        self.exclude_groups = exclude_groups
        self.exclude_groups_headers = exclude_groups_headers
        self.unique_sample_headers = unique_sample_headers
        self.image_number_heading = image_number_heading
        self.class_column_header =  class_column_header 
        self.intact_group_header = intact_group_header
        self.k_folds = int(k_folds)
        self.divide_on_header = divide_on_header
        self.make_train_valid = make_train_valid
        self.valid_fraction = valid_fraction
        self.leave_one_out = leave_one_out

    def main(self):
        print("Started make k-folds.")

        self.make_path_available()

        df_base = pd.read_csv(self.labels_path , delimiter= ",")
        df_base.dropna(subset = [self.class_column_header], inplace=True)
        df_base.drop_duplicates(inplace=True)

        all_groups_that_could_be_included = self.get_included_groups(df_base)
        df = df_base[df_base[self.include_header].isin(all_groups_that_could_be_included)]
        for index  in range(0,len(self.exclude_groups_headers)):
            exclude_groups_header = self.exclude_groups_headers[index]
            df = df[~df[exclude_groups_header].isin(self.exclude_groups[index])]

        if len(self.exclude_images_path) > 0:
            df = self.exclude_images(df)
            print("excluded images indicated with file")

        self.included_groups = self.get_included_groups(df)

        print("Included and exlcuded groups")
      
        if self.leave_one_out:
            k_folds_test = self.get_leave_one_out_test(df)
            print("Made " + str(len(k_folds_test)) + " test sets")
        else:
            k_folds_test = self.get_k_folds_test(df)
            df_test_statistics  = self.get_statistics(k_folds_test,df)
            df_test_statistics.to_csv(self.output_dir + "k_fold_test_statistics.csv", index = False)
            print("Test Statistics")
            print(df_test_statistics.to_latex())

        if self.make_train_valid:
            k_folds_train, k_folds_validation=  self.get_k_folds_tv(df, k_folds_test)
            df_validation_statistics  = self.get_statistics(k_folds_validation,df)
            df_train_statistics  = self.get_statistics(k_folds_train,df)
            
            df_validation_statistics.to_csv(self.output_dir + "k_fold_validation_statistics.csv", index = False)
            df_train_statistics.to_csv(self.output_dir + "k_fold_train_statistics.csv", index = False)

            print("Validation Statistics ")
            print(df_validation_statistics.to_latex())
            print("Train Statistics ")
            print(df_train_statistics.to_latex())

        for k_fold in range(0,self.k_folds):
            df_test_fold = k_folds_test[k_fold] 
            df_test_fold.to_csv(self.output_dir + "k_fold_test_"+ str(k_fold +1)+".csv", index = False)
            
            df_validation_fold = k_folds_validation[k_fold] 
            df_validation_fold.to_csv(self.output_dir + "k_fold_validation_"+ str(k_fold +1)+".csv", index = False)
            
            df_train_fold = k_folds_train[k_fold] 
            df_train_fold.to_csv(self.output_dir + "k_fold_train_"+ str(k_fold +1)+".csv", index = False)
           
        print("Finished. Find output in: " + self.output_dir)

    def get_included_groups(self, df):
        included_groups = self.included_groups
        if(len(included_groups) == 0):
            included_groups = df[self.include_header].unique()
            
        included_groups = [group for group in included_groups if group not in self.exclude_groups]
        return included_groups

    def exclude_images(self, df):
        df_bad_images = pd.read_csv(self.exclude_images_path , delimiter= ",")
        df_bad_images = df_bad_images.drop_duplicates(subset = self.unique_sample_headers)

        if 'Total' in df_bad_images.columns:
            bad_image_mask =  df_bad_images["Total"] == 1 
            df_bad_images = df_bad_images[bad_image_mask]
        elif 'total' in df_bad_images.columns:
            bad_image_mask =  df_bad_images["total"] == 1
            df_bad_images = df_bad_images[bad_image_mask]

        df_merged = pd.merge(df,df_bad_images[self.unique_sample_headers], on = self.unique_sample_headers, how = "outer", indicator = True ) 

        merge_both_mask = df_merged["_merge"] == "both" 
        merge_left_mask = df_merged["_merge"] == "left_only" 
        merge_right_mask = df_merged["_merge"] == "right_only" 

        print("Excluding images")
        print(str(df.shape[0]) +  " included images to start with")
        print(str(df_merged[merge_right_mask].shape[0]) + " images were not in the input labels")
        print(str(df_merged[merge_both_mask].shape[0]) + " images will be dropped")
        print(str(df_merged[merge_left_mask].shape[0]) + " images will be kept")

        df_new = df[df[self.image_number_heading].isin(df_merged[merge_left_mask][self.image_number_heading])]
        df = df_new.copy()

        return df


    def get_k_folds_test(self, df):
        number_of_folds = self.k_folds
        k_folds = [None]*number_of_folds
        group_n = {}
        df_unused = df.copy()

        for group in self.included_groups:
            df_group = df[df[self.include_header].isin([group])]
            group_n[group] = math.floor(df_group[self.divide_on_header].nunique()/number_of_folds)
            if group_n[group] < 1: 
                group_n[group] = 1
                print("A group did not have enough unique groupings to have 1 unique entry per fold. Group:" + str(group) )

        for k_fold in range(1,number_of_folds +1):
            df_fold= pd.DataFrame()
            for group in self.included_groups:
                df_group = df_unused[df_unused[self.include_header].isin([group])]
                if df_group.empty:
                    df_used = df[df[self.include_header].isin([group])]
                    df_unused = df_used.copy()
                    df_group = df_unused[df_unused[self.include_header].isin([group])]
                    print("Every unit from group had been used, re-using values for group " + str(group) + ".")
                df_group_coice = self.get_group_selection( df_group, group_n[group])
                df_fold = pd.concat([df_fold,df_group_coice],ignore_index = True)
            df_unused = pd.concat([df_unused, df_fold, df_fold]).drop_duplicates(keep=False)
            k_folds[k_fold-1] = df_fold

        # deal with the unused unique groups
        group_index = 0 
        for group in self.included_groups:
            df_group = df_unused[df_unused[self.include_header].isin([group])]
            
            if df_group.empty:
                continue
            
            k_fold_to_choose = [*range(1,number_of_folds +1)]
            while not df_group.empty:
                if len(k_fold_to_choose) == 0:
                    k_fold_to_choose = [*range(1,number_of_folds +1)]
                    print("Something is wrong with the way the number of groups are positioned, could put remainder into each fold") # TODO make this a warning
                
                chosen_k_fold = np.random.choice(k_fold_to_choose, size = 1, replace = False)[0]
                k_fold_to_choose.remove(chosen_k_fold)
                df_group_coice = self.get_group_selection(df_group, 1)

                k_folds[chosen_k_fold-1] = pd.concat([k_folds[chosen_k_fold-1],df_group_coice],ignore_index = True)
                df_unused = pd.concat([df_unused, df_group_coice, df_group_coice]).drop_duplicates(keep=False)
                df_group = pd.concat([df_group, df_group_coice, df_group_coice]).drop_duplicates(keep=False)
                
            group_index = group_index +1

        if not df_unused.empty:
            print("WARNING: Didn't put all avialiable compound into test k-folds! Left over:")
            print(df_unused)

        print("Made test sets for k-folds")    
        return k_folds

    def get_leave_one_out_test(self, df):
        number_of_folds = df[self.divide_on_header].nunique()
        self.k_folds = number_of_folds
        print("Leave one out will result in " + str(number_of_folds) + " folds.")
        k_folds = [None]*number_of_folds

        leave_out_list = df[self.divide_on_header].unique()

        k_fold = 1
        for entry in leave_out_list:
            df_fold= df[df[self.divide_on_header] == entry]
            k_folds[k_fold-1] = df_fold
            k_fold = k_fold +1

        print("Made test sets for leave one out. Made " + str(k_fold) + "folds")    
        return k_folds

    def get_k_folds_tv(self, df, k_fold_test):
        number_of_folds = self.k_folds
        validation_fraction = self.valid_fraction *(number_of_folds-1)/number_of_folds
        k_fold_train = [None]*number_of_folds
        k_fold_validation = [None]*number_of_folds
        group_n = {}

        for group in self.included_groups:
            df_group = df[df[self.include_header].isin([group])]
            group_n[group] = math.floor(df_group[self.divide_on_header].nunique()*validation_fraction)
            if group_n[group] < 1: 
                group_n[group] = 1
                print("A group did not have enough unique groupings to have 1 unique entry per validation fold. Using 1 unique entry anyway. Group:" + str(group) )

        for k_fold in range(1,number_of_folds +1):
            print("Starting k-fold: " + str(k_fold))
           
            df_fold_validation= pd.DataFrame()
            df_fold_train= pd.DataFrame()
            for group in self.included_groups:
                df_group = df_unused[df_unused[self.include_header].isin([group])]
                
                df_group_coice_validation = self.get_group_selection(df_group, group_n[group])
               
                df_fold_validation = pd.concat([df_fold_validation,df_group_coice_validation],ignore_index = True)
            df_fold_train = pd.concat([df_unused,df_fold_validation,df_fold_validation],ignore_index = True).drop_duplicates(keep=False)
            df_unused = pd.concat([df_unused, df_fold_validation, df_fold_train],ignore_index = True).drop_duplicates(keep=False)
            k_fold_validation[k_fold-1] = df_fold_validation
            k_fold_train[k_fold-1] = df_fold_train
            if not df_unused.empty:
                print("WARNING: Didn't put all available compound into train or valid k-folds! Left over:")
                print(df_unused)
 
        print("Made train and valid sets for k-folds")   

        return k_fold_train, k_fold_validation

    def get_group_selection(self, df_group, number_of_unique_entries):
        unique_entries = df_group[self.intact_group_header].unique()
        group_choice = np.random.choice(unique_entries, size = number_of_unique_entries, replace = False)
        df_group_choice = df_group[df_group[self.divide_on_header].isin(group_choice)]
        return df_group_choice   

    def get_statistics (self, k_folds, df):
        s_statistics = df.groupby(self.class_column_header)[self.divide_on_header].nunique()
        df_statistics = pd.DataFrame(s_statistics.index)
        df_statistics["Total"] = s_statistics.values

        number_of_folds = self.k_folds
        for k_fold in range(0,number_of_folds):
            df_fold = k_folds[k_fold] 
            df_grouped = df_fold.groupby(self.class_column_header)
            statistic_column_header = self.divide_on_header+"s_in_fold_"+str(k_fold +1)
            df_statistics[statistic_column_header] = df_grouped[self.divide_on_header].nunique().values
   
        return df_statistics        

    def make_path_available(self): 
        if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print("made the output dir") 


if __name__ == "__main__":
    MakeKFolds().main()