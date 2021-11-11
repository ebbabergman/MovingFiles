# Leave one compound out for testing/validation?
## DataFlowSort.py
## Moves files into a structure that works with data flow, 

##Backgroung from https://data.broadinstitute.org/bbbc/BBBC021/

# #NOTE: When evaluating accuracy of MOA classification, it is critical to ensure that the cross-validation is set up correctly. 
# #MOA classification is the task of classifying the MOA of an unseen compound. 
# #Therefore, the evaluation should be a leave-one-compound-out cross validation: 
# #in each iteration, hold out one compound (all replicates and at all concentrations), train on the remaining, and test on the held out compound.

import csv
import os
import shutil
import numpy as np
import random
import pandas as pd
import General

class LeaveOneOut:
    
    def __init__(self,
                labels_path = '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/Labels.csv',
                exclude_images_path = "/home/jovyan/Inputs/Kinase_Flagged_Sites/Kinase_Flags_CP_Strict.csv",
                output_dir = '/home/jovyan/Outputs/Kinase_Leave_One_Out_test/',
                save_labels_dir = '/home/jovyan/Outputs/Kinase_Leave_One_Out_test/',
                k_fold_dir = '/home/jovyan/Inputs/Kinase_Family_No_Compound_K_Fold/',
                k_fold_name = "k_fold_%s.csv",#where /%s is the k_fold number
                image_dir= '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/MiSyHo299/',
                image_name ='%s.png', #Where %s is the image number,
                validation_set_size  = 0.20, #Percentage written as decimal,
                include_groups = ['EGFR', 'PIKK','CDK'], #Empty for everything included,
                include_header = 'family',
                exclude_groups = ['P009063','P009083'], #Empty for everything included,
                exclude_header = 'plate',
                class_column_header = 'family',
                meta_data_header = ['plate', 'well', 'site'],
                well_index = 3,
                leave_out_index = 6,
                image_number_heading = "nr",
                name_to_leave_out = "" ,
                k_fold = "1",
                output_size = 1 # Percentage of original total size that should be used,
                ):
        self.labels_path = labels_path
        self.exclude_images_path = exclude_images_path
        self.output_dir = output_dir
        self.k_fold_dir = k_fold_dir
        self.k_fold_name = k_fold_name
        self.image_dir = image_dir
        self.image_name  = image_name
        self.validation_set_size = validation_set_size
        self.included_groups = include_groups
        self.include_header = include_header
        self.exclude_groups = exclude_groups
        self.exclude_header = exclude_header
        self.class_column_header =  class_column_header 
        self.meta_data_header = meta_data_header
        self.well_index =  well_index
        self.leave_out_index =  leave_out_index
        self.image_number_heading = image_number_heading
        self.name_to_leave_out =      name_to_leave_out
        self.output_size = output_size
        self.k_fold = int(k_fold)
        self.save_labels_dir = save_labels_dir
       
    def update_settings(self,
                labels_path = '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/Labels.csv',
                exclude_images_path = "/home/jovyan/Inputs/Kinase_Flagged_Sites/Kinase_Flags_CP_Strict.csv" ,
                output_dir = '/home/jovyan/Outputs/Kinase_Leave_One_Out_test/',
                save_labels_dir = '/home/jovyan/Outputs/Kinase_Leave_One_Out_test/',
                k_fold_dir = '/home/jovyan/Inputs/Kinase_compound_K_folds_one_by_one_Family/',
                k_fold_name = "k_fold_%s.csv",#where /%s is the k_fold number
                image_dir= '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/MiSyHo299/',
                image_name ='%s.png', #Where %s is the image number,
                validation_set_size  = 0.20, #Percentage written as decimal,
                include_groups = ['control', 'EGFR', 'PIKK','CDK'], #Empty for everything included,
                include_header = 'family',
                exclude_groups = ['P009063','P009083'], #Empty for everything included,
                exclude_header = 'plate',
                class_column_header = 'family',
                meta_data_header = ['plate', 'well', 'site'],
                well_index = 3,
                leave_out_index = 6,
                image_number_heading = "nr",
                name_to_leave_out = "" ,
                k_fold = "1",
                output_size = 1 # Percentage of original total size that should be used,
                ):
        self.labels_path = labels_path
        self.exclude_images_path = exclude_images_path
        self.output_dir = output_dir
        self.k_fold_dir = k_fold_dir
        self.k_fold_name = k_fold_name
        self.image_dir = image_dir
        self.image_name  = image_name
        self.validation_set_size = validation_set_size
        self.included_groups = include_groups
        self.include_header = include_header
        self.exclude_groups = exclude_groups
        self.exclude_header = exclude_header
        self.class_column_header =  class_column_header 
        self.meta_data_header = meta_data_header
        self.well_index =  well_index
        self.leave_out_index =  leave_out_index
        self.image_number_heading = image_number_heading
        self.name_to_leave_out =      name_to_leave_out
        self.output_size = output_size
        self.k_fold = int(k_fold)
        self.save_labels_dir = save_labels_dir

    def main(self):
        print("Starting leave one out")

        df = pd.read_csv(self.labels_path , delimiter= ",")
        
        if(len(self.included_groups) == 0):
            self.included_groups = df[self.include_header].unique()

        groups = self.included_groups
        df_used = self.get_usable_images(df,groups)
        
        k_fold_file = self.k_fold_dir + self.k_fold_name  % str(self.k_fold)
        df_test_all = pd.read_csv(k_fold_file)
        df_test =df_used[df_used[self.image_number_heading].isin(df_test_all[self.image_number_heading])] 
        
        if df_test.empty:
            df_test =df[df[self.image_number_heading].isin(df_test_all[self.image_number_heading])] 

        df_used = pd.concat([df_used, df_test, df_test]).drop_duplicates(keep=False)

        df_validation = df_used.groupby(self.class_column_header).sample(frac = self.validation_set_size)
        df_train = pd.concat([df_used, df_validation, df_validation]).drop_duplicates(keep=False)

        df["valid"] = df.index.isin(df_validation.index)
        df["train"] = df.index.isin(df_train.index)
        df["test"] = df.index.isin(df_test.index)

        ##Make some statistics 
        df_statistics_base = df[df[self.include_header].isin(groups) & ~df[self.exclude_header].isin(self.exclude_groups)]
        df_statistics_base = df_statistics_base[[self.class_column_header, "valid", "train", "test"]]
        
        df_used = self.get_usable_images(df,groups)
    
        df_statistics =pd.DataFrame(df_statistics_base.groupby(self.class_column_header).count()[["train"]].reset_index().values, columns=["group", "train"])
        df_statistics.rename(columns={"train": "total"}, inplace=True)
        df_statistics["used"] = df_used.groupby(self.class_column_header).count().reset_index()[[self.image_number_heading]]
        df_statistics["train"] = df_statistics_base[df_statistics_base["train"]==1].groupby(self.class_column_header).count().reset_index()[["train"]]
        df_statistics["percentage_train"] = df_statistics["train"] /df_statistics["used"] 
        df_statistics["valid"] = df_statistics_base[df_statistics_base["valid"]==1].groupby(self.class_column_header).count().reset_index()[["valid"]]
        df_statistics["percentage_valid"] = df_statistics["valid"] /df_statistics["used"] 
        df_statistics["test"] = df_statistics_base[df_statistics_base["test"]==1].groupby(self.class_column_header).count().reset_index()[["test"]]
        df_statistics["percentage_test"] = df_statistics["test"] /df_statistics["used"] 

        ## Save information 
        if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

        General.make_non_existing_path(self.output_dir)
        General.make_non_existing_path(self.save_labels_dir)

        df_save = df[df[self.include_header].isin(groups)& ~df[self.exclude_header].isin(self.exclude_groups)]
        df_save.to_csv(self.output_dir + "/Labels.csv", index = False)
        df_statistics.to_csv((self.output_dir + "/LabelStatistics.csv"), index = False)

        df_save.to_csv(self.save_labels_dir + "/Labels.csv", index = False)
        df_statistics.to_csv(self.save_labels_dir + "/LabelStatistics.csv", index = False)

        train_rows = df_train[[self.image_number_heading,self.class_column_header]].to_numpy()
        validation_rows = df_validation[[self.image_number_heading,self.class_column_header]].to_numpy()
        test_rows = df_test[[self.image_number_heading,self.class_column_header]].to_numpy()
        for row in train_rows:
            self.sort_into_class_folders(row[0],row[1], "Train")
        for row in validation_rows:
            self.sort_into_class_folders(row[0],row[1], "Validation")
        for row in test_rows:
            self.sort_into_test_folder(row[0], "Test")
            
        print("Finished leave one out")
    
    def run(self):
        self.main()
        
    def sort_into_class_folders(self, image_number, class_name, category): #where category is train, validation or test
        if(image_number == ''):
            return
        
        current_path = self.image_dir + self.image_name  % str(image_number)
        dir_path = self.output_dir+"/"  + category +"/" + str(class_name) 
        target_path = dir_path +"/" +str(image_number) + ".png"

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        shutil.copyfile(current_path, target_path)

    def sort_into_test_folder(self, image_number, category): #where category is train, validation or test
        if(image_number == ''):
            return
       
        current_path = self.image_dir + self.image_name  % str(image_number)
        dir_path = self.output_dir+"/"  + category + "/" +category #dataflow needs a subfolder, but test subfolder should not be class
        target_path = dir_path +"/" +str(image_number) + ".png"

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        shutil.copyfile(current_path, target_path)

    def sort_into_one_folder(self, image_number):
        current_path = self.image_dir + self.image_name  % str(image_number)
    
        dir_path = self.output_dir + "/images"
        target_path = dir_path +"/" +str(image_number) + ".png"

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        shutil.copyfile(current_path, target_path)

    def get_usable_images(self,df,groups):
        df_used = df[df[self.include_header].isin(groups) & ~df[self.exclude_header].isin(self.exclude_groups)]
        df_used = self.use_only_good_images(df_used)
        return df_used

    def use_only_good_images(self, df_used):
        if(self.exclude_images_path == ""):
            return df_used
        
        df_bad_images = pd.read_csv(self.exclude_images_path , delimiter= ";")
        df_bad_images.columns= df_bad_images.columns.str.lower()
        df_do_not_use = pd.merge(df_used,df_bad_images, on = self.meta_data_header, how = "left" )
        df_do_not_use = df_do_not_use[df_do_not_use["total"] == 1 ] ## total ==1 means at least one flag has been raised and the image should be excluded
        df_used = df_used[~df_used[self.image_number_heading].isin(df_do_not_use[self.image_number_heading])]
        return df_used


if __name__ == "__main__":
    LeaveOneOut().main()