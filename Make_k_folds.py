import csv
import os
import shutil
import numpy as np
import random
import pandas as pd

class MakeKFolds:
    ## TODO bryt ut funktioner
   
    def __init__(self,
                labels_path = '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/dataframe.csv',
                output_dir = '/home/jovyan/Outputs/Kinase_Leave_One_Out',
                image_dir= '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/MiSyHo299',
                image_name ='/%s.png', #Where %s is the image number,
                include_groups = ['control', 'TK','CMGC','AGC'], #Empty for everything included,
                include_header = 'group',
                class_column_header = 'group',
                other_header_with_numbers = 'well',
                k_folds = "10",
                ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.image_dir = image_dir
        self.image_name  = image_name
        self.included_groups = include_groups
        self.include_header = include_header
        self.class_column_header =  class_column_header 
        self.other_header_with_numbers = other_header_with_numbers
        self.k_folds = int(k_folds)

    def main(self):
        print("Started get info.")
        entries_list = []

        with open(self.labels_path, 'r') as read_obj:
            df = pd.read_csv(self.labels_path , delimiter= ";")
            ## Todo, remove test, or mark test somehow in labels
            groups = self.included_groups
            df_used = df[df[self.include_header].isin(groups)]

            k_folds = self.get_k_folds(df_used)

     
            ## todo: handle controls s o that not all of them end up in test - percentage value at top? Filter away from used and then just go?
             ##Make some statistics 
            df_statistics_base = df[df[self.include_header].isin(groups)]
            df_statistics_base = df_statistics_base[[self.class_column_header, self.other_header_with_numbers]]

            df_statistics =pd.DataFrame(df_statistics_base.groupby(self.class_column_header).count()[[self.other_header_with_numbers]].reset_index().values, columns=[self.class_column_header,self.other_header_with_numbers])
            df_statistics.rename(columns={self.other_header_with_numbers: "total"}, inplace=True)
            
            print (df_statistics)

            number_of_folds = self.k_folds
            for k_fold in range(0,number_of_folds):
                df_fold = k_folds[k_fold] 
                df_statistics[str(k_fold)] = df_fold.groupby(self.class_column_header).count().reset_index()[[self.other_header_with_numbers]]
           
            print (df_statistics)
            
            if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
                shutil.rmtree(self.output_dir)

            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)
                print("made the output dir")

        ## Write output 
        # file_object = open(self.output_dir + self. + ".txt", "w+")
        # for entry in entries_list:
        #     file_object.write("\"" output_file_name+ str(entry)+ "\"" +"\n") 

        print("Finished. Find output in: " + self.output_dir)


    def get_k_folds(self, df_used):
        number_of_folds = self.k_folds
        k_fold_frac = 1/number_of_folds
        k_folds = [None]*number_of_folds
        group_n = {}

        for group in self.included_groups:
            test = df_used[df_used[self.include_header].isin([group])]
            group_n[group] = int(test.count()[[self.other_header_with_numbers]]*k_fold_frac) 

        for k_fold in range(0,number_of_folds-1):
            df_fold= pd.DataFrame()
            for group in self.included_groups:
                df_group = df_used[df_used[self.include_header].isin([group])]
                df_group = df_group.sample(n = group_n[group])
                df_fold.append(df_group)
                df_fold.append(pd.DataFrame(df_group))
            df_used = pd.concat([df_used, df_fold, df_fold]).drop_duplicates(keep=False)
            k_folds[k_fold] = df_fold
        k_folds[number_of_folds-1] = df_used

        return k_folds

if __name__ == "__main__":
    MakeKFolds().main()
