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

class LeaveOneOut:
    
    def __init__(self,
                labels_path = '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/dataframe.csv',
                output_dir = '/home/jovyan/Outputs/Kinase_Leave_One_Out',
                image_dir= '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/MiSyHo299',
                image_name ='/%s.png', #Where %s is the image number,
                validation_set_size  = 0.20, #Percentage written as decimal,
                include_groups = ['control', 'TK','CMGC','AGC'], #Empty for everything included,
                include_header = 'group',
                class_column_header = 'group',
                well_index = 3,
                leave_out_index = 6,
                image_number_index = 1,
                name_to_leave_out = "" ,
                k_fold = "1",
                output_size = 1 # Percentage of original total size that should be used,
                ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.image_dir = image_dir
        self.image_name  = image_name
        self.validation_set_size = validation_set_size
        self.included_groups = include_groups
        self.include_header = include_header
        self.class_column_header =  class_column_header 
        self.well_index =  well_index
        self.leave_out_index =  leave_out_index
        self.image_number_index = image_number_index
        self.name_to_leave_out =      name_to_leave_out
        self.output_size = output_size
        self.k_fold = int(k_fold)
        self.k_folds = {1:['CBK289985', 'CBK278016', 'CBK290869', 'CBK290815', 'CBK278063', 'CBK277927', 'CBK013405', 'CBK293875', 'CBK288321G', 'CBK293873', 'CBK288278', 'CBK293881', 'CBK278132'],
2:['CBK293855', 'CBK277977', 'CBK290917', 'CBK293884', 'CBK288256', 'CBK277968', 'CBK040864', 'CBK293898', 'CBK278119G', 'CBK277922C', 'CBK293867', 'CBK290752', 'CBK290901'],
3:['CBK290262', 'CBK288270', 'CBK290915', 'CBK290847', 'CBK288351C', 'CBK288257', 'CBK293163', 'CBK290480', 'CBK278049G', 'CBK278128G', 'CBK288350', 'CBK293856', 'CBK293872'],
4:['CBK290836', 'CBK290877', 'CBK278057', 'CBK278050', 'CBK288271', 'CBK290830', 'CBK040892', 'CBK278131', 'CBK290809', 'CBK290269', 'CBK277936', 'CBK289904', 'CBK290556'],
5:['CBK293863', 'CBK290253', 'CBK290196C', 'CBK290589', 'CBK293908', 'CBK278067', 'CBK293174', 'CBK290852', 'CBK200518C', 'CBK288279C', 'CBK293911', 'CBK293897', 'CBK288314'],
6:['CBK278013', 'CBK041242C', 'CBK293876', 'CBK040896', 'CBK041143', 'CBK290224', 'CBK288316', 'CBK293852', 'CBK288330G', 'CBK293860G', 'CBK289902', 'CBK290859', 'CBK290457C'],
7:['CBK277987', 'CBK288268', 'CBK293871', 'CBK277983', 'CBK293861', 'CBK278056', 'CBK290768', 'CBK290822G', 'CBK289742', 'CBK290918', 'CBK041209', 'CBK288335', 'CBK277971C'],
8:['CBK293859', 'CBK277983G', 'CBK288323', 'CBK278127', 'CBK278084', 'CBK293879', 'CBK293904', 'CBK288321', 'CBK288266', 'CBK288280', 'CBK277959', 'CBK288311', 'CBK278019'],
9:['CBK290987', 'CBK288292', 'CBK277992', 'CBK288326C', 'CBK200981', 'CBK290853', 'CBK288253C', 'CBK290206', 'CBK288328', 'CBK041168', 'CBK290998', 'CBK290797', 'CBK290484'],
10:['CBK277962', 'CBK293880', 'CBK288272', 'CBK290649', 'CBK289937', 'CBK288255', 'CBK290855', 'CBK288269', 'CBK290794', 'CBK278081', 'CBK288310', 'CBK290671', 'CBK201383'],
}


    def update_settings(self,
                labels_path = '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/dataframe.csv',
                output_dir = '/home/jovyan/scratch-shared/Ebba/Kinase_Leave_One_Out',
                image_dir= '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/MiSyHo299',
                image_name ='/%s.png', #Where %s is the image number,
                validation_set_size  = 0.20, #Percentage written as decimal,
                include_groups = ['control', 'TK','CMGC','AGC'], #Empty for everything included,
                include_index = 10,
                class_index = 10,
                well_index = 3,
                leave_out_index = 6,
                image_number_index = 1,
                name_to_leave_out = "CBK013405" ,
                k_fold = 1,
                k_folds = {},
                output_size = 1 # Percentage of original total size that should be used,
                ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.image_dir = image_dir
        self.image_name  = image_name
        self.validation_set_size = validation_set_size
        self.included_groups = include_groups
        self.include_index = include_index
        self.class_index =  class_index 
        self.well_index =  well_index
        self.leave_out_index =  leave_out_index
        self.image_number_index = image_number_index
        self.name_to_leave_out =      name_to_leave_out
        self.output_size = output_size
        self.k_fold = int(k_fold)
        self.k_folds = {1:['CBK293871', 'CBK289742', 'CBK293861', 'CBK288271', 'CBK278067', 'CBK200518C', 'CBK293879', 'CBK288256', 'CBK288268', 'CBK290998', 'CBK293860G', 'CBK288311', 'CBK290756'],
2:['CBK293884', 'CBK277983', 'CBK290822G', 'CBK293876', 'CBK277987', 'CBK290852', 'CBK200981', 'CBK290815', 'CBK290649', 'CBK288297', 'CBK288310', 'CBK293911', 'CBK290956'],
3:['CBK277992', 'CBK290253', 'CBK278049G', 'CBK278131', 'CBK277927', 'CBK293904', 'CBK288292', 'CBK290830', 'CBK293875', 'CBK278001', 'CBK289902', 'CBK041209', 'CBK288314'],
4:['CBK288316', 'CBK290917', 'CBK288321', 'CBK288351C', 'CBK288330G', 'CBK290480', 'CBK277980', 'CBK040864', 'CBK288328', 'CBK290797', 'CBK277959', 'CBK293858C', 'CBK278129'],
5:['CBK288257', 'CBK293174', 'CBK277991', 'CBK290853', 'CBK040892', 'CBK288327', 'CBK290997C', 'CBK277977', 'CBK293908', 'CBK290859', 'CBK288279C', 'CBK277936', 'CBK290484'],
6:['CBK290196C', 'CBK293887', 'CBK041242C', 'CBK293852', 'CBK278057', 'CBK288269', 'CBK278013', 'CBK288270', 'CBK290224', 'CBK288350', 'CBK293891', 'CBK293867', 'CBK290561'],
7:['CBK290987', 'CBK288253C', 'CBK278119G', 'CBK013405', 'CBK290845', 'CBK293859', 'CBK278127', 'CBK277962', 'CBK040896', 'CBK278128G', 'CBK293881', 'CBK288278', 'CBK293865'],
8:['CBK293863', 'CBK278063', 'CBK289985', 'CBK290589', 'CBK288254', 'CBK290836', 'CBK293163', 'CBK278084', 'CBK293898', 'CBK290823', 'CBK290752', 'CBK290269', 'CBK290585'],
9:['CBK290869', 'CBK288326C', 'CBK290847', 'CBK278056', 'CBK290995', 'CBK290855', 'CBK288272', 'CBK288321G', 'CBK288323', 'CBK289904', 'CBK041168', 'CBK290022', 'CBK290454C'],
10:['CBK277968', 'CBK290267', 'CBK278016', 'CBK290877', 'CBK290206', 'CBK290915', 'CBK041143', 'CBK288255', 'CBK289937', 'CBK277922C', 'CBK288335', 'CBK293864', 'CBK201383'],
}


    def main(self):
        print("Starting leave one out")

        df = pd.read_csv(self.labels_path , delimiter= ";")
        ## Todo, remove test, or mark test somehow in labels
        groups = self.included_groups
        df_used = df[df[self.include_header].isin(groups)]
        
        #df_validation = df_used.group.sample(frac = self.validation_set_size) ##TODO Both lines same results, why ?
        df_validation = df_used.groupby(self.class_column_header).sample(frac=self.validation_set_size) ##TODO Both lines same results, why ?
        df_train = pd.concat([df_used, df_validation, df_validation]).drop_duplicates(keep=False)

        #df[test] = ???
        df["valid"] = df.index.isin(df_validation.index)
        df["train"] = df.index.isin(df_train.index)

        ## Get information of distribution between classes and groupings in the different set
        ## TODO
        ## Save information 
        if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print("made the output dir")

        df.to_csv(self.output_dir + "/Labels.csv")
        self.write_images_to_new_folders(df)
            
        print("Finished leave one out")
    
    def run(self):
        self.main()

    def write_images_to_new_folders(self, df_full):
        print("This is not implemented")
             # with open(self.output_dir + "/Labels.csv", 'w', newline = '') as new_labels_file:
        #     wr = csv.writer(new_labels_file, delimiter=",")
        #     wr.writerow(header)
        #     wr.writerows(train_rows)
        #     wr.writerows(validation_rows)
        #     wr.writerows(test_rows)

        #     for row in train_rows:
        #         if row != header:
        #             self.sort_into_class_folders(row, "Train")
        #     for row in validation_rows:
        #         if row != header:
        #             self.sort_into_class_folders(row, "Validation")
        #     for row in test_rows:
        #         if row != header:
        #             self.sort_into_test_folder(row, "Test")
        
    def sort_into_class_folders(self, row, category): #where category is train, validation or test
        current_path = self.image_dir + self.image_name  % str(row[self.image_number_index])
    
        dir_path = self.output_dir+"/"  + category +"/" + str(row[self.class_index]) 
        target_path = dir_path +"/" +str(row[self.image_number_index]) + ".png"

        if(str(row[self.image_number_index]) == ''):
            return

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print(str(row))
        shutil.copyfile(current_path, target_path)

    def sort_into_test_folder(self, row, category): #where category is train, validation or test
        current_path = self.image_dir + self.image_name  % str(row[self.image_number_index])
    
        dir_path = self.output_dir+"/"  + category + "/" +category #dataflow needs a subfolder, but test subfolder should not be class
        target_path = dir_path +"/" +str(row[self.image_number_index]) + ".png"

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print(str(row))
        shutil.copyfile(current_path, target_path)

    def sort_into_one_folder(self, row):
        current_path = self.image_dir + self.image_name  % str(row[self.image_number_index])
    
        dir_path = self.output_dir + "/images"
        target_path = dir_path +"/" +str(row[self.image_number_index]) + ".png"

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print(str(row))
        shutil.copyfile(current_path, target_path)

if __name__ == "__main__":
    LeaveOneOut().main()
