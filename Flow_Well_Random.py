# Leave one compound out for testing/validation?
## DataFlowSort.py
## Moves files into a structure that works with data flow, 

import csv
import os
import shutil
import numpy as np
import random

class LeaveOneOut:
    ## This version randomly assigns XX% to validation and the rest to training based on wells and nothing else
    def __init__(self,
                labels_path = '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/dataframe.csv',
                output_dir = '/home/jovyan/scratch-shared/Ebba/Kinase_Leave_One_Out',
                image_dir= '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/MiSyHo299',
                image_name ='/%s.png', #Where %s is the image number,
                validation_set_size  = 0.20, #Percentage written as decimal,
                included_classes = [], #Empty for all classes included,
                class_index = 11,
                well_index = 3,
                index_to_leave_out = 6,
                divide_by_index = 5,
                image_number_index = 1,
                name_to_leave_out = "" ,
                output_size = 1 # Percentage of original total size that should be used,
                ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.image_dir = image_dir
        self.image_name  = image_name
        self.validation_set_size = validation_set_size
        self.included_classes = included_classes
        self.class_index =  class_index 
        self.well_index =  well_index
        self.index_to_leave_out = index_to_leave_out
        self.divide_by_index =  divide_by_index
        self.image_number_index = image_number_index
        self.name_to_leave_out =      name_to_leave_out
        self.output_size = output_size

    def main(self):
        print("starting levave on out")

        with open(self.labels_path, 'r') as read_obj:
            csv_reader = csv.reader(read_obj, delimiter=";")
            csv_list = list(csv_reader)
            header = csv_list.pop(0) #remove header

            classes_to_include = self.included_classes
            train_rows, validation_rows, test_rows = self.get_randomized_sets_leave_one_out(csv_list, classes_to_include=classes_to_include )
            
            if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
                shutil.rmtree(self.output_dir)

            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)
                print("made the output dir")

            with open(self.output_dir + "/Labels.csv", 'w', newline = '') as new_labels_file:
                wr = csv.writer(new_labels_file, delimiter=",")
                wr.writerow(header)
                wr.writerows(train_rows)
                wr.writerows(validation_rows)
                wr.writerows(test_rows)

            for row in train_rows:
                if row != header:
                    self.sort_into_class_folders(row, "Train")
            for row in validation_rows:
                if row != header:
                    self.sort_into_class_folders(row, "Validation")
            for row in test_rows:
                if row != header:
                    self.sort_into_test_folder(row, "Test")
            
        print("Finished leave one out")
    
    def run(self):
        self.main()

    def update_settings(self,
                labels_path = '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/dataframe.csv',
                output_dir = '/home/jovyan/scratch-shared/Ebba/Kinase_Leave_One_Out',
                image_dir= '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/MiSyHo299',
                image_name ='/%s.png', #Where %s is the image number,
                validation_set_size  = 0.20, #Percentage written as decimal,
                included_classes = [], #Empty for all classes included,
                class_index = 11,
                well_index = 3,
                index_to_leave_out = 6,
                divide_by_index = 5,
                image_number_index = 1,
                name_to_leave_out = "" ,
                output_size = 1 # Percentage of original total size that should be used,
                ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.image_dir = image_dir
        self.image_name  = image_name
        self.validation_set_size = validation_set_size
        self.included_classes = included_classes
        self.class_index =  class_index 
        self.well_index =  well_index
        self.index_to_leave_out = index_to_leave_out
        self.divide_by_index =  divide_by_index
        self.image_number_index = image_number_index
        self.name_to_leave_out =      name_to_leave_out
        self.output_size = output_size

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

    def get_training_validation_rows(self,compound_dictionary):
        well_training_rows = []
        well_validation_rows = []

        well_keys = np.array(list(compound_dictionary.keys()))
        data_size =len(well_keys)
        if(data_size == 1):
            raise Exception("Note enough data to have both a validation and a training entry. Key: " + str(well_keys))

        validation_set_size = int(data_size * self.validation_set_size)
        if validation_set_size <= 0 :
            validation_set_size = 1

        training_set_size = int(data_size -validation_set_size)
        indices = np.arange(data_size)
        np.random.shuffle(indices)
        well_training_keys = well_keys[indices[:training_set_size]]
        well_validation_keys = well_keys[indices[training_set_size:training_set_size + validation_set_size]]

        for key in well_training_keys:
            well_training_rows.append(compound_dictionary[key])
        for key in well_validation_keys:
            well_validation_rows.append(compound_dictionary[key])
        
        well_training_rows = [item for sublist in well_training_rows for item in sublist]
        well_validation_rows = [item for sublist in well_validation_rows for item in sublist]
        return well_training_rows, well_validation_rows

    def get_randomized_sets_leave_one_out(self, csv_list, classes_to_include):

        val_train_dict = {}
        test_rows = []
        validation_rows = []
        train_rows = []
        included_rows = []
        leave_out = self.name_to_leave_out

        for entry in csv_list:
            class_for_row = entry[self.class_index] 
            leave_out_entry = entry[self.index_to_leave_out]
            divide_by = entry[self.divide_by_index]
            well = entry[self.well_index]
            if leave_out_entry == leave_out:
                test_rows.append([entry])
                continue
            if len(classes_to_include)==0 or class_for_row in classes_to_include :
                if divide_by not in val_train_dict:
                    val_train_dict[divide_by] = {}
                if class_for_row not in val_train_dict[divide_by]:
                    val_train_dict[divide_by][class_for_row] = {}
                if well not in val_train_dict[divide_by][class_for_row]:
                    val_train_dict[divide_by][class_for_row][well] = []
                val_train_dict[divide_by][class_for_row][well].append(entry)

        for key in val_train_dict.keys():
            for class_key in val_train_dict[key]:
                new_training_rows, new_validation_rows = self.get_training_validation_rows(val_train_dict[key][class_key])
                train_rows.append(new_training_rows) 
                validation_rows.append(new_validation_rows)

        train_rows = [item for sublist in train_rows for item in sublist]
        validation_rows = [item for sublist in validation_rows for item in sublist]
        test_rows = [item for sublist in test_rows for item in sublist]

        return train_rows, validation_rows,test_rows

if __name__ == "__main__":
    LeaveOneOut().main()
