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
                class_index = 10,
                well_index = 3,
                index_to_leave_out = 6,
                divide_by_index = 5,
                image_number_index = 1,
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
        self.output_size = output_size
        self.used_wells_dict = {}
        self.divisions_dict = {}
        self.runs = 0
        self.header = ""
        self.longest_class = 0

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
                runs = 0,
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
        self.output_size = output_size
        self.runs = runs
    
    def run(self):
        print("starting levave on out")

        if self.runs == 0:
            with open(self.labels_path, 'r') as read_obj:
                csv_reader = csv.reader(read_obj, delimiter=";")
                csv_list = list(csv_reader)
                self.header = csv_list.pop(0) #remove header

                self.set_divisions(csv_list)
               
        if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print("made the output dir")

        train_rows, validation_rows, test_rows = self.get_randomized_sets_leave_one_out()
                
        with open(self.output_dir + "/Labels.csv", 'w', newline = '') as new_labels_file:
            wr = csv.writer(new_labels_file, delimiter=",")
            wr.writerow(self.header)
            wr.writerows(train_rows)
            wr.writerows(validation_rows)
            wr.writerows(test_rows)

        for row in train_rows:
            if row != self.header:
                self.sort_into_class_folders(row, "Train")
        for row in validation_rows:
            if row != self.header:
                self.sort_into_class_folders(row, "Validation")
        for row in self.test_rows:
            if row != header:
                self.sort_into_test_folder(row, "Test")
        self.runs += 1
        print("Finished leave one out")
    
    def set_divisions(self, csv_list):

        division_dict = {}
        used_wells = {}

        for entry in csv_list:
            class_for_row = entry[self.class_index] 
            leave_out_entry = entry[self.index_to_leave_out]
            divide_by = entry[self.divide_by_index]
            well = entry[self.well_index]
            if class_for_row == '':
                continue
            if len(self.classes_to_include)==0 or class_for_row in self.classes_to_include :
                if divide_by not in division_dict:
                    division_dict[divide_by] = {}
                    used_wells[divide_by] = {}
                if class_for_row not in division_dict[divide_by]:
                    division_dict[divide_by][class_for_row] = {}
                    used_wells[divide_by][class_for_row] = []
                if well not in division_dict[divide_by][class_for_row]:
                    division_dict[divide_by][class_for_row][well] = []
                division_dict[divide_by][class_for_row][well].append(entry)

        self.divisions_dict = division_dict
        self.used_wells_dict = used_wells

        longest_class = 0
        for key in division_dict.keys():
            for class_key in division_dict[key]:
                if longest_class < len(division_dict[key][class_key]):
                    longest_class = len(division_dict[key][class_key])
        self.longest_class = longest_class

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

    def get_training_validation_rows(self,compound_dictionary, ignore_well = ""):
        well_training_rows = []
        well_validation_rows = []

        well_keys = np.array(list(compound_dictionary.keys()))
        available_wells = [w for w in well_keys if w not in ignore_well] 
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

    def get_randomized_sets_leave_one_out(self):

        val_train_dict = {}
        test_rows = []
        validation_rows = []
        train_rows = []
        included_rows = []
        leave_out = self.name_to_leave_out

        for divide_by_key in self.division_dict.keys():
            for class_key in self.division_dict[divide_by_key].keys():
                wells_for_class = self.division_dict[divide_by_key][class_key]
                used_wells_for_class = self.used_wells[divide_by_key][class_key]
                if len(wells_for_class) == len(used_wells_for_class):
                    used_wells_for_class = []
                available_wells = [w for w in self.division_dict[divide_by_key][class_key].keys() if w not in used_wells_for_class] 
                chosen_well = random.choice(available_wells)
                test_rows.append([self.division_dict[divide_by_key][class_key][chosen_well]])
                used_wells_for_class.append(chosen_well)
                used_wells[divide_by_key][class_key] = used_wells_for_class

                new_training_rows, new_validation_rows = self.get_training_validation_rows(val_train_dict[divide_by_key][class_key], chosen_well)
                train_rows.append(new_training_rows) 
                validation_rows.append(new_validation_rows)
                
        train_rows = [item for sublist in train_rows for item in sublist]
        validation_rows = [item for sublist in validation_rows for item in sublist]
        test_rows = [item for sublist in test_rows for item in sublist]

        return train_rows, validation_rows,test_rows

if __name__ == "__run__":
    LeaveOneOut().run()
