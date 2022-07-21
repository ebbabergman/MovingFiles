# Sort files based on a folder with images and csv files with corresponding names to be sorted

import os
import glob
import shutil
import csv
import numpy as np

import General_Moving



class SortFiles:

    def __init__(self, 
                input_dir = "/home/jovyan/Inputs/SPECS_Nuclei_Cutoff_CP_AUT_K_folds/",
                output_dir = "../Leave_One_Out/",
                image_name_header = "ImageNr",
                class_header = "selected_mechanism",
                current_k_fold = 1):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.image_name_header = image_name_header
        self.class_header = class_header
        self.current_k_fold = str(current_k_fold)
  
    def main(self):
        print("Starting Sort files")
        General_Moving.make_non_existing_path(self.output_dir)

        all_csv_files = glob.glob(self.input_dir+"*.csv")
        #all_files = [f for f in listdir(self.input_directory) if isfile(join(self.input_directory, f))]

        right_k_fold_paths = [path for path in all_csv_files if self.current_k_fold in path]
        train_file_path = [path for path in right_k_fold_paths if "train" in path][0]
        validation_file_path = [path for path in right_k_fold_paths if "valid" in path][0]
        test_file_path = [path for path in right_k_fold_paths if "test" in path][0]

        train_headers, train_data = self.read_data(train_file_path)
        validation_headers, validation_data = self.read_data(validation_file_path)
        test_headers, test_data = self.read_data(test_file_path)

        for row in train_data:
            image_nr_header = train_headers.index(self.image_name_header)
            class_header = train_headers.index(self.class_header)
            self.sort_into_sub_folders(row[image_nr_header],row[class_header], "Train")
        for row in validation_data:
            image_nr_header = validation_headers.index(self.image_name_header)
            class_header = validation_headers.index(self.class_header)
            self.sort_into_sub_folders(row[image_nr_header],row[class_header], "Validation")
        for row in test_data:
            image_nr_header = test_headers.index(self.image_name_header)
            #Tensorflow's dataflow needs a subfolder, but test subfolder should be "Test" not the actual class
            self.sort_into_sub_folders(row[image_nr_header],"Test", "Test")

    def read_data (self, path):
        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            headers = next(reader)
            data = np.array(list(reader))
        return headers, data

    def sort_into_sub_folders(self, image_number, sub_folder, folder,): 
        if(image_number == ''):
            return
        
        current_path = self.image_dir + self.image_name  % str(image_number)
        dir_path = self.output_dir+"/"  + folder +"/" + str(sub_folder) +"/" 
        target_path = dir_path +str(image_number) + ".png"

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        shutil.copyfile(current_path, target_path)

if __name__ == "__main__":
    SortFiles().main()