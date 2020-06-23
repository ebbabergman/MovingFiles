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

class LeaveOneOut:
    
    def __init__(self,
    labels_path = '/home/jovyan/kensert_CNN/bbbc021_labels.csv',
    output_dir = '/home/jovyan/scratch-shared/Ebba/Leave_One_Out',
    image_dir= '/home/jovyan/kensert_CNN/images_bbbc021',
    image_name ='/bbbc021_%s.png', #Where %s is the image number,
    validation_set_size = 0.20, #Percentage written as decimal,
    included_classes = [], #Empty for all classes included,
    #moa_to_leave_out = "DNA damage",
    #included_classes = ['Aurora kinase inhibitors', 'Eg5 inhibitors','DNA replication'] #Empty for all classes included,
    moa_to_leave_out = "DNA replication",
    compound_to_leave_out = "AZ-U",
    leave_out_moa = False,
    output_size = 1 # Percentage of original total size that should be used,
    ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.image_dir = image_dir
        self.image_name  = image_name
        self.validation_set_size = validation_set_size
        self.included_classes = included_classes
        self.moa_to_leave_out = moa_to_leave_out
        self.compound_to_leave_out = compound_to_leave_out
        self.leave_out_moa = leave_out_moa,
        self.output_size = output_size

    def update_settings(self,
        labels_path = '/home/jovyan/kensert_CNN/bbbc021_labels.csv',
        output_dir = '/home/jovyan/scratch-shared/Ebba/Leave_One_Out',
        image_dir= '/home/jovyan/kensert_CNN/images_bbbc021',
        image_name ='/bbbc021_%s.png', #Where %s is the image number,
        validation_set_size = 0.20, #Percentage written as decimal,
        included_classes = [], #Empty for all classes included,
        #moa_to_leave_out = "DNA damage",
        #included_classes = ['Aurora kinase inhibitors', 'Eg5 inhibitors','DNA replication'] #Empty for all classes included,
        moa_to_leave_out = "DNA replication",
        compound_to_leave_out = "AZ-U",
        leave_out_moa = False,
        output_size = 1 # Percentage of original total size that should be used,
        ):
            self.labels_path = labels_path
            self.output_dir = output_dir
            self.image_dir = image_dir
            self.image_name  = image_name
            self.validation_set_size = validation_set_size
            self.included_classes = included_classes
            self.moa_to_leave_out = moa_to_leave_out
            self.compound_to_leave_out = compound_to_leave_out
            self.leave_out_moa = leave_out_moa,
            self.output_size = output_size

    ##Assumes row structure is ['image_number', 'compound', 'concentration', 'moa', 'plate', 'well', 'replicate']
    def sort_into_class_folders(self,row, category): #where category is train, validation or test
        if(str(row[3]) == 'moa') :     #ignore header
            print("reached header!!!!")
            return
        current_path = self.image_dir + self.image_name  % str(row[0])
    
        dir_path = self.output_dir+"/"  + category +"/" + str(row[3]) 
        target_path = dir_path +"/" +str(row[0]) + ".png"

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print(str(row))
        shutil.copyfile(current_path, target_path)

    def sort_into_test_folder(self,row, category): #where category is train, validation or test
        if(str(row[3]) == 'moa') :     #ignore header
            print("reached header!!!!")
            return
        current_path = self.image_dir + self.image_name  % str(row[0])
    
        dir_path = self.output_dir+"/"  + category + "/" +category #dataflow needs a subfolder, but test subfolder should not be class
        target_path = dir_path +"/" +str(row[0]) + ".png"

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print(str(row))
        shutil.copyfile(current_path, target_path)

    def sort_into_one_folder(self,row):
        if(str(row[3]) == 'moa') :     #ignore header
            print("reached header!!!!")
            return
        current_path = self.image_dir + self.image_name  % str(row[0])
    
        dir_path = self.output_dir + "/images"
        target_path = dir_path +"/" +str(row[0]) + ".png"

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print(str(row))
        shutil.copyfile(current_path, target_path)

    def get_randomized_sets_leave_one_out(self,csv_list, classes_to_include):

        nested_dict = {}
        test_rows = []
        validation_rows = []
        train_rows = []
        included_rows = []

        for entry in csv_list:
            moa = entry[3] 
            compound = entry[1]
            if len(classes_to_include)==0 or moa in classes_to_include :
                if moa not in nested_dict:
                    nested_dict[moa] = {}
                if compound not in nested_dict[moa]:
                    nested_dict[moa][compound] = []
                nested_dict[moa][compound].append(entry)

        leave_out = self.compound_to_leave_out
        for moa in nested_dict:
            compound_dict = nested_dict[moa]
            if self.leave_out_moa and moa == self.moa_to_leave_out:
                leave_out = random.choice(list(compound_dict.keys()))
            for compound in compound_dict:
                if compound == leave_out:
                    if True:
                    # if random.random() < 0.5 :
                        test_rows = test_rows + compound_dict[compound]
                    # else: 
                        #validation_rows = validation_rows + compound_dict[compound]
                else:
                    included_rows =included_rows + compound_dict[compound]

        data_size =len(included_rows)
        validation_set_size = int(data_size * self.validation_set_size + 1)
        training_set_size = int(data_size -validation_set_size)
        indices = np.arange(data_size)
        np.random.shuffle(indices)

        train_rows =np.array(included_rows) [indices[:training_set_size]]
        validation_rows=np.array(included_rows) [indices[training_set_size:training_set_size + validation_set_size]]

        return train_rows, validation_rows,test_rows

    def run(self):
        print("starting levave on out")

        with open(self.labels_path, 'r') as read_obj:
            # pass the file object to reader() to get the reader object
            csv_reader = csv.reader(read_obj, delimiter=";")
            csv_list = list(csv_reader)
            header = ['image_number', 'compound', 'concentration', 'moa', 'plate', 'well', 'replicate']
            if(str(csv_list[0][3]) == 'moa'):
                csv_list.pop(0) #remove header

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
                self.sort_into_class_folders(row, "Train")
            for row in validation_rows:
                self.sort_into_class_folders(row, "Validation")
            for row in test_rows:
                self.sort_into_test_folder(row, "Test")
            
        print("Finished leave one out")