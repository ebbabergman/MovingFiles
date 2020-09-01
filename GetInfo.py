import csv
import os
import shutil
import numpy as np
import random

class GetInfo:
    ## TODO get 3 biggest classes
    ## print information

    def __init__(self,
                    labels_path = '/home/jovyan/scratch-shared/Ebba/KinaseInhibitorData/dataframe.csv',
                    output_dir = '/home/jovyan/scratch-shared/Ebba/Kinase_Leave_One_Out',
                    index_to_include = [],
                    class_index = 5,
                    well_index = 3,
                    index_to_leave_out = 6,
                    divide_by_index = 5,
                    loop_wells = 2
                    ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.classes_to_include = classes_to_include
        self.class_index =  class_index 
        self.well_index =  well_index
        self.index_to_leave_out = index_to_leave_out
        self.divide_by_index =  divide_by_index
        self.loop_wells = loop_wells

    def main(self):
        print("start get info")
        entries_list = []

        with open(self.labels_path, 'r') as read_obj:
            csv_reader = csv.reader(read_obj, delimiter=";")
            csv_list = list(csv_reader)
            header = csv_list.pop(0) #remove header

            entries_list = self.get_number(csv_list)
            
            if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
                shutil.rmtree(self.output_dir)

            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)
                print("made the output dir")

        ## Write output 
        file_object = open(self.output_dir + "/classes.txt", "w+")
        for entry in entries_list:
            file_object.write("\"" + str(entry)+ "\"" +"\n") 

        print("Finished. Find output in: " + self.output_dir)


    def get_number(self, csv_list):

        division_dict = {}
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
                if class_for_row not in division_dict[divide_by]:
                    division_dict[divide_by][class_for_row] = {}
                if well not in division_dict[divide_by][class_for_row]:
                    division_dict[divide_by][class_for_row][well] = []
                division_dict[divide_by][class_for_row][well].append(entry)

        longest_class = 0
        length_of_classes = {}
        for key in division_dict.keys():
            for class_key in division_dict[key]:
                length_of_classes[class_key] = len(division_dict[key][class_key]
                if longest_class < len(division_dict[key][class_key]):
                    longest_class = len(division_dict[key][class_key])        

        return length_of_classes


    def get_divisions(self, csv_list):

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

        longest_class = 0
        for key in division_dict.keys():
            for class_key in division_dict[key]:
                if longest_class < len(division_dict[key][class_key]):
                    longest_class = len(division_dict[key][class_key])
        
        iteration = 0 
        wells_to_run = [[] for x in range(int(longest_class * self.loop_wells))]
        while iteration < int(longest_class * self.loop_wells):
            for division in division_dict.keys():
                for class_key in division_dict[division].keys():
                    wells_for_class = division_dict[division][class_key]
                    used_wells_for_class = used_wells[division][class_key]
                    if len(wells_for_class) == len(used_wells_for_class):
                        used_wells_for_class = []
                    available_wells = [w for w in division_dict[division][class_key].keys() if w not in used_wells_for_class] 
                    well = random.choice(available_wells)
                    wells_to_run[iteration].append(str(well))
                    used_wells_for_class.append(well)
                    used_wells[division][class_key] = used_wells_for_class
            iteration += 1
        return wells_to_run

if __name__ == "__main__":
    GetInfo().main()
