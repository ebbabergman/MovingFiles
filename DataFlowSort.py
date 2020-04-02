## DataFlowSort.py
## Moves files into a structure that works with data flow, 

import csv
import os
import shutil
import numpy as np

DIR = '/Users/ebbbe288/Documents/TestData/'
FILE_NAME = '/Users/ebbbe288/Documents/TestData_bbc021v/bbbc021v1_labels.csv'
IMAGE_DIR= '/Users/ebbbe288/Documents/TestData_bbc021v/images_bbbc021'
IMAGE_NAME ='/bbbc021v1_%s.png' #Where %s is the image number

##Assumes row structure is ['image_number', 'compound', 'concentration', 'moa', 'plate', 'well', 'replicate']
def sort_into_class_folders(row, category): #Where category is train, validation or test
    if(str(row[3]) == 'moa') :     #Ignore header
        print("reached header!!!!")
        return
    current_path = IMAGE_DIR + IMAGE_NAME  % str(row[0])
   
    dir_path = DIR + category + str(row[3]) 
    target_path = dir_path +"/" +str(row[0]) + ".png"

    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print(str(row))
    shutil.copyfile(current_path, target_path)

def sort_into_test_folder(row, category): #Where category is train, validation or test
    if(str(row[3]) == 'moa') :     #Ignore header
        print("reached header!!!!")
        return
    current_path = IMAGE_DIR + IMAGE_NAME  % str(row[0])
   
    dir_path = DIR + category 
    target_path = dir_path +"/" +str(row[0]) + ".png"

    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print(str(row))
    shutil.copyfile(current_path, target_path)

def get_randomized_set_rows(read_obj):
    # Choose 3 random sets for training, validation and test
    data_size = sum(1 for row in read_obj)-1 # ignore header
    validation_set_size = int(data_size * 0.15 +1)
    training_set_size = int(data_size - 2* validation_set_size)

    print(data_size)
    indices = np.arange(data_size)
    print(indices[5])
    np.random.shuffle(indices)
    print("Random number at positoin 5: " + str(indices[5]))

    train_rows= indices[:training_set_size]
    validation_rows=indices[training_set_size:training_set_size + validation_set_size]
    test_rows = indices[training_set_size + validation_set_size:]
    return train_rows, validation_rows,test_rows

with open(FILE_NAME, 'r') as read_obj:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(read_obj, delimiter=";")
    train_row_numbers, validation_row_numbers,test_row_numbers = get_randomized_set_rows(read_obj)
    
    train_rows=[row for idx, row in enumerate(csv_reader) if idx in train_row_numbers]
    validation_rows=[row for idx, row in enumerate(csv_reader) if idx in validation_row_numbers]
    test_rows=[row for idx, row in enumerate(csv_reader) if idx in test_row_numbers]

    for row in train_rows:
        sort_into_class_folders(row, "Train")
    for row in validation_rows:
        sort_into_class_folders(row, "Validation")
    for row in test_rows:
        sort_into_test_folder(row, "Test")

