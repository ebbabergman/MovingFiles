## DataFlowSort.py
## Moves files into a structure that works with data flow, 

import csv
import os
import shutil
import numpy as np

OUTPUT_DIR = '/home/jovyan/scratch-shared/Ebba/Small_bbbc021'
LABELS_PATH = '/home/jovyan/kensert_CNN/bbbc021_labels.csv'
IMAGE_DIR= '/home/jovyan/kensert_CNN/images_bbbc021'
IMAGE_NAME ='/bbbc021_%s.png' #Where %s is the image number
OUTPUT_SIZE = 1 # Percentage of original total size that should be used
INCLUDED_CLASSES = ['Aurora kinase inhibitors', 'Eg5 inhibitors'] #Empty for all classes included, not recommended, this file only copies in that case

##Assumes row structure is ['image_number', 'compound', 'concentration', 'moa', 'plate', 'well', 'replicate']
def copy_image(row): #Where category is train, validation or test
    if(str(row[3]) == 'moa') :     #Ignore header
        print("reached header!!!!")
        return
    current_path = IMAGE_DIR + IMAGE_NAME  % str(row[0])
   
    dir_path = OUTPUT_DIR + "/Images"
    target_path = dir_path +"/" +str(row[0]) + ".png"

    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print(str(row))
    shutil.copyfile(current_path, target_path)

def get_randomized_set(csv_list):
    # Choose 3 random sets for training, validation and test
    included_rows = []

    if(len(INCLUDED_CLASSES)==0):
        included_rows = csv_list 
    else:
        for entry in csv_list:
            if(entry[3] in INCLUDED_CLASSES):
                included_rows.append(entry)

    data_size = len(included_rows )
    set_size = int(data_size * OUTPUT_SIZE)

    indices = np.arange(data_size)
    np.random.shuffle(indices)

    rows =np.array(included_rows) [indices[:set_size]]
    return rows

print("Starting program")

with open(LABELS_PATH, 'r') as read_obj:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(read_obj, delimiter=";")
    csv_list = list(csv_reader)
    if(str(csv_list[0][3]) == 'moa'):
        csv_list.pop(0) #Remove header
    
    rows = get_randomized_set(csv_list)

    with open(OUTPUT_DIR + "/Labels.csv", 'w', newline = '') as new_labels_file:
     wr = csv.writer(new_labels_file, quoting=csv.QUOTE_ALL)
     wr.writerow(rows)

    for row in rows:
       copy_image(row)


    
print("Finished program")