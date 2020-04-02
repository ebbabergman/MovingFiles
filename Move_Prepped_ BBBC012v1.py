##Moving of kensert structured files into folders defined by Methods Of Action (MOA)
## Data set first taken from bbbc021

import csv
import os


FILE_NAME = '/Users/ebbbe288/Documents/TestData_bbc021v/bbbc021v1_labels.csv'
# with open(FILE_NAME, newline='') as csvfile:
#     data = list(csv.reader(csvfile))
# print(data)

##Assumes row structure is ['image_number', 'compound', 'concentration', 'moa', 'plate', 'well', 'replicate']
def sort_into_class_folders(row):
     print(str(row))


with open(FILE_NAME, 'r') as read_obj:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(read_obj, delimiter=";")
    # Iterate over each row in the csv using reader object
    for row in csv_reader:
        sort_into_class_folders(row)

