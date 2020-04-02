##Moving of kensert structured files into folders defined by Methods Of Action (MOA)
## Data set first taken from bbbc021

import csv


FILE_NAME = '/Users/ebbbe288/Documents/TestData_bbc021v/bbbc021v1_labels.csv'
# with open(FILE_NAME, newline='') as csvfile:
#     data = list(csv.reader(csvfile))
# print(data)


with open(FILE_NAME, 'r') as read_obj:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(read_obj, delimiter=";")
    # Iterate over each row in the csv using reader object
    for row in csv_reader:
        # row variable is a list that represents a row in csv
        print(row[0])