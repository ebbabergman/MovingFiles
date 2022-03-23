import csv
import numpy as np


class GroupRows:

    staticmethod
    def group_rows(data, keep_complete_header, number_of_groups):
        groups = [None]*number_of_groups

        np_unique, np_index_unique, counts = np.unique(data[keep_complete_header], return_inverse = True, return_counts = True)

        min_per_fold = len(np_unique)//number_of_groups

        for i in number_of_groups:
            in_group = np_unique(size =(1,min_per_fold) ,replace = False)
            groups[i] = in_group
        
        #Should I return remaining np_unique or just divide them now?


        return groups, np_unique # return "left over" groups?
