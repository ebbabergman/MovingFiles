import csv
import numpy as np
from pymysql import NULL

import GroupRows

class MakeKFolds:

        def __init__(self,
                        labels_path = '/home/jovyan/Data/Specs/Specs_Labels.csv',
                        output_dir='/home/jovyan/Inputs/SPECS_QC_Automatic_Jordi_Controll_top3_K_folds/',
                        # include_groups = [], #Empty for everything included,
                        # include_groups = ["heat shock response signalling agonist", "phosphodiesterase inhibitor", "methyltransferase inhibitor","DILI","HDAC inhibitor","topoisomerase inhibitor", "mTOR inhibitor","NFkB pathway inhibitor","JAK inhibitor","pregnane x receptor agonist"], #Empty for everything included,
                        included_classes = ["negcon","DNA polymerase inhibitor", "mTOR inhibitor", "topoisomerase inhibitor"],
                        class_header = "selected_mechanism",
                        exclude_images_path = "/home/jovyan/Outputs/CellProfiler/QC_Specs/QC_Specs_OnlyFlaggedAut_AllPlates.csv",
                        intact_group_header = 'compound_id',
                        meta_data_header = ['plate', 'well', 'site'],
                        image_number_heading = "ImageNr",  
                        k_folds = "3",
                        divide_on_header = 'compound_id',
                        ):
                self.labels_path = labels_path
                self.output_dir = output_dir
                self.exclude_images_path = exclude_images_path
                self.included_classes = included_classes
                self.class_header = class_header
                self.meta_data_header = meta_data_header
                self.image_number_heading = image_number_heading
                self.intact_group_header = intact_group_header
                self.k_folds = int(k_folds)
                self.divide_on_header = divide_on_header

        def main(self):
                print("Started KFolds.")

                #Make numpy
                # df_base = pd.read_csv(self.labels_path , delimiter= ",")
                # df_base.dropna(subset = [self.class_column_header], inplace=True)
                # df_base.drop_duplicates(inplace=True)
              
                all_data = np.genfromtxt(self.labels_path, delimiter=',', names=True, dtype=None)
                self.intact_group_index = all_data.dtype.names.index(self.intact_group_header)

                test = all_data.dtype
                for i in range(0,len(test)):
                        print(i)

# # TODO: new_array = np.array(array, dtype = [("name", object), 
#                                      ("N1", int), 
#                                      ("N2", int),
#                                      ("N3", float)])

                all_data[self.class_header ]=  all_data[self.class_header].astype('str_')
                # Find out how to find all the data
                if len(self.included_classes) == 0:
                        # TODO set list to all unique values
                        print("todo - implement now")
                
                folds = np.empty(shape = (self.k_folds))
                remaining = []
                
                a = np.array(["hello", "world"])
                print(a == "world")
                # TODO Ebba change datatype so that it becomes string? Other way to solve it?
                # Make K-folds by GroupRows
                for group_name in self.included_classes:
                        #rows = np.where(all_data[self.class_header].equals(group_name) )
                        rows = np.where(all_data[self.class_header] == group_name)
                        group_data = all_data[rows]

                        new_k_folds, remaining_rows = GroupRows.GroupRows.group_rows(group_data, self.intact_group_index, self.k_folds)
                        print(remaining_rows)
                        remaining.append(remaining_rows)
                # Pitch out remaining unique values equally?

                # Save K-folds in output

        # self.included_groups = self.get_included_groups(df_base)
        #  df = df_base[df_base[self.include_header].isin(self.included_groups) & ~df_base[self.exclude_header].isin(self.exclude_groups)]
        
if __name__ == "__main__":
    MakeKFolds().main()