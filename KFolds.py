import csv
import numpy as np

import GroupRows

class MakeKFolds:

        def __init__(self,
                        labels_path = '/home/jovyan/Data/Specs/Labels.csv',
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
                # TODO Drop duplicates and NA or assume this is already done?
              
                all_data = np.genfromtxt(self.labels_path, delimiter=',', names = True, dtype = None, encoding = None)
                self.intact_group_index = all_data.dtype.names.index(self.intact_group_header)

                if len(self.included_classes) == 0:
                        # TODO set list to all unique values
                        print("todo - implement now")
                
                folds = [ [] for _ in range(self.k_folds) ]
                remaining = []
                
                classes =  all_data[self.class_header].astype('str_')   
                for group_name in self.included_classes:
                        rows = np.where(classes == group_name)
                        group_data = all_data[rows]
                        new_k_folds, remaining_rows = GroupRows.GroupRows.group_rows(group_data, self.intact_group_header, self.k_folds)
                        print(remaining_rows)
                        remaining.append(remaining_rows)
                        for fold in range(0,self.k_folds):
                               folds[fold].append(new_k_folds[fold]) 
                # TODO EBBA
                # Pitch out remaining unique values equally

                for fold in range(0,self.k_folds):
                        # TODO EBBA
                        # Get each k_folds unique value and the take *all of the rows* into a new csv file with the right name
                        
                        np.savetxt(self.output_dir + str(fold) + "_fold.csv", folds[fold], delimiter=",")


                print("Hello")

                # Save K-folds in output

if __name__ == "__main__":
    MakeKFolds().main()