import csv
import numpy as np

class MakeKFolds:

    def __init__(self,
                labels_path = '/home/jovyan/Data/Specs/Specs_Labels.csv',
                output_dir='/home/jovyan/Inputs/SPECS_QC_Automatic_Jordi_Controll_top3_K_folds/',
        #      include_groups = [], #Empty for everything included,
                #include_groups = ["heat shock response signalling agonist", "phosphodiesterase inhibitor", "methyltransferase inhibitor","DILI","HDAC inhibitor","topoisomerase inhibitor", "mTOR inhibitor","NFkB pathway inhibitor","JAK inhibitor","pregnane x receptor agonist"], #Empty for everything included,
                include_groups = ["negcon","DNA polymerase inhibitor", "mTOR inhibitor", "topoisomerase inhibitor"],
                include_header = "selected_mechanism",
                class_column_header = "selected_mechanism",
                exclude_groups = [],
                exclude_header = "selected_mechanism",
                exclude_images_path = "/home/jovyan/Outputs/CellProfiler/QC_Specs/QC_Specs_OnlyFlaggedAut_AllPlates.csv",
                intact_group_header = 'compound_id',
                intact_control_group_headers = ['plate', 'well'], # NOTE: hard coded for 2 headers to to troubles with dataframe
                meta_data_header = ['plate', 'well', 'site'],
                image_number_heading = "nr",   
                has_controls = False,
                frac_of_controls_to_use = 1,
                k_folds = "3",
                divide_on_header = 'compound_id',
                ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.exclude_images_path = exclude_images_path
        self.included_groups = include_groups
        self.include_header = include_header
        self.exclude_groups = exclude_groups
        self.exclude_header = exclude_header
        self.meta_data_header = meta_data_header
        self.image_number_heading = image_number_heading
        self.class_column_header =  class_column_header 
        self.intact_group_header = intact_group_header
        self.intact_control_group_headers = intact_control_group_headers
        self.k_folds = int(k_folds)
        self.divide_on_header = divide_on_header