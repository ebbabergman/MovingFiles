import pandas as pd
import numpy as np
import random
from make_k_folds_tvt import MakeTVTSets
from make_k_folds_tvt_specs import MakeKFoldsTVTSPECS
import General_Moving
# Make K-folds, including their train and validation parts, with csv files as outpus


class MakeKFoldsTVT_Tox21:

    def __init__(self,
                 labels_path="/scratch-shared/martin/003_SPECS1K_ML/001_data/Specs935_ImageMeans_AfterQC_AnnotatedWithMOA.csv",
                 output_dir='/home/jovyan/Inputs/tox_21_3_fold_tvt/',
                 include_groups=["heat shock response signalling agonist",
                                 "pregnane x receptor agonist",
                                 "estrogen receptor alpha agonist",
                                 "aryl hydrocarbon receptor agonist"],
                 include_header="selected_mechanism",
                 class_column_header="selected_mechanism",
                 excluded_groups=[],
                 excluded_groups_headers=[],
                 # exclude_images_path="/home/jovyan/Data/Specs/Flaggs/16Bit_images_nuclei_cut_above149.csv",
                 exclude_images_path="",
                 intact_group_header='Compound ID',
                 unique_sample_headers=['Plate', 'Well', 'Site'],
                 image_number_heading="ImageID",
                 k_folds="5",
                 divide_on_header='Compound ID',
                 # 1 = 100%,  Percentage of images remaining afte the test set has been excluded
                 valid_fraction=0.25,
                 # How many times more persumed negatives should we have?
                 proportion_persumed_negative=4,
                 leave_one_out=False,
                 make_unique_validation=True,
                 meta_data_headers=['ImageID', 'PlateID', 'Well', 'Site', 'Plate', 'Plate_Well', 'batch_id',
                                    'pertType', 'cmpd_conc', 'Flag', 'Count_nuclei', 'Batch nr', 'Compound ID', 'selected_mechanism']
                 ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.include_groups = include_groups
        self.include_header = include_header
        self.class_column_header = class_column_header
        self.excluded_groups = excluded_groups
        self.excluded_groups_headers = excluded_groups_headers
        self.exclude_images_path = exclude_images_path
        self.intact_group_header = intact_group_header
        self.unique_sample_headers = unique_sample_headers
        self.image_number_heading = image_number_heading
        self.k_folds = k_folds
        self.divide_on_header = divide_on_header
        self.valid_fraction = valid_fraction
        self.leave_one_out = leave_one_out
        self.make_unique_validation = make_unique_validation
        self.proportion_persumed_negative = proportion_persumed_negative
        self.meta_data_headers = meta_data_headers

    def main(self):
        print("Started making divisions for runs for tox21. Each included group will get their own k-fold.")
        General_Moving.make_non_existing_path(self.output_dir)

        self.make_lables()

        for predicted_group in self.include_groups:
            print("Making K-folds for " + str(predicted_group))
            group_name = predicted_group.replace(" ", "_")
            group_output_path = self.output_dir + "/" + group_name + "/"
            group_labels_path = self.output_dir + group_name + "_Labels.csv"

            makeKFoldsTVTSPECS = MakeKFoldsTVTSPECS(labels_path=group_labels_path,
                                                    output_dir=group_output_path,
                                                    include_groups=[],
                                                    include_header=group_name,
                                                    class_column_header=group_name,
                                                    excluded_groups=[],
                                                    excluded_groups_headers=[],
                                                    exclude_images_path="",
                                                    intact_group_header=self.intact_group_header,
                                                    unique_sample_headers=self.unique_sample_headers,
                                                    image_number_heading=self.image_number_heading,
                                                    k_folds=self.k_folds,
                                                    divide_on_header=self.divide_on_header,
                                                    valid_fraction=self.valid_fraction,
                                                    leave_one_out=self.leave_one_out,
                                                    make_unique_validation=self.make_unique_validation)

            makeKFoldsTVTSPECS.main()
        print("Finished making tox21 dataset. Find output in: " + self.output_dir)

    def make_lables(self):
        df_original = pd.read_csv(self.labels_path, sep=";")
        all_included_mask = df_original[self.include_header].isin(
            self.include_groups)
        df_available_neg = df_original[~all_included_mask]
        available_neg_intact_group = list(
            df_available_neg[self.intact_group_header].unique())

        for predicted_group in self.include_groups:
            self.make_labels(predicted_group, df_original,
                             available_neg_intact_group)

    def make_group_labels(self, predicted_group, df_original, available_neg_intact_group):
        print("Making labels for " + str(predicted_group))
        group_name = predicted_group.replace(" ", "_")
        group_mask = df_original[self.include_header] == predicted_group
        df_group = df_original[group_mask]
        size_of_true_group = len(
            df_group[self.intact_group_header].unique())  # number of rows

        negative_intact_group = random.sample(
            available_neg_intact_group, size_of_true_group * self.proportion_persumed_negative)
        negative_mask = df_original[self.intact_group_header].isin(
            negative_intact_group)

        df_new_group = df_original[self.meta_data_headers].copy()

        conditions = [group_mask, negative_mask]
        choices = [True, False]
        df_new_group[group_name] = np.select(
            conditions, choices, default=np.nan)

        df_new_group = df_new_group[df_new_group[group_name].notnull()]
        df_new_group.drop_duplicates(keep=False, inplace=True)

        group_labels_path = self.output_dir + group_name + "_Labels.csv"
        df_new_group.to_csv(group_labels_path,  index=False)

        print("Finished making labels for " + str(predicted_group))


if __name__ == "__main__":
    MakeKFoldsTVT_Tox21().main()
