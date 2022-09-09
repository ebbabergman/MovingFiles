import csv
import os
from pickle import FALSE
import shutil
from xxlimited import Str
from matplotlib import use
import numpy as np
import random
import pandas as pd
from pandas.core import indexing
import math

# Make K-folds, including their train and validation parts, with csv files as outpus


class MakeKFolds:

    def __init__(self,
                 labels_path="/home/jovyan/Data/Specs/Specs_Labels_First_MiPhHo.csv",
                 output_dir='/home/jovyan/Inputs/SPECS_Nuclei_Cutoff_no_sorbitol_no_test_negcon_plate_small_images_All_K_folds/',
                 #
                 include_groups = [], #Empty for everything included,
                 # top10 good images (nuclei cut off)
                 #include_groups=["heat shock response signalling agonist","phosphodiesterase inhibitor","DILI","methyltransferase inhibitor","estrogen receptor alpha modulator","cyclooxygenase inhibitor","pregnane x receptor agonist","PPAR receptor agonist","protein synthesis inhibitor","CC chemokine receptor antagonist"],
                 # include_groups = ["negcon", "heat shock response signalling agonist", "DILI", "estrogen receptor alpha modulator", "phosphodiesterase inhibitor",  "cyclooxygenase inhibitor"], #Empty for everything included,
                 # include_groups = ["negcon","heat shock response signalling agonist", "phosphodiesterase inhibitor", "methyltransferase inhibitor","DILI","HDAC inhibitor","topoisomerase inhibitor", "mTOR inhibitor","NFkB pathway inhibitor","JAK inhibitor","pregnane x receptor agonist"],                    #include_groups = ["negcon","DNA polymerase inhibitor", "mTOR inhibitor", "topoisomerase inhibitor"],
                 #include_groups = ['pregnane x receptor agonist', 'DILI', 'estrogen receptor alpha modulator', 'tubulin polymerization inhibitor', 'topoisomerase inhibitor', 'heat shock response signalling agonist', 'methyltransferase inhibitor', 'aryl hydrocarbon receptor agonist', 'estrogen receptor alpha agonist', 'mitochondrial toxicity  agonist', 'retinoid receptor agonist', 'protein synthesis inhibitor', 'phosphodiesterase inhibitor', 'DNA polymerase inhibitor', 'mTOR inhibitor', 'PPAR receptor agonist', 'glucocorticoid receptor agonist', 'ATPase inhibitor', 'cyclooxygenase inhibitor', 'NFkB pathway inhibitor', 'angiotensin converting enzyme inhibitor', 'adenosine receptor antagonist', 'PARP inhibitor', 'JAK inhibitor', 'HSP inhibitor', 'HDAC inhibitor',  'CC chemokine receptor antagonist', 'Aurora kinase inhibitor','negcon'],
                 include_header="selected_mechanism",
                 class_column_header="selected_mechanism",
                 excluded_groups=[["poscon", "empty"], ["P015085"], ['DNA synthesis inhibitor']],
                 excluded_groups_headers=["pert_type", "plate", "selected_mechanism"],
                 exclude_images_path="/home/jovyan/Data/Specs/Flaggs/old_MiPhHo_Labels_images_outside_nuclei_cut_82_149_no_sorbitol.csv",
                 intact_group_header='compound_id',
                 unique_sample_headers=['plate', 'well', 'site'],
                 image_number_heading="ImageNr",
                 k_folds="3",
                 divide_on_header='compound_id',
                 make_train_valid=True,
                 # 1 = 100%,  Percentage of images remaining afte the test set has been excluded
                 valid_fraction=0.25
                 ):
        self.labels_path = labels_path
        self.output_dir = output_dir
        self.exclude_images_path = exclude_images_path
        self.included_groups = include_groups
        self.include_header = include_header
        self.excluded_groups = excluded_groups
        self.excluded_groups_headers = excluded_groups_headers
        self.unique_sample_headers = unique_sample_headers
        self.image_number_heading = image_number_heading
        self.class_column_header = class_column_header
        self.intact_group_header = intact_group_header
        self.k_folds = int(k_folds)
        self.divide_on_header = divide_on_header
        self.make_train_valid = make_train_valid
        self.valid_fraction = valid_fraction

    def main(self):
        print("Started make k-folds.")

        self.make_path_available()

        df_base = pd.read_csv(self.labels_path, delimiter=",")
        df_base.dropna(subset=[self.class_column_header], inplace=True)
        df_base.drop_duplicates(inplace=True)

        df = self.exclude_groups(df_base)

        print("Included and exlcuded groups")
        df = self.exclude_images(df)
        print("excluded images indicated")

        df = self.include_groups(df)

        k_folds_test = self.get_k_folds_test(df)
        df_test_statistics = self.get_statistics(k_folds_test, df)
        df_test_statistics.to_csv(
            self.output_dir + "k_fold_test_statistics.csv", index=False)
        print("Test Statistics")
        print(df_test_statistics.to_latex())

        if self.make_train_valid:
            k_folds_train, k_folds_validation = self.get_k_folds_tv(
                df, k_folds_test)
            df_validation_statistics = self.get_statistics(
                k_folds_validation, df)
            df_train_statistics = self.get_statistics(k_folds_train, df)

            df_validation_statistics.to_csv(
                self.output_dir + "k_fold_validation_statistics.csv", index=False)
            df_train_statistics.to_csv(
                self.output_dir + "k_fold_train_statistics.csv", index=False)

            print("Validation Statistics ")
            print(df_validation_statistics.to_latex())
            print("Train Statistics ")
            print(df_train_statistics.to_latex())

        for k_fold in range(0, self.k_folds):
            df_test_fold = k_folds_test[k_fold]
            df_test_fold.to_csv(
                self.output_dir + "k_fold_test_" + str(k_fold + 1)+".csv", index=False)

            df_validation_fold = k_folds_validation[k_fold]
            df_validation_fold.to_csv(
                self.output_dir + "k_fold_validation_" + str(k_fold + 1)+".csv", index=False)

            df_train_fold = k_folds_train[k_fold]
            df_train_fold.to_csv(
                self.output_dir + "k_fold_train_" + str(k_fold + 1)+".csv", index=False)

        print("Finished. Find output in: " + self.output_dir)

    def exclude_groups(self,df):
        for index in range(0, len(self.excluded_groups_headers)):
            exclude_groups_header = self.excluded_groups_headers[index]
            df = df[~df[exclude_groups_header].isin(
                self.excluded_groups[index])]
        return df

    def exclude_images(self, df):
        df_bad_images = pd.read_csv(self.exclude_images_path, delimiter=",")
        df_bad_images = df_bad_images.drop_duplicates(
            subset=self.unique_sample_headers)

        if 'Total' in df_bad_images.columns:
            bad_image_mask = df_bad_images["Total"] == 1
            df_bad_images = df_bad_images[bad_image_mask]
        elif 'total' in df_bad_images.columns:
            bad_image_mask = df_bad_images["total"] == 1
            df_bad_images = df_bad_images[bad_image_mask]

        df_merged = pd.merge(df, df_bad_images[self.unique_sample_headers],
                             on=self.unique_sample_headers, how="outer", indicator=True)

        merge_both_mask = df_merged["_merge"] == "both"
        merge_left_mask = df_merged["_merge"] == "left_only"
        merge_right_mask = df_merged["_merge"] == "right_only"

        print("Excluding images")
        print(str(df.shape[0]) + " included images to start with")
        print(str(df_merged[merge_right_mask].shape[0]) +
              " images were not in the input labels")
        print(str(df_merged[merge_both_mask].shape[0]) +
              " images will be dropped")
        print(str(df_merged[merge_left_mask].shape[0]) +
              " images will be kept")

        df_new = df[df[self.image_number_heading].isin(
            df_merged[merge_left_mask][self.image_number_heading])]
        df = df_new.copy()

        return df

    def include_groups(self,df):
        included_groups = self.included_groups
        use_all_available_groups = len(included_groups) == 0
        available_groups = df[self.include_header].unique()

        if use_all_available_groups:
            print("Using all available groups")
            included_groups = available_groups
            self.included_groups = available_groups
        elif not all(elem in included_groups  for elem in available_groups):
            not_available = set(included_groups) - set(available_groups)
            raise GroupNotIncluded(not_available)

        print("Available groups: " + available_groups)
        df = df[df[self.include_header].isin(
            included_groups)]
        
        return df

    def get_k_folds_test(self, df):
        number_of_folds = self.k_folds
        k_folds = [None]*number_of_folds
        group_n = {}
        df_unused = df.copy()

        for group in self.included_groups:
            df_group = df[df[self.include_header].isin([group])]
            group_n[group] = math.floor(
                df_group[self.divide_on_header].nunique()/number_of_folds)
            if group_n[group] < 1:
                group_n[group] = 1
                print(
                    "A group did not have enough unique groupings to have 1 unique entry per fold. Group:" + str(group))

        for k_fold in range(1, number_of_folds + 1):
            df_fold = pd.DataFrame()
            for group in self.included_groups:
                df_group = df_unused[df_unused[self.include_header].isin([
                                                                         group])]
                if df_group.empty:
                    df_used = df[df[self.include_header].isin([group])]
                    df_unused = df_used.copy()
                    df_group = df_unused[df_unused[self.include_header].isin([
                                                                             group])]
                    print(
                        "Every unit from group had been used, re-using values for group " + str(group) + ".")
                df_group_coice = self.get_group_selection(
                    df_group, group_n[group])
                df_fold = pd.concat(
                    [df_fold, df_group_coice], ignore_index=True)
            df_unused = pd.concat(
                [df_unused, df_fold, df_fold]).drop_duplicates(keep=False)
            k_folds[k_fold-1] = df_fold

        # deal with the unused unique groups
        group_index = 0
        for group in self.included_groups:
            df_group = df_unused[df_unused[self.include_header].isin([group])]

            if df_group.empty:
                continue

            k_fold_to_choose = [*range(1, number_of_folds + 1)]
            while not df_group.empty:
                if len(k_fold_to_choose) == 0:
                    k_fold_to_choose = [*range(1, number_of_folds + 1)]
                    # TODO make this a warning
                    print(
                        "Something is wrong with the way the number of groups are positioned, could put remainder into each fold")

                chosen_k_fold = np.random.choice(
                    k_fold_to_choose, size=1, replace=False)[0]
                k_fold_to_choose.remove(chosen_k_fold)
                df_group_coice = self.get_group_selection(df_group, 1)

                k_folds[chosen_k_fold-1] = pd.concat(
                    [k_folds[chosen_k_fold-1], df_group_coice], ignore_index=True)
                df_unused = pd.concat(
                    [df_unused, df_group_coice, df_group_coice]).drop_duplicates(keep=False)
                df_group = pd.concat(
                    [df_group, df_group_coice, df_group_coice]).drop_duplicates(keep=False)

            group_index = group_index + 1

        if not df_unused.empty:
            print(
                "WARNING: Didn't put all avialiable compound into test k-folds! Left over:")
            print(df_unused)

        print("Made test sets for k-folds")
        return k_folds

    def get_k_folds_tv(self, df, k_fold_test):
        number_of_folds = self.k_folds
        #validation_fraction = (number_of_folds-1)/number_of_folds*self.valid_fraction
        k_fold_train = [None]*number_of_folds
        k_fold_validation = [None]*number_of_folds
        group_n = {}

        for group in self.included_groups:
            df_group = df[df[self.include_header].isin([group])]
            group_n[group] = math.floor(df_group[self.divide_on_header].nunique(
            )*(number_of_folds-1)/number_of_folds*self.valid_fraction)
            if group_n[group] < 1:
                group_n[group] = 1
                print("A group did not have enough unique groupings to have 1 unique entry per vakudation fold. Using 1 unique entry anyway. Group:" + str(group))

        for k_fold in range(1, number_of_folds + 1):
            df_unused = df.copy()
            df_k_fold_test = k_fold_test[k_fold-1]
            df_unused = pd.concat(
                [df_unused, df_k_fold_test, df_k_fold_test]).drop_duplicates(keep=False)

            df_fold_validation = pd.DataFrame()
            df_fold_train = pd.DataFrame()
            for group in self.included_groups:
                df_group = df_unused[df_unused[self.include_header].isin([
                                                                         group])]

                if k_fold == 1:
                    df_group_coice_validation = self.get_group_selection(
                        df_group, group_n[group])
                else:
                    df_unavailable_validation = pd.concat(
                        k_fold_validation, ignore_index=True)
                    df_group_unavailable_validation = df_unavailable_validation[df_unavailable_validation[self.include_header].isin([
                                                                                                                                    group])]
                    df_group_available_validation = pd.concat(
                        [df_group, df_group_unavailable_validation, df_group_unavailable_validation], ignore_index=True).drop_duplicates(keep=False)
                    unique_entries = df_group_available_validation[self.intact_group_header].unique(
                    )

                    if len(unique_entries) < group_n[group]:
                        number_of_added = group_n[group] - len(unique_entries)
                        add_to_validation = self.get_group_selection(
                            df_group_unavailable_validation, number_of_added)
                        df_group_coice_validation = pd.concat(
                            [df_group_available_validation, add_to_validation], ignore_index=True)
                        print("Reusing previous validation compounds for validation")
                    else:
                        df_group_coice_validation = self.get_group_selection(
                            df_group_available_validation, group_n[group])

                df_fold_validation = pd.concat(
                    [df_fold_validation, df_group_coice_validation], ignore_index=True)
            df_fold_train = pd.concat(
                [df_unused, df_fold_validation, df_fold_validation], ignore_index=True).drop_duplicates(keep=False)
            df_unused = pd.concat([df_unused, df_fold_validation, df_fold_train],
                                  ignore_index=True).drop_duplicates(keep=False)
            k_fold_validation[k_fold-1] = df_fold_validation
            k_fold_train[k_fold-1] = df_fold_train
            if not df_unused.empty:
                print(
                    "WARNING: Didn't put all avialiable compound into train or valid k-folds! Left over:")
                print(df_unused)

        print("Made train and valid sets for k-folds")

        return k_fold_train, k_fold_validation

    def get_group_selection(self, df_group, number_of_unique_entries):
        unique_entries = df_group[self.intact_group_header].unique()
        group_choice = np.random.choice(
            unique_entries, size=number_of_unique_entries, replace=False)
        df_group_choice = df_group[df_group[self.divide_on_header].isin(
            group_choice)]
        return df_group_choice

    def get_statistics(self, k_folds, df):
        # Make some statistics  # TODO make this a function and run for each type of data
        s_statistics = df.groupby(self.class_column_header)[
            self.divide_on_header].nunique()
        df_statistics = pd.DataFrame(s_statistics.index)
        df_statistics["Total"] = s_statistics.values

        number_of_folds = self.k_folds
        for k_fold in range(0, number_of_folds):
            df_fold = k_folds[k_fold]
            df_grouped = df_fold.groupby(self.class_column_header)
            statistic_column_header = self.divide_on_header + \
                "s_in_fold_"+str(k_fold)
            df_statistics[statistic_column_header] = df_grouped[self.divide_on_header].nunique(
            ).values

        return df_statistics

    def make_path_available(self):
        if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print("made the output dir")


class GroupNotIncluded(Exception):
    """Exception raised for errors in the group inclustion.

    Attributes:
        groups -- groups not available after exclusion
        message -- explanation of the error
    """

    def __init__(self, groups, message="A group could not be including after exclusion criteria was met."):
        self.groups = groups
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.groups} -> {self.message}'

if __name__ == "__main__":
    MakeKFolds().main()


