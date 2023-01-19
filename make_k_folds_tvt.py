import csv
from heapq import merge
import os
from pickle import FALSE
import shutil
from xxlimited import Str
import numpy as np
import random
import pandas as pd
from pandas.core import indexing
import math

# Make K-folds, including their train and validation parts, with csv files as outpus

class MakeTVTSets:

    def __init__(self,
                 labels_path="/home/jovyan/Data/BBBC021/BBBC021_Labels.csv",
                 output_dir='/home/jovyan/Inputs/test/',
                 include_groups=[],  # Empty for everything included,
                 include_header="moa",
                 class_column_header="moa",
                 excluded_groups=[
                     ["DMSO", "Cholesterol-lowering", "Eg5 inhibitors"]],
                 excluded_groups_headers=["moa"],
                 exclude_images_path="",
                 intact_group_header='compound',
                 unique_sample_headers=["ImageNumber"],
                 image_number_heading="ImageNumber",
                 k_folds="3",
                 divide_on_header='compound',
                 # 1 = 100%,  Percentage of images remaining afte the test set has been excluded
                 valid_fraction=0.25,
                 non_unique_divider=['concentration',
                                     'moa', 'compound', 'Replicate'],
                 make_unique_validation=True,
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
        self.valid_fraction = valid_fraction
        self.non_unique_divider = non_unique_divider
        self.make_unique_validation = make_unique_validation

    def make_k_folds(self):
        print("Started make k-folds.")

        self.make_path_available()

        df_base = self.get_base()
        print("Read in base for dataframe, starting inclusion and exclusion of rows")

        df = self.include_exclude_rows(df_base)
        print("Included and exlcuded groups")

        k_folds_test = self.get_k_folds_test(df)
        df_test_statistics = self.get_statistics(k_folds_test, df)
        df_test_statistics.to_csv(
            self.output_dir + "k_fold_test_statistics.csv", index=False)
        print("Test Statistics")
        print(df_test_statistics.to_latex())

        k_folds_train, k_folds_validation = self.get_k_folds_tv(
            df, k_folds_test)

        self.make_train_valid_statistics(k_folds_train, k_folds_validation, df)

        self.save_k_folds(k_folds_test, k_folds_validation, k_folds_train)
        print("Finished. Find output in: " + self.output_dir)

    def make_k_folds_train_test(self):
        print("Started make k-folds with only train and test sets, no validation.")

        self.make_path_available()

        df_base = self.get_base()
        print("Read in base for dataframe, starting inclusion and exclusion of rows")

        df = self.include_exclude_rows(df_base)
        print("Included and exlcuded groups")

        k_folds_test = self.get_k_folds_test(df)

        df_test_statistics = self.get_statistics(k_folds_test, df)
        df_test_statistics.to_csv(
            self.output_dir + "k_fold_test_statistics.csv", index=False)
        print("Test Statistics")
        print(df_test_statistics.to_latex())

        k_folds_train, k_folds_validation = self.get_k_folds_train_only(
            df, k_folds_test)

        self.make_train_valid_statistics(k_folds_train, k_folds_validation, df)

        self.save_k_folds(k_folds_test, k_folds_validation, k_folds_train)
        print("Finished. Find output in: " + self.output_dir)

    def make_leave_one_out(self):
        print("Started make leave one out.")

        self.make_path_available()

        df_base = self.get_base()
        print("Read in base for dataframe, starting inclusion and exclusion of rows")

        df = self.include_exclude_rows(df_base)
        print("Included and exlcuded groups")

        self.check_leave_one_out_validity(df)

        k_folds_test = self.get_leave_one_out_test(df)

        df_test_statistics = self.get_statistics_leave_one_out_test(
            k_folds_test, df)
        df_test_statistics.to_csv(
            self.output_dir + "k_fold_test_statistics.csv", index=False)
        print("Test Statistics")
        print(df_test_statistics.to_latex())

        k_folds_train, k_folds_validation = self.get_k_folds_tv(
            df, k_folds_test)

        self.check_no_overlap_between_tvt(
            self, k_folds_test, k_folds_validation, k_folds_train)

        self.make_train_valid_statistics(k_folds_train, k_folds_validation, df)

        self.save_k_folds(k_folds_test, k_folds_validation, k_folds_train)
        print("Finished. Find output in: " + self.output_dir)

    def check_leave_one_out_validity(self, df):
        df_sorted = df.sort_values(self.class_column_header, ascending=True)
        unique_combinations_per_group = df_sorted.groupby(self.class_column_header)[
            self.divide_on_header].nunique()

        unique_combinations_needed = 3  # train, validation, test
        if not self.make_unique_validation:
            unique_combinations_needed = 2  # train, test

        too_few_combinations_per_group = unique_combinations_per_group < unique_combinations_needed

        if too_few_combinations_per_group.any():
            raise Exception(
                "At least one grouping does not have enough unique combinations to proceed. The following do not hav enough unique combinations per group: " +  str(too_few_combinations_per_group
))

    def check_no_overlap_between_tvt(self, k_folds_test, k_folds_validation, k_folds_train):
        for fold in range(0, self.k_folds):
            df_train = k_folds_train[fold]
            df_test = k_folds_test [fold]

            overlaps, df_overlap =  self.get_merge_overlap(df_train, df_test)
            if overlaps:
                raise Exception(
                "Overlap between Test and Train set in fold:" + str(fold)+". Overlapping rows: " + print(df_overlap))
      
            if self.make_unique_validation:
                df_validation = k_folds_validation[fold]
                overlaps, df_overlap =  self.get_merge_overlap(df_validation, df_test)
                if overlaps:
                    raise Exception(
                    "Overlap between Test and Validation set in fold:" + str(fold)+". Overlapping rows: " + print(df_overlap))
        
                overlaps, df_overlap =  self.get_merge_overlap(df_validation, df_train)
                if overlaps:
                    raise Exception(
                    "Overlap between Train and Validation set in fold:" + str(fold)+". Overlapping rows: " + print(df_overlap))
        
    def get_merge_overlap(self, df1, df2, merge_on = []):
        if len(merge_on) == 0:
            merge_on = self.divide_on_header

        df3 = df2.merge(df1, on=merge_on, how='outer',indicator='present_in_both')
        
        df3['present_in_both'] = df3['present_in_both'].eq('both')
        if df3['present_in_both'].any():
            duplicates_unique_samples = df3[df3['present_in_both'] == True][self.unique_sample_headers]
            df1_duplicate_mask = df1[self.unique_sample_headers].is_in(duplicates_unique_samples)
            return True, df1[df1_duplicate_mask]
        
        return False, []
     

    def make_leave_one_out_train_test(self):
        print("Started make leave one out train only, no validation.")

        self.make_path_available()

        df_base = self.get_base()
        print("Read in base for dataframe, starting inclusion and exclusion of rows")

        df = self.include_exclude_rows(df_base)
        print("Included and exlcuded groups")

        k_folds_test = self.get_leave_one_out_test(df)
        df_test_statistics = self.get_statistics_leave_one_out_test(
            k_folds_test, df)
        df_test_statistics.to_csv(
            self.output_dir + "k_fold_test_statistics.csv", index=False)
        print("Test Statistics")
        print(df_test_statistics.to_latex())

        k_folds_train, k_folds_validation = self.get_k_folds_train_only(
            df, k_folds_test)

        self.make_train_valid_statistics(k_folds_train, k_folds_validation, df)

        self.save_k_folds(k_folds_test, k_folds_validation, k_folds_train)
        print("Finished. Find output in: " + self.output_dir)

    def get_base(self):
        df_base = pd.read_csv(self.labels_path, delimiter=",")

        # Drop any lingering index columns, i.e. columns containing the string "Unnamed"
        columns_to_drop = [s for s in df_base.columns if "Unnamed" in s]
        drop_filter = df_base.filter(columns_to_drop)
        df_base.drop(drop_filter, inplace=True, axis=1)

        df_base.dropna(subset=[self.class_column_header], inplace=True)
        df_base.drop_duplicates(
            subset=self.unique_sample_headers, inplace=True, ignore_index=True)
        return df_base

    def include_exclude_rows(self, df):
        all_groups_that_could_be_included = self.get_included_groups(df)
        df = df[df[self.include_header].isin(
            all_groups_that_could_be_included)]
        df = self.exclude_groups(df)
        if len(self.exclude_images_path) > 0:
            df = self.exclude_images(df)
            print("Excluded images indicated with file")
        self.included_groups = self.get_included_groups(df)
        return df

    def get_included_groups(self, df):
        included_groups = self.included_groups
        if (len(included_groups) == 0):
            included_groups = df[self.include_header].unique()

        included_groups = [
            group for group in included_groups if group not in self.excluded_groups]
        return included_groups

    def exclude_groups(self, df):
        for index in range(0, len(self.excluded_groups_headers)):
            excluded_groups_header = self.excluded_groups_headers[index]
            df = df[~df[excluded_groups_header].isin(
                self.excluded_groups[index])]
        return df

    def exclude_images(self, df):
        if self.exclude_images_path == "":
            return df
        df_bad_images = pd.read_csv(self.exclude_images_path, delimiter=",")
        df_bad_images = df_bad_images.drop_duplicates(
            subset=self.unique_sample_headers, ignore_index=True)

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

    def include_groups(self, df):
        included_groups = self.included_groups
        use_all_available_groups = len(included_groups) == 0
        available_groups = df[self.include_header].unique()

        if use_all_available_groups:
            print("Using all available groups")
            included_groups = available_groups
            self.included_groups = available_groups
        elif not set(included_groups).issubset(available_groups):
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
                [df_unused, df_fold, df_fold]).drop_duplicates(subset=self.unique_sample_headers, keep=False, ignore_index=True)
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
                    [df_unused, df_group_coice, df_group_coice]).drop_duplicates(subset=self.intact_group_header, keep=False, ignore_index=True)
                df_group = pd.concat(
                    [df_group, df_group_coice, df_group_coice]).drop_duplicates(subset=self.intact_group_header, keep=False, ignore_index=True)
            group_index = group_index + 1

        if not df_unused.empty:
            print(
                "WARNING: Didn't put all avialiable compound into test k-folds! Left over:")
            print(df_unused)

        print("Made test sets for k-folds")
        return k_folds

    def get_leave_one_out_test(self, df):
        number_of_folds = df[self.divide_on_header].nunique()
        self.k_folds = number_of_folds
        print("Leave one out will result in " +
              str(number_of_folds) + " folds.")
        k_folds = [None]*number_of_folds

        leave_out_list = df[self.divide_on_header].unique()

        k_fold = 1
        for entry in leave_out_list:
            df_fold = df[df[self.divide_on_header] == entry]
            k_folds[k_fold-1] = df_fold
            k_fold = k_fold + 1

        print("Made test sets for leave one out. Made " +
              str(number_of_folds) + "folds")
        return k_folds

    def get_k_folds_tv(self, df, k_fold_test):
        number_of_folds = self.k_folds
        validation_fraction = self.valid_fraction * \
            (number_of_folds-1)/number_of_folds
        k_fold_train = [None]*number_of_folds
        k_fold_validation = [None]*number_of_folds
        group_n = {}
        df_available_validation = df.copy()

        for group in self.included_groups:
            df_group = df[df[self.include_header].isin([group])]
            group_n[group] = math.floor(
                df_group[self.divide_on_header].nunique()*validation_fraction)
            if group_n[group] < 1:
                group_n[group] = 1
                print("A group did not have enough unique groupings to have 1 unique entry per validation fold. Using 1 unique entry anyway. Group:" + str(group))

        for k_fold in range(1, number_of_folds + 1):
            print("Starting k-fold: " + str(k_fold))
            df_unused = df.copy()
            df_fold_test = k_fold_test[k_fold-1]
            df_unused = pd.concat(
                [df_unused, df_fold_test, df_fold_test]).drop_duplicates(subset=self.unique_sample_headers, keep=False, ignore_index=True)

            df_fold_validation = pd.DataFrame()
            df_fold_train = pd.DataFrame()
            for group in self.included_groups:
                df_group = df_unused[df_unused[self.include_header].isin([
                                                                         group])]

                if k_fold == 1:
                    df_group_coice_validation = self.get_group_selection(
                        df_group, group_n[group])
                else:
                    df_available_group_validation = df_available_validation[df_available_validation[self.include_header].isin([
                                                                                                                              group])]
                    unique_entries = df_available_group_validation[self.intact_group_header].unique(
                    )

                    if  len(unique_entries) >= group_n[group]:
                        df_group_coice_validation = self.get_group_selection(
                        df_available_group_validation, group_n[group])
                    else:          
                        ## TODO if time it would probably be easier to have a column with "this has been sampled" rather than doing all this merging and dropping of dataframes 

                        df_group_coice_validation = df_available_group_validation.copy()

                        df_available_group_validation = self.reset_group_validation(df, df_available_group_validation, group)
                        df_available_validation = pd.concat(df_available_validation,df_available_group_validation )

                        number_of_added = group_n[group] - len(unique_entries)
                        add_to_validation = self.get_group_selection(
                            df_available_group_validation, number_of_added)
                        df_group_coice_validation = pd.concat(
                            [df_group_coice_validation, add_to_validation], ignore_index=True)
                        print(
                            "Reusing previous validation compounds for validation, for group: " + group)
                      
                df_fold_validation = pd.concat(
                    [df_fold_validation, df_group_coice_validation], ignore_index=True)

            df_available_validation = pd.concat(
                [df_available_validation, df_fold_validation, df_fold_validation], ignore_index=True).drop_duplicates(keep=False, ignore_index=True)
            df_fold_train = pd.concat(
                [df_unused, df_fold_validation, df_fold_validation], ignore_index=True).drop_duplicates(keep=False, ignore_index=True)
            df_unused = pd.concat([df_unused, df_fold_validation, df_fold_train],
                                  ignore_index=True).drop_duplicates(keep=False, ignore_index=True)
            k_fold_validation[k_fold-1] = df_fold_validation
            k_fold_train[k_fold-1] = df_fold_train
            if not df_unused.empty:
                print(
                    "WARNING: Didn't put all available compound into train or valid k-folds! Left over:")
                print(df_unused)

        print("Made train and valid sets for k-folds")

        return k_fold_train, k_fold_validation

    def get_k_folds_train_only(self, df, k_fold_test):
        """Get whatever's not in test into a df

        This function is meant to be used as an alternative for when we want to only have a train and a test set, but may need a nominal validation set to fit with programs.
        The train set for a fold will be whatever is not in the test set, and the validation set will be a copy of the train set.

        Keyword arguments:
        df -- the full dataset as a pandas dataframe
        k_fold_test -- a list of dataframe where each dataframe is the test set for a fold

        Returns:
        k_fold_train -- same structure as the k_fold_test but every row in df that is not in k_fold_test for that k_fold
        k_fold_valid -- a copy of k_fold_train
        """
        number_of_folds = self.k_folds
        k_fold_train = [None]*number_of_folds
        k_fold_validation = [None]*number_of_folds
        group_n = {}

        for k_fold in range(1, number_of_folds + 1):
            print("Starting k-fold: " + str(k_fold))
            df_unused = df.copy()
            df_k_fold_test = k_fold_test[k_fold-1]
            df_unused = pd.concat(
                [df_unused, df_k_fold_test, df_k_fold_test]).drop_duplicates(keep=False, ignore_index=True)

            k_fold_validation[k_fold-1] = df_unused
            k_fold_train[k_fold-1] = df_unused

        print("Made train and valid sets for k-folds")

        return k_fold_train, k_fold_validation

    def reset_sample_space_df(self, df, df_already_sampled, group):
        df_new_sample_space = df[df[self.include_header].isin([group])].copy()       
        df_new_sample_space = pd.concat(
            [df_new_sample_space, df_already_sampled, df_already_sampled], ignore_index=True).drop_duplicates(subset=self.intact_group_header, keep=False, ignore_index=True)

        df_new_sample_space.drop_duplicates(subset=self.unique_sample_headers, keep=False, ignore_index=True)
        return df_new_sample_space

    def get_group_selection(self, df_group, number_of_unique_entries):
        unique_entries = df_group[self.intact_group_header].unique()
        group_choice = np.random.choice(
            unique_entries, size=number_of_unique_entries, replace=False)
        df_group_choice = df_group[df_group[self.divide_on_header].isin(
            group_choice)]
        return df_group_choice

    def get_statistics(self, k_folds, df):
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

    def get_statistics_leave_one_out_test(self, k_folds, df):
        number_of_folds = self.k_folds
        for k_fold in range(0, number_of_folds):
            k_folds[k_fold]["k-fold"] = k_fold + 1

        non_unique_divider_test = self.non_unique_divider.copy()
        non_unique_divider_test.append("k-fold")
        df_statistics = pd.concat(k_folds, ignore_index=True)
        df_statistics = df_statistics.groupby(
            by=non_unique_divider_test).size().reset_index().rename(columns={0: 'records'})
        return df_statistics

    def make_path_available(self):
        if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print("made the output dir")

    def make_train_valid_statistics(self, k_folds_train, k_folds_validation, df):
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

    def save_k_folds(self, k_folds_test, k_folds_validation, k_folds_train):
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
