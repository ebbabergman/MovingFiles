
from email.mime import image
from make_k_folds_tvt import MakeTVTSets
# Make K-folds, including their train and validation parts, with csv files as outpus


class MakeKFoldsTVTBBBC021(MakeTVTSets):
                   
    def __init__(self,
                labels_path = "/home/jovyan/Data/BBBC021/BBBC021_Labels.csv",
                output_dir='/home/jovyan/Inputs/BBBC021_leave_one_out_No_Validation/',
                include_groups = [], #Empty for everything included,
                include_header = "moa",
                class_column_header = "moa",
                #exclude_groups = [["DMSO","Cholesterol-lowering","Eg5 inhibitors"]],
                excluded_groups = [["DMSO"]],
                excluded_groups_headers = ["moa"],
                exclude_images_path = "",
                intact_group_header = 'compound',
                unique_sample_headers = ["ImageNumber"],
                image_number_heading = "ImageNumber",   
                k_folds = "3",
                divide_on_header = 'compound',
                make_train_valid = True,
                valid_fraction = 0.25, # 1 = 100%,  Percentage of images remaining afte the test set has been excluded
                leave_one_out = True,
                make_unique_validation = False,
                ):
        super().__init__(labels_path = labels_path,
                output_dir =output_dir,
                include_groups = include_groups,
                include_header = include_header,
                class_column_header = class_column_header,
                excluded_groups = excluded_groups,
                excluded_groups_headers = excluded_groups_headers,
                exclude_images_path = exclude_images_path,
                intact_group_header = intact_group_header,
                unique_sample_headers = unique_sample_headers,
                image_number_heading = image_number_heading,   
                k_folds = k_folds,
                divide_on_header = divide_on_header,
                valid_fraction = valid_fraction)
        self.leave_one_out = leave_one_out
        self.make_unique_validation = make_unique_validation

    def main(self):
        print("Started making divisions for runs for BBBC021.")
        if self.leave_one_out :
            if self.make_unique_validation:
                self.make_leave_one_out()
            else:
                self.make_leave_one_out_train_test()
        else:
            if self.make_unique_validation:
                self.make_k_folds()
            else:
                self.make_k_folds_train_test()

        print("Finished. Find output in: " + self.output_dir)

if __name__ == "__main__":
    MakeKFoldsTVTBBBC021().main()