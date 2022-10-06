
from make_k_folds_tvt import MakeTVTSets
# Make K-folds, including their train and validation parts, with csv files as outpus

class MakeKFoldsTVTSPECS(MakeTVTSets):

    def __init__(self,
                 #labels_path="/home/jovyan/Data/Specs/Specs_Labels_First_MiPhHo.csv", #"smaller images"
                 labels_path="/home/jovyan/Data/Specs/Labels.csv",#bigger images
                 output_dir='/home/jovyan/Inputs/SPECS_16Bit_Top10_MoA_no_negcon_5_folds/',
                 #
                 #include_groups = [], #Empty for everything included,
                 #Jordis top 10
                 #include_groups = ["Aurora kinase inhibitor","tubulin polymerization inhibitor","JAK inhibitor","protein synthesis inhibitor","HDAC inhibitor","topoisomerase inhibitor","PARP inhibitor","ATPase inhibitor","retinoid receptor agonist","HSP inhibitor"],
                 # top10 good images (nuclei cut off)
                 #include_groups=["heat shock response signalling agonist","phosphodiesterase inhibitor","DILI","methyltransferase inhibitor","estrogen receptor alpha modulator","cyclooxygenase inhibitor","pregnane x receptor agonist","PPAR receptor agonist","protein synthesis inhibitor","CC chemokine receptor antagonist"],
                 # TOP 10 MoA
                 include_groups = ["cyclooxygenase inhibitor", "phosphodiesterase inhibitor", "EGFR inhibitor", "HDAC inhibitor", "topoisomerase inhibitor", "glucocorticoid receptor agonist", "CC chemokine receptor antagonist", "adenosine receptor antagonist",  "histone lysine methyltransferase inhibitor", "HSP inhibitor"],  
                 # include_groups = ["negcon", "heat shock response signalling agonist", "DILI", "estrogen receptor alpha modulator", "phosphodiesterase inhibitor",  "cyclooxygenase inhibitor"], #Empty for everything included,
                 # include_groups = ["negcon","heat shock response signalling agonist", "phosphodiesterase inhibitor", "methyltransferase inhibitor","DILI","HDAC inhibitor","topoisomerase inhibitor", "mTOR inhibitor","NFkB pathway inhibitor","JAK inhibitor","pregnane x receptor agonist"],                    #include_groups = ["negcon","DNA polymerase inhibitor", "mTOR inhibitor", "topoisomerase inhibitor"],
                 #include_groups = ['pregnane x receptor agonist', 'DILI', 'estrogen receptor alpha modulator', 'tubulin polymerization inhibitor', 'topoisomerase inhibitor', 'heat shock response signalling agonist', 'methyltransferase inhibitor', 'aryl hydrocarbon receptor agonist', 'estrogen receptor alpha agonist', 'mitochondrial toxicity  agonist', 'retinoid receptor agonist', 'protein synthesis inhibitor', 'phosphodiesterase inhibitor', 'DNA polymerase inhibitor', 'mTOR inhibitor', 'PPAR receptor agonist', 'glucocorticoid receptor agonist', 'ATPase inhibitor', 'cyclooxygenase inhibitor', 'NFkB pathway inhibitor', 'angiotensin converting enzyme inhibitor', 'adenosine receptor antagonist', 'PARP inhibitor', 'JAK inhibitor', 'HSP inhibitor', 'HDAC inhibitor',  'CC chemokine receptor antagonist', 'Aurora kinase inhibitor','negcon'],
                 include_header="moa",
                 class_column_header="moa",
                 excluded_groups=[["negcon","poscon", "empty"], ["P015085"]],
                 excluded_groups_headers=["pert_type", "plate"],
                 exclude_images_path="/home/jovyan/Data/Specs/Flaggs/16Bit_images_nuclei_cut_above149.csv",
                 #exclude_images_path="",
                 intact_group_header='compound_id',
                 unique_sample_headers=['plate', 'well', 'site'],
                 image_number_heading="ImageNr",
                 k_folds="5",
                 divide_on_header='compound_id',
                 make_train_valid=True,
                 # 1 = 100%,  Percentage of images remaining afte the test set has been excluded
                 valid_fraction=0.25,
                leave_one_out = False,
                make_unique_validation = True,
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
    MakeKFoldsTVTSPECS().main()


