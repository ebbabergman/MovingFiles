import os

def make_non_existing_path(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(path)