import os
import re
from pprint import pprint as pp
import glob


def list_files(path):

        folders = []
        

        for folder in glob.glob(path + "*/", recursive=True):
            files = []
            for file in glob.glob(folder + '*', recursive=True):
                files.append(file)
            folders.append([folder, files])

        return folders



