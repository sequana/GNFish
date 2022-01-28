import os
import re
from pprint import pprint as pp
import glob



class list_files():
    @staticmethod
    def list_files_method(self, path):

        folders = []
        

        for folder in glob.glob(path + "*/", recursive=True):
            files = []
            for file in glob.glob(folder + '*', recursive=True):
                files.append(file)
            folders.append([folder, files])

        return folders



