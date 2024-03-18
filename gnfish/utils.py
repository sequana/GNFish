import os
import re
import glob


def list_files(path):
    folders = []

    for folder in glob.glob(path + "*/", recursive=True):
        files = []
        for file in glob.glob(folder + "*", recursive=True):
            files.append(file)
        folders.append([folder, files])

    return folders
