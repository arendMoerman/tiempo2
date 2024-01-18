import shutil
import os
import json

class FileHandler(object):
    def __init__(self):
        self.savePaths = {}

    def addPath(self, path, path_key):
        self.savePaths[path_key] = path

    def saveInputs(self, interfaceObj, name, path_key):
        """!
        Save input dictionaries to JSON format
        """
        path = self.savePaths[path_key]

        os.makedirs(os.path.join(path, name), exist_ok=True)
        outfile = os.path.join(path, name+".json")

        with open(outfile, "w") as outfile: 
            json.dump(, outfile)

