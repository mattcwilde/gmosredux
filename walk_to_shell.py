import os
from glob import glob

commands = []
file = open("run_redux_on_project.sh", 'w')
path = "/Users/mwilde/Dropbox/COS-Gemini/gemini_data.GN-2014A-Q-1/reduced_data/"
pwd = os.path.abspath(os.path.curdir)
REDUX_FILE = "/reduce_science.py"
for root, dirs, files in os.walk(path, topdown=True):
    for name in dirs:
        if 'name' not in ['database']:
            if 'mask' in name:
                curdir = os.path.join(root, name)
                print(curdir)
                science_files = glob(curdir+'J*fits')
                if len(science_files < 2):
                    folder = os.path.join(root, name)
                    string = "cd {}".format(folder)
                    file.write(string)
                    file.write('\n')
                    file.write('python reduce_science.py')
                    file.write('\n')
file.write("cd {}".format(os.getcwd()))
file.close()
# with open("run_redux_on_project.sh", 'w') as f:
#     f.write(commands)