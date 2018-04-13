from __future__ import print_function

import os, sys, stat
import shutil
from glob import glob

def delete_file(name):
    """deletes the file in question if it exists

    """
    try:
        os.remove(name)
    except OSError:
        print('couldnt delete file:', name)
        pass


pid = str(sys.argv[1])
project_ids = ["GN-2014A-Q-1", "GN-2014B-LP-3", "GN-2015A-LP-3",
           "GS-2014A-Q-2", "GS-2014B-LP-4", "GS-2015A-LP-4"]

if pid in project_ids:
    print("running redux on project:", str(sys.argv[1]))
else:
    pass

path = "/Users/mwilde/Dropbox/COS-Gemini/gemini_data.{}/reduced_data/".format(pid)
raw_path  = "/Users/mwilde/Dropbox/COS-Gemini/gemini_data.{}/raw_data/".format(pid)
file = open("run_redux_on_{}.sh".format(pid), 'w')

pwd = os.path.abspath(os.path.curdir)
REDUX_FILE = "/reduce_science.py"
for root, dirs, files in os.walk(path, topdown=True):
    for name in dirs:
        if 'name' not in ['database'] and 'mask' in name:
            # make a shell script to run reduce_science.py
            target_dir = os.path.join(root, name)
            # print(target_dir)
            folder = os.path.join(root, name)
            string = "cd {}".format(folder)
            file.write(string)
            file.write('\n')
            file.write('python reduce_science.py')
            file.write('\n')

            # copy reduce_science.py symlink to that folder
            src = pwd+REDUX_FILE
            dst = target_dir+REDUX_FILE
            delete_file(dst)
            os.symlink(src, dst)

            # copy the raw fits file here
            files = glob(target_dir+"/[N|S]*fits")
            #print(target_dir + '/[N|S]*.fits')
            # print(" raw files:", files)
            # print(raw_path)

            file_name = os.path.basename(files[0])
            # print(files[0], file_name)
            # print(raw_path + files[0], "-->", target_dir + files[0])
            for f in files:
                file_name = os.path.basename(f)
                try:
                    if os.path.exists(os.path.join(raw_path, file_name)):
                        shutil.copy2(os.path.join(raw_path, file_name), os.path.join(target_dir, file_name))
                        # shutil.copy(raw_path+file_name, target_dir+file_name, follow_symlinks=False)
                        print(raw_path + file_name, "-->", target_dir + file_name)
                    else:
                        print('raw file doesnt exist:', os.path.join(raw_path, file_name))
                except (IOError, shutil.Error):
                    print('couldnt copy raw file:', file_name)
                    pass


file.write("cd {}".format(os.getcwd()))
file.close()

# make the shell script executable
os.chmod("run_redux_on_{}.sh".format(pid), stat.S_IRWXU)