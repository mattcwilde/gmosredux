from __future__ import print_function, division
import os
# import glob

# import reduce_science


def delete_file(name):
    """deletes the file in question if it exists 
    
    """
    try:
        os.remove(name)
    except OSError:
        pass

path = "/Users/mwilde/Dropbox/COS-Gemini/gemini_data.GN-2014A-Q-1/reduced_data/"
pwd = os.path.abspath(os.path.curdir)
REDUX_FILE = "/reduce_science.py"
for root, dirs, files in os.walk(path, topdown=True):
    for name in dirs:
        if 'name' not in ['database']:
            if 'mask' in name:

                print("walking through:",os.path.join(root, name))


                # src = pwd+REDUX_FILE
                # dst = os.path.join(root, name)+REDUX_FILE

                # delete the old version
                # delete_file(dst)

                # copy symlink to that folder
                # os.symlink(src, dst)

                # change to this directory
                os.chdir(os.path.join(root, name))
                print("changed directory to:",os.path.abspath(os.path.curdir))

os.chdir(pwd)
print("changed dir back to where we started:",os.path.abspath(os.path.curdir))