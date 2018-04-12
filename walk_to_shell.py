import os
import sys

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
file = open("run_redux_on_{}.sh".format(pid), 'w')

pwd = os.path.abspath(os.path.curdir)
REDUX_FILE = "/reduce_science.py"
for root, dirs, files in os.walk(path, topdown=True):
    for name in dirs:
        if 'name' not in ['database']:
            if 'mask' in name:
                curdir = os.path.join(root, name)
                print(curdir)
                folder = os.path.join(root, name)
                string = "cd {}".format(folder)
                file.write(string)
                file.write('\n')
                file.write('python reduce_science.py')
                file.write('\n')


                src = pwd+REDUX_FILE
                dst = os.path.join(root, name)+REDUX_FILE

                # delete the old version
                delete_file(dst)

                # copy symlink to that folder
                os.symlink(src, dst)
file.write("cd {}".format(os.getcwd()))
file.close()
# with open("run_redux_on_project.sh", 'w') as f:
#     f.write(commands)