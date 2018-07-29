import os

os.system("python Preprocess.py Project_Default_ModelA")
os.system("python CluesGW15.py Project_Default_ModelA")
os.system("python post_process.py Project_Default_ModelA")
