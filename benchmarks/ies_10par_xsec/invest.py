import os
import shutil
import numpy as np
import pandas as pd
import pyemu
base_d = "template"

def setup():
    if os.path.exists("master"):
        shutil.rmtree("master")
    shutil.copytree(base_d,"master")
    pst = pyemu.Pst(os.path.join("master","pest.pst"))
    pst.pestpp_options = {}
    pst.control_data.noptmax = 3
    pst.write(os.path.join("master","pest.pst"))
    pyemu.os_utils.run("pestpp-ies pest.pst",cwd="master")

if __name__ == "__main__":
    setup()