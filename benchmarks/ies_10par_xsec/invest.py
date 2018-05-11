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
    pst.write(os.path.join("template","pest_pareto.pst"))

    nreal_per = 5


    #pyemu.os_utils.run("pestpp-ies pest.pst /h :4004",cwd="master")
    pyemu.os_utils.start_slaves("template","pestpp-ies","pest_pareto.pst",num_slaves=20,master_dir="master")

if __name__ == "__main__":
    setup()