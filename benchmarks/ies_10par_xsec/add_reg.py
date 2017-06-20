import os
import pandas as pd
import pyemu

pst = pyemu.Pst(os.path.join("template_reg","pest.pst"))
pst.zero_order_tikhonov()
pst.control_data.pestmode = "regularization"
pst.pestpp_options["reg_frac"] = 0.1
pst.write(os.path.join("template_reg","pest.pst"))