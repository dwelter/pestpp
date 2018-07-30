import os
import pyemu

jco = pyemu.Jco.from_binary("test.jcb")

print(jco.row_names)