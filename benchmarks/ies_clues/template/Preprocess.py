##*********************************************************************************************************************
## Stand-alone script for running a steady-state surface-gw nitrogen leaching model, with linear exchanges
##*********************************************************************************************************************

# Imports

import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sparse
import scipy.sparse.linalg
from scipy import stats
import yaml  # for reading in parameter files
from inspect import getsourcefile
import time
from pprint import pprint # pretty print
# starttime = time.time()
##*********************************************************************************************************************
## Functions
##*********************************************************************************************************************

def check_fields(required_fields_list, checkframe,checkframename):
    # check all the list of required fields are present. Abort if they aren't
    test = set(required_fields_list)<= set(list(checkframe))
    if not test:
        missinglist = list(set(required_fields_list) - set(list(checkframe)))
        print('!!!!! Required fields not present in table: ',checkframename, '. Fields missing are: ',missinglist)
        sys.exit(0)

##*********************************************************************************************************************
##  Main program section
##*********************************************************************************************************************

pd.set_option('display.width', 150)

#  Set up project of folder name
projectdefault = 'Project_Default'

project = projectdefault # default project folder name, if none is passed to the script
if (len(sys.argv)) > 1: # if project folder name specified in command line
    if os.path.basename(sys.argv[0]) != 'pydevconsole.py': # to ignore the case where the code is called from pycharm console, where there are actually entries in argv
        project = (sys.argv[1]) # project folder name passed as command-line argument

# Set up folder names for inputs and outputs
path = os.path.abspath(getsourcefile(lambda:0))
if not(os.path.isdir(path)): # case where the file name is included
    path = os.path.dirname(path)

mydir = path + '/' + project  + '/'  # all files for a run are stored under this folder.

print ('Running project in folder ',mydir)
mydirInputs = mydir +'Inputs/'
mydirOutputs = mydir +'Inputs/'

# Read master file of file names. This is a dictionary of the filenames (without root directory)
with open(mydir + 'Filenames.yaml', 'r') as f:
    InputFiles = yaml.load(f)
#pprint(InputFiles)


######## 0. Read inputs

# Read parameters.
with open(mydirInputs + 'Parameters.yaml', 'r') as f:
    Parameters = yaml.load(f)

# Read reach file

# Read reaches (equivalent to subcatchments) in as a table.
ReachesFileName = InputFiles['ReachesFileName']
Reaches = pd.read_csv(mydirInputs + ReachesFileName,dtype={'ReachID': int, 'FromNode': int,'ToNode': int,'Length': float,'Area': float }) # length km, area km2
InitialReachesFieldList=list(Reaches)
checkfieldlist = str.split('ReachID,FromNode,ToNode,Length,Area,TotalRunoff,RechargeFraction,RechargeAgeY,StreamLossFraction,'
    'GWDischargeFraction,GWDeepDischargeFraction,GWShallowToDeepFraction,GWDeepToShallowFraction,GWExportFraction,GWDeepExportFraction,'
    'SourceYield,RechargeConcFactor,'
    'PointSourceLoadTperY,PointSourceLoadGWTperY,PointSourceFlow,PointSourceFlowGW,'
    'Abstraction,AbstractionGW,AbstractionGWDeep,'
    'AdditionalSourceLoadTperY,AdditionalSourceLoadGWTperY,AdditionalSourceLoadGWDeepTperY,'
    'AdditionalSourceFlow,AdditionalSourceFlowGW,AdditionalSourceFlowGWDeep,'
    'AdditionalSourceAgeY,AdditionalSourceAgeGWY,AdditionalSourceAgeGWDeepY,'
    'GWMixDepth,GWDeepMixDepth',sep=',')

if Parameters['IfGWDecayCoefficientFromAreaFractions'] == 1:
    checkfieldlist.extend(['GWFractionOxic','GWFractionMixed','GWFractionAnoxic','GWDeepFractionOxic','GWDeepFractionMixed','GWDeepFractionAnoxic'])
else:
    if Parameters['IfGWDecayFromCoefficients'] == 0:
        checkfieldlist.extend(['GWDecayFraction', 'GWDeepDecayFraction'])
    else:
        checkfieldlist.extend(['GWDecayK', 'GWDeepDecayK'])

check_fields(checkfieldlist, Reaches, 'Reaches')

Reaches.set_index('ReachID',inplace=True,drop=False) # set index of Reaches to the ReachID, for more efficient retrieval. Different from ReachIndex


nreaches = len(Reaches)
# Define a sequential ID (ReachIndex) for each reach. Starts at 0, ends at number of reaches minus 1.
Reaches["ReachIndex"] = list(range(nreaches))

# Read reach connectivity
ReachConnectivityFileName = InputFiles['ReachConnectivityFileName']  # filename for reaches connectivity override.
ReachConnectivityNew = pd.read_csv(mydirInputs + ReachConnectivityFileName,dtype={'FromReach': int, 'ToReach': int})
check_fields(["FromReach", "ToReach", "Fraction"], ReachConnectivityNew ,'ReachConnectivityNew')

# Read in groundwater lateral connectivity, for shallow then deep groundwater
# Don't take the discharge to streams or deep groundatwer into account at this stage

GWLateralConnectivityFileName = InputFiles['GWLateralConnectivityFileName']  # filename for groundwater connectivity.
GWConnectivity = pd.read_csv(mydirInputs + GWLateralConnectivityFileName,dtype={'FromReach':int,'ToReach':int})
check_fields(["FromReach", "ToReach", "FractionGW"], GWConnectivity,'GWConnectivity')

GWDeepLateralConnectivityFileName = InputFiles['GWDeepLateralConnectivityFileName']  # filename for groundwater connectivity.
GWDeepConnectivity = pd.read_csv(mydirInputs + GWDeepLateralConnectivityFileName,dtype={'FromReach':int,'ToReach':int})
check_fields(["FromReach", "ToReach", "FractionGWDeep"], GWDeepConnectivity,'GWDeepConnectivity')

#### Modify inputs

## First for shallow groundwater
# detect and correct for bidirectional exchanges

fromreaches= GWConnectivity['FromReach'].values
toreaches= GWConnectivity['ToReach'].values
fractiongw = GWConnectivity['FractionGW'].values
fractionreverse=[]
for i in range(0,len(GWConnectivity)):
    fromreach = fromreaches[i]
    toreach = toreaches[i]
    match = np.where((fromreaches == toreach) & (toreaches == fromreach))
    try:
        fractionreverse.append(fractiongw[match][0])
    except IndexError:
        fractionreverse.append(0)

fractiongw -= np.minimum(fractiongw,fractionreverse) #adjust the fraction groundwater to remove reverse flows
GWConnectivity['FractionGW'] = fractiongw

# Adjust GW lateral outflows so sum equals 1

#calculate the sum
SumOutGW  = GWConnectivity.pivot_table(index="FromReach", values=['FractionGW'], aggfunc=np.sum).reset_index()  # calculate the sums.  sum of lateral outflow first. Add in vertical later
SumOutGW.rename(columns={'FromReach': 'ReachID','FractionGW': 'FractionGWSum'}, inplace=True)
Reaches = pd.merge(Reaches,SumOutGW, how='left', on='ReachID')
Reaches.FractionGWSum = Reaches.FractionGWSum.fillna(0) # Replace sum of outflows where connectivity has not been defined (i.e. missing sum) with value of 0
Reaches.FractionGWSum += Reaches.GWExportFraction  #add vertical outflows and system loss

# make adjustments. Adjust all the outflows by the same fraction. (abstractions are not dealt with in this way, though, as they are fixed rates)
GWConnectivity = pd.merge(GWConnectivity,Reaches[['ReachID','FractionGWSum']],how='left',left_on='FromReach',right_on = 'ReachID')

sumgw = GWConnectivity.FractionGWSum.values
fractiongw = GWConnectivity.FractionGW.values
mask = sumgw > 1e-6
fractiongw[~mask]=0
fractiongw[mask] /= sumgw[mask]
GWConnectivity.FractionGW = fractiongw

sumgw =  Reaches.FractionGWSum.values
gwexportfraction = Reaches.GWExportFraction.values
mask = sumgw > 1e-6
gwexportfraction[~mask] = 0
gwexportfraction[mask] /= sumgw[mask]
Reaches.GWExportFraction = gwexportfraction

# modify vertical loss from gw to  ensure that sum is not greater than 1

Reaches['SumTemp'] = Reaches.GWDischargeFraction + Reaches.GWShallowToDeepFraction
Reaches.loc[Reaches.SumTemp > 1,'GWDischargeFraction'] /= Reaches.SumTemp
Reaches.loc[Reaches.SumTemp > 1,'GWShallowToDeepFraction'] /= Reaches.SumTemp
del Reaches['SumTemp']


## and for deep groundwater

# detect and correct for bidirectional exchanges

SumOutGWDeep  = GWDeepConnectivity.pivot_table(index="FromReach", values=['FractionGWDeep'], aggfunc=np.sum).reset_index()  # calculate the sums.  sum of lateral outflow first. Add in vertical later
SumOutGWDeep.rename(columns={'FromReach': 'ReachID','FractionGWDeep': 'FractionGWDeepSum'}, inplace=True)
Reaches = pd.merge(Reaches,SumOutGWDeep, how='left', on='ReachID')
Reaches.FractionGWDeepSum = Reaches.FractionGWDeepSum.fillna(0) # Replace sum of outflows where connectivity has not been defined (i.e. missing sum) with value of 0
Reaches.FractionGWDeepSum += Reaches.GWDeepExportFraction  #add vertical outflows and system loss

# make adjustments. Adjust all the outflows by the same fraction. (abstractions are not dealt with in this way, though, as they are fixed rates)
GWDeepConnectivity = pd.merge(GWDeepConnectivity,Reaches[['ReachID','FractionGWDeepSum']],how='left',left_on='FromReach',right_on = 'ReachID')

sumgwdeep = GWDeepConnectivity.FractionGWDeepSum.values
fractionwdeep = GWDeepConnectivity.FractionGWDeep.values
mask = sumgwdeep > 1e-6
fractionwdeep[~mask]=0
fractionwdeep[mask] /= sumgwdeep[mask]
GWDeepConnectivity.FractionGWDeep = fractionwdeep

sumgwdeep =  Reaches.FractionGWDeepSum.values
GWDeepexportfraction = Reaches.GWDeepExportFraction.values
mask = sumgwdeep > 1e-6
GWDeepexportfraction[~mask] = 0
GWDeepexportfraction[mask] /= sumgwdeep[mask]
Reaches.GWDeepExportFraction = GWDeepexportfraction

# Write the modified inputs
Reaches[InitialReachesFieldList].to_csv(mydirOutputs + 'Reaches_Modified.csv', index=False)
GWConnectivity[['FromReach','ToReach','FractionGW']].to_csv(mydirOutputs + 'GWConnectivity_Modified.csv', index=False)
GWDeepConnectivity[['FromReach','ToReach','FractionGWDeep']].to_csv(mydirOutputs + 'GWDeepConnectivity_Modified.csv', index=False)

#print(time.time() - starttime)