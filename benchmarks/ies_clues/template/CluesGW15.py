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
mydirOutputs = mydir +'Outputs/'

# Read master file of file names. This is a dictionary of the filenames (without root directory)
with open(mydir + 'Filenames.yaml', 'r') as f:
    InputFiles = yaml.load(f)
#pprint(InputFiles)


######## 0. Read inputs

# Read parameters.
with open(mydirInputs + 'Parameters.yaml', 'r') as f:
    Parameters = yaml.load(f)

# Read reach file and establish default connectivity

# Read reaches (equivalent to subcatchments) in as a table.
ReachesFileName = InputFiles['ReachesFileName']
Reaches = pd.read_csv(mydirInputs + ReachesFileName,dtype={'ReachID': int, 'FromNode': int,'ToNode': int,'Length': float,'Area': float }) # length km, area km2
checkfieldlist = str.split('ReachID,FromNode,ToNode,Length,Area,TotalRunoff,RechargeFraction,RechargeAgeY,StreamLossFraction,'
    'GWDischargeFraction,GWDeepDischargeFraction,GWShallowToDeepFraction,GWDeepToShallowFraction,GWExportFraction,GWDeepExportFraction,'
    'SourceYield,RechargeConcFactor,'
    'PointSourceLoadTperY,PointSourceLoadGWTperY,PointSourceFlow,PointSourceFlowGW,'
    'Abstraction,AbstractionGW,AbstractionGWDeep,'
    'AdditionalSourceLoadTperY,AdditionalSourceLoadGWTperY,AdditionalSourceLoadGWDeepTperY,'
    'AdditionalSourceFlow,AdditionalSourceFlowGW,AdditionalSourceFlowGWDeep,'
    'AdditionalSourceAgeY,AdditionalSourceAgeGWY,AdditionalSourceAgeGWDeepY,'
    'GWMixDepth,GWDeepMixDepth,ConcFactorStream',sep=',')


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

##### 1. Calculate connectivity

# Derive default connectivity from stream reaches

# Connectivity for each reach is a table of from reach, to reach, and the fractions of outflow that goes to each of those reaches (excluding stream losses to GW). Can have multiple entries for each from reach.

NodeConnDict = {}  # create dictionary of lists. Key is from node, item is list of reaches flowing from that node.
ReachConnList = [] # list of dictionaries, Each item of list has information about one connection, stored as a dictionary.
Reaches['NReachesDS'] = 0

for reachid in Reaches.ReachID.tolist():  # For each reach, find the from node and add to the list of reaches flowing from that node
    fromnode = Reaches.at[reachid, 'FromNode']
    if fromnode in NodeConnDict: #item already in the dictionary
        NodeConnDict[fromnode].append(reachid) #append the reach to the list for that node
    else: # new item in the dictionary
        NodeConnDict[fromnode] = [reachid] #add new list of reaches to the dictionary
for reachid in Reaches.ReachID.tolist():
    tonode = Reaches.at[reachid, 'ToNode'] # for each reach, find the number of reaches downstream and create a new connectivity file
    if tonode in NodeConnDict:
        nreachesds = len(NodeConnDict[tonode])
        Reaches.at[reachid, 'NReachesDS'] = nreachesds
        for reachds in NodeConnDict[tonode]:
            itemadd = {"FromReach": reachid, "ToReach": reachds,"Fraction": 1 / nreachesds}  # Add reach to the combined list of downstream reaches. Split flow amongst downstream reaches
        ReachConnList.append(itemadd)
ReachConnectivity = pd.DataFrame(ReachConnList)

# Override connectivity reach flow allocation read in from a connectivity file
# This would generally only be needed if there is a flow split and the flows are not split evenly.

for i in range(0, len(ReachConnectivityNew)):
    reach = ReachConnectivityNew.iloc[i]
    if len(ReachConnectivity.loc[(ReachConnectivity.FromReach == reach.FromReach) & (ReachConnectivity.ToReach == reach.ToReach),"Fraction"]) == 0: # a new connection
        ReachConnectivity.append({'FromReach':reach.FromReach,'ToReach':reach.FromReach,'Fraction':reach.Fraction})
    else: # replace existing fraction
        ReachConnectivity.loc[(ReachConnectivity.FromReach == reach.FromReach) & (ReachConnectivity.ToReach == reach.ToReach),"Fraction"] = reach.Fraction
ReachConnectivity.to_csv(mydirOutputs +'SurfaceConnectivity.csv',index=False) #Write connectivity. Handy for setting up the groundwater connectivity defaults.

# Check or ensure that stream connectivity sums to 1. This is before stream losses are considered. # Currently not taking any action on this
# Note that for terminal reaches, diversions from that reach are not permitted.
SumToDS = ReachConnectivity.pivot_table(index="FromReach", values=['Fraction'], aggfunc=np.sum).reset_index()  # calculate the sums
SumToDS.rename(columns={'FromReach': 'ReachID'}, inplace=True)
Temp = pd.merge(Reaches[['ReachID', 'NReachesDS']], SumToDS, how='left', on='ReachID')
Temp.Fraction = Temp.Fraction.fillna(0)
IncompleteList = Temp[(np.abs(Temp.Fraction -1)  > 1e-3) & (Temp.NReachesDS != 0)].ReachID.tolist()

# modify stream lateral connectivity to take account of StreamLoss. (don't correct for abstrations, because they are a fixed volumetric rate rather than proportion)

ReachConnectivity = ReachConnectivity.merge(Reaches[['ReachID','StreamLossFraction']], how='left', left_on='FromReach',right_on='ReachID') # firt merge the stream loss into the reach connectivity frame. This allows for multiple records with the same fromreach
ReachConnectivity.Fraction *= (1-ReachConnectivity.StreamLossFraction)

# Ensure GW lateral outflows sum to 1

SumOutGW  = GWConnectivity.pivot_table(index="FromReach", values=['FractionGW'], aggfunc=np.sum).reset_index()  # calculate the sums.  sum of lateral outflow first. Add in vertical later
SumOutGW.rename(columns={'FromReach': 'ReachID','FractionGW': 'FractionGWSum'}, inplace=True)
Reaches = pd.merge(Reaches,SumOutGW, how='left', on='ReachID')

Reaches.FractionGWSum = Reaches.FractionGWSum.fillna(0) # Replace sum of outflows where connectivity has not been defined (i.e. missing sum) with value of 0
Reaches.FractionGWSum += Reaches.GWExportFraction  #add vertical outflows and system loss
IncompleteList = Reaches[(np.abs(Reaches.FractionGWSum - 1) > 1e-3) ].ReachID.tolist()

if Parameters['IfAdjustGWOutflowsForContinuity'] == 0: #Don't make adjustments. Outflows must sum to 1, unless the reach is flagged.
    if IncompleteList:  # there are some items in the incomplete list
        print("Error!!. There are shallow groundwater reservoirs which do not have sum of outflows equal to 1. Reach identifiers are: ", IncompleteList)
        sys.exit(0)
else: # make adjustments. Adjust all the lateral outflows by the same fraction. (abstractions are not dealt with in this way, though, as they are fixed rates)
    GWConnectivity = pd.merge(GWConnectivity,Reaches[['ReachID','FractionGWSum']],how='left',left_on='FromReach',right_on = 'ReachID')
    ZeroList = Reaches[(np.abs(Reaches.FractionGWSum)  < 1e-3) ].ReachID.tolist()
    if ZeroList:
        print("Error!!. There are groundwater reservoirs that have a sum of outflows equal to 0. ReachID: ",ZeroList)
        sys.exit(0)
    GWConnectivity.FractionGW /= GWConnectivity.FractionGWSum
    Reaches.GWExportFraction /= Reaches.FractionGWSum

#...and for deep groundwater
SumOutGWDeep  = GWDeepConnectivity.pivot_table(index="FromReach", values=['FractionGWDeep'], aggfunc=np.sum).reset_index()  # calculate the sums
SumOutGWDeep.rename(columns={'FromReach': 'ReachID','FractionGWDeep': 'FractionGWDeepSum'}, inplace=True)
Reaches = pd.merge(Reaches,SumOutGWDeep, how='left', on='ReachID')
Reaches.FractionGWDeepSum = Reaches.FractionGWDeepSum.fillna(0) # Replace sum of outflows where connectivity has not be defined (i.e. missing sum) with value of 0
Reaches.FractionGWDeepSum += Reaches.GWDeepExportFraction
IncompleteList = Reaches[(np.abs(Reaches.FractionGWDeepSum - 1)  > 1e-3) ].ReachID.tolist()

if Parameters['IfAdjustGWOutflowsForContinuity'] == 0: #Don't make adjustments. Outflows must sum to 1
    if IncompleteList:  # there are some items in the incomplete list
        print("Error!!. There are deep groundwater reservoirs which do not have sum of outflows equal to 1. Reach identifiers are: ", IncompleteList)
        sys.exit(0)
else:  # make adjustments. Adjust all the outflows by the same fraction. (abstractions are not dealt with in this way, though, as they are fixed rates)
    GWDeepConnectivity = pd.merge(GWDeepConnectivity, Reaches[['ReachID', 'FractionGWDeepSum']], how='left', left_on='FromReach',right_on = 'ReachID')
    ZeroList = ZeroList = Reaches[(np.abs(Reaches.FractionGWDeepSum)  < 1e-3) ].ReachID.tolist()
    if ZeroList:
        print("Error!!. There are deep groundwater reservoirs that have a sum of outflows equal to 0. ReachID: ",ZeroList)
        sys.exit(0)
    GWDeepConnectivity.FractionGWDeep /= GWDeepConnectivity.FractionGWDeepSum
    Reaches.GWDeepExportFraction /= Reaches.FractionGWDeepSum

# modify GW connectivity and lateral export fraction to take account of discharge to stream and deep groundwater. Formula assumes lateral outflows add to 1.
GWConnectivity = GWConnectivity.merge(Reaches[['ReachID','GWDischargeFraction','GWShallowToDeepFraction']], how='left', left_on='FromReach',right_on='ReachID')
GWConnectivity.FractionGW *= (1-GWConnectivity.GWDischargeFraction-GWConnectivity.GWShallowToDeepFraction)
Reaches.GWExportFraction *= (1-Reaches.GWDischargeFraction-Reaches.GWShallowToDeepFraction) #also adjust the lateral system loss

GWDeepConnectivity = GWDeepConnectivity.merge(Reaches[['ReachID','GWDeepDischargeFraction','GWDeepToShallowFraction']], how='left', left_on='FromReach',right_on='ReachID')
GWDeepConnectivity.FractionGWDeep *= (1-GWDeepConnectivity.GWDeepDischargeFraction-GWDeepConnectivity.GWDeepToShallowFraction)
Reaches.GWDeepExportFraction *= (1-Reaches.GWDeepDischargeFraction-Reaches.GWDeepToShallowFraction) #also adjust the lateral system loss

Reaches.set_index('ReachID',inplace=True,drop=False) # reset the index, which was lost in the merge. This is used later for speed

######## 2. Set up flow calculations

nelements = 3 * nreaches

A = sparse.lil_matrix((nelements, nelements), dtype= float)  #Empty lil matrix. First nreaches are for stream reaches, second are for groundwater

# add the surface water connectivity
for r in ReachConnectivity.itertuples():
    LocFromReach = Reaches.at[np.long(r.FromReach),"ReachIndex"] #Index of the from reach
    LocToReach = Reaches.at[np.long(r.ToReach),"ReachIndex"]
    A[LocFromReach,LocToReach]=r.Fraction

# add surface water losses to groundwater.
for s in Reaches.itertuples():
    A[s.ReachIndex, s.ReachIndex+nreaches] = s.StreamLossFraction

# add the shallow groundwater lateral connectivity
for g in GWConnectivity.itertuples():
    LocFromReach = Reaches.at[np.long(g.FromReach),"ReachIndex"]
    LocToReach = Reaches.at[np.long(g.ToReach),"ReachIndex"]
    A[LocFromReach+nreaches,LocToReach+nreaches]=g.FractionGW

# add shallow groundwater vertical connectivity
for v in Reaches.itertuples():
    A[nreaches + v.ReachIndex, v.ReachIndex] = v.GWDischargeFraction
    A[nreaches + v.ReachIndex, 2 * nreaches + v.ReachIndex] = v.GWShallowToDeepFraction

# add the deep groundwater lateral connectivity
for g in GWDeepConnectivity.itertuples():
    LocFromReach = Reaches.at[np.long(g.FromReach),"ReachIndex"]
    LocToReach = Reaches.at[np.long(g.ToReach),"ReachIndex"]
    A[LocFromReach+2*nreaches,LocToReach+2*nreaches]=g.FractionGWDeep

# add deep groundwater vertical connectivity
for v in Reaches.itertuples():
    A[2* nreaches + v.ReachIndex, v.ReachIndex] = v.GWDeepDischargeFraction
    A[2* nreaches + v.ReachIndex, nreaches + v.ReachIndex] = v.GWDeepToShallowFraction

# AA = A.todense() #for viewing

A = A.tocsr() # convert to csr format ready for linear algebra calculations

# set up B matrix from A matrix
I = sparse.identity(nelements, format="csr")  # sparse identity matrix, csr format
B = I - A.transpose()

# Set up flow sources

# runoff and recharge
brunoff = Reaches.TotalRunoff * (1 - Reaches.RechargeFraction)
brecharge = Reaches.TotalRunoff *  Reaches.RechargeFraction

b = np.concatenate((brunoff.values,brecharge.values,np.full(nreaches, 0, dtype= float)))
#Point source flows
p  = np.concatenate((Reaches['PointSourceFlow'].values,Reaches['PointSourceFlowGW'].values,np.full(nreaches, 0, dtype= float)))
#Abstractions
x  = np.concatenate((Reaches['Abstraction'].values,Reaches['AbstractionGW'].values,Reaches['AbstractionGWDeep'].values))

#Additional source flows
s  = np.concatenate((Reaches['AdditionalSourceFlow'].values,Reaches['AdditionalSourceFlowGW'].values,Reaches['AdditionalSourceFlowGW'].values))

########## 3.Run flow calculations
i=0
while True: #while loop to deal with potential negative flows associated with abstractions
    #x=np.full(nelements,0)
    FlowOut = sparse.linalg.spsolve(B, b + p + s - x)
    #stats.describe(FlowOut[0:nreaches])

    Reaches['FlowOut']= FlowOut[0:nreaches]
    Reaches['FlowOutGW'] =  FlowOut[nreaches:2*nreaches]
    Reaches['FlowOutGWDeep'] =  FlowOut[2*nreaches:3*nreaches]

    FlowOutTotal = FlowOut + x # Total flow out including abstractions
    Reaches['FlowOutTotal'] =  FlowOutTotal[0:nreaches]
    Reaches['FlowOutGWTotal'] =  FlowOutTotal[nreaches:2*nreaches]
    Reaches['FlowOutGWDeepTotal'] =  FlowOutTotal[2*nreaches:3*nreaches]

    # Check for negative flows in reaches with abstraction, and make modifications to avoid negative flows
    ReachMask = ((Reaches.Abstraction > 0) & (Reaches['FlowOut']<-1e-10)) #allow a small negative flow out for numerical tolerance
    ReachMaskGW = ((Reaches.AbstractionGW > 0) & (Reaches['FlowOutGW'] < -1e-10))
    ReachMaskGWDeep = ((Reaches.AbstractionGWDeep > 0) & (Reaches['FlowOutGWDeep'] < -1e-10))
    ReachesException = Reaches[ReachMask  | ReachMaskGW | ReachMaskGWDeep]
    if i == 0:
        ReachesException.to_csv(mydirOutputs + 'AbstractionsTooLarge.csv', index=False)
    i += 1
    if len(ReachesException) > 0: #there were negative flows associated with abstraction
        print('!!! Abstractions caused negative flow. Reaches with negative flow and abstractions are listed in file AbstractionsTooLarge.csv. Abstractions reduced to give zero total flow out')
        # adjust abstraction
        Reaches.loc[ReachMask,'Abstraction'] =  np.maximum(Reaches.loc[ReachMask,'Abstraction'] + Reaches.loc[ReachMask,'FlowOut'],0)#adjust surface water abstraction
        Reaches.loc[ReachMaskGW, 'AbstractionGW'] =  np.maximum(Reaches.loc[ReachMaskGW,'AbstractionGW'] + Reaches.loc[ReachMaskGW,'FlowOutGW'],0)
        Reaches.loc[ReachMaskGWDeep, 'AbstractionGWDeep'] =  np.maximum(Reaches.loc[ReachMaskGWDeep,'AbstractionGWDeep'] + Reaches.loc[ReachMaskGWDeep,'FlowOutGWDeep'],0)
        x  = np.concatenate((Reaches['Abstraction'].values,Reaches['AbstractionGW'].values,Reaches['AbstractionGWDeep'].values)) # update the extractions vector
    else:
        break # there were no negative flows associated with abstraction
    if i > 10:
        print('!!! More than 10 iterations to get resolve negative flows by adjusting abstractions. Stopping programme')
        sys.exit(0)


# Calculate the flows out of reaches reach after losses (subtract stream losses from the flow through the reach, or gw dischage from flow through the reservoir)
FlowOutLateral = np.full(nelements, 0, dtype= float)
FlowOutLateral[0:nreaches] = FlowOut[0:nreaches]*(1 - Reaches.StreamLossFraction)  #OK to multiply these like this, because both are sorted by ReachIndex at this stage.
FlowOutLateral[nreaches:2 * nreaches] = FlowOut[nreaches:2 * nreaches] * (1 - Reaches.GWDischargeFraction - Reaches.GWShallowToDeepFraction)  # get the dischargefraction from the reach files to ensure correct sort order. Originally from GWVerticalConnectivity
FlowOutLateral[2*nreaches:nelements] = FlowOut[2*nreaches:nelements] * (1 - Reaches.GWDeepDischargeFraction - Reaches.GWDeepToShallowFraction)  # get the dischargefraction from the reach files to ensure correct sort order. Originally from GWVerticalConnectivity
Reaches['FlowOutLateral'] = FlowOutLateral[0:nreaches]

######## 4. Set up water quality calculations

# Derive load recharge and direct runoff
LoadDiffuse = Reaches.Area * Reaches.SourceYield*1000000/86400/365 #g/s
Reaches['CDirect'] = LoadDiffuse/((1-Reaches.RechargeFraction*(1-Reaches.RechargeConcFactor))*Reaches.TotalRunoff) #(g/m3)

#  Calculate surface water decay for each element (this could be flow-dependent, so do it now that flows have been calculated).

Qtemp = FlowOutTotal[0:nreaches].copy()
Qtemp[Qtemp < 0.001] = 0.001 # replace very small flows with 1 L/s. This mirrors the sparrow code. Avoids unrealistically large decay or difficulties with negative flows
Reaches['StreamDecayK'] = Parameters['StreamDecayCoefficient'] * Qtemp ** Parameters['StreamDecayExponent'] # per km
Reaches['StreamTransferFraction'] = np.exp(-Reaches.StreamDecayK*Reaches.Length) # this is the factor remaining after decay. Contrast with the groundwater loss fractions reach in.

# Groundwater decay

Reaches['VolumeGWShallow'] = Reaches.GWMixDepth.values * Reaches.Area.values * 1000000
Reaches['VolumeGWDeep'] = Reaches.GWDeepMixDepth.values * Reaches.Area.values * 1000000

if Parameters['IfGWDecayCoefficientFromAreaFractions'] == 1:
    Reaches['GWDecayK'] = (Reaches['GWFractionOxic']*Parameters['GWDecayCoefficientOxic'] +
        Reaches['GWFractionMixed'] * Parameters['GWDecayCoefficientMixed'] +
        Reaches['GWFractionAnoxic'] * Parameters['GWDecayCoefficientAnoxic'])
    Reaches['GWDeepDecayK'] = (Reaches['GWDeepFractionOxic'] * Parameters['GWDeepDecayCoefficientOxic'] +
        Reaches['GWDeepFractionMixed'] * Parameters['GWDeepDecayCoefficientMixed'] +
        Reaches['GWDeepFractionAnoxic'] * Parameters['GWDeepDecayCoefficientAnoxic'])

if Parameters['IfGWDecayFromCoefficients'] == 1: #otherwise decay fractions are in input value for each subcatchment
    Tempmask = Reaches.FlowOutGWTotal>1e-4
    Reaches['GWDecayFraction']=1 #default value for small flows, all decayed
    Reaches.loc[Tempmask,'GWDecayFraction'] = 1- np.exp(-Reaches.loc[Tempmask,'GWDecayK']  * Reaches.loc[Tempmask,'VolumeGWShallow'] / Reaches.loc[Tempmask,'FlowOutGWTotal']/86400 / 365)
    Tempmask = Reaches.FlowOutGWDeepTotal>1e-4
    Reaches['GWDeepDecayFraction']=1
    Reaches.loc[Tempmask,'GWDeepDecayFraction'] = 1- np.exp(-Reaches.loc[Tempmask,'GWDeepDecayK']  * Reaches.loc[Tempmask,'VolumeGWDeep'] / Reaches.loc[Tempmask,'FlowOutGWDeepTotal']/86400 / 365)

TransferVector = np.concatenate((Reaches['StreamTransferFraction'].values,1-Reaches['GWDecayFraction'].values,1-Reaches['GWDeepDecayFraction'].values)) # Stacked vector of surface and groundwater transfer fractions. 1-delta in the documentation. where delta is the decay fraction.
TransferArray = sparse.diags(TransferVector,0, format = 'csr')

# set up B matrix for contamination

AbstractionFractionVector = np.full(nelements,0,dtype = float) # default value of 0 will apply if total flow out is less than zero

AbstractionFractionVector[FlowOutTotal > 0] = x[FlowOutTotal > 0]/FlowOutTotal[FlowOutTotal > 0]
NonAbstractionFractionArray =  sparse.diags(1-AbstractionFractionVector,0, format = 'csr')
Bcontam = I - TransferArray * A.transpose()* NonAbstractionFractionArray

# Set up contaminant sources (combined diffuse and Additional)

Reaches['PointSourceLoad'] = Reaches.PointSourceLoadTperY /86400/365 * 1000000 # convert to g/s
Reaches['PointSourceLoadGW'] = Reaches.PointSourceLoadGWTperY /86400/365 * 1000000 # convert to g/s
Reaches['AdditionalSourceLoad'] = Reaches.AdditionalSourceLoadTperY /86400/365 * 1000000 # convert to g/s
Reaches['AdditionalSourceLoadGW'] = Reaches.AdditionalSourceLoadGWTperY /86400/365 * 1000000 # convert to g/s
Reaches['AdditionalSourceLoadGWDeep'] = Reaches.AdditionalSourceLoadGWDeepTperY /86400/365 * 1000000 # convert to g/s

d = np.concatenate((Reaches.CDirect.values*b[0:nreaches],
                    Reaches.CDirect.values*b[nreaches:2 * nreaches]*Reaches.RechargeConcFactor.values,
                    np.full(nreaches, 0, dtype=float))) # diffuse sources g/s
ps  = np.concatenate((Reaches['PointSourceLoad'].values,Reaches['PointSourceLoadGW'].values,np.full(nreaches, 0, dtype=float)))
v  = np.concatenate((Reaches['AdditionalSourceLoad'].values,Reaches['AdditionalSourceLoadGW'].values,Reaches['AdditionalSourceLoadGWDeep'].values))

bcontam = TransferVector *(d + ps + v)

######## 5. Run water quality calculations

FluxOutTotal = sparse.linalg.spsolve(Bcontam, bcontam)  # Solve for flux out including that associated with abstractions. Units g/s
FluxOut = FluxOutTotal * (1-AbstractionFractionVector) # flux out to other elements

# Calculate the flux out of reaches reach after losses (subtract losses from the flow through the reach)
FluxOutLateral = np.full(nelements, 0, dtype= float)
FluxOutLateral[0:nreaches] = FluxOut[0:nreaches]*(1 - Reaches.StreamLossFraction)
FluxOutLateral[nreaches:2*nreaches] = FluxOut[nreaches:2*nreaches]*(1 - Reaches.GWDischargeFraction - Reaches.GWShallowToDeepFraction) # get the dischargefraction from the reach files to ensure correct sort order. Originally from GWVerticalConnectivity
FluxOutLateral[2*nreaches:nelements] = FluxOut[2*nreaches:nelements]*(1 - Reaches.GWDischargeFraction - Reaches.GWDeepToShallowFraction) # get the dischargefraction from the reach files to ensure correct sort order. Originally from GWVerticalConnectivity

Concentration = np.full(nelements,0,dtype = float)
Concentration[FlowOutTotal > 0] = FluxOutTotal[FlowOutTotal > 0]/FlowOutTotal[FlowOutTotal > 0] #g/m3
Reaches['ConcentrationFW'] = Concentration[0:nreaches]
Reaches['ConcentrationMedian'] = Reaches['ConcentrationFW'] * Reaches['ConcFactorStream']


Reaches['FluxOutTperY']= FluxOut[0:nreaches] * 86400 * 365 / 1e6  # Converted from g/s to t/year
Reaches['FluxOutTotalTperY']= FluxOutTotal[0:nreaches] * 86400 * 365 / 1e6
Reaches['FluxOutLateralTperY']= FluxOutLateral[0:nreaches]  * 86400 * 365 / 1e6

######## 6. Set up age calculations


Reaches['StreamVolume'] = Reaches.Length *1000 * Parameters['StreamAreaCoefficient'] * FlowOutTotal[0:nreaches]** Parameters['StreamAreaExponent'] # Calculate Stream volume m3

AgeFluxGWShallow = np.where(FluxOutTotal[nreaches:2*nreaches]== 0, 0, Reaches.VolumeGWShallow)#units of m3 #correct the age flux to zero if there is no flow out
AgeFluxGWDeep = np.where(FluxOutTotal[2*nreaches:nelements]== 0, 0, Reaches.VolumeGWDeep)

Volume = np.concatenate((Reaches.StreamVolume.values,AgeFluxGWShallow,AgeFluxGWDeep)) # Stream and groundwater volume vector

AgeFluxAdditional = np.concatenate((Reaches.AdditionalSourceAgeY.values,Reaches.AdditionalSourceAgeGWY.values,Reaches.AdditionalSourceAgeGWDeepY.values))*365*86400*s
AgeFluxRecharge = np.concatenate((np.full(nreaches,0,dtype=float), Reaches.RechargeAgeY.values,np.full(nreaches,0,dtype=float)))*b


######## 7. Run age calculations
AgeFluxTotal = sparse.linalg.spsolve(B, AgeFluxAdditional + AgeFluxRecharge + Volume )  # Solve for age flux, m3

Age = np.full(nelements,0,dtype = float) # default value of 0 will apply if total flow out is less than zero
Age[FlowOutTotal > 0] = AgeFluxTotal[FlowOutTotal > 0]/FlowOutTotal[FlowOutTotal > 0]
AgeY = Age/86400/365

Reaches['AgeY'] = AgeY[0:nreaches]
Reaches['AgeGWY'] = AgeY[nreaches:2*nreaches]
Reaches['AgeGWDeepY'] = AgeY[2*nreaches:nelements]

####### 6. Write outputs

OutFields = [
'ReachID','FromNode','ToNode','Length','Area','TotalRunoff','RechargeFraction','StreamLossFraction','GWDischargeFraction',
'GWDeepDischargeFraction','GWShallowToDeepFraction','GWDeepToShallowFraction','PointSourceFlow','PointSourceFlowGW',
'AdditionalSourceFlow','AdditionalSourceFlowGW','AdditionalSourceFlowGWDeep','Abstraction','AbstractionGW','AbstractionGWDeep',
'SourceYield','RechargeConcFactor','PointSourceLoadTperY','PointSourceLoadGWTperY','AdditionalSourceLoadTperY','AdditionalSourceLoadGWTperY',
'AdditionalSourceLoadGWDeepTperY','RechargeAgeY','AdditionalSourceAgeY','AdditionalSourceAgeGWY','AdditionalSourceAgeGWDeepY',
'GWDecayK','GWDeepDecayK','GWDecayFraction','GWDeepDecayFraction','GWMixDepth','GWDeepMixDepth',
'FlowOutTotal','FlowOutLateral','FluxOutTotalTperY',
'FluxOutLateralTperY','ConcentrationFW','ConcentrationMedian','AgeY','AgeGWY','AgeGWDeepY'
]

Reaches[OutFields].to_csv(mydirOutputs + 'Reaches_Output_Full.csv', index=False)

OutFields2 = ['ReachID','TotalRunoff','RechargeFraction','StreamLossFraction','GWDischargeFraction','SourceYield','RechargeConcFactor',
            'GWDecayK','GWDeepDecayK','GWDecayFraction','GWDeepDecayFraction','GWMixDepth',
            'FlowOutLateral','FluxOutLateralTperY','ConcentrationMedian','AgeY','AgeGWY','GWExportFraction']
Reaches[OutFields].to_csv(mydirOutputs + 'Reaches_Output_Short.csv', index=False)

print ('Finished run')