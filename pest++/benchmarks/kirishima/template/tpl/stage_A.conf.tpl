ptf ~
VENT_NORTHING 3538377
VENT_EASTING 110616
VENT_ELEVATION 1700

PLUME_HEIGHT ~PLUME_HEA           ~ 
ALPHA ~ALPHAA              ~ 
BETA ~BETAA               ~ 
ERUPTION_MASS ~ERUPTIONA           ~ 
MAX_GRAINSIZE -7 
MIN_GRAINSIZE 7   
MEDIAN_GRAINSIZE ~MEDIAN_GA           ~
STD_GRAINSIZE ~STD_GRAIA           ~ 

EDDY_CONST ~EDDY_CONA           ~

# diffusion coeff for large particles (m2/s)
DIFFUSION_COEFFICIENT ~DIFFUSIOA           ~ 

# threshold for change in diffusion (seconds fall time)
FALL_TIME_THRESHOLD ~FALL_TIMA           ~ 

# density model for the pyroclasts
LITHIC_DENSITY 	~   lithdena    ~
PUMICE_DENSITY 	~   pumdena     ~

#define column integration steps
COL_STEPS 400 
PART_STEPS 140 

# Note: 
# 0 = uniform distribution using threshold at PLUME_RATIO (no longer used)
# 1 = log-normal distribution using beta (no longer used)
# 2 = beta distribution using parameters alpha and beta (set below)
PLUME_MODEL 2
