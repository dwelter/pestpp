ptf #
// This parameter file corresponds to the 'species 1' column of 
// the sheet '3PG_Parameters' in the Excel / VB version of 3PG. 

// Allometric relationships & partitioning
"Foliage:stem partitioning ratio @ D=2 cm"          #fsp0_2  #
"Foliage:stem partitioning ratio @ D=20 cm"         #fsp1_2  #
"Constant in the stem mass v. diam. relationship"   #smdc_2  #
"Power in the stem mass v. diam. relationship"      #smdp_2  #
"Maximum fraction of NPP to roots"                  #npp1_2  #
"Minimum fraction of NPP to roots"                  #npp0_2  #

// Temperature modifier (fT)
"Minimum temperature for growth"                             5
"Optimum temperature for growth"                            25
"Maximum temperature for growth"                            40

// Frost modifier (fFRost)
"Days production lost per frost day"                         0

// Litterfall & root turnover
"Maximum litterfall rate"                          #lrm_2    #
"Litterfall rate at t = 0"                         #lr0_2    #
"Age at which litterfall rate has median value"             24
"Average monthly root turnover rate"                     0.015

// Conductance
"Maximum canopy conductance"                              0.02
"Maximum stomatal conductance"                     #msc_2    #
"Defines stomatal response to VPD"                        0.05
"Canopy boundary layer conductance"                        0.2

// Fertitlity effects
"Value of 'm' when FR = 0"                                   0
"Value of 'fNutr' when FR = 0"                               1

// Soil water modifier (fSW)
"Moisture ratio deficit for fq = 0.5"                      0.7
"Power of moisture ratio deficit"                            9

// Stem numbers
"Max. stem mass per tree @ 1000 trees/hectare"             300

// Age modifier (fAge)
"Maximum stand age used in age modifier"                    50
"Power of relative age in function for fAge"                 4
"Relative age to give fAge = 0.5"                         0.95

// Canopy structure and processes
"Specific leaf area at age 0"                      #sla0_2   #
"Specific leaf area for mature leaves"             #sla1_2   #
"Age at which specific leaf area = (SLA0+SLA1)/2"         0.92
"Extinction coefficient for absorption of PAR by canopy"   0.5
"Age at canopy cover"                                        1
"Proportion of intercepted rainfall evaporated from canopy" 0.15
"Canopy quantum efficiency"                                0.055

// Branch and bark fraction (fracBB)
"Branch and bark fraction at age 0"                       0.75
"Branch and bark fraction for mature stands"              0.15
"Age at which fracBB = (fracBB0+fracBB1)/2"                1.5

// Various
"Ratio NPP/GPP"                                           0.47
"Basic density"                                            0.5
