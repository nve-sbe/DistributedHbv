
#Define model domain for distributed HBV hydrological model
stationMask control_mask.txt

#Preprocessing for distributed HBV hydrological model
prehbv control_pre.txt

#Set environment variable for SeNorge meteorological data
source_set_env

#Run distributed HBV hydrological model
hbv control_hbv.txt

