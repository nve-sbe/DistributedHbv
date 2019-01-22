#!/usr/bin/tcsh -x

#./run_all.sh 

set RCM = OBS
set RCMVer = OBS
set hbvmodel = hbv

set Initialize = 0
set Intermezzo = 0
set RunHbv = 1
set ExtractProjections = 0
set CalcEnergyPotential = 0
set Gzip = 0

if( $Initialize == 1 ) then
#cp /data/HBVSim/CommonDataDistHbv/* .
mkdir timeseries
#mkdir energypotential
#mkdir results
#ln -s /home/shh/project/Code/hbv_2.11_clima_control/$hbvmodel .
#cp /data02/Ican/hbv/test/MuDyAero/code/$hbvmodel .
set Start = 1982
set End = 2012
while ( $Start <= $End ) 
mkdir $Start
@ Start++
end 
else
endif

if( $Intermezzo == 1 ) then
set Start = 1960
set End = 2005
while ( $Start <= $End ) 
cd $Start
gunzip hgw*
gunzip hsm*
cd ../
@ Start++
end 
else
endif

if( $RunHbv == 1 ) then
setenv METDATA /data02/Ican/vic_sim/past_1km/bil
#setenv METDATA /hdata/grid/metdata
./$hbvmodel  control_hbv_pm.txt
#mv  *var results
#R CMD BATCH ./nse_all.r
else
endif

if( $ExtractProjections == 1 ) then 
/home/shh/project/energy/run_all_2014.sh 1971 2000 $RCM $RCM $RCMVera 1
else
endif

if( $CalcEnergyPotential == 1 ) then #base: always 1981-2010 -> gjør bare for rcp45 og rcp84
else
endif


if( $Gzip == 1 ) then
set Start = 1960
set End = 2005
while ( $Start <= $End ) 
cd $Start
#gzip scf*
gzip hgw*
gzip hsm*
cd ../
@ Start++
end 
else
endif
