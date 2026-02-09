#pragma once
/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments.          *
 *                                                                                      *
 *  Author: Stein Beldring, Norway                                                      *
 *                                                                                      *
  ****************************************************************************************/

#include "DateTime.h"
#include "Date.h"
#include "CTime.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

//enum LANDSURFACE { OPEN, BOG, FOREST, ALPINE, HEATHER, ROCK, GLACIER };
//enum SOIL { OPEN_SOIL, PEAT, FOREST_SOIL, ALPINE_SOIL, HEATHER_SOIL, BEDROCK, GLACIER_BED };
enum LANDSURFACE { SURF0, SURF1, SURF2, SURF3, SURF4, SURF5, SURF6, SURF7, SURF8, SURF9,  
                   SURF10, SURF11, SURF12, SURF13, SURF14, SURF15, SURF16, SURF17, SURF18, SURF19, 
                   GLACIER };
enum SOIL { SOIL0, SOIL1, SOIL2, SOIL3, SOIL4, SOIL5, GLACIER_BED }; //NB! Navn i hbv_soil_parameters.dta må stemme overens med disse
#define _USE_MATH_DEFINES   //pi
#define ELEMENT(a,b) (((a)*nCols)+(b))
const int numberPotentialEvaporationValuesPerYear=12;      // Number of long-term mean potential evaporation values per year
const int maximumNumberLandClasses=3;        // Maximum number of land/soil classes in use for each computational element
//const int numberLandSurfaceClasses=7;        // All possible land surface types including glaciers, excluding lakes
const int numberLandSurfaceClasses=21;       // All possible land surface types including glaciers, excluding lakes
const int numberSoilClasses=7;              // All possible soil/subsurface types including glaciers, excluding lakes
const int maximumCorrectionCatchments=1000;
const int numberInputSeries=7;
const int numberSnowClasses=9;
const int minimumTimeStep=3600;
const int finalYear=2100;
const int finalMonth=12;
const int finalDay=31;
const int finalHour=23;
const int finalMinute=59;
const double numberSecondsDay=86400.0;
const double missingData=-9999.0;
const double probNorm[9] = {0.01,0.04,0.1,0.2,0.3,0.2,0.1,0.04,0.01};
const double epsilon=1.0e-4;

class DistributedHbv;


//class ParametersGeneral
class ParametersGeneral
{
 public:
  ParametersGeneral();
  ~ParametersGeneral();
  void SetSECONDS_TIMESTEP(int value) { SECONDS_TIMESTEP = value; }
  int GetSECONDS_TIMESTEP() const { return SECONDS_TIMESTEP; }
  void SetNUM_PREC_SERIES(int value) { NUM_PREC_SERIES = value; }
  int GetNUM_PREC_SERIES() const { return NUM_PREC_SERIES; }
  void SetNUM_TEMP_SERIES(int value) { NUM_TEMP_SERIES = value; }
  int GetNUM_TEMP_SERIES() const { return NUM_TEMP_SERIES; }
  void SetPREC_GRAD_LOW(double value) { PREC_GRAD_LOW = value; }
  double GetPREC_GRAD_LOW() const { return PREC_GRAD_LOW; }
  void SetPREC_GRAD_HIGH(double value) { PREC_GRAD_HIGH = value; }
  double GetPREC_GRAD_HIGH() const { return PREC_GRAD_HIGH; }
  void SetGRAD_CHANGE_ALT(double value) { GRAD_CHANGE_ALT = value; }
  double GetGRAD_CHANGE_ALT() const { return GRAD_CHANGE_ALT; }
  void SetPREC_CORR_RAIN(double value) { PREC_CORR_RAIN = value; }
  double GetPREC_CORR_RAIN() const { return PREC_CORR_RAIN; }
  void SetPREC_CORR_SNOW(double value) { PREC_CORR_SNOW = value; }
  double GetPREC_CORR_SNOW() const { return PREC_CORR_SNOW; }
  void SetLAPSE_DRY(double value) { LAPSE_DRY = value; }
  double GetLAPSE_DRY() const { return LAPSE_DRY; }
  void SetLAPSE_WET(double value) { LAPSE_WET = value; }
  double GetLAPSE_WET() const { return LAPSE_WET; }
  void SetDAY_TEMP_MEMORY(double value) { DAY_TEMP_MEMORY = value; }
  double GetDAY_TEMP_MEMORY() const { return DAY_TEMP_MEMORY; }
  void SetLAKE_EPOT_PAR(double value) { LAKE_EPOT_PAR = value; }
  double GetLAKE_EPOT_PAR() const { return LAKE_EPOT_PAR; }
  void SetKLAKE(double value) { KLAKE = value; }
  double GetKLAKE() const { return KLAKE; }
  void SetDELTA_LEVEL(double value) { DELTA_LEVEL = value; }
  double GetDELTA_LEVEL() const { return DELTA_LEVEL; }
  void SetNLAKE(double value) { NLAKE = value; }
  double GetNLAKE() const { return NLAKE; }
  void SetINITIAL_SOIL_MOISTURE(double value) { INITIAL_SOIL_MOISTURE = value; }
  double GetINITIAL_SOIL_MOISTURE() const { return INITIAL_SOIL_MOISTURE; }
  void SetINITIAL_UPPER_ZONE(double value) { INITIAL_UPPER_ZONE = value; }
  double GetINITIAL_UPPER_ZONE() const { return INITIAL_UPPER_ZONE; }
  void SetINITIAL_LOWER_ZONE(double value) { INITIAL_LOWER_ZONE = value; }
  double GetINITIAL_LOWER_ZONE() const { return INITIAL_LOWER_ZONE; }
  void SetINITIAL_LAKE_TEMP(double value) { INITIAL_LAKE_TEMP = value; }
  double GetINITIAL_LAKE_TEMP() const { return INITIAL_LAKE_TEMP; }
  void SetINITIAL_LAKE_LEVEL(double value) { INITIAL_LAKE_LEVEL = value; }
  double GetINITIAL_LAKE_LEVEL() const { return INITIAL_LAKE_LEVEL; }
  void SetINITIAL_SNOW(double value) { INITIAL_SNOW = value; }
  double GetINITIAL_SNOW() const { return INITIAL_SNOW; }
  void SetDAY_SNOW_ZERO(int value) { DAY_SNOW_ZERO = value; }
  int GetDAY_SNOW_ZERO() const { return DAY_SNOW_ZERO; }
  void SetHeight_WIND_INSTRUMENT(double  value) { Height_WIND_INSTRUMENT = value; }
  double  GetHeight_WIND_INSTRUMENT() const { return Height_WIND_INSTRUMENT; }
  void SetHeight_HUMI_INSTRUMENT(double  value) { Height_HUMI_INSTRUMENT = value; }
  double  GetHeight_HUMI_INSTRUMENT() const { return Height_HUMI_INSTRUMENT; }
  
 private:
  int SECONDS_TIMESTEP;
  int NUM_PREC_SERIES;
  int NUM_TEMP_SERIES;
  int DAY_SNOW_ZERO;
  double PREC_GRAD_LOW;
  double PREC_GRAD_HIGH;
  double GRAD_CHANGE_ALT;
  double PREC_CORR_RAIN;
  double PREC_CORR_SNOW;
  double LAPSE_DRY;
  double LAPSE_WET;
  double DAY_TEMP_MEMORY;
  double LAKE_EPOT_PAR;
  double KLAKE;
  double DELTA_LEVEL;
  double NLAKE;
  double INITIAL_SOIL_MOISTURE;
  double INITIAL_UPPER_ZONE;
  double INITIAL_LOWER_ZONE;
  double INITIAL_LAKE_TEMP;
  double INITIAL_LAKE_LEVEL;
  double INITIAL_SNOW;
  double   Height_WIND_INSTRUMENT;
  double  Height_HUMI_INSTRUMENT;
};


//class ParametersLandSurface
class ParametersLandSurface
{
 public:
  ParametersLandSurface();
  ~ParametersLandSurface();
  void SetINTER_MAX(double value) { INTER_MAX = value; }
  double GetINTER_MAX() const { return INTER_MAX; }
  void SetEPOT_PAR(double value) { EPOT_PAR = value; }
  double GetEPOT_PAR() const { return EPOT_PAR; }
  void SetWET_PER_CORR(double value) { WET_PER_CORR = value; }
  double GetWET_PER_CORR() const { return WET_PER_CORR; }
  void SetACC_TEMP(double value) { ACC_TEMP = value; }
  double GetACC_TEMP() const { return ACC_TEMP; }
  void SetMELT_TEMP(double value) { MELT_TEMP = value; }
  double GetMELT_TEMP() const { return MELT_TEMP; }
  void SetSNOW_MELT_RATE(double value) { SNOW_MELT_RATE = value; }
  double GetSNOW_MELT_RATE() const { return SNOW_MELT_RATE; }
  void SetICE_MELT_RATE(double value) { ICE_MELT_RATE = value; }
  double GetICE_MELT_RATE() const { return ICE_MELT_RATE; }
  void SetFREEZE_EFF(double value) { FREEZE_EFF = value; }
  double GetFREEZE_EFF() const { return FREEZE_EFF; }
  void SetMAX_REL(double value) { MAX_REL = value; }
  double GetMAX_REL() const { return MAX_REL; }
  void SetALBEDO(double value) { ALBEDO = value; }
  double GetALBEDO() const { return ALBEDO; }
  void SetCV_SNOW(double value) { CV_SNOW = value; }
  double GetCV_SNOW() const { return CV_SNOW; }
  void SetSNOW_WEIGHT(int k, double value) { SNOW_WEIGHT[k]= value; }
  double GetSNOW_WEIGHT(int k) const { return SNOW_WEIGHT[k]; }
  void SetTREE_HEIGHT(double value) { TREE_HEIGHT = value; }
  double GetTREE_HEIGHT() const { return TREE_HEIGHT; }
  void SetTREE_LAI(double value) { TREE_LAI = value; }
  double GetTREE_LAI() const { return TREE_LAI; }
  void SetTREE_LAI_CORR(double value) { TREE_LAI_CORR = value; }
  double GetTREE_LAI_CORR() const { return TREE_LAI_CORR; }
  void SetBULK_RESISTANCE(double value) { BULK_RESISTANCE = value; }
  double GetBULK_RESISTANCE() const { return BULK_RESISTANCE; }
  void SetALBEDO_SNOW(double value) { ALBEDO_SNOW = value; }
  double GetALBEDO_SNOW() const { return ALBEDO_SNOW; }
  void SetDECIDUOUS_SHARE(double value) { DECIDUOUS_SHARE = value; }
  double GetDECIDUOUS_SHARE() const { return DECIDUOUS_SHARE; }
  void SetWind_H(double value) { WIND_H = value; }
  double GetWind_H() const { return WIND_H; }
  void SetTopen_min(double value) { Topen_min = value; }
  double GetTopen_min() const { return Topen_min; }
  void SetTclose_min(double value) { Tclose_min = value; }
  double GetTclose_min() const { return Tclose_min; }
  void SetVPDclose(double value) { VPDclose = value; }
  double GetVPDclose() const { return VPDclose; }
  void SetVPDopen(double value) { VPDopen = value; }
  double GetVPDopen() const { return VPDopen; }
  void SetGH(double value) { gh = value; }
  double GetGH() const { return gh; }
  void SetCL(double value) { cl = value; }
  double GetCL() const { return cl; }
  void SetZ0G(double value) { z0g = value; }
  double GetZ0G() const { return z0g; }
  void SetGSMAX(double value) { gsmax = value; }
  double GetGSMAX() const { return gsmax; }
  void SetCR(double value) { CR = value; }
  double GetCR() const { return CR; }
  void SetD50(double value) { D50 = value; }
  double GetD50() const { return D50; }
  void SetQ50(double value) { Q50 = value; }
  double GetQ50() const { return Q50; }
  
 private:
  double INTER_MAX;
  double EPOT_PAR;
  double WET_PER_CORR;
  double ACC_TEMP;
  double MELT_TEMP;  
  double SNOW_MELT_RATE;
  double ICE_MELT_RATE;
  double FREEZE_EFF;
  double MAX_REL;
  double ALBEDO;    
  double CV_SNOW;
  double SNOW_WEIGHT[numberSnowClasses];
  double TREE_HEIGHT;
  double TREE_LAI;
  double TREE_LAI_CORR;
  double BULK_RESISTANCE;
  double ALBEDO_SNOW;
  double DECIDUOUS_SHARE;
  double WIND_H;
  double Topen_min;
  double Tclose_min;
  double VPDclose;
  double VPDopen;
  double gh;
  double cl;
  double z0g;
  double gsmax;
  double CR;
  double D50;
  double Q50;
};


//class ParametersSubSurfaceHbv
class ParametersSubSurfaceHbv
{
 public:
  ParametersSubSurfaceHbv();
  ~ParametersSubSurfaceHbv();
  void SetFC(double value) { FC = value; }
  double GetFC() const { return FC; }
  void SetFCDEL(double value) { FCDEL = value; }
  double GetFCDEL() const { return FCDEL; }
  void SetBETA(double value) { BETA = value; }
  double GetBETA() const { return BETA; }
  void SetINFMAX(double value) { INFMAX = value; }
  double GetINFMAX() const { return INFMAX; }
  void SetKUZ(double value) { KUZ = value; }
  double GetKUZ() const { return KUZ; }
  void SetALFA(double value) { ALFA = value; }
  double GetALFA() const { return ALFA; }
  void SetPERC(double value) { PERC = value; }
  double GetPERC() const { return PERC; }
  void SetKLZ(double value) { KLZ = value; }
  double GetKLZ() const { return KLZ; }
  void SetDRAW(double value) { DRAW = value; }
  double GetDRAW() const { return DRAW; }
  
 private:
  double FC;
  double FCDEL;
  double BETA;  
  double INFMAX;         
  double KUZ;         
  double ALFA;  
  double PERC;
  double KLZ;               
  double DRAW;              
};


//class SelectedTimeSeriesElements
class SelectedTimeSeriesElements
{
 public:
  SelectedTimeSeriesElements();
  ~SelectedTimeSeriesElements();
  void SetNumberElements(int value) { numberElements = value; }
  int GetNumberElements() const { return numberElements; }
  void SetTimeSeriesElement(int index, int value) { TimeSeriesElements[index] = value; }
  int GetTimeSeriesElement(int index) const { return TimeSeriesElements[index]; }
  void SelectedTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout);

 private:
  int numberElements;
  int * TimeSeriesElements;
};


// class MeteorologicalStations
class MeteorologicalStations
{
 public:
  MeteorologicalStations();
  ~MeteorologicalStations();
  void SetNumPrecStations(int value) { numberPrecStations = value; }
  int GetNumPrecStations() { return numberPrecStations; }
  void SetNumTempStations(int value) { numberTempStations = value; }
  int GetNumTempStations() { return numberTempStations; }
  int GetStationNumber(int k)  { return stationNumber[k]; }
  double GetStationCoordX(int k)  { return stationCoordX[k]; }
  double GetStationCoordY(int k)  { return stationCoordY[k]; }
  double GetStationAltitude(int k)  { return stationAltitude[k]; }
  void SetMeteorologicalStations(ifstream &fileControl, ofstream &fout);

 private:
  int numberPrecStations;
  int numberTempStations;
  int * stationNumber;
  double * stationCoordX;
  double * stationCoordY;
  double * stationAltitude;
};


// class InputTimeSeries
class InputTimeSeries
{
 public:
  InputTimeSeries(int numberRows, int numberColums, DateTime firstTime, DateTime lastTime, int secondsPerTimeStep);
  ~InputTimeSeries();
  int GetNumberTimeSteps() const { return timeSteps; }
  int GetNumberInputSeries() const { return numberSeries; }
  void SetGeneralPar(ParametersGeneral *parObj) { commonPar = parObj; }
  ParametersGeneral *GetGeneralPar() const { return commonPar; }
  DateTime GetDateTime(int i) const { return datetime[i]; }
  void SetInput(ifstream &fin);
  double GetInput(int i, int j) const { return inputArray[i*numberSeries+j]; }
  void WriteInput();

 private:
  ParametersGeneral *commonPar;
  int timeSteps;
  int numberSeries;
  DateTime * datetime;
  double * inputArray;
};


// class InputElement
class InputElement
{
 public:
  InputElement(int value);
  ~InputElement();
  int GetNumberValues() const { return numberValues; }
  void SetInput(int i, double value) { inputArray[i] = value; }
  double GetInput(int i) const { return inputArray[i]; }

 private:
  int numberValues;
  double * inputArray;
};


//class EvaporationControl
class EvaporationControl
{
public:
  void SetEvaporationModellingControl(char value);
  char GetEvaporationModellingControl() const;
  double GetEvaporationArray(int i) const { return evaporationArray[i]; }
  EvaporationControl();
  ~EvaporationControl();

private:
  char evaporationModellingControl;
  int numberValues;
  double * evaporationArray;
};


// class LakeWaterBalance
class LakeWaterBalance
{
 public:
  LakeWaterBalance();
  ~LakeWaterBalance();
  void WaterBalance(int timeStep, DateTime datetime, int dayofyear_PM2);
  void SetLakeValues(double temperature, double waterlevel);
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetLakeEvap() const;
  double GetRunoff() const;
  double GetWaterLevel() const;
  double GetWaterLevelChange() const;
  void SetGeneralPar(ParametersGeneral *parObj) { commonPar = parObj; }
  ParametersGeneral *GetGeneralPar() const { return commonPar; }
  //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
  //  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
  void SetInputElement(InputElement *inElementObj) { inElement = inElementObj; }
  InputElement *GetInputElement() const { return inElement; }
  void SetLandScapeElement(DistributedHbv *theElement) { landScapeElement = theElement; }
  DistributedHbv *GetLandScapeElement() const { return landScapeElement; }

 private:
  ParametersGeneral *commonPar;
  //  InputTimeSeries *inTimeSeries;
  InputElement *inElement;
  DistributedHbv *landScapeElement;
  double precipitation;                       /*  Precipitation (m/timestep)  */
  double temp;                                /*  Air temperature (deg. C)   */
  double lakeTemp;                            /*  Lake temperature (deg. C)  */
  double lakeEvaporation;                     /*  Evaporation (m/timestep)  */
  double waterLevel;                          /*  Lake water level (m)  */
  double waterLevelChange;                    /*  Lake water level change (m)  */
  double runoff;                              /*  Runoff (m/timestep)  */
  double discharge;                           /*  Lake outflow (m3/s)  */
  double tempMax;                          // maximum temperature (deg. C)
  double tempMin;                          // minimum temperature (deg. C)
  double wind1;                            // wind speed (m/s)
  double radiationS;                      // solar radiation (MJ/m2/day)
  double vp;                     // actual vapor pressure (Pa)
};


// class Vegetation
class Vegetation
{
 public:
  Vegetation();
  ~Vegetation();
  void WaterBalance(int timeStep, DateTime datetime, int dayofyear_PM2, double snowCoverFraction);
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetInterceptionStore() const;
  double GetInterceptionLoss() const;
  double GetThroughFall() const;
  void SetGeneralPar(ParametersGeneral *parObj) { commonPar = parObj; }
  ParametersGeneral *GetGeneralPar() const { return commonPar; }
  void SetLandSurfacePar(ParametersLandSurface *parObj) { landSurfacePar = parObj; }
  ParametersLandSurface *GetLandSurfacePar() { return landSurfacePar; }
  //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
  //  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
  void SetInputElement(InputElement *inElementObj) { inElement = inElementObj; }
  InputElement *GetInputElement() const { return inElement; }
  void SetLandScapeElement(DistributedHbv *theElement) { landScapeElement = theElement; }
  DistributedHbv *GetLandScapeElement() const { return landScapeElement; }
  double GetDryPeriod() { return dryPeriod; }
  double GetPET_PM() { return potev; }

 private:
  ParametersGeneral *commonPar;
  ParametersLandSurface *landSurfacePar;
  //  InputTimeSeries *inTimeSeries;
  InputElement *inElement;
  DistributedHbv *landScapeElement;
  double precipitation;                    /*  Precipitation (m/timestep)  */
  double temp;                             /*  Air temperature (deg. C)  */
  double tempMax;                          // maximum temperature (deg. C)
  double tempMin;                          // minimum temperature (deg. C)
  double wind1;                            // wind speed (m/s)
  double radiationS;                      // solar radiation (MJ/m2/day)
  double vp;                              // actual vapor pressure (Pa)
  double potev;                            /*  Potential evapotranspiration (m/timestep)  */
  double prevInterception;                 /*  Interception from previous time step (m)  */
  double interceptionStore;                /*  Interception store (m)  */
  double interceptionLoss;                 /*  Interception loss from vegetation (m)  */
  double throughFall;                      /*  Throughfall (m/timestep)  */
  double wetPeriod, dryPeriod;             /*  Length of wet and dry period during evapotranspiration (fraction of timestep)  */
  int    budburst;                          /* account the days after budburst*/
  double sft;                              /*Growing Degree Day (GDD) sum*/
  double calc_lai;                          /* calcualted LAI */
  double g;                                /* heat unit index*/
  double dm ;                              /* biomass */
  double huharv ;                          /* head unit for harvest*/
  double olai ;                            /* = calc_lai */
};


// class Snow
class Snow
{
 public:
  Snow();
  ~Snow();
  void WaterBalance(int timeStep, DateTime datetime, double waterInput);
  void SetSnowStore(double value);
  double GetSnowCoverFraction() const;
  double GetSnowStore() const;
  double GetSnowWaterEquivalentChange() const;
  double GetMeltWater() const;
  double GetWaterOutput() const;
  double GetDistSnowStore(int index) const { return distSnowStore[index]; }
  double GetDistMeltWater(int index) const { return distMeltWater[index]; }
  void SetSnowValues(double snowstore, double meltwater);
  void SeDistSnowStore(int index, double value) { distSnowStore[index] = value; }
  void SeDistMeltWater(int index, double value) { distMeltWater[index] = value; }
  void SetGeneralPar(ParametersGeneral *parObj) { commonPar = parObj; }
  ParametersGeneral *GetGeneralPar() const { return commonPar; }
  void SetLandSurfacePar(ParametersLandSurface *parObj) { landSurfacePar = parObj; }
  ParametersLandSurface *GetLandSurfacePar() { return landSurfacePar; }
  //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
  //  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
  void SetInputElement(InputElement *inElementObj) { inElement = inElementObj; }
  InputElement *GetInputElement() const { return inElement; }
  void SetLandScapeElement(DistributedHbv *theElement) { landScapeElement = theElement; }
  DistributedHbv *GetLandScapeElement() const { return landScapeElement; }

 private:
  ParametersGeneral *commonPar;
  ParametersLandSurface *landSurfacePar;
  //  InputTimeSeries *inTimeSeries;
  InputElement *inElement;
  DistributedHbv *landScapeElement;
  double temp;                               /*  Air temperature (deg. C)  */
  double snowStore;                          /*  Snow store (m)  */
  double meltWater;                          /*  Meltwater in snow (m)  */
  double waterOutput;                        /*  Output of meltwater from snow store (m/timestep)  */
  double snowWaterEquivalentChange;          /*  Change of snow water equivalent (m)  */
  double snowCoverFraction;                  /*  Fraction of area covered by snow  */
  double distSnowStore[numberSnowClasses];   /*  Distributed snow store (m)  */
  double distMeltWater[numberSnowClasses];   /*  Distributed meltwater in snow (m)  */
  double distWaterOutput[numberSnowClasses]; /*  Output of meltwater from snow store (m/timestep)  */
};


// class GlacierSurface
class GlacierSurface
{
 public:
  GlacierSurface();
  ~GlacierSurface();
  void WaterBalance(int timeStep, DateTime datetime);
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetIceMelt() const;
  void SetGeneralPar(ParametersGeneral *parObj) { commonPar = parObj; }
  ParametersGeneral *GetGeneralPar() const { return commonPar; }
  void SetLandSurfacePar(ParametersLandSurface *parObj) { landSurfacePar = parObj; }
  ParametersLandSurface *GetLandSurfacePar() { return landSurfacePar; }
  //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
  //  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
  void SetInputElement(InputElement *inElementObj) { inElement = inElementObj; }
  InputElement *GetInputElement() const { return inElement; }
  void SetLandScapeElement(DistributedHbv *theElement) { landScapeElement = theElement; }
  DistributedHbv *GetLandScapeElement() const { return landScapeElement; }

 private:
  ParametersGeneral *commonPar;
  ParametersLandSurface *landSurfacePar;
  //  InputTimeSeries *inTimeSeries;
  InputElement *inElement;
  DistributedHbv *landScapeElement;
  double precipitation;                    /*  Precipitation (m/timestep)  */
  double temp;                             /*  Air temperature (deg. C)  */
  double iceMelt;                          /*  Meltwater from ice (m/timestep)  */
};


// class HBV
class HBV
{
 public:
  HBV();
  ~HBV();
  void SetInitialHbvValues();
  void SetSubSurfaceHbvStore(double sm, double uz, double lz);
  void WaterBalance(int timeStep, DateTime datetime, double waterInput, double snowCoverFraction, double dryPeriod, double potev);
  double GetSoilMoisture() const;
  double GetSoilMoistureDeficit() const;
  double GetPercSoilUpper() const;
  double GetUpperZone() const;
  double GetLowerZone() const;
  double GetTranspSoilEvap() const;
  double GetRunoff() const;
  void SetGeneralPar(ParametersGeneral *parObj) { commonPar = parObj; }
  ParametersGeneral *GetGeneralPar() const { return commonPar; }
  void SetLandSurfacePar(ParametersLandSurface *parObj) { landSurfacePar = parObj; }
  ParametersLandSurface *GetLandSurfacePar() { return landSurfacePar; }
  void SetSubSurfaceHbvPar(ParametersSubSurfaceHbv *parObj) { subSurfacePar = parObj; }
  ParametersSubSurfaceHbv *GetSubSurfaceHbvPar() const { return subSurfacePar; }
  //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
  //  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
  void SetInputElement(InputElement *inElementObj) { inElement = inElementObj; }
  InputElement *GetInputElement() const { return inElement; }
  void SetLandScapeElement(DistributedHbv *theElement) { landScapeElement = theElement; }
  DistributedHbv *GetLandScapeElement() const { return landScapeElement; }

 private:
  ParametersGeneral *commonPar;
  ParametersSubSurfaceHbv *subSurfacePar;
  ParametersLandSurface *landSurfacePar;
  //  InputTimeSeries *inTimeSeries;
  InputElement *inElement;
  DistributedHbv *landScapeElement;
  double temp;                                /*  Temperature (deg. C)  */
  double soilMoisture;                        /*  Soil moisture content (m)  */ 
  double percSoilUpper;                       /*  Percolation from soil moisture zone to upper zone (m/timestep)  */ 
  double upperZone;                           /*  Upper groundwater zone water content (m)  */
  double lowerZone;                           /*  Lower groundwater zone water content (m)  */
  double lowerRunoff;                         /*  Runoff from lower layer (m/timestep)  */
  double upperRunoff;                         /*  Runoff from upper layer (m/timestep)  */
  double transpSoilEvap;                      /*  Water lost from subsurface by evapotranspiration (m)  */
  double runoff;                              /*  Runoff (m/timestep)  */
};


// class Lake
class Lake
{
 public:
  Lake();
  ~Lake();
  void SetLakeWaterBalance(LakeWaterBalance *theLakeWaterBalance) { ptrLakeWaterBalance = theLakeWaterBalance; }
  LakeWaterBalance *GetLakeWaterBalance() const { return ptrLakeWaterBalance; }
  void SetAreaFraction(double value) { areaFraction = value; }
  double GetAreaFraction() const { return areaFraction; }
  void WaterBalance(int timeStep, DateTime datetime, int dayofyear_PM2);
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetLakeEvap() const;
  double GetRunoff() const;
  double GetLakeStorage() const;
  double GetLakeStorageChange() const;

 private:
  LakeWaterBalance *ptrLakeWaterBalance;
  double areaFraction;
};


// class Glacier
class Glacier
{
 public:
  Glacier();
  ~Glacier();
  void SetGlacierSurface(GlacierSurface *theGlacierSurface) { ptrSurface = theGlacierSurface; }
  GlacierSurface *GetGlacierSurface() const { return ptrSurface; }
  void SetSnow(Snow *theSnow) { ptrSnow = theSnow; }
  Snow *GetSnow() const { return ptrSnow; }
  void SetHBV(HBV *theHBV) { ptrHbv = theHBV; }
  HBV *GetHBV() const { return ptrHbv; }
  void SetAreaFraction(double value) { areaFraction = value; }
  double GetAreaFraction() const { return areaFraction; }
  void SetSnowStore(double value);
  void SetSubSurfaceHbvStore(double sm, double uz, double lz);
  void WaterBalance(int timeStep, DateTime datetime, int initialTimeSteps, int numberTimeSteps);
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetSnowCoverFraction() const;
  double GetSnowStore() const;
  double GetMeltWater() const;
  double GetWaterOutput() const;
  double GetGlacierMassBalance() const;
  double GetSurfaceMassBalance() const;
  void SetSurfaceMassBalance(double value);
  double GetGlacierIceMelt() const;
  double GetSoilMoisture() const;
  double GetSoilMoistureDeficit() const;
  double GetPercSoilUpper() const;
  double GetUpperZone() const;
  double GetLowerZone() const;
  double GetRunoff() const;

 private:
  GlacierSurface *ptrSurface;
  Snow *ptrSnow;
  HBV *ptrHbv;
  double areaFraction;
  double glacierMassBalance;
};


// class HbvAquifer
class HbvAquifer
{
 public:
  HbvAquifer();
  ~HbvAquifer();
  void SetNextHbvAquifer(HbvAquifer *theHbvAquifer) { nextHbvAquifer = theHbvAquifer; }
  HbvAquifer *GetNextHbvAquifer() const { return nextHbvAquifer; }
  void SetVegetation(Vegetation *theVegetation) { ptrVeg = theVegetation; }
  Vegetation *GetVegetation() const { return ptrVeg; }
  void SetSnow(Snow *theSnow) { ptrSnow = theSnow; }
  Snow *GetSnow() const { return ptrSnow; }
  void SetHBV(HBV *theHBV) { ptrHbv = theHBV; }
  HBV *GetHBV() const { return ptrHbv; }
  void SetAreaFraction(double value) { areaFraction = value; }
  double GetAreaFraction() const { return areaFraction; }
  void SetSnowStore(double value);
  void SetSubSurfaceHbvStore(double sm, double uz, double lz);
  void WaterBalance(int timeStep, DateTime datetime, int dayofyear_PM2);
  double GetTotalHbvAreaFraction() const;
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetSnowCoverFraction() const;
  double GetSnowStore() const;
  double GetMeltWater() const;
  double GetWaterOutput() const;
  double GetSoilMoisture() const;
  double GetSoilMoistureDeficit() const;
  double GetPercSoilUpper() const;
  double GetUpperZone() const;
  double GetLowerZone() const;
  double GetInterceptionLoss() const;
  double GetTranspSoilEvap() const;
  double GetRunoff() const;
  double GetPET() const;

 private:
  Vegetation *ptrVeg;
  Snow *ptrSnow;
  HBV *ptrHbv;
  HbvAquifer *nextHbvAquifer;
  double areaFraction;
};


// class SubCatchment
class SubCatchment
{
public:
  SubCatchment();
  ~SubCatchment(); 
  void SetSubCatchmentIndex(int value) { subCatchmentIndex = value; }
  int GetSubCatchmentIndex() const { return subCatchmentIndex; }
  void SetIdentifier(int value) { identifier = value; }
  int GetIdentifier() const { return identifier; }
  void SetNumUpStream(int value);
  int GetNumUpStream() const { return numUpStream; }
  void SetUpStream(int k, SubCatchment *theSubCatchment) { UpStream[k] = theSubCatchment; }
  SubCatchment *GetUpStream(int k) const { return UpStream[k]; }
  void SetNumLandScape(int value) { numLandScape = value; }
  int GetNumLandScape() const { return numLandScape; }
  void SetLandScapeElement(DistributedHbv *theElement) { landScapeElement = theElement; }
  DistributedHbv *GetLandScapeElement() const { return landScapeElement; }
  void SetCorrection(double value) { correction = value; }
  double GetCorrection() const { return correction; }
  void AllocateAccumulatedDischarge(int numberTimeSteps); 
  void SetAccumulatedDischarge(int index, double value) { accumulatedDischarge[index] = value; }
  double GetAccumulatedDischarge(int index) const { return accumulatedDischarge[index]; }
  void AllocateAccumulatedWaterBalance(int numberTimeSteps); 
  void AllocateWaterBalance(int numberTimeSteps); 
  void SetAccumulatedPrecipitation(int index, double value) { accumulatedPrecipitation[index] = value; }
  double GetAccumulatedPrecipitation(int index) const { return accumulatedPrecipitation[index]; }
  void SetAccumulatedTemperature(int index, double value) { accumulatedTemperature[index] = value; }
  double GetAccumulatedTemperature(int index) const { return accumulatedTemperature[index]; }
  void SetAccumulatedLakeStorage(int index, double value) { accumulatedLakeStorage[index] = value; }
  double GetAccumulatedLakeStorage(int index) const { return accumulatedLakeStorage[index]; }
  void SetAccumulatedSnowStore(int index, double value) { accumulatedSnowStore[index] = value; }
  double GetAccumulatedSnowStore(int index) const { return accumulatedSnowStore[index]; }
  void SetAccumulatedMeltWater(int index, double value) { accumulatedMeltWater[index] = value; }
  double GetAccumulatedMeltWater(int index) const { return accumulatedMeltWater[index]; }
  void SetAccumulatedWaterOutput(int index, double value) { accumulatedWaterOutput[index] = value; }
  double GetAccumulatedWaterOutput(int index) const { return accumulatedWaterOutput[index]; }
  void SetAccumulatedSnowCoverFraction(int index, double value) { accumulatedSnowCoverFraction[index] = value; }
  double GetAccumulatedSnowCoverFraction(int index) const { return accumulatedSnowCoverFraction[index]; }
  void SetAccumulatedGlacierMassBalance(int index, double value) { accumulatedGlacierMassBalance[index] = value; }
  double GetAccumulatedGlacierMassBalance(int index) const { return accumulatedGlacierMassBalance[index]; }
  void SetAccumulatedGlacierIceMelt(int index, double value) { accumulatedGlacierIceMelt[index] = value; }
  double GetAccumulatedGlacierIceMelt(int index) const { return accumulatedGlacierIceMelt[index]; }
  void SetAccumulatedEvapotranspiration(int index, double value) { accumulatedEvapotranspiration[index] = value; }
  double GetAccumulatedEvapotranspiration(int index) const { return accumulatedEvapotranspiration[index]; }
  void SetAccumulatedRunoff(int index, double value) { accumulatedRunoff[index] = value; }
  double GetAccumulatedRunoff(int index) const { return accumulatedRunoff[index]; }
  void SetAccumulatedHbvSoilMoisture(int index, double value) { accumulatedHbvSoilMoisture[index] = value; }
  double GetAccumulatedHbvSoilMoisture(int index) const { return accumulatedHbvSoilMoisture[index]; }
  void SetAccumulatedHbvSoilMoistureDeficit(int index, double value) { accumulatedHbvSoilMoistureDeficit[index] = value; }
  double GetAccumulatedHbvSoilMoistureDeficit(int index) const { return accumulatedHbvSoilMoistureDeficit[index]; }
  void SetAccumulatedHbvPercSoilUpper(int index, double value) { accumulatedHbvPercSoilUpper[index] = value; }
  double GetAccumulatedHbvPercSoilUpper(int index) const { return accumulatedHbvPercSoilUpper[index]; }
  void SetAccumulatedHbvUpperZone(int index, double value) { accumulatedHbvUpperZone[index] = value; }
  double GetAccumulatedHbvUpperZone(int index) const { return accumulatedHbvUpperZone[index]; }
  void SetAccumulatedHbvLowerZone(int index, double value) { accumulatedHbvLowerZone[index] = value; }
  double GetAccumulatedHbvLowerZone(int index) const { return accumulatedHbvLowerZone[index]; }
  void SetAccumulatedSum(int index, double value) { accumulatedSum[index] = value; }
  double GetAccumulatedSum(int index) const { return accumulatedSum[index]; }
  void SetAccumulatedSumLake(int index, double value) { accumulatedSumLake[index] = value; }
  double GetAccumulatedSumLake(int index) const { return accumulatedSumLake[index]; }
  void SetAccumulatedSumSnow(int index, double value) { accumulatedSumSnow[index] = value; }
  double GetAccumulatedSumSnow(int index) const { return accumulatedSumSnow[index]; }
  void SetAccumulatedSumGlacier(int index, double value) { accumulatedSumGlacier[index] = value; }
  double GetAccumulatedSumGlacier(int index) const { return accumulatedSumGlacier[index]; }
  void SetAccumulatedSumHbv(int index, double value) { accumulatedSumHbv[index] = value; }
  double GetAccumulatedSumHbv(int index) const { return accumulatedSumHbv[index]; }
  void SetSubCatchmentPrecipitation(int index, double value) { subCatchmentPrecipitation[index] = value; }
  double GetSubCatchmentPrecipitation(int index) const { return subCatchmentPrecipitation[index]; }
  void SetSubCatchmentTemperature(int index, double value) { subCatchmentTemperature[index] = value; }
  double GetSubCatchmentTemperature(int index) const { return subCatchmentTemperature[index]; }
  void SetSubCatchmentLakeStorage(int index, double value) { subCatchmentLakeStorage[index] = value; }
  double GetSubCatchmentLakeStorage(int index) const { return subCatchmentLakeStorage[index]; }
  void SetSubCatchmentSnowStore(int index, double value) { subCatchmentSnowStore[index] = value; }
  double GetSubCatchmentSnowStore(int index) const { return subCatchmentSnowStore[index]; }
  void SetSubCatchmentMeltWater(int index, double value) { subCatchmentMeltWater[index] = value; }
  double GetSubCatchmentMeltWater(int index) const { return subCatchmentMeltWater[index]; }
  void SetSubCatchmentWaterOutput(int index, double value) { subCatchmentWaterOutput[index] = value; }
  double GetSubCatchmentWaterOutput(int index) const { return subCatchmentWaterOutput[index]; }
  void SetSubCatchmentSnowCoverFraction(int index, double value) { subCatchmentSnowCoverFraction[index] = value; }
  double GetSubCatchmentSnowCoverFraction(int index) const { return subCatchmentSnowCoverFraction[index]; }
  void SetSubCatchmentGlacierMassBalance(int index, double value) { subCatchmentGlacierMassBalance[index] = value; }
  double GetSubCatchmentGlacierMassBalance(int index) const { return subCatchmentGlacierMassBalance[index]; }
  void SetSubCatchmentGlacierIceMelt(int index, double value) { subCatchmentGlacierIceMelt[index] = value; }
  double GetSubCatchmentGlacierIceMelt(int index) const { return subCatchmentGlacierIceMelt[index]; }
  void SetSubCatchmentEvapotranspiration(int index, double value) { subCatchmentEvapotranspiration[index] = value; }
  double GetSubCatchmentEvapotranspiration(int index) const { return subCatchmentEvapotranspiration[index]; }
  void SetSubCatchmentRunoff(int index, double value) { subCatchmentRunoff[index] = value; }
  double GetSubCatchmentRunoff(int index) const { return subCatchmentRunoff[index]; }
  void SetSubCatchmentHbvSoilMoisture(int index, double value) { subCatchmentHbvSoilMoisture[index] = value; }
  double GetSubCatchmentHbvSoilMoisture(int index) const { return subCatchmentHbvSoilMoisture[index]; }
  void SetSubCatchmentHbvSoilMoistureDeficit(int index, double value) 
    { subCatchmentHbvSoilMoistureDeficit[index] = value; }
  double GetSubCatchmentHbvSoilMoistureDeficit(int index) const 
    { return subCatchmentHbvSoilMoistureDeficit[index]; }
  void SetSubCatchmentHbvPercSoilUpper(int index, double value) { subCatchmentHbvPercSoilUpper[index] = value; }
  double GetSubCatchmentHbvPercSoilUpper(int index) const { return subCatchmentHbvPercSoilUpper[index]; }
  void SetSubCatchmentHbvUpperZone(int index, double value) { subCatchmentHbvUpperZone[index] = value; }
  double GetSubCatchmentHbvUpperZone(int index) const { return subCatchmentHbvUpperZone[index]; }
  void SetSubCatchmentHbvLowerZone(int index, double value) { subCatchmentHbvLowerZone[index] = value; }
  double GetSubCatchmentHbvLowerZone(int index) const { return subCatchmentHbvLowerZone[index]; }
  void ObsDataInput(char * fileObsStreamflow, DateTime firstTime, DateTime lastTime, int numberTimeSteps, int secondsPerTimeStep);
  double GetObsData(int index) const { return observedData[index]; }

private:
  int subCatchmentIndex;
  int identifier;
  int numUpStream;
  int numLandScape;
  double correction;
  double *accumulatedDischarge;
  double *accumulatedPrecipitation;
  double *accumulatedTemperature;
  double *accumulatedLakeStorage;
  double *accumulatedSnowStore;
  double *accumulatedMeltWater;
  double *accumulatedWaterOutput;
  double *accumulatedSnowCoverFraction;
  double *accumulatedGlacierMassBalance;
  double *accumulatedGlacierIceMelt;
  double *accumulatedEvapotranspiration;
  double *accumulatedRunoff;
  double *accumulatedHbvSoilMoisture;
  double *accumulatedHbvSoilMoistureDeficit;
  double *accumulatedHbvPercSoilUpper;
  double *accumulatedHbvUpperZone;
  double *accumulatedHbvLowerZone;
  double *accumulatedSum;
  double *accumulatedSumLake;
  double *accumulatedSumSnow;
  double *accumulatedSumGlacier;
  double *accumulatedSumHbv;
  double *subCatchmentPrecipitation;
  double *subCatchmentTemperature;
  double *subCatchmentLakeStorage;           
  double *subCatchmentSnowStore;           
  double *subCatchmentMeltWater;           
  double *subCatchmentWaterOutput;           
  double *subCatchmentSnowCoverFraction;           
  double *subCatchmentGlacierMassBalance;           
  double *subCatchmentGlacierIceMelt;           
  double *subCatchmentEvapotranspiration;  
  double *subCatchmentRunoff;              
  double *subCatchmentHbvSoilMoisture;     
  double *subCatchmentHbvSoilMoistureDeficit;     
  double *subCatchmentHbvPercSoilUpper;     
  double *subCatchmentHbvUpperZone;        
  double *subCatchmentHbvLowerZone;  
  double *observedData;
  SubCatchment **UpStream;
  DistributedHbv *landScapeElement;
};


// class DistributedHbv
class DistributedHbv
{
 public:
  DistributedHbv();
  ~DistributedHbv(); 
  void SetGeoIndex(int value) { geoIndex = value; }
  int GetGeoIndex() const { return geoIndex; }
  void SetLandIndex(int value) { landIndex = value; }
  int GetLandIndex() const { return landIndex; }
  void SetArea(double value) { area = value; }
  double GetArea() const { return area; }
  void SetLatitude(double value) { latitude = value; }
  double GetLatitude() const { return latitude; }
  void SetElevation(double value) { elevation = value; }
  double GetElevation() const { return elevation; }
  void SetTimeFalling(double value) { timefalling = value; }
  double GetTimeFalling() const { return timefalling; }
  void SetSlopeAngle(double value) { slopeAngle = value; }
  double GetSlopeAngle() const { return slopeAngle; }
  void SetAspect(double value) { aspect = value; }
  double GetAspect() const { return aspect; }
  void SetPcorr(double value) { Pcorr = value; }
  double GetPcorr() const { return Pcorr; }
  void SetSelectedTimeSeriesElements(SelectedTimeSeriesElements *object) { selectedElements = object; }
  SelectedTimeSeriesElements *GetSelectedTimeSeriesElements() const { return selectedElements; }
  void SetGeneralPar(ParametersGeneral *parObj) { commonPar = parObj; }
  ParametersGeneral *GetGeneralPar() const { return commonPar; }
  void SetEvaporationControlObj(EvaporationControl *controlObj) { evaporationControlObject = controlObj; }
  EvaporationControl *GetEvaporationControlObj() const { return evaporationControlObject; }
  void SetPrecipitationCorrection(double value) { precipitationCorrection = value; }
  double GetPrecipitationCorrection() const { return precipitationCorrection; }
  void SetTemperatureCorrection(double value) { temperatureCorrection = value; }
  double GetTemperatureCorrection() const { return temperatureCorrection; }
  void SetNextElement(DistributedHbv *theElement) { nextElement = theElement; }
  DistributedHbv *GetNextElement() const { return nextElement; }
  void SetLake(Lake *theLake) { lake = theLake; }
  Lake *GetLake() const { return lake; }
  void SetGlacier(Glacier *theGlacier) { glacier = theGlacier; }
  Glacier *GetGlacier() const { return glacier; }
  void SetHbvAquifer(HbvAquifer *theHbvAquifer) { hbvAquifer = theHbvAquifer; }
  HbvAquifer *GetHbvAquifer() const { return hbvAquifer; }
  void SetPrecStationsWeightedElevation(double value) { precStationsWeightedElevation = value; }
  double GetPrecStationsWeightedElevation() { return precStationsWeightedElevation; }
  void SetTempStationsWeightedElevation(double value) { tempStationsWeightedElevation = value; }
  double GetTempStationsWeightedElevation() { return tempStationsWeightedElevation; }
  void AllocateMetSeries(int numberPrecSeries, int numberTempSeries);
  void SetMetSeriesNumber(int k, int number) { metSeriesNumber[k] = number; }
  void SetMetSeriesWeight(int k, double weight) { metSeriesWeight[k] = weight; }
  int GetMetSeriesNumber(int k) { return metSeriesNumber[k]; }
  double GetMetSeriesWeight(int k) { return metSeriesWeight[k]; }
  void AllocateWaterBalance(int numberTimeSteps); 
  void SetDistributedHbvPrecipitation(int index, double value) { distributedHbvPrecipitation[index] = value; }
  double GetDistributedHbvPrecipitation(int index) const { return distributedHbvPrecipitation[index]; }
  void SetDistributedHbvTemperature(int index, double value) { distributedHbvTemperature[index] = value; }
  double GetDistributedHbvTemperature(int index) const { return distributedHbvTemperature[index]; }
  void SetDistributedHbvLakeStorage(int index, double value) { distributedElementLakeStorage[index] = value; }
  double GetDistributedHbvLakeStorage(int index) const { return distributedElementLakeStorage[index]; }
  void SetDistributedHbvSnowStore(int index, double value) { distributedHbvSnowStore[index] = value; }
  double GetDistributedHbvSnowStore(int index) const { return distributedHbvSnowStore[index]; }
  void SetDistributedHbvSnowCoverFraction(int index, double value) { distributedHbvSnowCoverFraction[index] = value; }
  double GetDistributedHbvSnowCoverFraction(int index) const { return distributedHbvSnowCoverFraction[index]; }
  void SetDistributedHbvMeltWater(int index, double value) { distributedHbvMeltWater[index] = value; }
  double GetDistributedHbvMeltWater(int index) const { return distributedHbvMeltWater[index]; }
  void SetDistributedHbvWaterOutput(int index, double value) { distributedHbvWaterOutput[index] = value; }
  double GetDistributedHbvWaterOutput(int index) const { return distributedHbvWaterOutput[index]; }
  void SetDistributedHbvGlacierMassBalance(int index, double value) { distributedHbvGlacierMassBalance[index] = value; }
  double GetDistributedHbvGlacierMassBalance(int index) const { return distributedHbvGlacierMassBalance[index]; }
  void SetDistributedHbvGlacierIceMelt(int index, double value) { distributedHbvGlacierIceMelt[index] = value; }
  double GetDistributedHbvGlacierIceMelt(int index) const { return distributedHbvGlacierIceMelt[index]; }
  void SetDistributedHbvEvapotranspiration(int index, double value) { distributedHbvEvapotranspiration[index] = value; }
  double GetDistributedHbvEvapotranspiration(int index) const { return distributedHbvEvapotranspiration[index]; }
  void SetDistributedHbvRunoff(int index, double value) { distributedHbvRunoff[index] = value; }
  double GetDistributedHbvRunoff(int index) const { return distributedHbvRunoff[index]; }
  void SetDistributedHbvSoilMoisture(int index, double value) { distributedHbvSoilMoisture[index] = value; }
  double GetDistributedHbvSoilMoisture(int index) const { return distributedHbvSoilMoisture[index]; }
  void SetDistributedHbvSoilMoistureDeficit(int index, double value)  
    { distributedHbvSoilMoistureDeficit[index] = value; }
  double GetDistributedHbvSoilMoistureDeficit(int index) const { return distributedHbvSoilMoistureDeficit[index]; }
  void SetDistributedHbvPercSoilUpper(int index, double value) { distributedHbvPercSoilUpper[index] = value; }
  double GetDistributedHbvPercSoilUpper(int index) const { return distributedHbvPercSoilUpper[index]; }
  void SetDistributedHbvUpperZone(int index, double value) { distributedHbvUpperZone[index] = value; }
  double GetDistributedHbvUpperZone(int index) const { return distributedHbvUpperZone[index]; }
  void SetDistributedHbvLowerZone(int index, double value) { distributedHbvLowerZone[index] = value; }
  double GetDistributedHbvLowerZone(int index) const { return distributedHbvLowerZone[index]; }
  void SetAccumulatedDischarge(double localValue);
  double GetAccumulatedDischarge() const { return accumulatedDischarge; }
  void SetAccumulatedSum(double value) { accumulatedSum = value; }
  double GetAccumulatedSum() const { return accumulatedSum; }
  void SetAccumulatedSumLake(double value) { accumulatedSumLake = value; }
  double GetAccumulatedSumLake() const { return accumulatedSumLake; }
  void SetAccumulatedSumSnow(double value) { accumulatedSumSnow = value; }
  double GetAccumulatedSumSnow() const { return accumulatedSumSnow; }
  void SetAccumulatedSumGlacier(double value) { accumulatedSumGlacier = value; }
  double GetAccumulatedSumGlacier() const { return accumulatedSumGlacier; }
  void SetAccumulatedSumHbv(double value) { accumulatedSumHbv = value; }
  double GetAccumulatedSumHbv() const { return accumulatedSumHbv; }
  void SetAccumulatedPrecipitation(double value) { accumulatedPrecipitation = value; }
  double GetAccumulatedPrecipitation() const { return accumulatedPrecipitation; }
  void SetAccumulatedTemperature(double value) { accumulatedTemperature = value; }
  double GetAccumulatedTemperature() const { return accumulatedTemperature; }
  void SetAccumulatedLakeStorage(double value) { accumulatedLakeStorage = value; }
  double GetAccumulatedLakeStorage() const { return accumulatedLakeStorage; }
  void SetAccumulatedSnowStore(double value) { accumulatedSnowStore = value; }
  double GetAccumulatedSnowStore() const { return accumulatedSnowStore; }
  void SetAccumulatedMeltWater(double value) { accumulatedMeltWater = value; }
  double GetAccumulatedMeltWater() const { return accumulatedMeltWater; }
  void SetAccumulatedWaterOutput(double value) { accumulatedWaterOutput = value; }
  double GetAccumulatedWaterOutput() const { return accumulatedWaterOutput; }
  void SetAccumulatedSnowCoverFraction(double value) { accumulatedSnowCoverFraction = value; }
  double GetAccumulatedSnowCoverFraction() const { return accumulatedSnowCoverFraction; }
  void SetAccumulatedGlacierMassBalance(double value) { accumulatedGlacierMassBalance = value; }
  double GetAccumulatedGlacierMassBalance() const { return accumulatedGlacierMassBalance; }
  void SetAccumulatedGlacierIceMelt(double value) { accumulatedGlacierIceMelt = value; }
  double GetAccumulatedGlacierIceMelt() const { return accumulatedGlacierIceMelt; }
  void SetAccumulatedEvapotranspiration(double value) { accumulatedEvapotranspiration = value; }
  double GetAccumulatedEvapotranspiration() const { return accumulatedEvapotranspiration; }
  void SetAccumulatedRunoff(double value) { accumulatedRunoff = value; }
  double GetAccumulatedRunoff() const { return accumulatedRunoff; }
  void SetAccumulatedHbvSoilMoisture(double value) { accumulatedHbvSoilMoisture = value; }
  double GetAccumulatedHbvSoilMoisture() const { return accumulatedHbvSoilMoisture; }
  void SetAccumulatedHbvSoilMoistureDeficit(double value) { accumulatedHbvSoilMoistureDeficit = value; }
  double GetAccumulatedHbvSoilMoistureDeficit() const { return accumulatedHbvSoilMoistureDeficit; }
  void SetAccumulatedHbvPercSoilUpper(double value) { accumulatedHbvPercSoilUpper = value; }
  double GetAccumulatedHbvPercSoilUpper() const { return accumulatedHbvPercSoilUpper; }
  void SetAccumulatedHbvUpperZone(double value) { accumulatedHbvUpperZone = value; }
  double GetAccumulatedHbvUpperZone() const { return accumulatedHbvUpperZone; }
  void SetAccumulatedHbvLowerZone(double value) { accumulatedHbvLowerZone = value; }
  double GetAccumulatedHbvLowerZone() const { return accumulatedHbvLowerZone; }
  void SetSnowStore(double value);
  void SetSubSurfaceHbvStore(double sm, double uz, double lz);
  int GetNumberSum() const;
  void SetSumWaterBalance();
  void SetInitialStorage();
  double GetInitialStorage() const;
  void SetFinalStorage();
  double GetFinalStorage() const;
  double GetSumPrecipitation() const;
  double GetSumEvapotranspiration() const;
  double GetSumRunoff() const;
  double GetLandArea() const;
  double GetLakeArea() const;
  double GetGlacierArea() const;
  double GetHbvArea() const;
  void WaterBalance(int timeStep, DateTime datetime, int initialTimeSteps, int numberTimeSteps, int dayofyear_PM2) const;
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetLakeStorage() const;
  double GetLakeStorageChange() const;
  double GetSnowCoverFraction() const;
  double GetSnowStore() const;
  double GetMeltWater() const;
  double GetWaterOutput() const;
  double GetGlacierMassBalance() const;
  double GetGlacierIceMelt() const;
  double GetHbvSoilMoisture() const;
  double GetHbvSoilMoistureDeficit() const;
  double GetHbvPercSoilUpper() const;
  double GetHbvUpperZone() const;
  double GetHbvLowerZone() const;
  double GetEvapotranspiration() const;
  double GetInterceptionLoss() const;
  double GetTranspSoilEvap() const;
  double GetLakeEvap() const;
  double GetRunoff() const;
  double GetDischarge() const;
  double GetPET() const;

private:
  int geoIndex;
  int landIndex;
  int numberSum;
  double initialStorage;
  double finalStorage;
  double sumPrecipitation;
  double sumEvapotranspiration;
  double sumRunoff;
  double accumulatedSum;
  double accumulatedSumLake;
  double accumulatedSumSnow;
  double accumulatedSumGlacier;
  double accumulatedSumHbv;
  double area;
  double latitude;
  double elevation;
  double timefalling;
  double slopeAngle;
  double aspect;
  double Pcorr;
  double precipitationCorrection;
  double temperatureCorrection;
  double accumulatedDischarge;
  double accumulatedPrecipitation;
  double accumulatedTemperature;
  double accumulatedLakeStorage;
  double accumulatedSnowStore;
  double accumulatedMeltWater;
  double accumulatedWaterOutput;
  double accumulatedSnowCoverFraction;
  double accumulatedGlacierMassBalance;
  double accumulatedGlacierIceMelt;
  double accumulatedEvapotranspiration;
  double accumulatedRunoff;
  double accumulatedHbvSoilMoisture;
  double accumulatedHbvSoilMoistureDeficit;
  double accumulatedHbvPercSoilUpper;
  double accumulatedHbvUpperZone;
  double accumulatedHbvLowerZone;
  double precStationsWeightedElevation;
  double tempStationsWeightedElevation;
  int * metSeriesNumber;
  double * metSeriesWeight;
  double *distributedHbvPrecipitation;
  double *distributedHbvTemperature;
  double *distributedElementLakeStorage;           
  double *distributedHbvSnowStore;           
  double *distributedHbvSnowCoverFraction;           
  double *distributedHbvMeltWater;           
  double *distributedHbvWaterOutput;           
  double *distributedHbvGlacierMassBalance;           
  double *distributedHbvGlacierIceMelt;           
  double *distributedHbvEvapotranspiration;  
  double *distributedHbvRunoff;              
  double *distributedHbvSoilMoisture;     
  double *distributedHbvSoilMoistureDeficit;     
  double *distributedHbvPercSoilUpper;     
  double *distributedHbvUpperZone;        
  double *distributedHbvLowerZone;        
  DistributedHbv *nextElement;
  Lake *lake;
  Glacier *glacier;
  HbvAquifer *hbvAquifer;
  SelectedTimeSeriesElements *selectedElements;
  ParametersGeneral *commonPar;
  EvaporationControl * evaporationControlObject;
};

