#pragma once
/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments.          *
 *  Preprocessing.                                                                      *
 *                                                                                      *
 *  Author: Stein Beldring, Norway                                                      *
 *                                                                                      *
 ****************************************************************************************/

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
enum SOIL { SOIL0, SOIL1, SOIL2, SOIL3, SOIL4, SOIL5, SOIL6, SOIL7, SOIL8, SOIL9, 
            SOIL10, SOIL11, SOIL12, SOIL13, SOIL14, SOIL15, SOIL16, SOIL17, SOIL18, SOIL19, 
            GLACIER_BED };
#define ELEMENT(a,b) (((a)*nCols)+(b))
const int maximumNumberLandClasses=2;        // Maximum number of land/soil classes in use for each computational element
//const int numberLandSurfaceClasses=7;        // All possible land surface types including glaciers, excluding lakes
//const int numberSoilClasses=7;               // All possible soil/subsurface types including glaciers, excluding lakes
const int numberLandSurfaceClasses=21;       // All possible land surface types including glaciers, excluding lakes
const int numberSoilClasses=21;              // All possible soil/subsurface types including glaciers, excluding lakes
const int finalYear=2099;
const int finalMonth=12;
const int finalDay=31;
const int finalHour=23;
const int finalMinute=59;
const double missingData=-9999.0;

class ParametersGeneral;
class MeteorologicalStations;
class DistributedHbv;


//class SubCatchment
class SubCatchment
{
public:
  SubCatchment();
  ~SubCatchment(); 
  void SetSubCatchmentIndex(int value) { subCatchmentIndex = value; }
  int GetSubCatchmentIndex() const { return subCatchmentIndex; }
  void SetIdentifier(int value) { identifier = value; }
  int GetIdentifier() const { return identifier; }
  void SetNumLandScape(int value) { numLandScape = value; }
  int GetNumLandScape() const { return numLandScape; }
  void SetCorrection(double value) { correction = value; }
  double GetCorrection() const { return correction; }
  void SetNumUpStream(int value);
  int GetNumUpStream() const { return numUpStream; }
  void SetUpStream(int k, SubCatchment *theSubCatchment) { UpStream[k] = theSubCatchment; }
  SubCatchment *GetUpStream(int k) const { return UpStream[k]; }
  void SetLandScapeElement(DistributedHbv *theElement) { landScapeElement = theElement; }
  DistributedHbv *GetLandScapeElement() const { return landScapeElement; }

private:
  int subCatchmentIndex, identifier;
  int numLandScape;
  int numUpStream;
  double correction;
  SubCatchment **UpStream;
  DistributedHbv *landScapeElement;
};


//class DistributedHbv
class DistributedHbv
{
 public:
  DistributedHbv();
  ~DistributedHbv(); 
  void SetGeoIndex(int value) { geoIndex = value; }
  int GetGeoIndex() const { return geoIndex; }
  void SetLandIndex(int value) { landIndex = value; }
  int GetLandIndex() const { return landIndex; }
  void SetSubCatchmentValue(int value) { subCatchmentValue = value; }
  int GetSubCatchmentValue() const { return subCatchmentValue; }
  void SetArea(double value) { area = value; }
  double GetArea() const { return area; }
  void SetElevation(double value) { elevation = value; }
  double GetElevation() const { return elevation; }
  void SetSlopeAngle(double value) { slopeAngle = value; }
  double GetSlopeAngle() const { return slopeAngle; }
  void SetAspect(double value) { aspect = value; }
  double GetAspect() const { return aspect; }
  void SetLakePercent(double value) { lakePercent = value; }
  double GetLakePercent() const { return lakePercent; }
  void SetForestPercent(double value) { forestPercent = value; }
  double GetForestPercent() const { return forestPercent; }
  void SetBogPercent(double value) { bogPercent = value; }
  double GetBogPercent() const { return bogPercent; }
  void SetGlacierPercent(double value) { glacierPercent = value; }
  double GetGlacierPercent() const { return glacierPercent; }
  void SetSoil(int value) { Soiltype = value; }
  int GetSoil() const { return Soiltype; }
  void SetOpenLandPercent(double value) { openLandPercent = value; }
  double GetOpenLandPercent() const { return openLandPercent; }
  void SetAlpinePercent(double value) { alpinePercent = value; }
  double GetAlpinePercent() const { return alpinePercent; }
  void SetHeatherPercent(double value) { heatherPercent = value; }
  double GetHeatherPercent() const { return heatherPercent; }
  void SetBedrockPercent(double value) { bedrockPercent = value; }
  double GetBedrockPercent() const { return bedrockPercent; }
  void SetTreeLevel(double value) { treeLevel = value; }
  double GetTreeLevel() const { return treeLevel; }
  void SetXCoord(double value) { xCoord = value; }
  double GetXCoord() const { return xCoord; }
  void SetYCoord(double value) { yCoord = value; }
  double GetYCoord() const { return yCoord; }
  void SetLandSurfacePercent(int index, double value) { landSurfacePercent[index] = value; }
  double GetLandSurfacePercent(int index) { return landSurfacePercent[index]; }
  void SetNextElement(DistributedHbv *theElement) { nextElement = theElement; }
  DistributedHbv *GetNextElement() const { return nextElement; }

private:
  int geoIndex;
  int landIndex;
  int subCatchmentValue;
  int Soiltype;
  double area;
  double elevation;
  double slopeAngle;
  double aspect;
  double lakePercent;
  double forestPercent;
  double bogPercent;
  double glacierPercent;
  double openLandPercent;
  double alpinePercent;
  double heatherPercent;
  double bedrockPercent;
  double treeLevel;
  double xCoord;
  double yCoord;
  double landSurfacePercent[numberLandSurfaceClasses-1];
  DistributedHbv *nextElement;
};


//class ParametersGeneral
class ParametersGeneral
{
 public:
  ParametersGeneral();
  ~ParametersGeneral();
  void SetNUM_PREC_SERIES(int value) { NUM_PREC_SERIES = value; }
  int GetNUM_PREC_SERIES() const { return NUM_PREC_SERIES; }
  void SetNUM_TEMP_SERIES(int value) { NUM_TEMP_SERIES = value; }
  int GetNUM_TEMP_SERIES() const { return NUM_TEMP_SERIES; }
  
 private:
  int NUM_PREC_SERIES;
  int NUM_TEMP_SERIES;
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
