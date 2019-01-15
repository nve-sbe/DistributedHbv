/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments.          *
 *  Preprocessing.                                                                      *
 *                                                                                      *
 *  Author: Stein Beldring, Norway                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "classPreHbv.h"


//class SubCatchment
SubCatchment::SubCatchment():
  numLandScape(0),
  numUpStream(0),
  correction(0.0)
{
  UpStream = 0;
  SetLandScapeElement(0); 
}

SubCatchment::~SubCatchment()
{ 
}

void SubCatchment::SetNumUpStream(int value) 
{
  numUpStream = value;
  SubCatchment **upStr = new SubCatchment * [value];
  UpStream = upStr;
}


//class DistributedHbv
DistributedHbv::DistributedHbv():
  lakePercent(0.0),
  forestPercent(0.0),
  bogPercent(0.0),
  glacierPercent(0.0),
  openLandPercent(0.0),
  alpinePercent(0.0),
  heatherPercent(0.0),
  bedrockPercent(0.0),
  xCoord(0.0),
  yCoord(0.0)
{
  SetNextElement(0);
  for (int i=0; i<numberLandSurfaceClasses-1; i++) 
    landSurfacePercent[i]=0.0;
}

DistributedHbv::~DistributedHbv()
{ 
}


//class ParametersGeneral
ParametersGeneral::ParametersGeneral():
  NUM_PREC_SERIES(0),
  NUM_TEMP_SERIES(0)
{
}
     
ParametersGeneral::~ParametersGeneral()
{     
}


// class MeteorologicalStations
MeteorologicalStations::MeteorologicalStations():
  numberPrecStations(0),
  numberTempStations(0)
{
}

MeteorologicalStations::~MeteorologicalStations()
{
}

void MeteorologicalStations::SetMeteorologicalStations(ifstream &fileControl, ofstream &fout) 
{
  int i, numberPrecStations, numberTempStations;
  char stationType;
  char fileName[80];

 /*  cout << " File with meteorological stations: ";
      cin >> fileName;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  //  cout << fileName << endl;
  ifstream finMetSta(fileName);  // Open for reading
  if (!finMetSta.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finMetSta.ignore(100,':'); finMetSta >> numberPrecStations;
  finMetSta.ignore(100,':'); finMetSta >> numberTempStations;
  SetNumPrecStations(numberPrecStations);
  SetNumTempStations(numberTempStations);

  stationNumber = new int [numberPrecStations+numberTempStations];
  stationCoordX = new double [numberPrecStations+numberTempStations];
  stationCoordY = new double [numberPrecStations+numberTempStations];
  stationAltitude = new double [numberPrecStations+numberTempStations];
  for (i=0; i<numberPrecStations; i++) {
    finMetSta >> stationType >> stationNumber[i] >> stationCoordX[i] >> stationCoordY[i] >> stationAltitude[i];
    finMetSta.ignore(256,'\n');
    if (stationType != 'P' && stationType != 'p') {
      cout << endl << " Not legal station type : " << stationType << " , Station number = " << i << endl << endl;
      exit(1);
    }
  }
  for (i=0; i<numberTempStations; i++) {
    finMetSta >> stationType >> stationNumber[numberPrecStations+i] >> stationCoordX[numberPrecStations+i] 
        >> stationCoordY[numberPrecStations+i] >> stationAltitude[numberPrecStations+i];
    finMetSta.ignore(256,'\n');
    if (stationType != 'T' && stationType != 't') {
      cout << endl << " Not legal station type : " << stationType  
           << " station number = " << numberPrecStations+i << endl << endl;
      exit(1);
    }
  }

  finMetSta.close();

  fout << endl << "Meteorological stations: \n";
  fout << GetNumPrecStations() << endl;
  fout << GetNumTempStations() << endl;
  fout << "Prec.    " << endl; 
  for (i=0; i<GetNumPrecStations(); i++) {
    fout << GetStationNumber(i) << "  ";
    fout << GetStationCoordX(i) << "  ";
    fout << GetStationCoordY(i) << "  ";
    fout << GetStationAltitude(i) << endl;
  }
  fout << "Temp.    " << endl;
  for (i=0; i<GetNumTempStations(); i++) {
    fout << GetStationNumber(GetNumPrecStations()+i) << "  ";
    fout << GetStationCoordX(GetNumPrecStations()+i) << "  ";
    fout << GetStationCoordY(GetNumPrecStations()+i) << "  ";
    fout << GetStationAltitude(GetNumPrecStations()+i) << endl;
  }
  fout << endl;

}

