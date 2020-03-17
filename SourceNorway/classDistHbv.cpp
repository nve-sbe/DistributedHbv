/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments.          *
 *                                                                                      *
 *  Author: Stein Beldring, Norway                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "classDistHbv.h"


//class ParametersGeneral
ParametersGeneral::ParametersGeneral():
  SECONDS_TIMESTEP(0),
  NUM_PREC_SERIES(0),
  NUM_TEMP_SERIES(0),
  DAY_SNOW_ZERO(0),
  PREC_GRAD_LOW(0.0),
  PREC_GRAD_HIGH(0.0),
  GRAD_CHANGE_ALT(0.0),          
  PREC_CORR_RAIN(0.0),
  PREC_CORR_SNOW(0.0),
  LAPSE_DRY(0.0),
  LAPSE_WET(0.0),
  DAY_TEMP_MEMORY(0.0),
  LAKE_EPOT_PAR(0.0),
  KLAKE(0.0),
  DELTA_LEVEL(0.0),
  NLAKE(0.0),
  INITIAL_SOIL_MOISTURE(0.0),
  INITIAL_UPPER_ZONE(0.0),
  INITIAL_LOWER_ZONE(0.0),
  INITIAL_LAKE_TEMP(0.0),
  INITIAL_LAKE_LEVEL(0.0),
  INITIAL_SNOW(0.0)
{
}
     
ParametersGeneral::~ParametersGeneral()
{     
}


//class ParametersLandSurface
ParametersLandSurface::ParametersLandSurface():
  INTER_MAX(0.0),
  EPOT_PAR(0.0),
  WET_PER_CORR(0.0),
  ACC_TEMP(0.0),        
  MELT_TEMP(0.0),       
  SNOW_MELT_RATE(0.0),  
  ICE_MELT_RATE(0.0),  
  FREEZE_EFF(0.0),      
  MAX_REL(0.0),         
  ALBEDO(0.0),
  CV_SNOW(0.0)
{
  int i;
  for (i=0; i<numberSnowClasses; i++) SNOW_WEIGHT[i]=0.0;
}
     
ParametersLandSurface::~ParametersLandSurface()
{     
}


//class ParametersSubSurfaceHbv
ParametersSubSurfaceHbv::ParametersSubSurfaceHbv():
  FC(0.0),
  FCDEL(0.0),
  BETA(0.0),  
  INFMAX(0.0),         
  KUZ(0.0),         
  ALFA(0.0),  
  PERC(0.0),
  KLZ(0.0),               
  DRAW(0.0)
{
}
     
ParametersSubSurfaceHbv::~ParametersSubSurfaceHbv()
{     
}


//class SelectedTimeSeriesElements
SelectedTimeSeriesElements::SelectedTimeSeriesElements():
  numberElements(0)
{
}
     
SelectedTimeSeriesElements::~SelectedTimeSeriesElements()
{     
}

void SelectedTimeSeriesElements::SelectedTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  int i, numElements, elementId;

  /*  cout << " File with landscape elements selected for time series output: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream fileTimeSeries(fileName);  // Open for reading
  if (!fileTimeSeries.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  fileTimeSeries.ignore(100,':');
  fileTimeSeries >> numElements;
  SetNumberElements(numElements);
  TimeSeriesElements = new int [numElements];
  for (i=0; i<numElements; i++) {
    fileTimeSeries >> elementId;
    SetTimeSeriesElement(i, elementId);
    fileTimeSeries.ignore(256,'\n');
  }
  fout << endl << "Number of landscape elements selected for time series output: " << GetNumberElements() << endl;
  for (i=0; i<GetNumberElements(); i++) {
    fout << GetTimeSeriesElement(i) << endl;
  }
  fout << endl;
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


// class InputTimeSeries
InputTimeSeries::InputTimeSeries(int numberRows, int numberColums, DateTime firstTime, DateTime lastTime, int secondsPerTimeStep):
  timeSteps(numberRows),
  numberSeries(numberColums)
{
  int i;
  SetGeneralPar(0);
  datetime = new DateTime [timeSteps];
  inputArray = new double [timeSteps*numberSeries];
  datetime[0] = firstTime;
  for (i=1; i<timeSteps; i++) {
    datetime[i] = datetime[i-1] + secondsPerTimeStep; 
  }
  /*  for (i=0; i<timeSteps; i++) {
      datetime[i] = firstTime + i*secondsPerTimeStep;
      }*/
  for (i=0; i<timeSteps*numberSeries; i++) {
    inputArray[i] = missingData;
    //    cout << i << " " << inputArray[i] << "  ";
  }
  if (datetime[0] != firstTime || datetime[timeSteps-1] != lastTime) {
    cout << endl << " DateTime error during initialisation of InputTimeSeries array: " <<
      datetime[0] << "  " << firstTime  << "  " << datetime[timeSteps-1]  << "  " << lastTime << endl << endl;
    exit(1);
  }
  //  cout << datetime[0] << "  " << firstTime  << "  " << datetime[timeSteps-1]  << "  " << lastTime << endl << endl;
}

InputTimeSeries::~InputTimeSeries()
{
}

void InputTimeSeries::SetInput(ifstream &fin)
{
  char ch;
  char buffer[10240];
  int i, j;
  int date, time, year, mth, day, hour, min;
  double * tempArray = new double [numberSeries];
  fin.getline(buffer, 10240);
  fin.getline(buffer, 10240);
  i = 0;
  while (fin >> date >> ch >> time) {
    for (j=0;j<numberSeries;j++) {
      fin >> tempArray[j];
    }
    year = date/10000;
    mth = (date%10000)/100;
    day = date%100;
    hour = time/100;
    min = time%100;
    //    cout << year << " " << mth << " " << day << " " << hour << " " << min << " " << endl;
    while ((datetime[i].getYear()<year ||
       (datetime[i].getYear()==year && datetime[i].getMonth()<mth) ||
       (datetime[i].getYear()==year && datetime[i].getMonth()==mth && datetime[i].getDay()<day) ||
       (datetime[i].getYear()==year && datetime[i].getMonth()==mth && datetime[i].getDay()==day && datetime[i].getHour()<hour) ||
       (datetime[i].getYear()==year && datetime[i].getMonth()==mth && datetime[i].getDay()==day && 
        datetime[i].getHour()==hour && datetime[i].getMinute()<min)) &&
        i<timeSteps-1) {
      i++;
      //      cout << " i = " << i << endl;
    }
    //    cout << datetime[i].getYear() << " " << datetime[i].getMonth() << " " << datetime[i].getDay() << " " << 
    //      datetime[i].getHour() << " " << datetime[i].getMinute() << endl;
    if (datetime[i].getYear()==year && datetime[i].getMonth()==mth && datetime[i].getDay()==day && 
        datetime[i].getHour()==hour && datetime[i].getMinute()==min) {
      //      cout << " data found  " << i << " " << year << " " << mth << " " << day << " " << hour << " " << min << " " << endl;
      for (j=0;j<numberSeries;j++) {
        inputArray[i*numberSeries+j] = tempArray[j];
      }
    }
  }
  delete [] tempArray;
}

void InputTimeSeries::WriteInput()
{
  FILE *fp_out;
  int i, j;
  if ((fp_out = fopen("input_out.txt", "w")) == NULL ) {
    printf("\n File input_out.txt not found!\n\n");
    exit(1);
  }
  for (i=0; i<timeSteps; i++) {
    fprintf(fp_out,"%04d%02d%02d/%02d%02d",GetDateTime(i).getYear(),GetDateTime(i).getMonth(),GetDateTime(i).getDay(),
            GetDateTime(i).getHour(),GetDateTime(i).getMinute());
    for (j=0;j<numberSeries;j++){
      fprintf(fp_out,"%15.5f",GetInput(i,j));
    }
    fprintf(fp_out,"\n");
  }
  fclose(fp_out);
}


// class InputElement
InputElement::InputElement(int value):
  numberValues(value)
{
  int i;
  inputArray = new double [numberValues];
  for (i=0; i<numberValues; i++) inputArray[i] = missingData;
}

InputElement::~InputElement()
{
}


//class EvaporationControl
EvaporationControl::EvaporationControl():
  numberValues(numberPotentialEvaporationValuesPerYear),
  evaporationModellingControl('E')
{
  int i;
  evaporationArray = new double [numberValues];
  for (i=0; i<numberValues; i++) {
    evaporationArray[i] = missingData;
  }
}

EvaporationControl::~EvaporationControl()
{
}

void EvaporationControl::SetEvaporationModellingControl(char value)
{
  char fileNameInput[100];
  char buffer[256];
  int i;
  evaporationModellingControl = value;
  if (evaporationModellingControl == 'M' || evaporationModellingControl == 'm') {
    strcpy(fileNameInput,"monthlyEvaporation.txt");
    cout << " " << fileNameInput << "\n";
    ifstream finInput(fileNameInput);  // Open for reading
    if (!finInput.is_open()) {
      cout << endl << " Error opening file " << fileNameInput << endl << endl;
      exit(1);
    }
    finInput.ignore(256,'\n');
    for (i=0; i<numberValues; i++) {
      finInput.ignore(100,':'); finInput >> evaporationArray[i];
    }
    finInput.close();
    for (i=0; i<numberValues; i++) {
      cout << " " << i << "    " << GetEvaporationArray(i) << "\n";
    }
  }
}

char EvaporationControl::GetEvaporationModellingControl() const
{
    return evaporationModellingControl;
}


// class LakeWaterBalance
LakeWaterBalance::LakeWaterBalance():
  precipitation(0.0),
  temp(0.0),
  lakeEvaporation(0.0),
  waterLevelChange(0.0),
  runoff(0.0),
  discharge(0.0)
{
  SetGeneralPar(0);
  //  SetInputTimeSeries(0);
  SetInputElement(0);
  SetLandScapeElement(0);
}

LakeWaterBalance::~LakeWaterBalance()
{
}

void LakeWaterBalance::SetLakeValues(double temperature, double waterlevel)
{
  lakeTemp = temperature;
  waterLevel = waterlevel;
}

double LakeWaterBalance::GetPrecipitation() const
{
  return precipitation;
}

double LakeWaterBalance::GetTemperature() const
{
  return temp;
}

double LakeWaterBalance::GetLakeEvap() const
{
  return lakeEvaporation;
}

double LakeWaterBalance::GetRunoff() const
{
  return runoff;
}

double LakeWaterBalance::GetWaterLevel() const
{
  return waterLevel;
}

double LakeWaterBalance::GetWaterLevelChange() const
{
  return waterLevelChange;
}


// class Vegetation
Vegetation::Vegetation():
  precipitation(0.0),
  temp(0.0),
  potev(0.0),
  prevInterception(0.0),
  interceptionStore(0.0),
  interceptionLoss(0.0),
  throughFall(0.0),
  wetPeriod(0.0),
  dryPeriod(0.0)
{
  SetGeneralPar(0);
  SetLandSurfacePar(0);
  //  SetInputTimeSeries(0);
  SetInputElement(0);
  SetLandScapeElement(0);
}

Vegetation::~Vegetation()
{
}

double Vegetation::GetPrecipitation() const
{
  return precipitation;
}

double Vegetation::GetTemperature() const
{
  return temp;
}

double Vegetation::GetInterceptionStore() const
{
  return interceptionStore;
}

double Vegetation::GetInterceptionLoss() const
{
  return interceptionLoss;
}

double Vegetation::GetThroughFall() const
{
  return throughFall;
}


// class Snow
Snow::Snow():
  temp(0.0),
  snowStore(0.0),
  meltWater(0.0),
  snowWaterEquivalentChange(0.0),
  waterOutput(0.0),
  snowCoverFraction(0.0)
{
  int i;
  for (i=0; i<numberSnowClasses; i++) {
    distSnowStore[i]=0.0;
    distMeltWater[i]=0.0;
    distWaterOutput[i]=0.0;
  }
  SetGeneralPar(0);
  SetLandSurfacePar(0);
  //  SetInputTimeSeries(0);
  SetInputElement(0);
  SetLandScapeElement(0);
}

Snow::~Snow()
{
}

void Snow::SetSnowStore(double value)
{
  meltWater = value*landSurfacePar->GetMAX_REL();
  snowStore = value - meltWater;
  if (snowStore < 0.0) snowStore = 0.0;
  for (int i=0; i<numberSnowClasses; i++) {
    distMeltWater[i] = value*landSurfacePar->GetMAX_REL();
    distSnowStore[i] = value - distMeltWater[i];
    if (distSnowStore[i] < 0.0) distSnowStore[i] = 0.0;
  }
}

double Snow::GetSnowCoverFraction() const
{ 
return snowCoverFraction; 
}

double Snow::GetSnowStore() const
{
  return snowStore;
}

double Snow::GetMeltWater() const
{
  return meltWater;
}

double Snow::GetSnowWaterEquivalentChange() const
{
  return snowWaterEquivalentChange;
}

double Snow::GetWaterOutput() const
{
  return waterOutput;
}

void Snow::SetSnowValues(double snowstore, double meltwater)
{
  snowStore = snowstore;
  meltWater = meltwater;
}


// class GlacierSurface
GlacierSurface::GlacierSurface():
  precipitation(0.0),
  temp(0.0),
  iceMelt(0.0)
{
  SetGeneralPar(0);
  SetLandSurfacePar(0);
  //  SetInputTimeSeries(0);
  SetInputElement(0);
  SetLandScapeElement(0);
}

GlacierSurface::~GlacierSurface()
{
}

double GlacierSurface::GetPrecipitation() const
{
  return precipitation;
}

double GlacierSurface::GetTemperature() const
{
  return temp;
}

double GlacierSurface::GetIceMelt() const
{
  return iceMelt;
}


// class HBV
HBV::HBV():
  temp(0.0),
  soilMoisture(0.0),
  percSoilUpper(0.0),
  upperZone(0.0),
  lowerZone(0.0),
  lowerRunoff(0.0),
  transpSoilEvap(0.0),
  upperRunoff(0.0),
  runoff(0.0)
{
  SetGeneralPar(0);
  SetSubSurfaceHbvPar(0);
  SetLandSurfacePar(0);
  //  SetInputTimeSeries(0);
  SetInputElement(0);
  SetLandScapeElement(0);
}

HBV::~HBV()
{
}

void HBV::SetInitialHbvValues()
{
  soilMoisture = commonPar->GetINITIAL_SOIL_MOISTURE();
  upperZone = commonPar->GetINITIAL_UPPER_ZONE();
  lowerZone = commonPar->GetINITIAL_LOWER_ZONE();
}

void HBV::SetSubSurfaceHbvStore(double sm, double uz, double lz)
{
  soilMoisture = sm;
  upperZone = uz;
  lowerZone = lz;
}

double HBV::GetSoilMoisture() const
{
  return soilMoisture;
}

double HBV::GetSoilMoistureDeficit() const
{
  return subSurfacePar->GetFC()-soilMoisture;
}

double HBV::GetPercSoilUpper() const
{
  return percSoilUpper;
}

double HBV::GetUpperZone() const
{
  return upperZone;
}

double HBV::GetLowerZone() const
{
  return lowerZone;
}

double HBV::GetTranspSoilEvap() const
{
  return transpSoilEvap;
}

double HBV::GetRunoff() const
{
  return runoff;
}


// class Lake
Lake::Lake()
{
  SetAreaFraction(0.0);
}

Lake::~Lake()
{
}

void Lake::WaterBalance(int timeStep, DateTime datetime)
{
  GetLakeWaterBalance()->WaterBalance(timeStep, datetime);
}

double Lake::GetPrecipitation() const
{
  double precipitation=GetLakeWaterBalance()->GetPrecipitation();
  return precipitation*GetAreaFraction()/100.0;
}

double Lake::GetTemperature() const
{
  double temperature=GetLakeWaterBalance()->GetTemperature();
  return temperature*GetAreaFraction()/100.0;
}

double Lake::GetLakeEvap() const
{
  double lakeEvap=GetLakeWaterBalance()->GetLakeEvap();
  return lakeEvap*GetAreaFraction()/100.0;
}

double Lake::GetRunoff() const
{
  double runoff=GetLakeWaterBalance()->GetRunoff();
  return runoff*GetAreaFraction()/100.0;
}

double Lake::GetLakeStorage() const
{
  double lakeStorage=GetLakeWaterBalance()->GetWaterLevel();
  return lakeStorage*GetAreaFraction()/100.0;
}

double Lake::GetLakeStorageChange() const
{
  double storageChange=GetLakeWaterBalance()->GetWaterLevelChange();
  return storageChange*GetAreaFraction()/100.0;
}


// class Glacier
Glacier::Glacier():
  glacierMassBalance(0.0)
{
  SetAreaFraction(0.0);
}

Glacier::~Glacier()
{
}

void Glacier::SetSnowStore(double value) 
{
  GetSnow()->SetSnowStore(value);
}

void Glacier::SetSubSurfaceHbvStore(double sm, double uz, double lz)
{
  GetHBV()->SetSubSurfaceHbvStore(sm, uz, lz);
}

void Glacier::WaterBalance(int timeStep, DateTime datetime, int initialTimeSteps, int numberTimeSteps)
{
  GetGlacierSurface()->WaterBalance(timeStep, datetime);
  GetSnow()->WaterBalance(timeStep, datetime, GetGlacierSurface()->GetPrecipitation());
  //  GetGlacierIce()->WaterBalance(timeStep, datetime, GetGlacierSurface()->GetWaterOutput());
  if (timeStep > initialTimeSteps) {
    glacierMassBalance = glacierMassBalance + GetSnow()->GetSnowWaterEquivalentChange() - GetGlacierSurface()->GetIceMelt()*(1.0-GetSnow()->GetSnowCoverFraction());
  }
  else if (initialTimeSteps > 0) {
    glacierMassBalance = 0.0;
  }
  GetHBV()->WaterBalance(timeStep, datetime, GetGlacierSurface()->GetIceMelt()*(1.0-GetSnow()->GetSnowCoverFraction()) 
                         + GetSnow()->GetWaterOutput(), GetSnow()->GetSnowCoverFraction(), missingData);
}

double Glacier::GetPrecipitation() const
{
  double precipitation=GetGlacierSurface()->GetPrecipitation();
  return precipitation*GetAreaFraction()/100.0;
}

double Glacier::GetTemperature() const
{
  double temperature=GetGlacierSurface()->GetTemperature();
  return temperature*GetAreaFraction()/100.0;
}

double Glacier::GetSnowCoverFraction() const
{
  double snowCoverFraction=GetSnow()->GetSnowCoverFraction();
  return snowCoverFraction*GetAreaFraction()/100.0;
}

double Glacier::GetSnowStore() const
{
  double snowStore=GetSnow()->GetSnowStore();
  return snowStore*GetAreaFraction()/100.0;
}

double Glacier::GetMeltWater() const
{
  double meltWater=GetSnow()->GetMeltWater();
  return meltWater*GetAreaFraction()/100.0;
}

double Glacier::GetWaterOutput() const
{
  double waterOutput=GetGlacierSurface()->GetIceMelt()*(1.0-GetSnow()->GetSnowCoverFraction())+GetSnow()->GetWaterOutput();
  return waterOutput*GetAreaFraction()/100.0;
}

double Glacier::GetGlacierMassBalance() const
{
  return glacierMassBalance*GetAreaFraction()/100.0;
}

double Glacier::GetSurfaceMassBalance() const
{
  return glacierMassBalance;
}

void Glacier::SetSurfaceMassBalance(double value)
{
  glacierMassBalance = value;
}

double Glacier::GetGlacierIceMelt() const
{
  double iceMelt=GetGlacierSurface()->GetIceMelt()*(1.0-GetSnow()->GetSnowCoverFraction());
  return iceMelt*GetAreaFraction()/100.0;
}


double Glacier::GetSoilMoisture() const
{
  double soilMoisture=GetHBV()->GetSoilMoisture();
  return soilMoisture*GetAreaFraction()/100.0;
}

double Glacier::GetSoilMoistureDeficit() const
{
  double soilMoistureDeficit=GetHBV()->GetSoilMoistureDeficit();
  return soilMoistureDeficit*GetAreaFraction()/100.0;
}

double Glacier::GetPercSoilUpper() const
{
  double percSoilUpper=GetHBV()->GetPercSoilUpper();
  return percSoilUpper*GetAreaFraction()/100.0;
}

double Glacier::GetUpperZone() const
{
  double upperZone=GetHBV()->GetUpperZone();
  return upperZone*GetAreaFraction()/100.0;
}

double Glacier::GetLowerZone() const
{
  double lowerZone=GetHBV()->GetLowerZone();
  return lowerZone*GetAreaFraction()/100.0;
}

double Glacier::GetRunoff() const
{
  double runoff=GetHBV()->GetRunoff();
  //  cout << "glacier runoff " << runoff << "  " << GetGlacierSurface()->GetIceMelt() + GetSnow()->GetWaterOutput() << endl;
  return runoff*GetAreaFraction()/100.0;
}


// class HbvAquifer
HbvAquifer::HbvAquifer()
{
  SetAreaFraction(0.0);
  SetNextHbvAquifer(0);
}

HbvAquifer::~HbvAquifer()
{
}

void HbvAquifer::SetSnowStore(double value)
{
  if (GetNextHbvAquifer()) 
    GetNextHbvAquifer()->SetSnowStore(value);
  GetSnow()->SetSnowStore(value);
}

void HbvAquifer::SetSubSurfaceHbvStore(double sm, double uz, double lz)
{
  if (GetNextHbvAquifer()) 
    GetNextHbvAquifer()->SetSubSurfaceHbvStore(sm, uz, lz);
  GetHBV()->SetSubSurfaceHbvStore(sm, uz, lz);
}

void HbvAquifer::WaterBalance(int timeStep, DateTime datetime)
{
  GetVegetation()->WaterBalance(timeStep, datetime);
  GetSnow()->WaterBalance(timeStep, datetime, GetVegetation()->GetThroughFall());
  GetHBV()->WaterBalance(timeStep, datetime, GetSnow()->GetWaterOutput(), GetSnow()->GetSnowCoverFraction(), GetVegetation()->GetDryPeriod());
}

double HbvAquifer::GetTotalHbvAreaFraction() const
{
  double totalHbvAreaFraction=0.0;
  if (GetNextHbvAquifer()) 
    totalHbvAreaFraction=totalHbvAreaFraction+GetNextHbvAquifer()->GetTotalHbvAreaFraction();
  totalHbvAreaFraction=totalHbvAreaFraction+GetAreaFraction();
  return totalHbvAreaFraction;
}

double HbvAquifer::GetPrecipitation() const
{
  double precipitation=0.0;
  if (GetNextHbvAquifer()) 
    precipitation=precipitation+GetNextHbvAquifer()->GetPrecipitation();
  precipitation=precipitation+(GetVegetation()->GetPrecipitation()*GetAreaFraction()/100.0);
  return precipitation;
}

double HbvAquifer::GetTemperature() const
{
  double temperature=0.0;
  if (GetNextHbvAquifer()) 
    temperature=temperature+GetNextHbvAquifer()->GetTemperature();
  temperature=temperature+(GetVegetation()->GetTemperature()*GetAreaFraction()/100.0);
  return temperature;
}

double HbvAquifer::GetSnowCoverFraction() const
{
  double snowCoverFraction=0.0;
  if (GetNextHbvAquifer()) 
    snowCoverFraction=snowCoverFraction+GetNextHbvAquifer()->GetSnowCoverFraction();
  snowCoverFraction=snowCoverFraction+(GetSnow()->GetSnowCoverFraction()*GetAreaFraction()/100.0);
  return snowCoverFraction;
}

double HbvAquifer::GetSnowStore() const
{
  double snowStore=0.0;
  if (GetNextHbvAquifer()) 
    snowStore=snowStore+GetNextHbvAquifer()->GetSnowStore();
  snowStore=snowStore+(GetSnow()->GetSnowStore()*GetAreaFraction()/100.0);
  return snowStore;
}

double HbvAquifer::GetMeltWater() const
{
  double meltWater=0.0;
  if (GetNextHbvAquifer()) 
    meltWater=meltWater+GetNextHbvAquifer()->GetMeltWater();
  meltWater=meltWater+(GetSnow()->GetMeltWater()*GetAreaFraction()/100.0);
  return meltWater;
}

double HbvAquifer::GetWaterOutput() const
{
  double waterOutput=0.0;
  if (GetNextHbvAquifer()) 
    waterOutput=waterOutput+GetNextHbvAquifer()->GetWaterOutput();
  waterOutput=waterOutput+(GetSnow()->GetWaterOutput()*GetAreaFraction()/100.0);
  return waterOutput;
}

double HbvAquifer::GetSoilMoisture() const
{
  double soilMoisture=0.0;
  if (GetNextHbvAquifer()) 
    soilMoisture=soilMoisture+GetNextHbvAquifer()->GetSoilMoisture();
  soilMoisture=soilMoisture+(GetHBV()->GetSoilMoisture()*GetAreaFraction()/100.0);
  return soilMoisture;
}

double HbvAquifer::GetSoilMoistureDeficit() const
{
  double soilMoistureDeficit=0.0;
  if (GetNextHbvAquifer()) 
    soilMoistureDeficit=soilMoistureDeficit+GetNextHbvAquifer()->GetSoilMoistureDeficit();
  soilMoistureDeficit=soilMoistureDeficit+(GetHBV()->GetSoilMoistureDeficit()*GetAreaFraction()/100.0);
  return soilMoistureDeficit;
}

double HbvAquifer::GetPercSoilUpper() const
{
  double percSoilUpper=0.0;
  if (GetNextHbvAquifer()) 
    percSoilUpper=percSoilUpper+GetNextHbvAquifer()->GetPercSoilUpper();
  percSoilUpper=percSoilUpper+(GetHBV()->GetPercSoilUpper()*GetAreaFraction()/100.0);
  return percSoilUpper;
}

double HbvAquifer::GetUpperZone() const
{
  double upperZone=0.0;
  if (GetNextHbvAquifer()) 
    upperZone=upperZone+GetNextHbvAquifer()->GetUpperZone();
  upperZone=upperZone+(GetHBV()->GetUpperZone()*GetAreaFraction()/100.0);
  return upperZone;
}

double HbvAquifer::GetLowerZone() const
{
  double lowerZone=0.0;
  if (GetNextHbvAquifer()) 
    lowerZone=lowerZone+GetNextHbvAquifer()->GetLowerZone();
  lowerZone=lowerZone+(GetHBV()->GetLowerZone()*GetAreaFraction()/100.0);
  return lowerZone;
}

double HbvAquifer::GetInterceptionLoss() const
{
  double interceptionLoss=0.0;
  if (GetNextHbvAquifer()) 
    interceptionLoss=interceptionLoss+GetNextHbvAquifer()->GetInterceptionLoss();
  interceptionLoss=interceptionLoss+(GetVegetation()->GetInterceptionLoss()*GetAreaFraction()/100.0);
  return interceptionLoss;
}

double HbvAquifer::GetTranspSoilEvap() const
{
  double transpSoilEvap=0.0;
  if (GetNextHbvAquifer()) 
    transpSoilEvap=transpSoilEvap+GetNextHbvAquifer()->GetTranspSoilEvap();
  transpSoilEvap=transpSoilEvap+(GetHBV()->GetTranspSoilEvap()*GetAreaFraction()/100.0);
  return transpSoilEvap;
}

double HbvAquifer::GetRunoff() const
{
  double runoff=0.0;
  if (GetNextHbvAquifer()) 
    runoff=runoff+GetNextHbvAquifer()->GetRunoff();
  runoff=runoff+(GetHBV()->GetRunoff()*GetAreaFraction()/100.0);
  return runoff;
}


// class SubCatchment
SubCatchment::SubCatchment():
  numUpStream(0),
  numLandScape(0),
  maxBas(1.0),
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

void SubCatchment::AllocateAccumulatedDischarge(int numberTimeSteps) 
{
  accumulatedDischarge = new double [numberTimeSteps];
}

void SubCatchment::AllocateAccumulatedInFlow(int numberTimeSteps)
{
  accumulatedInFlow = new double[numberTimeSteps];
}

void SubCatchment::AllocateAccumulatedWaterBalance(int numberTimeSteps) 
{
  accumulatedPrecipitation = new double [numberTimeSteps];
  accumulatedTemperature = new double [numberTimeSteps];
  accumulatedLakeStorage = new double [numberTimeSteps];
  accumulatedSnowStore = new double [numberTimeSteps];
  accumulatedMeltWater = new double [numberTimeSteps];
  accumulatedWaterOutput = new double [numberTimeSteps];
  accumulatedSnowCoverFraction = new double [numberTimeSteps];
  accumulatedGlacierMassBalance = new double [numberTimeSteps];
  accumulatedGlacierIceMelt = new double [numberTimeSteps];
  accumulatedEvapotranspiration = new double [numberTimeSteps];
  accumulatedRunoff = new double [numberTimeSteps];
  accumulatedHbvSoilMoisture = new double [numberTimeSteps];
  accumulatedHbvSoilMoistureDeficit = new double [numberTimeSteps];
  accumulatedHbvPercSoilUpper = new double [numberTimeSteps];
  accumulatedHbvUpperZone = new double [numberTimeSteps];
  accumulatedHbvLowerZone = new double [numberTimeSteps];
  accumulatedSum = new double [numberTimeSteps];
  accumulatedSumLake = new double [numberTimeSteps];
  accumulatedSumSnow = new double [numberTimeSteps];
  accumulatedSumGlacier = new double [numberTimeSteps];
  accumulatedSumHbv = new double [numberTimeSteps];
}

void SubCatchment::AllocateWaterBalance(int numberTimeSteps) 
{
  subCatchmentPrecipitation = new double [numberTimeSteps];
  subCatchmentTemperature = new double [numberTimeSteps];
  subCatchmentLakeStorage = new double [numberTimeSteps];
  subCatchmentSnowStore = new double [numberTimeSteps];
  subCatchmentMeltWater = new double [numberTimeSteps];
  subCatchmentWaterOutput = new double [numberTimeSteps];
  subCatchmentSnowCoverFraction = new double [numberTimeSteps];
  subCatchmentGlacierMassBalance = new double [numberTimeSteps];
  subCatchmentGlacierIceMelt = new double [numberTimeSteps];
  subCatchmentEvapotranspiration = new double [numberTimeSteps];
  subCatchmentRunoff = new double [numberTimeSteps];
  subCatchmentHbvSoilMoisture = new double [numberTimeSteps];
  subCatchmentHbvSoilMoistureDeficit = new double [numberTimeSteps];
  subCatchmentHbvPercSoilUpper = new double [numberTimeSteps];
  subCatchmentHbvUpperZone = new double [numberTimeSteps];
  subCatchmentHbvLowerZone = new double [numberTimeSteps];
}

void SubCatchment::ObsDataInput(char * fileObsStreamflow, DateTime firstTime, DateTime lastTime, int numberTimeSteps, int secondsPerTimeStep)
{
  char c, buffer[256];
  int i, subCatchmentId, year, month, day, hour, min;
  double value;
  bool catchmentFound=false;
  DateTime * datetime = new DateTime [numberTimeSteps];
  datetime[0] = firstTime;
  for (i=1; i<numberTimeSteps; i++) {
    datetime[i] = datetime[i-1] + secondsPerTimeStep; 
  }
  /*  for (i=0; i<numberTimeSteps; i++) {
      datetime[i] = firstTime + i*secondsPerTimeStep;
      }*/
  if (datetime[0] != firstTime || datetime[numberTimeSteps-1] != lastTime) {
    cout << endl << " DateTime error during initialisation of observed data array: " <<
      datetime[0] << "  " << "  " << firstTime  << "  " << datetime[numberTimeSteps-1]  << "  " << lastTime << endl << endl;
    exit(1);
  }
  observedData = new double [numberTimeSteps];
  for (i=0; i<numberTimeSteps; i++) observedData[i] = missingData;
  ifstream fileIn(fileObsStreamflow);
  if (!fileIn.is_open()) {
    cout << endl << " Error opening file " << fileObsStreamflow << endl << endl;
    //    exit (1);
  }
  while (fileIn.getline(buffer, 256) && !catchmentFound) {
    if (buffer[0]=='#') {
      sscanf(buffer,"%c %c %c %c %c %d ",&c,&c,&c,&c,&c,&subCatchmentId);
      if (subCatchmentId==GetIdentifier()) {
        i=0;
        while (fileIn.peek() != '#' && fileIn.peek() != EOF) {
          fileIn >> buffer >> value;
          sscanf(buffer,"%4d %2d %2d %c %2d %2d ", &year, &month, &day, &c, &hour, &min);
          while ((datetime[i].getYear()<year ||
                  (datetime[i].getYear()==year && datetime[i].getMonth()<month) ||
                  (datetime[i].getYear()==year && datetime[i].getMonth()==month && datetime[i].getDay()<day) ||
                  (datetime[i].getYear()==year && datetime[i].getMonth()==month && datetime[i].getDay()==day && 
                   datetime[i].getHour()<hour) ||
                  (datetime[i].getYear()==year && datetime[i].getMonth()==month && datetime[i].getDay()==day && 
                   datetime[i].getHour()==hour && datetime[i].getMinute()<min)) &&
                 i<numberTimeSteps-1) {
            i++;
            //      cout << " i = " << i << endl;
          }
          //    cout << datetime[i].getYear() << " " << datetime[i].getMonth() << " " << datetime[i].getDay() << " " << 
          //      datetime[i].getHour() << " " << datetime[i].getMinute() << endl;
          if (datetime[i].getYear()==year && datetime[i].getMonth()==month && datetime[i].getDay()==day && 
              datetime[i].getHour()==hour && datetime[i].getMinute()==min) {
            //      cout << " data found  " << i << " " << year << " " << month << " " << day << " " << hour << " " << min << " " << endl;
            if (value <= -9.0) 
              observedData[i] = missingData;
            else if (observedData[i] > -9.0 && observedData[i] < 0.0) 
              observedData[i] = 0.0;
            else
              observedData[i] = value;
            //    cout << subCatchmentId << " " << i << " " << year << " " << month << " ";
            //    cout << day << " " << hour << " " << min << " " << " " << value << endl;
          }
          fileIn.ignore(256,'\n');
        }
      } else {
        while (fileIn.peek() != '#' && fileIn.peek() != EOF) {
          //      fileIn.getline(buffer, 256);
          fileIn.ignore(256,'\n');
        }
      }
    }
  }
  fileIn.close();
}  


// class DistributedHbv
DistributedHbv::DistributedHbv():
  accumulatedSum(0.0),
  accumulatedSumLake(0.0),
  accumulatedSumSnow(0.0),
  accumulatedSumGlacier(0.0),
  accumulatedSumHbv(0.0),
  area(0.0),
  elevation(0.0),
  slopeAngle(0.0),
  aspect(0.0),
  precipitationCorrection(1.0),
  temperatureCorrection(0.0),
  accumulatedDischarge(0.0),
  accumulatedPrecipitation(0.0),
  accumulatedTemperature(0.0),
  accumulatedLakeStorage(0.0),
  accumulatedSnowStore(0.0),
  accumulatedMeltWater(0.0),
  accumulatedWaterOutput(0.0),
  accumulatedSnowCoverFraction(0.0),
  accumulatedGlacierMassBalance(0.0),
  accumulatedGlacierIceMelt(0.0),
  accumulatedEvapotranspiration(0.0),
  accumulatedRunoff(0.0),
  accumulatedHbvSoilMoisture(0.0),
  accumulatedHbvSoilMoistureDeficit(0.0),
  accumulatedHbvPercSoilUpper(0.0),
  accumulatedHbvUpperZone(0.0),
  accumulatedHbvLowerZone(0.0),
  precStationsWeightedElevation(0.0),
  tempStationsWeightedElevation(0.0),
  numberSum(0),
  initialStorage(0),
  finalStorage(0),
  sumPrecipitation(0.0),
  sumEvapotranspiration(0.0),
  sumRunoff(0.0)
{
  SetNextElement(0);
  SetLake(0);
  SetGlacier(0);
  SetHbvAquifer(0);
  SetSelectedTimeSeriesElements(0);
  SetGeneralPar(0);
  SetEvaporationControlObj(0);
}

DistributedHbv::~DistributedHbv()
{ 
}

void DistributedHbv::AllocateMetSeries(int numberPrecSeries, int numberTempSeries) 
{
  metSeriesNumber = new int [numberPrecSeries+numberTempSeries];
  metSeriesWeight = new double [numberPrecSeries+numberTempSeries];
}

void DistributedHbv::AllocateWaterBalance(int numberTimeSteps) 
{
  distributedHbvPrecipitation = new double [numberTimeSteps];
  distributedHbvTemperature = new double [numberTimeSteps];
  distributedElementLakeStorage = new double [numberTimeSteps];
  distributedHbvSnowStore = new double [numberTimeSteps];
  distributedHbvSnowCoverFraction = new double [numberTimeSteps];
  distributedHbvMeltWater = new double [numberTimeSteps];
  distributedHbvWaterOutput = new double [numberTimeSteps];
  distributedHbvGlacierMassBalance = new double [numberTimeSteps];
  distributedHbvGlacierIceMelt = new double [numberTimeSteps];
  distributedHbvEvapotranspiration = new double [numberTimeSteps];
  distributedHbvRunoff = new double [numberTimeSteps];
  distributedHbvSoilMoisture = new double [numberTimeSteps];
  distributedHbvSoilMoistureDeficit = new double [numberTimeSteps];
  distributedHbvPercSoilUpper = new double [numberTimeSteps];
  distributedHbvUpperZone = new double [numberTimeSteps];
  distributedHbvLowerZone = new double [numberTimeSteps];
}

void DistributedHbv::SetAccumulatedDischarge(double localValue) 
{
   accumulatedDischarge = localValue;
}

void DistributedHbv::SetSnowStore(double value)
{
  if (GetGlacier()) GetGlacier()->SetSnowStore(value);
  if (GetHbvAquifer()) GetHbvAquifer()->SetSnowStore(value);
}

void DistributedHbv::SetSubSurfaceHbvStore(double sm, double uz, double lz)
{
  if (GetGlacier())
    if (GetGlacier()->GetHBV()) GetGlacier()->SetSubSurfaceHbvStore(sm, uz, lz);
  if (GetHbvAquifer()) GetHbvAquifer()->SetSubSurfaceHbvStore(sm, uz, lz);
}

double DistributedHbv::GetLandArea() const
{
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetLandArea " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return sumAreaFraction*GetArea();
  else return missingData;
}

double DistributedHbv::GetLakeArea() const
{
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetLake()) {
    sumAreaFraction=sumAreaFraction+GetLake()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetLakeArea " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return sumAreaFraction*GetArea();
  else return missingData;
}

double DistributedHbv::GetGlacierArea() const
{
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetGlacierArea " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return sumAreaFraction*GetArea();
  else return missingData;
}

double DistributedHbv::GetHbvArea() const
{
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    if (GetGlacier()->GetHBV()) {
      sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
      elementFound=true;
    }
  }
  if (GetHbvAquifer()) {
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetHbvArea " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return sumAreaFraction*GetArea();
  else return missingData;
}

void DistributedHbv::WaterBalance(int timeStep, DateTime datetime, int initialTimeSteps, int numberTimeSteps) const
{
  HbvAquifer *lastHbvAquifer;
  if (GetLake()) {
    GetLake()->WaterBalance(timeStep, datetime);
  }
  if (GetGlacier()) {
    GetGlacier()->WaterBalance(timeStep, datetime, initialTimeSteps, numberTimeSteps);
  }
  if (GetHbvAquifer()) {
    lastHbvAquifer=GetHbvAquifer();
    while (lastHbvAquifer) {
      lastHbvAquifer->WaterBalance(timeStep, datetime);
      lastHbvAquifer=lastHbvAquifer->GetNextHbvAquifer();
    }
  }
}

double DistributedHbv::GetPrecipitation() const
{
  double precipitation=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetLake()) { 
    precipitation=GetLake()->GetPrecipitation();
    sumAreaFraction=sumAreaFraction+GetLake()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetGlacier()) { 
    precipitation=precipitation+GetGlacier()->GetPrecipitation();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) { 
    precipitation=precipitation+GetHbvAquifer()->GetPrecipitation();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetPrecipitation " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return precipitation/sumAreaFraction;
  //  if (elementFound) return precipitation;
  else return missingData;
}

double DistributedHbv::GetTemperature() const
{
  double temperature=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetLake()) {
    temperature=GetLake()->GetTemperature();
    sumAreaFraction=sumAreaFraction+GetLake()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetGlacier()) {
    temperature=temperature+GetGlacier()->GetTemperature();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) { 
    temperature=temperature+GetHbvAquifer()->GetTemperature();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetTemperature " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return temperature/sumAreaFraction;
  //  if (elementFound) return temperature;
  else return missingData;
}

double DistributedHbv::GetLakeStorage() const
{
  double lakeStorage=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetLake()) {
    lakeStorage=GetLake()->GetLakeStorage();
    sumAreaFraction=sumAreaFraction+GetLake()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetLakeStorage " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return lakeStorage/sumAreaFraction;
  //  if (elementFound) return lakeStorage;
  else return missingData;
}

double DistributedHbv::GetLakeStorageChange() const
{
  double lakeStorageChange=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetLake()) {
    lakeStorageChange=GetLake()->GetLakeStorageChange();
    sumAreaFraction=sumAreaFraction+GetLake()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetLakeStorageChange " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return lakeStorageChange/sumAreaFraction;
  //  if (elementFound) return lakeStorageChange;
  else return missingData;
}

double DistributedHbv::GetSnowStore() const
{
  double snowStore=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    snowStore=GetGlacier()->GetSnowStore();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    snowStore=snowStore+GetHbvAquifer()->GetSnowStore();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetSnowStore " << elementFound << "  " << snowStore;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetSnowStore " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return snowStore/sumAreaFraction;
  //  if (elementFound) return snowStore;
  else return missingData;
}

double DistributedHbv::GetSnowCoverFraction() const
{
  double snowCoverFraction=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    snowCoverFraction=GetGlacier()->GetSnowCoverFraction();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    snowCoverFraction=snowCoverFraction+GetHbvAquifer()->GetSnowCoverFraction();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetSnowCoverFraction " << elementFound << "  " << snowCoverFraction;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetSnowCoverFraction " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return snowCoverFraction/sumAreaFraction;
  //  if (elementFound) return snowCoverFraction;
  else return missingData;
}

double DistributedHbv::GetMeltWater() const
{
  double meltWater=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    meltWater=GetGlacier()->GetMeltWater();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    meltWater=meltWater+GetHbvAquifer()->GetMeltWater();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetMeltWater " << elementFound << "  " << meltWater;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetMeltWater " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return meltWater/sumAreaFraction;
  //  if (elementFound) return meltWater;
  else return missingData;
}

double DistributedHbv::GetWaterOutput() const
{
  double waterOutput=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    waterOutput=GetGlacier()->GetWaterOutput();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    waterOutput=waterOutput+GetHbvAquifer()->GetWaterOutput();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetWaterOutput " << elementFound << "  " << waterOutput;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetWaterOutput " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return waterOutput/sumAreaFraction;
  //  if (elementFound) return waterOutput;
  else return missingData;
}

double DistributedHbv::GetGlacierMassBalance() const
{
  double glacierMassBalance=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  //  cout << "start GetGlacierMassBalance\n";
  if (GetGlacier()) {
    glacierMassBalance=GetGlacier()->GetGlacierMassBalance();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetGlacierMassBalance " << elementFound << "  " << glacierMassBalance;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetGlacierMassBalance " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  //  cout << "  end GetGlacierMassBalance\n";
  if (elementFound) return glacierMassBalance/sumAreaFraction;
  //  if (elementFound) return glacierMassBalance;
  else return missingData;
}

double DistributedHbv::GetGlacierIceMelt() const
{
  double glacierIceMelt=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  //  cout << "start GetGlacierIceMelt\n";
  if (GetGlacier()) {
    glacierIceMelt=GetGlacier()->GetGlacierIceMelt();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetGlacierIceMelt " << elementFound << "  " << glacierIceMelt;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetGlacierIceMelt " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  //  cout << "  end GetGlacierIceMelt\n";
  if (elementFound) return glacierIceMelt/sumAreaFraction;
  //  if (elementFound) return glacierIceMelt;
  else return missingData;
}

double DistributedHbv::GetHbvSoilMoisture() const
{
  double soilMoisture=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    if (GetGlacier()->GetHBV()) {
      soilMoisture=GetGlacier()->GetSoilMoisture();
      sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
      elementFound=true;
    }
  }
  if (GetHbvAquifer()) {
    soilMoisture=soilMoisture+GetHbvAquifer()->GetSoilMoisture();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetSoilMoisture " << elementFound << "  " << soilMoisture;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetSoilMoisture " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return soilMoisture/sumAreaFraction;
  //  if (elementFound) return soilMoisture;
  else return missingData;
}

double DistributedHbv::GetHbvSoilMoistureDeficit() const
{
  double soilMoistureDeficit=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    if (GetGlacier()->GetHBV()) {
      soilMoistureDeficit=GetGlacier()->GetSoilMoistureDeficit();
      sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
      elementFound=true;
    }
  }
  if (GetHbvAquifer()) {
    soilMoistureDeficit=soilMoistureDeficit+GetHbvAquifer()->GetSoilMoistureDeficit();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetSoilMoistureDeficit " << elementFound << "  " << soilMoistureDeficit;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetSoilMoistureDeficit " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return soilMoistureDeficit/sumAreaFraction;
  //  if (elementFound) return soilMoistureDeficit;
  else return missingData;
}

double DistributedHbv::GetHbvPercSoilUpper() const
{
  double percSoilUpper=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    if (GetGlacier()->GetHBV()) {
      percSoilUpper=GetGlacier()->GetPercSoilUpper();
      sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
      elementFound=true;
    }
  }
  if (GetHbvAquifer()) {
    percSoilUpper=percSoilUpper+GetHbvAquifer()->GetPercSoilUpper();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetPercSoilUpper " << elementFound << "  " << percSoilUpper;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetPercSoilUpper " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return percSoilUpper/sumAreaFraction;
  //  if (elementFound) return percSoilUpper;
  else return missingData;
}

double DistributedHbv::GetHbvUpperZone() const
{
  double upperZone=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    if (GetGlacier()->GetHBV()) {
      upperZone=GetGlacier()->GetUpperZone();
      sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
      elementFound=true;
    }
  }
  if (GetHbvAquifer()) {
    upperZone=upperZone+GetHbvAquifer()->GetUpperZone();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetHbvUpperZone " << elementFound << "  " << upperZone;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetHbvUpperZone " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return upperZone/sumAreaFraction;
  //  if (elementFound) return upperZone;
  else return missingData;
}

double DistributedHbv::GetHbvLowerZone() const
{
  double lowerZone=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    if (GetGlacier()->GetHBV()) {
      lowerZone=GetGlacier()->GetLowerZone();
      sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
      elementFound=true;
    }
  }
  if (GetHbvAquifer()) {
    lowerZone=lowerZone+GetHbvAquifer()->GetLowerZone();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetHbvLowerZone " << elementFound << "  " << lowerZone;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetHbvLowerZone " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return lowerZone/sumAreaFraction;
  //  if (elementFound) return lowerZone;
  else return missingData;
}

double DistributedHbv::GetEvapotranspiration() const
{
  double evapotranspiration=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    evapotranspiration=evapotranspiration+GetHbvAquifer()->GetInterceptionLoss()+GetHbvAquifer()->GetTranspSoilEvap();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetLake()) {
    evapotranspiration=evapotranspiration+GetLake()->GetLakeEvap();
    sumAreaFraction=sumAreaFraction+GetLake()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetInterceptionLoss " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return evapotranspiration/sumAreaFraction;
  //  if (elementFound) return evapotranspiration;
  else return missingData;
}

double DistributedHbv::GetInterceptionLoss() const
{
  double interceptionLoss=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    interceptionLoss=interceptionLoss+GetHbvAquifer()->GetInterceptionLoss();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetInterceptionLoss " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return interceptionLoss/sumAreaFraction;
  //  if (elementFound) return interceptionLoss;
  else return missingData;
}

double DistributedHbv::GetTranspSoilEvap() const
{
  double transpSoilEvap=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    transpSoilEvap=transpSoilEvap+GetHbvAquifer()->GetTranspSoilEvap();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetTranspSoilEvap " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return transpSoilEvap/sumAreaFraction;
  //  if (elementFound) return transpSoilEvap;
  else return missingData;
}

double DistributedHbv::GetLakeEvap() const
{
  double lakeEvap=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetLake()) {
    lakeEvap=GetLake()->GetLakeEvap();
    sumAreaFraction=sumAreaFraction+GetLake()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetLakeEvap " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return lakeEvap/sumAreaFraction;
  //  if (elementFound) return lakeEvap;
  else return missingData;
}

double DistributedHbv::GetRunoff() const
{
  double runoff=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetLake()) {
    runoff=GetLake()->GetRunoff();
    sumAreaFraction=sumAreaFraction+GetLake()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetGlacier()) {
    runoff=runoff+GetGlacier()->GetRunoff();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    runoff=runoff+GetHbvAquifer()->GetRunoff();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0+epsilon)) {
    cout << " sumAreaFraction GetRunoff " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return runoff/sumAreaFraction;
  //  if (elementFound) return runoff;
  else return missingData;
}

double DistributedHbv::GetDischarge() const
{
  if (GetRunoff() != missingData) return GetRunoff()*area/commonPar->GetSECONDS_TIMESTEP();         /*  runoff (m) -> discharge (m3/s)  */
  else return missingData;
}

int DistributedHbv::GetNumberSum() const
{
    return numberSum;
}

void DistributedHbv::SetSumWaterBalance()
{
  if (GetRunoff() != missingData) {
    sumPrecipitation = sumPrecipitation + GetPrecipitation();
    sumEvapotranspiration = sumEvapotranspiration + GetEvapotranspiration();
    sumRunoff = sumRunoff + GetRunoff();
    numberSum++;
  }
}

void DistributedHbv::SetInitialStorage()
{
  if (GetHbvAquifer()) {
    initialStorage = initialStorage + GetSnowStore() + GetMeltWater() + GetHbvSoilMoisture() + GetHbvUpperZone() + GetHbvLowerZone() + GetLakeStorage();
    //    initialStorage = initialStorage + GetSnowStore() + GetMeltWater() + GetHbvSoilMoisture() + GetHbvUpperZone() + GetHbvLowerZone() + GetLakeStorage() + GetGlacierIceThickness()*GetGeneralPar()->GetDENSITY_ICE();
  }
}

double DistributedHbv::GetInitialStorage() const
{
    return initialStorage;
}

void DistributedHbv::SetFinalStorage()
{
  if (GetHbvAquifer()) {
    finalStorage = finalStorage + GetSnowStore() + GetMeltWater() + GetHbvSoilMoisture() + GetHbvUpperZone() + GetHbvLowerZone() + GetLakeStorage();
    //    finalStorage = finalStorage + GetSnowStore() + GetMeltWater() + GetHbvSoilMoisture() + GetHbvUpperZone() + GetHbvLowerZone() + GetLakeStorage() + GetGlacierIceThickness()*GetGeneralPar()->GetDENSITY_ICE();
  }
}

double DistributedHbv::GetFinalStorage() const
{
    return finalStorage;
}

double DistributedHbv::GetSumPrecipitation() const
{
    return sumPrecipitation;
}

double DistributedHbv::GetSumEvapotranspiration() const
{
    return sumEvapotranspiration;
}

double DistributedHbv::GetSumRunoff() const
{
    return sumRunoff;
}

