/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments.          *
 *                                                                                      *
 *  Author: Stein Beldring, Norway                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "classDistHbv.h"
#include "parameters.h"
#include "utilities.h"
#include "netcdf.h"

#ifndef FALSE
#define FALSE 0
#define TRUE !FALSE
#endif

#define NODATA 1.e+20f  //nodata in .nc files made for KSS. (HySN5, land masks and climate projections)
#define NPARAM 4
#define S_PR_DAY 86400
#define KELVIN 273.15

void ReadSubCatchmentIdentifier(DistributedHbv * const DistHbv, SubCatchment * const CatchmentElement,
                                int numWatc, ifstream &fileControl, ofstream &fout);
void SnowGlacierIceReDistribution(SubCatchment ** const Outlet, DistributedHbv * const DistHbv, ParametersGeneral * ParGeneralStore, int initialTimeSteps, int numberTimeSteps,
                                  int numLand, int numWatcOut, int timeStep, DateTime datetime, 
                                  int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, bool modelCalibration, ofstream &fout);
void WaterBalanceTimeSeries(DistributedHbv * const DistHbv, ParametersGeneral * const ParGeneralStore,
                            MeteorologicalStations * const MetStations, int initialTimeSteps, int numberTimeSteps, 
                            InputTimeSeries * InputTimeSeriesStore, InputElement * InputElementStore, 
                            int numLand, int timeStep, DateTime datetime, bool * inputDataFound);
void WaterBalanceGrid(DistributedHbv * const DistHbv,  ParametersGeneral * ParGeneralStore, 
                      InputElement * InputElementStore, int initialTimeSteps, int numberTimeSteps,
                      int numLand, int timeStep, int nRows, int nCols, DateTime datetime, char * precPath,
                      unsigned short int * precip10, unsigned short int * temp10K, unsigned short int * tmax10K, 
		      unsigned short int * tmin10K, unsigned short int * wind10, unsigned short int * solar10, 
	              unsigned short int * vp10,  bool * inputDataFound,
                      int * indexStore, int numberIndexStore, ofstream &fout);
void WaterBalanceGridNetcdf(DistributedHbv * const DistHbv,  ParametersGeneral * ParGeneralStore, 
			    InputElement * InputElementStore, int initialTimeSteps, int numberTimeSteps,
			    int numLand, int timeStep, int nRows, int nCols, DateTime datetime,
			    float *prec_in, float *temp_in, float *tmax_in, 
			    float *tmin_in, float *wind_in, float *srad_in, 
			    float *vp_in, bool * inputDataFound,
			    int * indexStore, int numberIndexStore, ofstream &fout);
void ReadNetcdf(int initialTimeSteps, int numberTimeSteps,
		int numLand, int timeStep, int nRows, int nCols, 
		DateTime datetime,char ETscheme,
		char * precPath,char * tmeanPath, char *tmaxPath,
		char * tminPath,char *windPath,char * metPath3,
		float *prec_in, float *temp_in, float *tmax_in, 
		float *tmin_in, float *wind_in, float *srad_in, 
		float *vp_in,float *hurs_in);
void ReadNetcdfClimateProj(int initialTimeSteps, int numberTimeSteps,
			   int numLand, int timeStep, int nRows, int nCols, 
			   DateTime datetime,char ETscheme,char *precPath,
			   char *tmeanPath, char *tmaxPath,char *tminPath,
			   char *windPath,char *rsdsPath,char *vpPath,
			   char *forcingName,char *rcpName,char *biasName,			   
			   float *prec_in, float *temp_in, float *tmax_in, 
			   float *tmin_in, float *wind_in, float *srad_in, 
			   float *vp_in,float *hurs_in);
void HandleError(int);
void TraverseCorrectionSubCatchment(SubCatchment * const thisSubCatchment, int numberCorrectionCatchments,
                                   int * correctionCatchments, double * correctionPrecipitation, 
                                   double * correctionTemperature, ofstream &fout);
void TraverseCorrectionLandScape(DistributedHbv * const thisElement, double precCorr, double tempCorr);
void TraverseSubCatchment(SubCatchment * const thisSubCatchment, int timeStep, ofstream &fout);
void TraverseLandScape(DistributedHbv * const thisElement, int timeStep, ofstream &fout);
void TraverseMissingDataSubCatchment(SubCatchment * const thisSubCatchment, int timeStep, ofstream &fout);
void TraverseMissingDataLandScape(DistributedHbv * const thisElement, int timeStep, ofstream &fout);
void WriteSubCatchmentIdentifier(SubCatchment * const CatchmentElement, int numWatc, ofstream &fout);
void WriteSubCatchmentDischarge(SubCatchment * CatchmentElement, int numWatc, DateTime startSimulationTime, 
                                DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                                int secondsPerTimeStep, bool modelCalibration, ofstream &fout);
void WriteSubCatchmentWaterBalance(SubCatchment * CatchmentElement, int numWatc, DateTime startSimulationTime, 
                                   DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                                   int secondsPerTimeStep);
void WriteBinaryGrid(DistributedHbv * const DistHbv, DateTime datetime, int numLand, int timeStep, int nRows,int nCols, 
                    int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout, int writePET);
void WriteAsciiGrid(DistributedHbv * const DistHbv, DateTime datetime, int numLand, int timeStep, int nRows,int nCols, 
                    int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void WriteAsciiGridWaterBalance(DistributedHbv * const DistHbv, DateTime startSimulationTime, DateTime endSimulationTime, int numLand, int nRows, int nCols, 
                    int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void ObjectiveCriteria(int numberTimeSteps, double *obs, double *sim, 
                       double *ns, double *rmse, double *bias, double *pears);
void WriteDistributedHbvTimeSeries(DistributedHbv * const DistHbv, int numLand, 
                                   DateTime startSimulationTime, DateTime endSimulationTime, 
                                   int initialTimeSteps, int numberTimeSteps, int secondsPerTimeStep);
void ReadModelStateVariables(DistributedHbv * const DistHbv, int numLand);
void WriteModelStateVariables(DistributedHbv * const DistHbv, int numLand, DateTime datetime, int timeStep, int secondsTimestep);

int main(int argc, char *argv[])
{
  std::cout << "\n\n Distributed HBV water balance model \n\n";

  bool inputDataFound=false;
  bool modelCalibration=false;
  bool readModelStates=false;
  bool writeModelStates=false;
  char fileName[100];
  char metPathFileName[100];
  char metPath2FileName[100];
  char metPath3FileName[100];
  char fileNameInput[100];
  char fileObsStreamflow[100];
  char buffer[256];
  char modelRun = 'R';
  char modelStates = 'S';
  char inputFormat = 'F';
  char inputFileFormat = 'N'; //see control file
  char evaporationModellingControl = 'T'; //T: temperature index, P: penman monteith, see control file
  char forcingType = 'S'; //S: senorge obs based, C: climate projections, see control file
  char forcingName [50];
  char rcpName [50];  
  char biasName [50];
  char ch;
  char precPath[200];
  char tmeanPath[200];
  char tmaxPath[200];
  char tminPath[200];
  char windPath[200];
  char rsdsPath[200];
  char vpPath[200];
  char hursPath[200];
  int i, j, k;
  int numLand;
  int geoIndex;
  int landIndex;
  int nRows, nCols, noData;
  int landSurf;
  int soil;
  int numberTimeSteps, initialTimeSteps, timeStep;
  int numWatc, numWatcUp, numWatcOut;
  int subCatchmentId;
  int day, mth, year, hour, minute;
  int writePET;
  //  int time, startModelTime, startSimulationTime, endSimulationTime;
  int numberCorrectionCatchments;
  int numberIndexStore;
  int seriesNumber;
  int ndays,oldyear;
  double seriesWeight;
  double sumWeight;
  double sumArea;
  double precStationsWeightedElevation, tempStationsWeightedElevation;
  double correction;
  double xllCorner, yllCorner, cellSize;
  double elementArea, elementLatitude, elementElevation, elementTimefalling, elementSlopeAngle, elementAspect, lakePercent, glacierPercent;
  double areaFraction[maximumNumberLandClasses];
  DateTime datetime;
  DateTime datetime_nc;
  Lake *ptrLake;
  Glacier *ptrGlacier;
  HbvAquifer *ptrHbvAquifer;
  HbvAquifer *lastHbvAquifer;
  LakeWaterBalance *ptrLakeWaterBalance;
  Vegetation *ptrVegetation;
  Snow *ptrSnow;
  GlacierSurface * ptrSurface;
  //  GlacierIce *ptrIce;
  HBV *ptrHbv;
  LANDSURFACE landSurfType[numberLandSurfaceClasses];
  SOIL soilType[numberSoilClasses];
  float *prec_in,*temp_in,*tmax_in,*tmin_in,*srad_in,*vp_in,*hurs_in,*wind_in;
  //FORCESTRUCT Forcings[6];

  DateTime nowDate;
  nowDate.now();
  DateTime finalDate(finalYear,finalMonth,finalDay,finalHour,finalMinute,0);
  if (nowDate > finalDate) {
    cout << "\n *\n *  Program is terminated !  \n *\n\n\n";
    exit(1);
  }

  if (argc != 2) {
    cout << " " << argv[0] << "  <control file name>\n\n";
    exit(1);
  }
  ifstream fileControl(argv[1]);
  if (!fileControl.is_open()) {
    cout << " Error opening file " << argv[1] << endl << endl;
    exit (1);
  }

  fileControl.ignore(100,':');
  fileControl >> modelRun;
  fileControl.ignore(256,'\n');
  if (modelRun != 'S' && modelRun != 'C' && modelRun != 's' && modelRun != 'c') {
    cout << "\n Type of model run, simulation(S) or calibration(C) \n\n";
    exit (1);
  }
  if (modelRun == 'C' || modelRun == 'c') modelCalibration = true;

  fileControl.ignore(100,':');
  fileControl >> inputFormat;
  fileControl.ignore(256,'\n');
  if (inputFormat != 'G' && inputFormat != 'T' && inputFormat != 'g' && inputFormat != 't') {
    cout << "\n Input data format, grid files(G) or time series file(T) \n\n";
    exit (1);
  }

  fileControl.ignore(100,':');
  fileControl >> inputFileFormat;
  fileControl.ignore(256,'\n');
  fprintf(stderr,"fileformat %c\n",inputFileFormat);

  if (inputFileFormat != 'N' && inputFileFormat != 'B' && inputFileFormat != 'n' && inputFileFormat != 'b') {
    cout << "\n Input file format, binary grid files(B) or netcdf(N) \n\n";
    exit (1);
  }

  fileControl.ignore(100,':');
  fileControl >> evaporationModellingControl;
  fileControl.ignore(256,'\n');
  fprintf(stderr,"Evapotranspiration scheme: %c\n",evaporationModellingControl);

  if (evaporationModellingControl != 'T' && evaporationModellingControl != 'P') {
    cout << "\n Evapotranspiration scheme, temperature index (T) or Penman Monteith (P) \n\n";
    exit (1);
  }

  // Object for storing evaporation modelling scheme information. 
  EvaporationControl * EvaporationControlObject = new EvaporationControl;
  EvaporationControlObject->SetEvaporationModellingControl(evaporationModellingControl); 

  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ofstream fout(fileName);  // Open for writing
  fprintf(stderr,"output file %s\n",fileName);
 
  fileControl.ignore(100, ':');
  fileControl >> writePET;
  fileControl.ignore(256, '\n');
  fprintf(stderr,"writePET %d\n\n",writePET);

  // Read dates for start model, start simulation, end simulation
  /*  cout << " Start model date and time (day, month, year, hour, minute): ";
      cin >> day >> mth >> year >> hour >> minute;*/
  fileControl.ignore(100,':');
  fileControl >> day >> mth >> year >> hour >> minute;
  fileControl.ignore(256,'\n');
  DateTime startModelTime(year,mth,day,hour,minute,0);
  if (startModelTime.legal() != 1) {
    cout << endl << " Not legal date and time " << endl << endl;
    exit(1);
  }
  
  //cout << " Start simulation date and time (day, month, year, hour, minute): ";
  //cin >> day >> mth >> year >> hour >> minute;
      
  fileControl.ignore(100,':');
  fileControl >> day >> mth >> year >> hour >> minute;
  fileControl.ignore(256,'\n');
  DateTime startSimulationTime(year,mth,day,hour,minute,0);
  if (startSimulationTime.legal() != 1) {
    cout << endl << " Not legal date and time " << endl << endl;
    exit(1);
  }
  /*  cout << endl;
      cout << " End simulation date and time (day, month, year, hour, minute): ";
      cin >> day >> mth >> year >> hour >> minute;*/
  fileControl.ignore(100,':');
  fileControl >> day >> mth >> year >> hour >> minute;
  fileControl.ignore(256,'\n');
  DateTime endSimulationTime(year,mth,day,hour,minute,0);
  if (endSimulationTime.legal() != 1) {
    cout << endl << " Not legal date and time " << endl << endl;
    exit(1);
  }
  cout << endl;
  if (startModelTime > startSimulationTime || startSimulationTime > endSimulationTime) {
    cout << endl << " DateTime error " << endl << endl;
    exit(1);
  }

  cout << startSimulationTime << endl << endSimulationTime << endl;

  fileControl.ignore(100, ':');
  fileControl >> forcingType;
  fileControl.ignore(256, '\n');
  fprintf(stderr,"forcingtype %c\n\n",forcingType);
  if (forcingType != 'S' && forcingType != 'C') {
    cout << "\n Input forcing file type, senorge(S) or climateprojections(C) \n\n";
    exit (1);
  }

  // Read information on paths to input forcings
  // (meteorological data, 7 paths, in order to accomodate seNorge, klinogrid, HySN5 and climateprojections)
  // if evaporation scheme = T: will read only P and T forcing files
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  strcpy(precPath,fileName);
  printf(" precpath (prec) %s\n",fileName);
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  strcpy(tmeanPath,fileName);
  printf(" tmeanpath (tmean) %s\n",fileName);
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  strcpy(tmaxPath,fileName);
  printf(" tmaxpath (tmax) %s\n",fileName);
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  strcpy(tminPath,fileName);
  printf(" tminpath (tmin) %s\n",fileName);  
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  strcpy(windPath,fileName);
  printf(" windpath (wind) %s\n",fileName);
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  strcpy(rsdsPath,fileName);
  printf(" rsdspath (rsds) %s\n",fileName);
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  strcpy(hursPath,fileName);
  printf(" hurspath (hurs) %s\n",fileName);

  // Read information on forcing name (used only when reading climate projections)
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  strcpy(forcingName,fileName);
  printf(" \n forcingname: %s\n",fileName);
  // Read information on rcp name (used only when reading climate projections)
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  strcpy(rcpName,fileName);
  printf(" rcpname: %s\n",fileName);
  // Read information on bias name (used only when reading climate projections)
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  strcpy(biasName,fileName);
  printf(" biasname: %s\n",fileName);

  // Object for storing meteorological station information
  MeteorologicalStations * MetStations = new MeteorologicalStations;
  MetStations->SetMeteorologicalStations(fileControl, fout);

  // Read common parameters file and set parameter values
  ParametersGeneral * ParGeneralStore = new ParametersGeneral;
  SetGeneralParameters(ParGeneralStore, fileControl, fout);
  if ((ParGeneralStore->GetSECONDS_TIMESTEP() % minimumTimeStep) != 0) {
    cout << endl << " Temporal resolution must be a multiple of " << minimumTimeStep << " seconds : " << ParGeneralStore->GetSECONDS_TIMESTEP() << endl << endl;
    exit(1);
  }
  cout << endl << " Common parameters file read " << endl;
  // End read common parameters

  // Read landsurface parameters file and set parameter values
  ParametersLandSurface * ParLandSurfaceStore = new ParametersLandSurface [numberLandSurfaceClasses];
  SetLandSurfaceParameters(ParLandSurfaceStore, fileControl, fout);
  // End read landsurface parameters

  // Read subsurface parameters file and set parameter values
  ParametersSubSurfaceHbv * ParSubSurfaceHbvStore = new ParametersSubSurfaceHbv [numberSoilClasses];
  SetSubSurfaceHbvParameters(ParSubSurfaceHbvStore, fileControl, fout);
  // End read subsurface parameters

  // Object for storing landscape elements selected for time series output
  SelectedTimeSeriesElements * SelectedTimeSeriesElementsStore = new SelectedTimeSeriesElements;
  SelectedTimeSeriesElementsStore->SelectedTimeSeriesElementsInput(fileControl, fout);

  // Object for storing input data for each landscape element
  InputElement * InputElementStore = new InputElement(numberInputSeries);

  // Calculate no. of initial and simulation time steps
  initialTimeSteps = (startSimulationTime.date2jday() - startModelTime.date2jday()) * (int)(numberSecondsDay/ParGeneralStore->GetSECONDS_TIMESTEP());
  numberTimeSteps = (1 + endSimulationTime.date2jday() - startSimulationTime.date2jday()) * (int)(numberSecondsDay/ParGeneralStore->GetSECONDS_TIMESTEP());

  // Read landscape element file and generate landscape element objects  
  /*  cout << " File with landscape element information: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finLandScape(fileName);  // Open for reading
  if (!finLandScape.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  else  {
    cout << endl << " File with landscape element information opened: " << fileName << endl;
  }
  finLandScape >> buffer >> nCols;
  finLandScape >> buffer >> nRows;
  finLandScape >> buffer >> xllCorner;
  finLandScape >> buffer >> yllCorner;
  finLandScape >> buffer >> cellSize;
  finLandScape >> buffer >> noData;
  finLandScape.ignore(100,':');
  finLandScape >> numLand;
  DistributedHbv * DistHbv = new DistributedHbv [numLand];
  if (!DistHbv) {
    cout << endl << " Error during dynamic memory allocating for DistHbv[" << numLand << "]" << endl << endl;
    exit(1);
  } 
  for (i=0; i<numLand; i++) {
    finLandScape >> landIndex >> geoIndex >> elementArea >> elementLatitude >> elementElevation >> elementTimefalling;
    finLandScape >> elementSlopeAngle >> elementAspect >> lakePercent >> glacierPercent;
    sumArea = lakePercent + glacierPercent;
    for (j=0; j<maximumNumberLandClasses; j++) {
      finLandScape >> landSurf >> soil >> areaFraction[j]; 
      landSurfType[j]=LANDSURFACE(landSurf);
      soilType[j]=SOIL(soil);
      sumArea = sumArea + areaFraction[j];
    }

    if (sumArea != 100.0) {
      lakePercent = lakePercent*100.0/sumArea;
      glacierPercent = glacierPercent*100.0/sumArea;
      for (j=0; j<maximumNumberLandClasses; j++) 
        areaFraction[j] = areaFraction[j]*100.0/sumArea;
    }

    DistHbv[landIndex].SetEvaporationControlObj(EvaporationControlObject);
    DistHbv[landIndex].SetGeoIndex(geoIndex);
    DistHbv[landIndex].SetLandIndex(landIndex);
    DistHbv[landIndex].SetArea(elementArea);
    DistHbv[landIndex].SetLatitude(elementLatitude);
    DistHbv[landIndex].SetElevation(elementElevation);
    DistHbv[landIndex].SetTimeFalling(elementTimefalling);
    DistHbv[landIndex].SetSlopeAngle(elementSlopeAngle);
    DistHbv[landIndex].SetAspect(elementAspect);
    DistHbv[landIndex].SetSelectedTimeSeriesElements(SelectedTimeSeriesElementsStore);
    DistHbv[landIndex].SetGeneralPar(ParGeneralStore);

    // Allocate space for water balance time series for landscape elements 
    for (j=0; j<DistHbv[landIndex].GetSelectedTimeSeriesElements()->GetNumberElements(); j++) {
      if (DistHbv[landIndex].GetLandIndex() == DistHbv[landIndex].GetSelectedTimeSeriesElements()->GetTimeSeriesElement(j)) {
        DistHbv[landIndex].AllocateWaterBalance(initialTimeSteps+numberTimeSteps);
      }
    }
    if (inputFormat == 'T' || inputFormat == 't') {
      sumWeight=0.0;
      precStationsWeightedElevation=0.0;
      DistHbv[landIndex].AllocateMetSeries(ParGeneralStore->GetNUM_PREC_SERIES(),ParGeneralStore->GetNUM_TEMP_SERIES());
      for (j=0; j<ParGeneralStore->GetNUM_PREC_SERIES(); j++) {
        finLandScape >> seriesNumber >> seriesWeight; 
	//	cout << " P  " << seriesNumber << "  " << seriesWeight << "    ";
        DistHbv[landIndex].SetMetSeriesNumber(j, seriesNumber);
        DistHbv[landIndex].SetMetSeriesWeight(j, seriesWeight);
        precStationsWeightedElevation = precStationsWeightedElevation+MetStations->GetStationAltitude(seriesNumber)*seriesWeight;
        sumWeight = sumWeight+seriesWeight;
      }
      if (sumWeight != 1.0) {
        for (j=0; j<ParGeneralStore->GetNUM_PREC_SERIES(); j++) {
          seriesWeight = DistHbv[landIndex].GetMetSeriesWeight(j)/sumWeight;
          DistHbv[landIndex].SetMetSeriesWeight(j, seriesWeight);
          precStationsWeightedElevation = precStationsWeightedElevation/sumWeight;
        }
      }
      DistHbv[landIndex].SetPrecStationsWeightedElevation(precStationsWeightedElevation);
      //      cout << " precStationsWeightedElevation = " << precStationsWeightedElevation << endl;
      sumWeight=0.0;
      tempStationsWeightedElevation=0.0;
      for (j=0; j<ParGeneralStore->GetNUM_TEMP_SERIES(); j++) {
        finLandScape >> seriesNumber >> seriesWeight; 
	//	cout << " T  " << seriesNumber << "  " << seriesWeight << "    ";
        DistHbv[landIndex].SetMetSeriesNumber(ParGeneralStore->GetNUM_PREC_SERIES()+j,seriesNumber);
        DistHbv[landIndex].SetMetSeriesWeight(ParGeneralStore->GetNUM_PREC_SERIES()+j, seriesWeight);
        tempStationsWeightedElevation = tempStationsWeightedElevation+MetStations->GetStationAltitude(seriesNumber)*seriesWeight;
        sumWeight = sumWeight+seriesWeight;
      }
      //    cout << endl;
      if (sumWeight != 1.0) {
        for (j=0; j<ParGeneralStore->GetNUM_TEMP_SERIES(); j++) {
          seriesWeight = DistHbv[landIndex].GetMetSeriesWeight(ParGeneralStore->GetNUM_PREC_SERIES()+j)/sumWeight;
          DistHbv[landIndex].SetMetSeriesWeight(ParGeneralStore->GetNUM_PREC_SERIES()+j, seriesWeight);
          tempStationsWeightedElevation = tempStationsWeightedElevation/sumWeight;
        }
      }
      DistHbv[landIndex].SetTempStationsWeightedElevation(tempStationsWeightedElevation);
    }
    else {    
      finLandScape.ignore(256,'\n');
    }
    if (lakePercent>0.0) {
      ptrLake = new Lake;
      ptrLake->SetAreaFraction(lakePercent);
      ptrLakeWaterBalance = new LakeWaterBalance;
      ptrLakeWaterBalance->SetGeneralPar(ParGeneralStore);
      //      ptrLakeWaterBalance->SetInputTimeSeries(InputTimeSeriesStore);
      ptrLakeWaterBalance->SetInputElement(InputElementStore);
      ptrLakeWaterBalance->SetLandScapeElement(&DistHbv[landIndex]);
      ptrLakeWaterBalance->SetLakeValues(ParGeneralStore->GetINITIAL_LAKE_TEMP(), ParGeneralStore->GetINITIAL_LAKE_LEVEL());
      ptrLake->SetLakeWaterBalance(ptrLakeWaterBalance);
      DistHbv[landIndex].SetLake(ptrLake);
    }
    if (glacierPercent>0.0) {
      ptrGlacier = new Glacier;
      ptrGlacier->SetAreaFraction(glacierPercent);
      ptrSurface = new GlacierSurface;
      ptrSurface->SetGeneralPar(ParGeneralStore);
      ptrSurface->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]); 
      //      ptrSurface->SetInputTimeSeries(InputTimeSeriesStore);
      ptrSurface->SetInputElement(InputElementStore);
      ptrSurface->SetLandScapeElement(&DistHbv[landIndex]);
      ptrGlacier->SetGlacierSurface(ptrSurface);
      ptrSnow = new Snow;
      ptrSnow->SetGeneralPar(ParGeneralStore);
      ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]); 
      //      ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
      ptrSnow->SetInputElement(InputElementStore);
      ptrSnow->SetLandScapeElement(&DistHbv[landIndex]);
      ptrSnow->SetSnowStore(ParGeneralStore->GetINITIAL_SNOW());
      ptrGlacier->SetSnow(ptrSnow);
      //      ptrIce = new GlacierIce;
      //      ptrIce->SetGeneralPar(ParGeneralStore);
      //      ptrIce->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]); 
      //      ptrIce->SetInputTimeSeries(InputTimeSeriesStore);
      //      ptrIce->SetInputElement(InputElementStore);
      //      ptrIce->SetLandScapeElement(&DistHbv[landIndex]);
      //      ptrGlacier->SetGlacierIce(ptrIce);
      ptrHbv = new HBV;
      ptrHbv->SetGeneralPar(ParGeneralStore);
      ptrHbv->SetSubSurfaceHbvPar(&ParSubSurfaceHbvStore[GLACIER_BED]);
      ptrHbv->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]); 
      //      ptrHbv->SetInputTimeSeries(InputTimeSeriesStore);
      ptrHbv->SetInputElement(InputElementStore);
      ptrHbv->SetLandScapeElement(&DistHbv[landIndex]);
      ptrHbv->SetSubSurfaceHbvStore(ParGeneralStore->GetINITIAL_SOIL_MOISTURE(), 
				    ParGeneralStore->GetINITIAL_UPPER_ZONE(), ParGeneralStore->GetINITIAL_LOWER_ZONE());
      ptrGlacier->SetHBV(ptrHbv);
      DistHbv[landIndex].SetGlacier(ptrGlacier);
    }
    if (areaFraction[0]>0.0) {
      ptrHbvAquifer = new HbvAquifer;
      ptrHbvAquifer->SetAreaFraction(areaFraction[0]);
      ptrVegetation = new Vegetation;
      ptrVegetation->SetGeneralPar(ParGeneralStore);
      ptrVegetation->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]); 
      //      ptrVegetation->SetInputTimeSeries(InputTimeSeriesStore);
      ptrVegetation->SetInputElement(InputElementStore);
      ptrVegetation->SetLandScapeElement(&DistHbv[landIndex]);
      ptrHbvAquifer->SetVegetation(ptrVegetation);
      ptrSnow = new Snow;
      ptrSnow->SetGeneralPar(ParGeneralStore);
      ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]); 
      //      ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
      ptrSnow->SetInputElement(InputElementStore);
      ptrSnow->SetLandScapeElement(&DistHbv[landIndex]);
      ptrSnow->SetSnowStore(ParGeneralStore->GetINITIAL_SNOW());
      ptrHbvAquifer->SetSnow(ptrSnow);
      ptrHbv = new HBV;
      ptrHbv->SetGeneralPar(ParGeneralStore);
      ptrHbv->SetSubSurfaceHbvPar(&ParSubSurfaceHbvStore[soilType[0]]);
      ptrHbv->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]); 
      //      ptrHbv->SetInputTimeSeries(InputTimeSeriesStore);
      ptrHbv->SetInputElement(InputElementStore);
      ptrHbv->SetLandScapeElement(&DistHbv[landIndex]);
      ptrHbv->SetSubSurfaceHbvStore(ParGeneralStore->GetINITIAL_SOIL_MOISTURE(), 
                                    ParGeneralStore->GetINITIAL_UPPER_ZONE(), ParGeneralStore->GetINITIAL_LOWER_ZONE());
      ptrHbvAquifer->SetHBV(ptrHbv);
      DistHbv[landIndex].SetHbvAquifer(ptrHbvAquifer);
      j=1;
      while (j<maximumNumberLandClasses && areaFraction[j]>0.0) {
        lastHbvAquifer=ptrHbvAquifer;
        ptrHbvAquifer = new HbvAquifer;
        ptrHbvAquifer->SetAreaFraction(areaFraction[j]);
        ptrVegetation = new Vegetation;
        ptrVegetation->SetGeneralPar(ParGeneralStore);
        ptrVegetation->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]); 
        //        ptrVegetation->SetInputTimeSeries(InputTimeSeriesStore);
        ptrVegetation->SetInputElement(InputElementStore);
        ptrVegetation->SetLandScapeElement(&DistHbv[landIndex]);
        ptrHbvAquifer->SetVegetation(ptrVegetation);
        ptrSnow = new Snow;
        ptrSnow->SetGeneralPar(ParGeneralStore);
        ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]); 
        //        ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
        ptrSnow->SetInputElement(InputElementStore);
        ptrSnow->SetLandScapeElement(&DistHbv[landIndex]);
        ptrSnow->SetSnowStore(ParGeneralStore->GetINITIAL_SNOW());
        ptrHbvAquifer->SetSnow(ptrSnow);
        ptrHbv = new HBV;
        ptrHbv->SetGeneralPar(ParGeneralStore);
        ptrHbv->SetSubSurfaceHbvPar(&ParSubSurfaceHbvStore[soilType[j]]); 
        ptrHbv->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]); 
        //        ptrHbv->SetInputTimeSeries(InputTimeSeriesStore);
        ptrHbv->SetInputElement(InputElementStore);
        ptrHbv->SetLandScapeElement(&DistHbv[landIndex]);
        ptrHbv->SetSubSurfaceHbvStore(ParGeneralStore->GetINITIAL_SOIL_MOISTURE(), 
                                      ParGeneralStore->GetINITIAL_UPPER_ZONE(), ParGeneralStore->GetINITIAL_LOWER_ZONE());
        ptrHbvAquifer->SetHBV(ptrHbv);
        lastHbvAquifer->SetNextHbvAquifer(ptrHbvAquifer);
        j++;
      }
    }
  }

  ptrLake=0;
  ptrGlacier=0;
  ptrHbvAquifer=0;
  lastHbvAquifer=0;
  ptrVegetation=0;
  ptrSnow=0;
  //  ptrIce=0;
  ptrHbv=0;
  finLandScape.close();
  // End read file with landscape elements 

  // Read name of file with observed streamflow data 
  /*  cout << " File with observed streamflow data: ";
  cin >> fileObsStreamflow;
  cout << endl;*/
  fileControl.ignore(100, ':');
  fileControl >> fileObsStreamflow;
  fileControl.ignore(1024, '\n');
  // End read name of file with observed streamflow data

  // Read sub-catchment information file and generate sub-catchment objects  
  /*  cout << " File with sub-catchment hierarchy: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream fileWCo(fileName);
  if (!fileWCo.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  cout << endl << " Sub catchment file opened " << fileName << endl << endl; 
  // Sub-catchment elements
  fileWCo.ignore(100,':');
  fileWCo >> numWatc;
  cout << "\n # Number of sub-catchments " << numWatc << endl;
  SubCatchment *CatchmentElement = new SubCatchment [numWatc];
  for (i=0; i<numWatc; i++) {
    fileWCo >> j >> ch >> subCatchmentId >> correction;
    if (j != i) {
      cout << endl << "Error reading file " << fileName << "\t" << i << "\t" << j << endl;
      exit (1);
    }
    fileWCo.ignore(256,'\n');
    CatchmentElement[i].SetSubCatchmentIndex(i);
    CatchmentElement[i].SetIdentifier(subCatchmentId);
    CatchmentElement[i].SetCorrection(correction);
    CatchmentElement[i].AllocateAccumulatedDischarge(initialTimeSteps+numberTimeSteps);
    CatchmentElement[i].AllocateAccumulatedWaterBalance(initialTimeSteps+numberTimeSteps);
    CatchmentElement[i].AllocateWaterBalance(initialTimeSteps+numberTimeSteps);
    CatchmentElement[i].ObsDataInput(fileObsStreamflow, startSimulationTime, endSimulationTime, numberTimeSteps, 
                                     ParGeneralStore->GetSECONDS_TIMESTEP());
    cout << " Sub-catchment index " << i << "  " << "Sub-catchment identifier " << CatchmentElement[i].GetIdentifier() << "  " << "Sub-catchment correction " << CatchmentElement[i].GetCorrection() << endl;
  }
  // Watercourse outlets
  fileWCo.ignore(100,':');
  fileWCo >> numWatcOut;
  cout << "\n # Number of watercourse outlets " << numWatcOut << endl;
  SubCatchment ** Outlet = new SubCatchment * [numWatcOut];
  for (i=0; i<numWatcOut; i++) {
    fileWCo >> j;
    Outlet[i] = &CatchmentElement[j];
    cout << " Outlet no. " << i << "\t" << " Sub-catchment no. " << j << "\t" << endl;
    fileWCo.ignore(256,'\n');
  }
  // Hierarchy of sub-catchments
  fileWCo.getline(buffer, 256);
  cout << "\n " << buffer << endl;
  while (fileWCo >> i) {
    fileWCo >> numWatcUp;
    CatchmentElement[i].SetNumUpStream(numWatcUp);
    fileWCo.ignore(100,':');
    cout << " Downstream, sub-catchment no.  " << i << "    Identifier  " << CatchmentElement[i].GetIdentifier() << endl;
    cout << " No. of upstream sub-catchments " << numWatcUp << endl;
    k = 0;
    while (fileWCo.peek() != '\n') {
      fileWCo >> j;
      cout << "\t" << "Upstream, sub-catchment no. " << j ;
      CatchmentElement[i].SetUpStream(k, &CatchmentElement[j]);
      cout  << "\t" << "UpStream[" << k << "]" << "    Identifier  " << CatchmentElement[i].GetUpStream(k)->GetIdentifier() << endl;
      while (fileWCo.peek() == ' ') fileWCo.ignore(1,' ');
      k++;
    }
    fileWCo.ignore(256,'\n');
    if (numWatcUp!=k) {
      cout << endl << " Error in number of upstream pointers for sub-catchment no. " << i << endl << endl;
      exit (1);
    } 
  }
  fileWCo.close();
  // End read file with sub-catchment information

  // Read information about sub-catchment elements and landscape elements
  ReadSubCatchmentIdentifier(DistHbv, CatchmentElement, numWatc, fileControl, fout);

  // Write information about sub-catchment outlets and landscape elements to file test_waterland.txt
  //  WriteSubCatchmentIdentifier(CatchmentElement, numWatc, fout);

  // Precipitation and temperature correction for catchments
  int * correctionCatchments = new int[maximumCorrectionCatchments];
  double * correctionPrecipitation = new double[maximumCorrectionCatchments];
  double * correctionTemperature = new double[maximumCorrectionCatchments];
  /*  cout << " File with precipitation and temperature correction for catchments: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  //  cout << fileName << endl;
  ifstream finCorrection(fileName);  // Open for reading
  if (!finCorrection.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finCorrection.getline(buffer, 256);
  i=0;
  fout << "Precipitation and temperature correction for catchments: \n";
  while (finCorrection >> correctionCatchments[i]) {
    finCorrection >> correctionPrecipitation[i] >> correctionTemperature[i];
    fout << i << "  " << correctionCatchments[i] << "  " <<  correctionPrecipitation[i] << "  " <<  correctionTemperature[i] << endl;
    i++;
  }
  numberCorrectionCatchments = i;
  for (i=0; i<numWatcOut; i++) { 
    TraverseCorrectionSubCatchment(Outlet[i], numberCorrectionCatchments, correctionCatchments, 
                                   correctionPrecipitation, correctionTemperature, fout);
  }
  finCorrection.close();
  delete [] correctionCatchments;
  delete [] correctionPrecipitation;
  delete [] correctionTemperature;
  // End precipitation and temperature correction for catchments

  // File with grid cell numbers
  /* cout << " File with grid cell numbers: ";
  cin >> fileName;
  cout << endl; */
  fileControl.ignore(100, ':');
  fileControl >> fileName;
  fileControl.ignore(256, '\n');
  FILE * finGridNumbers;
  if ((finGridNumbers = fopen(fileName, "r")) == NULL) {  // Open for reading
	  printf("\n File %s not found!\n\n", fileName);
	  exit(1);
  } 
  fgets(buffer, 256, finGridNumbers);
  i = 0;
  while (fgets(buffer, 256, finGridNumbers)) {
	  i++;
  }
  numberIndexStore = i;
  fgets(buffer, 256, finGridNumbers);
  cout << endl << " Antall element i indexStore = " << numberIndexStore << endl;
  int * indexStore = new int[numberIndexStore];
  rewind(finGridNumbers);
  fgets(buffer, 256, finGridNumbers);
  i = 0;
  while (fgets(buffer, 256, finGridNumbers)) {
	  sscanf(buffer, "%c %d", &ch, &indexStore[i]);
	  i++;
  }
  if (i != numberIndexStore) {
	  printf("\n Index error file %s, i = %d numberIndexStore = %d !\n\n", fileName, i, numberIndexStore);
	  exit(1);
  }
  fclose(finGridNumbers);
  
  // Read model state variables
  if (readModelStates) { 
    ReadModelStateVariables(DistHbv, numLand);
  }

  // Time series format input data (e.g. from met stations) 
  if (inputFormat == 'T' || inputFormat == 't') {
    // Object for storing input data for time series
    InputTimeSeries * InputTimeSeriesStore = new
      InputTimeSeries(initialTimeSteps+numberTimeSteps, MetStations->GetNumPrecStations()+MetStations->GetNumTempStations(),
                      startModelTime, endSimulationTime, ParGeneralStore->GetSECONDS_TIMESTEP());
    InputTimeSeriesStore->SetGeneralPar(ParGeneralStore);
    //cout << " File with meteorological input data: ";
    //cin >> fileName;
    //cout << endl;
    fileControl.ignore(100,':');
    fileControl >> fileNameInput;
    fileControl.ignore(256,'\n');
    cout << fileNameInput << endl;
    ifstream finInput(fileNameInput);  // Open for reading
    if (!finInput.is_open()) {
      cout << endl << " Error opening file " << fileNameInput << endl << endl;
      exit(1);
    }
    InputTimeSeriesStore->SetInput(finInput);
    finInput.close();
    // Write input data to test file
    //    InputTimeSeriesStore->WriteInput();

    // Water balance for all elements and time steps for for spin-up period and simulation period
    timeStep=0;
    for (datetime=startModelTime; datetime<=endSimulationTime; datetime+=ParGeneralStore->GetSECONDS_TIMESTEP()) {
      //    cout << "  " << datetime.getYear() << "  " << datetime.getMonth() << "  " << datetime.getDay() << "  " 
      //         << datetime.getHour() << "  " << datetime.getMinute() << "  " << datetime.getSecond() << endl;
      inputDataFound=true;
      WaterBalanceTimeSeries(DistHbv, ParGeneralStore, MetStations, initialTimeSteps, numberTimeSteps, 
                             InputTimeSeriesStore, InputElementStore, 
                             numLand, timeStep, datetime, &inputDataFound);
      // Traverse sub-catchments and landscape elements
      if (inputDataFound) {
        for (i=0; i<numWatcOut; i++) { 
          TraverseSubCatchment(Outlet[i], timeStep, fout);
        }
      }
      else {
        for (i=0; i<numWatcOut; i++) { 
          TraverseMissingDataSubCatchment(Outlet[i], timeStep, fout);
        }
      }
      // Write state variables for all landscape elements 
      //      if (timeStep == (int)(initialTimeSteps+numberTimeSteps/2) || timeStep == initialTimeSteps+numberTimeSteps-1) {
      WriteBinaryGrid(DistHbv, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout, writePET);
      //	WriteAsciiGrid(DistHbv, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
      //      }
      SnowGlacierIceReDistribution(Outlet, DistHbv, ParGeneralStore, initialTimeSteps, numberTimeSteps, numLand, numWatcOut, timeStep, datetime, 
				   nRows, nCols, noData, xllCorner, yllCorner, cellSize, modelCalibration, fout); 
      timeStep++;
    }
    if (timeStep != initialTimeSteps+numberTimeSteps) {
      cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
      exit(1);
    }
    // Write discharge from all sub-catchment elements in sub-catchment hierarchy to output files
    WriteSubCatchmentDischarge(CatchmentElement, numWatc, startSimulationTime, endSimulationTime, initialTimeSteps, 
                               numberTimeSteps, ParGeneralStore->GetSECONDS_TIMESTEP(), modelCalibration, fout);
    WriteSubCatchmentWaterBalance(CatchmentElement, numWatc, startSimulationTime, endSimulationTime, 
                                  initialTimeSteps, numberTimeSteps, ParGeneralStore->GetSECONDS_TIMESTEP());
    delete InputTimeSeriesStore;
  }

  // Grid file format input data (e.g. senorge data)
  else { 
    // allocate memory for input daily bil files (should only be necessary  if inputFileFormat == 'B' || inputFileFormat == 'b'  ) 
      unsigned short int * precip10 = new unsigned short int [nRows*nCols];
      unsigned short int * temp10K = new unsigned short int [nRows*nCols];
      unsigned short int * tmax10K = new unsigned short int[nRows*nCols];
      unsigned short int * tmin10K = new unsigned short int[nRows*nCols];
      unsigned short int * wind10 = new unsigned short int[nRows*nCols];
      unsigned short int * solar10 = new unsigned short int[nRows*nCols];
      unsigned short int * vp10 = new unsigned short int[nRows*nCols];
      unsigned short int * hurs = new unsigned short int[nRows*nCols];
      
      // allocate memory for input annual .nc files  (should only be necessary  if inputFileFormat == 'N' || inputFileFormat == 'n' ) 
      prec_in = (float*)calloc(366*nRows*nCols,sizeof(float));
      temp_in = (float*)calloc(366*nRows*nCols,sizeof(float));
      tmax_in = (float*)calloc(366*nRows*nCols,sizeof(float));
      tmin_in = (float*)calloc(366*nRows*nCols,sizeof(float));
      srad_in = (float*)calloc(366*nRows*nCols,sizeof(float));
      vp_in = (float*)calloc(366*nRows*nCols,sizeof(float));
      hurs_in = (float*)calloc(366*nRows*nCols,sizeof(float));
      wind_in = (float*)calloc(366*nRows*nCols,sizeof(float));
    
    // Water balance for all elements and time steps for for spin-up period and simulation period
    timeStep=0;
    oldyear=datetime.getYear();
    for (datetime=startModelTime; datetime<=endSimulationTime; datetime+=ParGeneralStore->GetSECONDS_TIMESTEP()) {
      inputDataFound=true;

      if (inputFileFormat == 'N' || inputFileFormat == 'n') {
	if(datetime==startModelTime ||  datetime.getYear()>oldyear) {
	  //readnetcdf, returner forcings for inneværende år (KSS netcdf filer har data for et år pr fil).
	  if (forcingType == 'S')  //i.e. historiske, senorge .nc filer
	    ReadNetcdf(initialTimeSteps, numberTimeSteps, 
		       numLand,timeStep,nRows,nCols,datetime,evaporationModellingControl,
		       precPath,tmeanPath,tmaxPath,tminPath,windPath,rsdsPath,
		       prec_in,temp_in,tmax_in,tmin_in,wind_in,srad_in,vp_in,hurs_in);
           
           // cout << "  " << datetime.getYear() << "  " << datetime.getMonth() << "  " << datetime.getDay() << "  " 
           //    << datetime.getHour() << "  " << datetime.getMinute() << "  " << datetime.getSecond() << endl;  
	  if (forcingType == 'C')   //i.e. klimadata, KSS-style .nc filer
	    ReadNetcdfClimateProj(initialTimeSteps, numberTimeSteps, 
				  numLand, timeStep, nRows, nCols, datetime,evaporationModellingControl,
				  precPath,tmeanPath,tmaxPath,tminPath,windPath,rsdsPath,hursPath,
				  forcingName,rcpName,biasName,
				  prec_in,temp_in,tmax_in,tmin_in,wind_in,srad_in,vp_in,hurs_in);
	  WaterBalanceGridNetcdf(DistHbv, ParGeneralStore,InputElementStore,initialTimeSteps,numberTimeSteps, 
				 numLand,timeStep,nRows,nCols,datetime,
				 prec_in,temp_in,tmax_in,tmin_in,wind_in,srad_in,vp_in, 
				 &inputDataFound,indexStore,numberIndexStore,fout);
	  oldyear=datetime.getYear();	  
	}
	else { //forcings for this time step already read, go straight to waterbalancegridnetcdf()
	  WaterBalanceGridNetcdf(DistHbv, ParGeneralStore, InputElementStore, initialTimeSteps, numberTimeSteps, 
				 numLand, timeStep, nRows, nCols, datetime,
				 prec_in, temp_in, tmax_in, tmin_in, wind_in, srad_in, vp_in, 
				 &inputDataFound, indexStore, numberIndexStore, fout);

	}
      }
      else   //met input = bil files 
	WaterBalanceGrid(DistHbv, ParGeneralStore, InputElementStore, initialTimeSteps, numberTimeSteps, 
			 numLand, timeStep, nRows, nCols, 
			 datetime, precPath, precip10, temp10K, tmax10K, tmin10K, wind10, solar10, vp10, 
			 &inputDataFound, indexStore, numberIndexStore, fout);
        
       
      // Traverse sub-catchments and landscape elements
      if (inputDataFound) {
	for (i=0; i<numWatcOut; i++) { 
	  TraverseSubCatchment(Outlet[i], timeStep, fout);
	}
      }
      else {
	for (i=0; i<numWatcOut; i++) { 
	  TraverseMissingDataSubCatchment(Outlet[i], timeStep, fout);
	}
      }
      // Write state variables for all landscape elements
      //if (timeStep == (int)(initialTimeSteps+numberTimeSteps/2) || timeStep == initialTimeSteps+numberTimeSteps-1) {
      if (timeStep >= initialTimeSteps) {
	//printf("%d %d\n",timeStep,initialTimeSteps);
	WriteBinaryGrid(DistHbv,datetime,numLand,timeStep,nRows,nCols,noData,xllCorner,yllCorner,cellSize,fout,writePET);
	//WriteAsciiGrid(DistHbv, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
      }
      SnowGlacierIceReDistribution(Outlet, DistHbv, ParGeneralStore, initialTimeSteps, numberTimeSteps, numLand, numWatcOut, timeStep, datetime, 
				   nRows, nCols, noData, xllCorner, yllCorner, cellSize, modelCalibration, fout); 
      timeStep++;
    }
    if (timeStep != initialTimeSteps+numberTimeSteps) {
      cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
      exit(1);
    }
    // Write discharge from all sub-catchment elements in sub-catchment hierarchy to output files
    WriteSubCatchmentDischarge(CatchmentElement, numWatc, startSimulationTime, endSimulationTime, initialTimeSteps, 
			       numberTimeSteps, ParGeneralStore->GetSECONDS_TIMESTEP(), modelCalibration, fout);
    WriteSubCatchmentWaterBalance(CatchmentElement, numWatc, startSimulationTime, endSimulationTime, 
				  initialTimeSteps, numberTimeSteps, ParGeneralStore->GetSECONDS_TIMESTEP());

    //if (inputFileFormat == 'B' || inputFileFormat == 'b') {
      delete [] precip10;
      delete [] temp10K;
      delete [] tmax10K;
      delete [] tmin10K;
      delete [] wind10;
      delete [] solar10;
      delete [] vp10;
      //}
      //if (inputFileFormat == 'N' || inputFileFormat == 'n') {
      free(prec_in);
      free(temp_in);
      free(tmax_in);
      free(tmin_in);
      free(srad_in);
      free(vp_in);
      free(wind_in);
      free(hurs_in);
  }
  
  // Write model state variables
  if (writeModelStates) { 
    WriteModelStateVariables(DistHbv, numLand, datetime, timeStep-1, ParGeneralStore->GetSECONDS_TIMESTEP());
  }

  // Write water balance grid
   WriteAsciiGridWaterBalance(DistHbv, startSimulationTime, endSimulationTime, numLand, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
 
  // Write state variable time series for landscape elements selected for output
  WriteDistributedHbvTimeSeries(DistHbv, numLand, startSimulationTime, endSimulationTime, 
                                initialTimeSteps, numberTimeSteps, ParGeneralStore->GetSECONDS_TIMESTEP());


  fout << endl;
  fout.close();
  fileControl.close();

  delete [] CatchmentElement;
  delete [] Outlet;
  delete [] DistHbv;
  delete ParGeneralStore;
  delete [] ParLandSurfaceStore;
  delete [] ParSubSurfaceHbvStore;
  delete InputElementStore;
  delete [] indexStore;

  return 0;
}


void ReadSubCatchmentIdentifier(DistributedHbv * const DistHbv, SubCatchment * const CatchmentElement, int numWatc, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char ch;
  int i,j;
  int subCatchmentId, numLandScape;
  int landIndex, geoIndex;
  DistributedHbv * thisElement;

  // Read information about sub-catchment elements and landscape elements
  /*  cout << "\n File with information about sub-catchment elements and landscape elements: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finSubCatchment(fileName);  // Open for reading
  if (!finSubCatchment.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  } 
  // Connect landscape elements to sub-catchment elements
  for (i=0; i<numWatc; i++) {
    finSubCatchment >> ch >> subCatchmentId >> ch >> numLandScape;
    finSubCatchment.ignore(256,'\n');
    CatchmentElement[i].SetNumLandScape(numLandScape);
    if (subCatchmentId!=CatchmentElement[i].GetIdentifier()) {
      cout << endl << " Error reading file " << fileName << " for sub-catchment " << i << "\t" 
           << subCatchmentId << endl << endl;
      exit (1);
    }
    if (numLandScape > 0) {
      finSubCatchment >> landIndex >> geoIndex;
      thisElement=&DistHbv[landIndex];
      CatchmentElement[i].SetLandScapeElement(thisElement);
      //    cout << thisElement->GetLandIndex() << "  " << thisElement->GetGeoIndex() << endl;
      for (j=1; j<numLandScape; j++) {
        finSubCatchment >> landIndex >> geoIndex;
        thisElement->SetNextElement(&DistHbv[landIndex]);
        thisElement=&DistHbv[landIndex];
        //      cout << thisElement->GetLandIndex() << "  " << thisElement->GetGeoIndex() << endl;
      }
    }
  }
}


void SnowGlacierIceReDistribution(SubCatchment ** const Outlet, DistributedHbv * const DistHbv, ParametersGeneral * ParGeneralStore, int initialTimeSteps, int numberTimeSteps,
                                  int numLand, int numWatcOut, int timeStep, DateTime datetime, 
                                  int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, bool modelCalibration, ofstream &fout)
{
    int i;
    // Snow store is removed at day no. DAY_SNOW_ZERO
    if (ParGeneralStore->GetDAY_SNOW_ZERO() > 0 &&
            dayNumber(datetime.getYear(), datetime.getMonth(), datetime.getDay()) ==
            ParGeneralStore->GetDAY_SNOW_ZERO() + leapYear(datetime.getYear()))
    {
      //        cout << "            Snow " << datetime.getYear() << " " << datetime.getMonth() << " " << datetime.getDay() << " " << endl;
        for (i = 0; i < numLand; i++)
        {
	  //cout << "  ingjerd snow set to zero.... " << ParGeneralStore->GetDAY_SNOW_ZERO() << endl;
	  DistHbv[i].SetSnowStore(0.0);
        }
    }
}


void WaterBalanceTimeSeries(DistributedHbv * const DistHbv, ParametersGeneral * const ParGeneralStore, 
                            MeteorologicalStations * const MetStations, int initialTimeSteps, int numberTimeSteps, 
                            InputTimeSeries * InputTimeSeriesStore, InputElement * InputElementStore,
                            int numLand, int timeStep, DateTime datetime, bool * inputDataFound)
{
  int i,indSta,indMet,indWeight,indexPrec,indexTemp, dayofyear_PM2;
  double precipitation, temperature, weight, elementPrecipitation, elementTemperature;
  double precGradLow, precGradHigh, lapseRate, lapseRateDry, lapseRateWet;
  double elementElevation, gradElevation, metStationElevation, metStationWeight;
  double rainSnowTemperature;
  gradElevation = ParGeneralStore->GetGRAD_CHANGE_ALT();
  precGradLow = ParGeneralStore->GetPREC_GRAD_LOW();
  precGradHigh = ParGeneralStore->GetPREC_GRAD_HIGH();
  lapseRateDry = ParGeneralStore->GetLAPSE_DRY();
  lapseRateWet = ParGeneralStore->GetLAPSE_WET();
  for (i=0; i<numLand; i++) {
    elementElevation = DistHbv[i].GetElevation();  
    // Precipitation data
    elementPrecipitation = 0.0;
    precipitation = 0.0;
    weight = 0.0;
    for (indexPrec=0; indexPrec<ParGeneralStore->GetNUM_PREC_SERIES(); indexPrec++) {
      indSta = DistHbv[i].GetMetSeriesNumber(indexPrec);
      indWeight = indexPrec;
      indMet = indSta;
      metStationElevation = MetStations->GetStationAltitude(indMet);
      metStationWeight = DistHbv[i].GetMetSeriesWeight(indexPrec);
      precipitation = InputTimeSeriesStore->GetInput(timeStep,indMet);
      if (precipitation >= 0.0) { 
        precipitation = precipitation/1000.0;
        if (gradElevation==0.0 || (elementElevation<=gradElevation && metStationElevation<=gradElevation)) 
          precipitation = (precipitation)*pow(precGradLow,(elementElevation-metStationElevation)/100.0)*metStationWeight;
        else if (elementElevation>gradElevation && metStationElevation<gradElevation)
          precipitation = (precipitation)*pow(precGradLow,(gradElevation-metStationElevation)/100.0)*
            pow(precGradHigh,(elementElevation-gradElevation)/100.0)*metStationWeight;
        else if (elementElevation<gradElevation && metStationElevation>gradElevation)
          precipitation = (precipitation)*pow(precGradLow,(elementElevation-gradElevation)/100.0)*
            pow(precGradHigh,(gradElevation-metStationElevation)/100.0)*metStationWeight;
        else
          precipitation = (precipitation)*pow(precGradHigh,(elementElevation-metStationElevation)/100.0)*metStationWeight;
        elementPrecipitation = elementPrecipitation + precipitation;
        weight = weight + metStationWeight;
      }
      /*      cout << indexPrec << "  " << indSta << "  " << indMet << "  " << metStationWeight << "  " << weight << endl;
              if (indexPrec==ParGeneralStore->GetNUM_PREC_SERIES()-1) cout << endl;*/
    }
    if (weight > 0.0)
      elementPrecipitation = elementPrecipitation/weight;
    else
      elementPrecipitation = missingData;
    //    cout << indexPrec << "  " << metStationElevation << "  " << weight << "  " << elementPrecipitation << endl;
    /*    for (j=0; j<ParGeneralStore->GetNUM_PREC_SERIES(); j++) {
          cout << "P " << DistHbv[i].GetMetSeriesNumber(j) << "  " << DistHbv[i].GetMetSeriesWeight(j) << endl;
          }
          for (j=0; j<ParGeneralStore->GetNUM_TEMP_SERIES(); j++) {
          cout << "T " << DistHbv[i].GetMetSeriesNumber(ParGeneralStore->GetNUM_PREC_SERIES()+j) << 
          "  " << DistHbv[i].GetMetSeriesWeight(ParGeneralStore->GetNUM_PREC_SERIES()+j) << endl;
          }*/
    // Temperature data    
    if (elementPrecipitation == 0.0) 
      lapseRate = lapseRateDry;
    else
      lapseRate = lapseRateWet;
    elementTemperature = 0.0;
    temperature = 0.0;
    weight = 0.0;
    for (indexTemp=0; indexTemp<ParGeneralStore->GetNUM_TEMP_SERIES(); indexTemp++) {
      indSta = DistHbv[i].GetMetSeriesNumber(ParGeneralStore->GetNUM_PREC_SERIES() + indexTemp);
      indWeight = ParGeneralStore->GetNUM_PREC_SERIES() + indexTemp;
      indMet = MetStations->GetNumPrecStations() + indSta;
      metStationElevation = MetStations->GetStationAltitude(indMet);
      metStationWeight = DistHbv[i].GetMetSeriesWeight(indWeight);
      temperature = InputTimeSeriesStore->GetInput(timeStep,indMet); 
      if (temperature > -99.0) { 
        temperature = (temperature/1.0 + lapseRate*(elementElevation-metStationElevation)/100.0)*
          metStationWeight;
        elementTemperature = elementTemperature + temperature; 
        weight = weight + metStationWeight;
      }
      /*      cout << indexTemp << "  " << indSta << "  " << indMet << "  " << metStationWeight << "  " << weight << endl;
              if (indexTemp==ParGeneralStore->GetNUM_TEMP_SERIES()-1) cout << endl;*/
    }
    if (weight > 0.0)
      elementTemperature = elementTemperature/weight;
    else
      elementTemperature = missingData;
    if (elementPrecipitation > 0.0) {
      rainSnowTemperature = elementTemperature + 
        lapseRate*(elementElevation-DistHbv[i].GetPrecStationsWeightedElevation())/100.0;
      //        lapseRate*(DistHbv[i].GetTempStationsWeightedElevation()-DistHbv[i].GetPrecStationsWeightedElevation())/100.0;
      if (rainSnowTemperature >= 0.0) 
        elementPrecipitation = elementPrecipitation*ParGeneralStore->GetPREC_CORR_RAIN();
      else 
        elementPrecipitation = elementPrecipitation*ParGeneralStore->GetPREC_CORR_RAIN()*ParGeneralStore->GetPREC_CORR_SNOW();
    }
    if (elementPrecipitation > missingData && elementTemperature > missingData) {
      InputElementStore->SetInput(0,elementPrecipitation*DistHbv[i].GetPrecipitationCorrection());
      InputElementStore->SetInput(1,elementTemperature+DistHbv[i].GetTemperatureCorrection());

      dayofyear_PM2 = dayNumber(InputTimeSeriesStore->GetDateTime(timeStep).getYear(),
				InputTimeSeriesStore->GetDateTime(timeStep).getMonth(),
				InputTimeSeriesStore->GetDateTime(timeStep).getDay());
      DistHbv[i].WaterBalance(timeStep,datetime,initialTimeSteps,numberTimeSteps,dayofyear_PM2);

    }
    else {
      *inputDataFound=false;
    }
  }
}

void WaterBalanceGrid(DistributedHbv * DistHbv,  ParametersGeneral * ParGeneralStore, InputElement * InputElementStore,
                      int initialTimeSteps, int numberTimeSteps, 
                      int numLand, int timeStep, int nRows, int nCols, DateTime datetime, char * precPath,
                      unsigned short int * precip10, unsigned short int * temp10K, unsigned short int * tmax10K,
		      unsigned short int * tmin10K, unsigned short int * wind10, unsigned short int * solar10,
		      unsigned short int * vp10, bool * inputDataFound,
                      int * indexStore, int numberIndexStore, ofstream &fout)
{
  ifstream filePrec, fileTemp, fileTmax, fileTmin, fileSolar, fileWind, fileVP;
  int i,j,k, jNew, dayofyear_PM2;
  unsigned short int metMissing = 10000;
  double precipitation, preci_corr, temperature, Tmax, Tmin, wind, solarRadiation, vp;
  char fileName[100];
  char precFileName[100];
  char tempFileName[100];
  char TmaxFileName[100];
  char TminFileName[100];
  char windFileName[100];
  char solarFileName[100];
  char VPFileName[100];
  char hydYear[5];

  strcpy(precFileName,precPath); //same path for all variables when using .bil input files.
  strcpy(tempFileName,precPath);
  strcpy(TmaxFileName,precPath);
  strcpy(TminFileName,precPath);
  strcpy(windFileName,precPath);
  strcpy(solarFileName,precPath);
  strcpy(VPFileName,precPath);

  strcat(precFileName,"/rr/");
  strcat(tempFileName,"/tm/");
  strcat(TmaxFileName, "/tmax/");
  strcat(TminFileName, "/tmin/");
  strcat(windFileName, "/wind/");
  strcat(solarFileName, "/srad/");
  strcat(VPFileName, "/vp/");
  
  sprintf(hydYear,"%04d",datetime.getYear());
  strcat(precFileName,hydYear);
  strcat(tempFileName,hydYear);
  strcat(TmaxFileName, hydYear);
  strcat(TminFileName, hydYear);
  strcat(windFileName, hydYear);
  strcat(solarFileName, hydYear);
  strcat(VPFileName, hydYear);

  sprintf(fileName,"/tm_%04d_%02d_%02d.bil",datetime.getYear(),datetime.getMonth(),datetime.getDay());
  strcat(tempFileName,fileName);
  sprintf(fileName,"/rr_%04d_%02d_%02d.bil",datetime.getYear(),datetime.getMonth(),datetime.getDay());
  strcat(precFileName,fileName);

  sprintf(fileName, "/tmax_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());
  strcat(TmaxFileName, fileName);

  sprintf(fileName, "/tmin_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());
  strcat(TminFileName, fileName);

  sprintf(fileName, "/wind_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());
  strcat(windFileName, fileName);

  sprintf(fileName, "/srad_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());
  strcat(solarFileName, fileName);

  sprintf(fileName, "/vp_%04d_%02d_%02d.bil", datetime.getYear(), datetime.getMonth(), datetime.getDay());
  strcat(VPFileName, fileName);

  filePrec.open(precFileName,ios::in | ios::binary);
  if (!filePrec.is_open()) {
    cout << endl << "Error opening file " << precFileName << endl << endl;
    exit (1);
  }

  fileTemp.open(tempFileName,ios::in | ios::binary);
  if (!fileTemp.is_open()) {
    cout << endl << "Error opening file " << tempFileName << endl << endl;
    exit (1);
  }

  fileTmax.open(TmaxFileName, ios::in | ios::binary);
  if (!fileTmax.is_open()) {
	  cout << endl << "Error opening file " << TmaxFileName << endl << endl;
	  exit(1);
  }

  fileTmin.open(TminFileName, ios::in | ios::binary);
  if (!fileTmin.is_open()) {
	  cout << endl << "Error opening file " << TminFileName << endl << endl;
	  exit(1);
  }

  fileWind.open(windFileName, ios::in | ios::binary);
  if (!fileWind.is_open()) {
	  cout << endl << "Error opening file " << windFileName << endl << endl;
	  exit(1);
  }

  fileSolar.open(solarFileName, ios::in | ios::binary);
  if (!fileSolar.is_open()) {
	  cout << endl << "Error opening file " << solarFileName << endl << endl;
	  exit(1);
  }

  fileVP.open(VPFileName, ios::in | ios::binary);
  if (!fileVP.is_open()) {
	  cout << endl << "Error opening file " << VPFileName << endl << endl;
	  exit(1);
  }
  
  dayofyear_PM2 = dayNumber(datetime.getYear(), datetime.getMonth(), datetime.getDay());

  streamoff newPosition;
  k=0;
  for (i = 0; i<numberIndexStore; i++) {
     if (k<numLand) {
       //  if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
        while (DistHbv[k].GetGeoIndex()<indexStore[i]) {
	    cout << " element not found " << i << " " << indexStore[i] << "    " << k << " " << DistHbv[k].GetGeoIndex() << endl;
	    k++;
	}
	if(dayofyear_PM2==244) {
	  cout << "ingjerd bil i " << i << " indexStore(i) " << indexStore[i] << endl;
	  cout << "ingjerd bil k "  << k << " getgeo(k) " << DistHbv[k].GetGeoIndex() << endl;
	}

	if (DistHbv[k].GetGeoIndex() == indexStore[i]) {
          //newPosition = ELEMENT(i,j)*(sizeof(unsigned short int));
	  newPosition = i*(sizeof(unsigned short int));
          filePrec.seekg(newPosition, ios::beg);
          fileTemp.seekg(newPosition, ios::beg);
	  fileTmax.seekg(newPosition, ios::beg);
	  fileTmin.seekg(newPosition, ios::beg);
	  fileWind.seekg(newPosition, ios::beg);
	  fileSolar.seekg(newPosition, ios::beg);
	  fileVP.seekg(newPosition, ios::beg);
          //          filePrec.read((unsigned short int *) &(precip10[ELEMENT(i,j)]), sizeof (unsigned short int));
          //          fileTemp.read((unsigned short int *) &(temp10K[ELEMENT(i,j)]), sizeof (unsigned short int));
	  /*          filePrec.read(reinterpret_cast<char *> (&(precip10[ELEMENT(i,j)])), sizeof (unsigned short int));
          fileTemp.read(reinterpret_cast<char *> (&(temp10K[ELEMENT(i,j)])), sizeof (unsigned short int));
		  fileTmax.read(reinterpret_cast<char *> (&(tmax10K[ELEMENT(i, j)])), sizeof(unsigned short int));
		  fileTmin.read(reinterpret_cast<char *> (&(tmin10K[ELEMENT(i, j)])), sizeof(unsigned short int));
		  fileWind.read(reinterpret_cast<char *> (&(wind10[ELEMENT(i, j)])), sizeof(unsigned short int));
		  fileSolar.read(reinterpret_cast<char *> (&(solar10[ELEMENT(i, j)])), sizeof(unsigned short int));
		  fileMaxRh.read(reinterpret_cast<char *> (&(maxrh10[ELEMENT(i, j)])), sizeof(unsigned short int));
		  fileMinRh.read(reinterpret_cast<char *> (&(minrh10[ELEMENT(i, j)])), sizeof(unsigned short int));*/

		  filePrec.read(reinterpret_cast<char *> (&(precip10[indexStore[i]])), sizeof(unsigned short int));
		  fileTemp.read(reinterpret_cast<char *> (&(temp10K[indexStore[i]])), sizeof(unsigned short int));
		  fileTmax.read(reinterpret_cast<char *> (&(tmax10K[indexStore[i]])), sizeof(unsigned short int));
		  fileTmin.read(reinterpret_cast<char *> (&(tmin10K[indexStore[i]])), sizeof(unsigned short int));
		  fileWind.read(reinterpret_cast<char *> (&(wind10[indexStore[i]])), sizeof(unsigned short int));
		  fileSolar.read(reinterpret_cast<char *> (&(solar10[indexStore[i]])), sizeof(unsigned short int));
		  fileVP.read(reinterpret_cast<char *> (&(vp10[indexStore[i]])), sizeof(unsigned short int));
          //      printf("%d  %hu  %hu\n",ELEMENT(i,j),precip10[ELEMENT(i,j)],temp10K[ELEMENT(i,j)]);
          /*if (precip10[ELEMENT(i,j)]<metMissing && temp10K[ELEMENT(i,j)]<metMissing) {
            *inputDataFound=true;
            precipitation = DistHbv[k].GetPrecipitationCorrection()*(double)precip10[ELEMENT(i,j)]/10000.0;
            temperature = DistHbv[k].GetTemperatureCorrection()+(double)(temp10K[ELEMENT(i,j)]-2731)/10.0;
			Tmax = DistHbv[k].GetTemperatureCorrection() + (double)(tmax10K[ELEMENT(i, j)] - 2731) / 10.0;
			Tmin = DistHbv[k].GetTemperatureCorrection() + (double)(tmin10K[ELEMENT(i, j)] - 2731) / 10.0;
			wind = (double)wind10[ELEMENT(i, j)]/ 10.0;
			solarRadiation = (double)solar10[ELEMENT(i, j)]/ 10.0;
			max_humidity = (double)maxrh10[ELEMENT(i, j)]/ 10.0;
			min_humidity =(double)minrh10[ELEMENT(i, j)]/ 10.0;*/

		  if (precip10[indexStore[i]]<metMissing && temp10K[indexStore[i]]<metMissing) {
			  *inputDataFound = true;
			  precipitation = DistHbv[k].GetPrecipitationCorrection()*(double)precip10[indexStore[i]] / 10000.0;
			  temperature = DistHbv[k].GetTemperatureCorrection() + (double)(temp10K[indexStore[i]] - 2731) / 10.0;
			  Tmax = DistHbv[k].GetTemperatureCorrection() + (double)(tmax10K[indexStore[i]] - 2731) / 10.0;
			  Tmin = DistHbv[k].GetTemperatureCorrection() + (double)(tmin10K[indexStore[i]] - 2731) / 10.0;
			  wind = (double)wind10[indexStore[i]] / 10.0;
			  solarRadiation = (double)solar10[indexStore[i]] / 10.0;
			  vp = (double)vp10[indexStore[i]] / 1000.0/10.0;   // unit: K Pa
			  //      printf("%d  %f  %f\n",ELEMENT(i,j),precipitation,temperature);

			  //preci_corr = precipitation/((0.82-(0.81*exp((temperature-0.66)/1.07))/(1+exp((temperature-0.66)/1.07)))*exp(-pow((wind/4.24),1.81))
                          //        +(0.81*exp((temperature-0.66)/1.07))/(1+exp((temperature-0.66)/1.07))+0.18);

			  if(temperature < 0.5) {
			    preci_corr = ParGeneralStore->GetPREC_CORR_SNOW()*precipitation ;
			  } else {
			    preci_corr = ParGeneralStore->GetPREC_CORR_RAIN()*precipitation ;
			  }


			  InputElementStore->SetInput(0, preci_corr);
			  InputElementStore->SetInput(1, temperature);
			  InputElementStore->SetInput(2, Tmax);
			  InputElementStore->SetInput(3, Tmin);
			  InputElementStore->SetInput(4, wind);
			  InputElementStore->SetInput(5, solarRadiation);
			  InputElementStore->SetInput(6, vp);
			  //            InputElementStore->SetInput(0,precipitation);
			  //            InputElementStore->SetInput(1,temperature);
		  }
		  /*else {
            jNew=j;
            while ((precip10[ELEMENT(i,jNew)]>=metMissing || temp10K[ELEMENT(i,jNew)]>=metMissing) && jNew>0) {
              jNew--;
              newPosition = ELEMENT(i,jNew)*(sizeof(unsigned short int));
              filePrec.seekg(newPosition, ios::beg);
              fileTemp.seekg(newPosition, ios::beg);
              filePrec.read(reinterpret_cast<char *> (&(precip10[ELEMENT(i,jNew)])), sizeof (unsigned short int));
              fileTemp.read(reinterpret_cast<char *> (&(temp10K[ELEMENT(i,jNew)])), sizeof (unsigned short int));
            }
            if (precip10[ELEMENT(i,jNew)]<metMissing && temp10K[ELEMENT(i,jNew)]<metMissing) {
              *inputDataFound=true;
              precipitation = DistHbv[k].GetPrecipitationCorrection()*(double)precip10[ELEMENT(i,jNew)]/10000.0;
              temperature = DistHbv[k].GetTemperatureCorrection()+(double)(temp10K[ELEMENT(i,jNew)]-2731)/10.0;
              //      printf("%d  %f  %f\n",ELEMENT(i,jNew),precipitation,temperature);
              InputElementStore->SetInput(0,precipitation);
              InputElementStore->SetInput(1,temperature);
	      }*/
		  else {
              //      cout << endl << " Missing meterological data for: " << endl;
              //      cout << precFileName << " or " << tempFileName << endl;
              //      cout << " row = " << i << "  col = " << j << "  element no. = " << ELEMENT(i,j) << endl;
              //      printf("  Precipitation %hu  Temperature %hu\n",precip10,temp10K);
              //      cout << endl << endl;
              //      exit (1);
              //      *inputDataFound=false;
		    *inputDataFound=true;
              //InputElementStore->SetInput(0,0.0);
              //InputElementStore->SetInput(1,5.0);
		    if (InputElementStore->GetInput(0) == missingData) InputElementStore->SetInput(0, 0.0);
		    if (InputElementStore->GetInput(1) == missingData) InputElementStore->SetInput(1, 5.0);
		    if (InputElementStore->GetInput(2) == missingData) InputElementStore->SetInput(2, 10.0);
		    if (InputElementStore->GetInput(3) == missingData) InputElementStore->SetInput(3, 0.0);
		    if (InputElementStore->GetInput(4) == missingData) InputElementStore->SetInput(4, 2.0);
		    if (InputElementStore->GetInput(5) == missingData) InputElementStore->SetInput(5, 20.0);
		    if (InputElementStore->GetInput(6) == missingData) InputElementStore->SetInput(6, 800.0);
            // }
		  }
          // Snow store is removed at day no. DAY_SNOW_ZERO //iha removed, was commented out anyway
		  
		  DistHbv[k].WaterBalance(timeStep,datetime,initialTimeSteps,numberTimeSteps,dayofyear_PM2);
		  k++;
        }
      //}
     }
  }
  filePrec.close();
  fileTemp.close();
  fileTmax.close();
  fileTmin.close();
  fileSolar.close();
  fileWind.close();
  fileVP.close();
}

void WaterBalanceGridNetcdf(DistributedHbv * DistHbv,  ParametersGeneral * ParGeneralStore,
			    InputElement * InputElementStore,
			    int initialTimeSteps, int numberTimeSteps, 
			    int numLand, int timeStep, int nRows, int nCols, DateTime datetime,
			    float *prec_in, float *temp_in, float *tmax_in, 
			    float *tmin_in, float *wind_in, float *srad_in, 
			    float *vp_in, bool * inputDataFound,
			    int * indexStore, int numberIndexStore, ofstream &fout)
{
  int i,j,k,dayofyear,iday,icell;
  int metMissing = -998;
  double precipitation, temperature, Tmax, Tmin, wind, solarRadiation, vp;
  char fileName[100];
  char hydYear[5];

  iday = dayNumber(datetime.getYear(), datetime.getMonth(), datetime.getDay()-1); //iday to use when finding todays forcings
  dayofyear=iday;
  //if(iday==244) 
  //  printf("ingjerd waterbalancegridnetcdf %d %d %d prec today %.2f %.2f %.2f %.2f %.2f\n",
  //	   datetime.getYear(),datetime.getMonth(),datetime.getDay(),prec_in[iday*nRows*nCols+72748],temp_in[iday*nRows*nCols+72748],tmax_in[iday*nRows*nCols+72748],tmin_in[iday*nRows*nCols+72748],srad_in[iday*nRows*nCols+72748]);
 
  k=0;
  for (i = 0; i<numberIndexStore; i++) {
    if (k<numLand) {
      while (DistHbv[k].GetGeoIndex()<indexStore[i]) {
	cout << " 2 element not found i " << i << " indexStore " << indexStore[i] << " k " << k << " getgeo " << DistHbv[k].GetGeoIndex() << endl;
	k++;
       }
       if (DistHbv[k].GetGeoIndex() == indexStore[i]) {
  	 //cout << " 2 element found i " << i << " indexStore " << indexStore[i] << " k " << k << " getgeo " << DistHbv[k].GetGeoIndex() << endl;
	 icell=iday*nRows*nCols+indexStore[i];
	 precipitation = DistHbv[k].GetPrecipitationCorrection()*prec_in[icell]/1000.; // prec from mm to m
	 temperature = DistHbv[k].GetTemperatureCorrection()+temp_in[icell]; // temperatures in C
	 Tmax = DistHbv[k].GetTemperatureCorrection()+tmax_in[icell]; 
	 Tmin = DistHbv[k].GetTemperatureCorrection()+tmin_in[icell];
	 solarRadiation = srad_in[icell]*ParGeneralStore->GetSECONDS_TIMESTEP()/1e6; // W/m2-> MJ/m2/timestep   
	 vp = vp_in[icell]/1000.; // input: Pa. HBV: uses kPa.
	 wind = wind_in[icell]; // m/s

	 //sanity check
	 if (precipitation < 0 ) precipitation=0.0;
	 if (temperature < metMissing) temperature=5.0;
	 if (Tmax < metMissing) {
	   //cout << " Tmax " << i << " " << indexStore[i] << " tmax  " << Tmax << " iday " << iday << endl;
	  Tmax=temperature+5;
	  //cout << " Tmax " << i << " " << indexStore[i] << " tmax  " << Tmax << endl;
	 }
	 if (Tmin < metMissing) Tmin=temperature-5;
	 if (wind < 0 || wind > 50) {
	   wind=1.;
	 }
	 if (solarRadiation < 0) solarRadiation = 2.0; //var 20
	 if (vp < 0) vp=0.8; // ingjerd: denne var opprinnelig 800, men skal vel være i kPa?

	 *inputDataFound = true;
	 InputElementStore->SetInput(0, precipitation);
	 InputElementStore->SetInput(1, temperature);
	 InputElementStore->SetInput(2, Tmax);
	 InputElementStore->SetInput(3, Tmin);
	 InputElementStore->SetInput(4, wind);
	 InputElementStore->SetInput(5, solarRadiation);
	 InputElementStore->SetInput(6, vp);

	 DistHbv[k].WaterBalance(timeStep,datetime,initialTimeSteps,numberTimeSteps,dayofyear);
	 k++;
       }
     }
  }
}

void ReadNetcdf(int initialTimeSteps, int numberTimeSteps, 
		int numLand, int timeStep, int nRows, int nCols, DateTime datetime, 
		char ETscheme,char * precPath,char *tmeanPath,char *tmaxPath,
		char *tminPath,char *windPath,char *rsdsPath,
		float *prec_in, float *temp_in, float *tmax_in,
		float *tmin_in, float *wind_in, float *srad_in,
		float *vp_in,float *hurs_in)
{
  ifstream filePrec, fileTemp, fileTmax, fileTmin, fileSolar, fileWind, fileVP;
  float dummy;
  int i,j,k,dayofyear;
  int ncprecin,nctminin,nctmeanin,nctmaxin;
  int ncsradin,nclradin,ncvpin,ncwindin,nchursin;
  int status;
  int precidin,tempidin,tmaxidin,tminidin,windidin,sradidin,lradidin,vpidin,hursidin;
  int flag,cellid;
  int ndays;
  unsigned short int metMissing = 10000;
  double precipitation, temperature, Tmax, Tmin, wind, solarRadiation, vp, hurs;
  char fileName[100];
  char hydYear[5];
  char metPathFileName[100];
  char metPath2FileName[100];
  char metPath3FileName[100];
  char precFileName[200];
  char tmeanFileName[200];
  char tmaxFileName[200];
  char tminFileName[200];
  char sradFileName[200];
  char vpFileName[200];
  char windFileName[200];  
  char hursFileName[200];
  char VarName[NPARAM][20] = { "rr", "tg", "tx", "tn" }; //senorge variable names (nc files)
  static size_t start[]={0,0,0};
  static size_t count[]={366,1550,1195}; //hm. får advarsel når skriver nrows og ncols her. typen stemmer ikke?

  ndays=0;
  //fprintf(stderr,"%d %d\n",nRows,nCols);
  sprintf(hydYear,"%04d",datetime.getYear());
  strcpy(precFileName,precPath);
  strcat(precFileName,hydYear);
  strcat(precFileName,".nc");

  sprintf(hydYear,"%04d",datetime.getYear());
  strcpy(tmeanFileName,tmeanPath);
  strcat(tmeanFileName,hydYear);
  strcat(tmeanFileName,".nc");

  sprintf(hydYear,"%04d",datetime.getYear());
  strcpy(tmaxFileName,tmaxPath);
  strcat(tmaxFileName,hydYear);
  strcat(tmaxFileName,".nc");

  sprintf(hydYear,"%04d",datetime.getYear());
  strcpy(tminFileName,tminPath);
  strcat(tminFileName,hydYear);
  strcat(tminFileName,".nc");

  sprintf(hydYear,"%04d",datetime.getYear());
  strcpy(windFileName,windPath);
  strcat(windFileName,hydYear);
  strcat(windFileName,".nc4");
 
  sprintf(hydYear,"%04d",datetime.getYear());
  strcpy(sradFileName,rsdsPath);
  strcat(sradFileName,"rsds_daily_");
  strcat(sradFileName,hydYear);
  strcat(sradFileName,".nc4");

  sprintf(hydYear,"%04d",datetime.getYear());
  strcpy(vpFileName,rsdsPath);
  strcat(vpFileName,"vp_daily_");
  strcat(vpFileName,hydYear);
  strcat(vpFileName,".nc4");

  sprintf(hydYear,"%04d",datetime.getYear());
  strcpy(hursFileName,rsdsPath);
  strcat(hursFileName,"hurs_daily_");
  strcat(hursFileName,hydYear);
  strcat(hursFileName,".nc4");

  dayofyear = dayNumber(datetime.getYear(), datetime.getMonth(), datetime.getDay());
  ndays=dayNumber(datetime.getYear(),12,31);
  printf("ingjerd readnetcdf %d %d %d DoY %d, ET scheme: %c\n",
	 datetime.getYear(),datetime.getMonth(),datetime.getDay(),dayofyear,ETscheme);

  count[0]=ndays;

  /* Open and read metfiles */
  printf("readnetcdf prec file %s\n",precFileName);
  status = nc_open(precFileName,NC_NOWRITE,&ncprecin); HandleError(status);
  status = nc_inq_varid(ncprecin,VarName[0],&precidin); HandleError(status);
  status = nc_get_vara_float(ncprecin,precidin,start,count,prec_in); HandleError(status);
  nc_close(ncprecin);
  
  printf("readnetcdf tmean file %s\n",tmeanFileName);
  status = nc_open(tmeanFileName,NC_NOWRITE,&nctmeanin); HandleError(status);
  status = nc_inq_varid(nctmeanin,VarName[1],&tempidin); HandleError(status);
  status = nc_get_vara_float(nctmeanin,tempidin,start,count,temp_in); HandleError(status);
  status = nc_close(nctmeanin);
  

 
  if (ETscheme == 'P') { //etscheme = penman monteith
    printf("readnetcdf tmax file %s\n",tmaxFileName);
    status = nc_open(tminFileName,NC_NOWRITE,&nctmaxin); HandleError(status);
    status = nc_inq_varid(nctmaxin,VarName[2],&tmaxidin); HandleError(status);
    status = nc_get_vara_float(nctmaxin,tmaxidin,start,count,tmax_in); HandleError(status);
    status = nc_close(nctmaxin); HandleError(status);

    printf("readnetcdf tmin file %s\n",tminFileName);
    status = nc_open(tminFileName,NC_NOWRITE,&nctminin); HandleError(status);
    status = nc_inq_varid(nctminin,VarName[3],&tminidin); HandleError(status);
    status = nc_get_vara_float(nctminin,tminidin,start,count,tmin_in); HandleError(status);
    status = nc_close(nctminin); HandleError(status);

    printf("readnetcdf sradfile %s\n",sradFileName);
    status = nc_open(sradFileName,NC_NOWRITE,&ncsradin); HandleError(status);
    status = nc_inq_varid(ncsradin,"rsds",&sradidin); HandleError(status);
    status = nc_get_vara_float(ncsradin,sradidin,start,count,srad_in); HandleError(status);
    status = nc_close(ncsradin); HandleError(status);

  /*  printf("readnetcdf %s\n",vpFileName);
    status = nc_open(vpFileName,NC_NOWRITE,&ncvpin); HandleError(status);
    status = nc_inq_varid(ncvpin,"vp",&vpidin); HandleError(status);
    status = nc_get_vara_float(ncvpin,vpidin,start,count,vp_in); HandleError(status);
    status = nc_close(ncvpin); HandleError(status);   */

    printf("readnetcdf windfile %s\n",windFileName);
    status = nc_open(windFileName,NC_NOWRITE,&ncwindin); HandleError(status);
    status = nc_inq_varid(ncwindin,"windspeed_10m",&windidin); HandleError(status);
    status = nc_get_vara_float(ncwindin,windidin,start,count,wind_in); HandleError(status);
    status = nc_close(ncwindin); HandleError(status);

    printf("readnetcdf hursfile %s\n",hursFileName);
    status = nc_open(hursFileName,NC_NOWRITE,&nchursin); HandleError(status);
    status = nc_inq_varid(nchursin,"hurs",&hursidin); HandleError(status);
    status = nc_get_vara_float(nchursin,hursidin,start,count,hurs_in); HandleError(status);
    status = nc_close(nchursin); HandleError(status);   
  }

  if (ETscheme == 'P') {
    for(k=0;k<ndays*(nRows*nCols);k++) {
      vp_in[k]=hurs_in[k]*(6.11*exp((17.27*temp_in[k])/(237.3+temp_in[k]))); //hurs=percent,->vp_in=Pa (nb! later divided by 1000!)
    }
  }
  
   /* flag=0;
  for(k=0;k<ndays;k++) {
    cellid=0; //cellid counter
    for(i=0;i<nRows;i++) {
      for(j=0;j<nCols;j++) {
	if(k<=180 && cellid==1573912) {
	  dummy=hurs_in[flag]*(6.11*exp((17.27*temp_in[flag])/(237.3+temp_in[flag]))); //Pa
	  printf("ingjerd readnetcdf day %d row %d col %d cellnr %d \t p %.1f t %.2f srad %.1f vp %.1f vp-est %.3f hurs %.1f\n",
		 k,i,j,i*nCols+j,prec_in[flag],temp_in[flag],srad_in[flag],vp_in[flag],dummy,hurs_in[flag]);
	}
	flag+=1;
	cellid+=1;
      }
    }
    }*/
    //printf("readnetcdf file ready\n");
}

void ReadNetcdfClimateProj(int initialTimeSteps, int numberTimeSteps, 
			   int numLand, int timeStep, int nRows, int nCols, DateTime datetime, 
			   char ETscheme,char *precPath,char *tmeanPath,char *tmaxPath,
			   char *tminPath,char *windPath,char *rsdsPath,char *hursPath,
			   char *forcingName,char *rcpName,char *biasName,
			   float *prec_in, float *temp_in, float *tmax_in,
			   float *tmin_in, float *wind_in, float *srad_in,
			   float *vp_in,float *hurs_in)
{
  ifstream filePrec, fileTemp, fileTmax, fileTmin, fileSolar, fileWind, fileHurs;
  int i,j,k,dayofyear;
  int ncprecin,nctminin,nctmeanin,nctmaxin;
  int ncsradin,nclradin,nchursin,ncpressin,ncwindin;
  int status;
  int precidin,tempidin,tmaxidin,tminidin,windidin,sradidin,lradidin,hursidin,pressidin;
  int flag,cellid;
  int ndays;
  double precipitation, temperature, Tmax, Tmin, wind,solarRadiation,hurs;
  char fileName[100];
  char hydYear[5];
  char precFileName[200];
  char tmeanFileName[200];
  char tmaxFileName[200];
  char tminFileName[200];
  char sradFileName[200];
  char hursFileName[200];
  char windFileName[200];
  char projfix[100];
  int ThisYear;
  static size_t start[]={0,0,0};
  static size_t count[]={366,1550,1195}; //hm. får advarsel når skriver nrows og ncols her. typen stemmer ikke?

  ndays=0; 
  sprintf(hydYear,"%04d",datetime.getYear());
  ThisYear=(int)datetime.getYear();
  
  if(ThisYear<=2014) {
    strcpy(projfix,"hist/");
    strcat(projfix,forcingName);
    strcat(projfix,"_hist_");
    strcat(projfix,biasName);
    printf("\nprojfix: %s\n",projfix);
  }
  else  {
    strcpy(projfix,rcpName);
    strcat(projfix,"/");
    strcat(projfix,forcingName);
    strcat(projfix,"_");
    strcat(projfix,rcpName);
    strcat(projfix,"_");
    strcat(projfix,biasName);
    printf("\nprojfix: %s\n\n",projfix);
  }
    
  strcpy(precFileName,precPath);
  strcat(precFileName,forcingName);
  strcat(precFileName,"/pr/");
  strcat(precFileName,projfix);
  strcat(precFileName,"-sn2018v2005_rawbc_norway_1km_pr_daily_");
  strcat(precFileName,hydYear);
  strcat(precFileName,".nc4");
  
  strcpy(tmeanFileName,tmeanPath);
  strcat(tmeanFileName,forcingName);
  strcat(tmeanFileName,"/tas/");
  strcat(tmeanFileName,projfix);
  strcat(tmeanFileName,"-sn2018v2005_rawbc_norway_1km_tas_daily_");  
  strcat(tmeanFileName,hydYear);
  strcat(tmeanFileName,".nc4");

  strcpy(tmaxFileName,tmaxPath);
  strcat(tmaxFileName,forcingName);
  strcat(tmaxFileName,"/tasmax/");
  strcat(tmaxFileName,projfix);
  strcat(tmaxFileName,"-sn2018v2005_rawbc_norway_1km_tasmax_daily_"); 
  strcat(tmaxFileName,hydYear);
  strcat(tmaxFileName,".nc4");

  strcpy(tminFileName,tminPath);
  strcat(tminFileName,forcingName);
  strcat(tminFileName,"/tasmin/");
  strcat(tminFileName,projfix);
  strcat(tminFileName,"-sn2018v2005_rawbc_norway_1km_tasmin_daily_"); 
  strcat(tminFileName,hydYear);
  strcat(tminFileName,".nc4");

  strcpy(windFileName,windPath);
  strcat(windFileName,forcingName);
  strcat(windFileName,"/sfcWind/");
  strcat(windFileName,projfix);
  strcat(windFileName,"-klinogrid1612_rawbc_norway_1km_sfcWind_daily_"); 
  strcat(windFileName,hydYear);
  strcat(windFileName,".nc4");
 
  strcpy(sradFileName,rsdsPath);
  strcat(sradFileName,forcingName);
  strcat(sradFileName,"/rsds/");
  strcat(sradFileName,projfix);
  strcat(sradFileName,"-hysn2018v2005era5_rawbc_norway_1km_rsds_daily_"); 
  strcat(sradFileName,hydYear);
  strcat(sradFileName,".nc4");

  strcpy(hursFileName,hursPath);
  strcat(hursFileName,forcingName);
  strcat(hursFileName,"/hurs/");
  strcat(hursFileName,projfix);
  strcat(hursFileName,"-hysn2018v2005era5_rawbc_norway_1km_hurs_daily_"); 
  strcat(hursFileName,hydYear);
  strcat(hursFileName,".nc4");

  dayofyear = dayNumber(datetime.getYear(), datetime.getMonth(), datetime.getDay());
  ndays=dayNumber(datetime.getYear(),12,31);
  printf("ingjerd readnetcdf %d %d %d DoY %d, ET scheme: %c\n",
	 datetime.getYear(),datetime.getMonth(),datetime.getDay(),dayofyear,ETscheme);

  count[0]=ndays;

  /* Open and read metfiles */
  printf("readnetcdf prec file %s\n",precFileName);
  status = nc_open(precFileName,NC_NOWRITE,&ncprecin); HandleError(status);
  status = nc_inq_varid(ncprecin,"pr",&precidin); HandleError(status);
  status = nc_get_vara_float(ncprecin,precidin,start,count,prec_in); HandleError(status);
  nc_close(ncprecin);
  
  printf("readnetcdf tmean file %s\n",tmeanFileName);
  status = nc_open(tmeanFileName,NC_NOWRITE,&nctmeanin); HandleError(status);
  status = nc_inq_varid(nctmeanin,"tas",&tempidin); HandleError(status);
  status = nc_get_vara_float(nctmeanin,tempidin,start,count,temp_in); HandleError(status);
  status = nc_close(nctmeanin); HandleError(status);

  if (ETscheme == 'P') { //etscheme = penman-monteith
    printf("readnetcdf tmax file %s\n",tmaxFileName);
    status = nc_open(tmaxFileName,NC_NOWRITE,&nctmaxin); HandleError(status);
    status = nc_inq_varid(nctmaxin,"tasmax",&tmaxidin); HandleError(status);
    status = nc_get_vara_float(nctmaxin,tmaxidin,start,count,tmax_in); HandleError(status);
    status = nc_close(nctmaxin); HandleError(status);

    printf("readnetcdf tmin file %s\n",tminFileName);
    status = nc_open(tminFileName,NC_NOWRITE,&nctminin); HandleError(status);
    status = nc_inq_varid(nctminin,"tasmin",&tminidin); HandleError(status);
    status = nc_get_vara_float(nctminin,tminidin,start,count,tmin_in); HandleError(status);
    status = nc_close(nctminin); HandleError(status);

    printf("readnetcdf sradfile %s\n",sradFileName);
    status = nc_open(sradFileName,NC_NOWRITE,&ncsradin); HandleError(status);
    status = nc_inq_varid(ncsradin,"rsds",&sradidin); HandleError(status);
    status = nc_get_vara_float(ncsradin,sradidin,start,count,srad_in); HandleError(status);
    status = nc_close(ncsradin); HandleError(status);

    printf("readnetcdf hursfile %s\n",hursFileName);
    status = nc_open(hursFileName,NC_NOWRITE,&nchursin); HandleError(status);
    status = nc_inq_varid(nchursin,"hurs",&hursidin); HandleError(status);
    status = nc_get_vara_float(nchursin,hursidin,start,count,hurs_in); HandleError(status);
    status = nc_close(nchursin); HandleError(status);

    printf("readnetcdf windfile %s\n",windFileName);
    status = nc_open(windFileName,NC_NOWRITE,&ncwindin); HandleError(status);
    status = nc_inq_varid(ncwindin,"sfcWind",&windidin); HandleError(status);
    status = nc_get_vara_float(ncwindin,windidin,start,count,wind_in); HandleError(status);
    status = nc_close(ncwindin); HandleError(status);
  }

  if (ETscheme == 'T') {
    for(k=0;k<ndays*(nRows*nCols);k++) {
      if(prec_in[k]<NODATA) {
	prec_in[k]*=S_PR_DAY; 	
        temp_in[k]-=KELVIN;
      }
    }
  }
  else {
    for(k=0;k<ndays*(nRows*nCols);k++) {
      if(prec_in[k]<NODATA) {
	prec_in[k]*=S_PR_DAY; 	
 	temp_in[k]-=KELVIN;
	tmax_in[k]-=KELVIN;
	tmin_in[k]-=KELVIN;
	vp_in[k]=hurs_in[k]*(6.11*exp((17.27*temp_in[k])/(237.3+temp_in[k]))); //hurs=percent,->vp_in=Pa (nb! later divided by 1000!)
      }
    }
  }
  
  /*flag=0; 
  for(k=0;k<ndays;k++) {
    cellid=0; 
    for(i=0;i<nRows;i++) {
      for(j=0;j<nCols;j++) {
	if(k<=150 && cellid==1573912) 
	  printf("readnetcdfclimate day %d row %d col %d cellnr %d \t p %.1f t %.2f tmin %.2f tmax %.2f srad %.1f wind %.1f hurs %.1f vp %.1f\n",
	  	 k,i,j,i*nCols+j,prec_in[flag],temp_in[flag],tmin_in[flag],tmax_in[flag],srad_in[flag],wind_in[flag],hurs_in[flag],vp_in[flag]);
	flag+=1;
	cellid+=1;
      }
    }
    }*/
}

void HandleError(int status) 
{
  if (status != NC_NOERR) {
    fprintf(stderr, "%s\n", nc_strerror(status));
    exit(-1);
  }
}

void TraverseCorrectionSubCatchment(SubCatchment * const thisSubCatchment, int numberCorrectionCatchments,
     int * correctionCatchments, double * correctionPrecipitation, double * correctionTemperature, ofstream &fout)
{
  int i;
  double precCorr, tempCorr;
  DistributedHbv * thisElement;

  for (i=0; i<thisSubCatchment->GetNumUpStream(); i++) {
    TraverseCorrectionSubCatchment(thisSubCatchment->GetUpStream(i), numberCorrectionCatchments, correctionCatchments, 
                                   correctionPrecipitation, correctionTemperature, fout);
  }
  //  fout << "numberCorrectionCatchments  " << numberCorrectionCatchments << endl;
  if (numberCorrectionCatchments > 0) {
    precCorr = correctionPrecipitation[0];
    tempCorr = correctionTemperature[0];
  }
  else {
    precCorr = 1.0;
    tempCorr = 0.0;
  }
  //  fout << "initial precCorr, tempCorr  " << precCorr << "  " << tempCorr << endl;
  for (i=0;i<numberCorrectionCatchments;i++) {
    if (thisSubCatchment->GetIdentifier() == correctionCatchments[i]) {
      precCorr = correctionPrecipitation[i];
      tempCorr = correctionTemperature[i];
    }
  }
  thisElement=thisSubCatchment->GetLandScapeElement();
  fout << thisSubCatchment->GetIdentifier() << "   " << "precCorr  " << precCorr << "    tempCorr  " << tempCorr << endl;
  while (thisElement) {
    TraverseCorrectionLandScape(thisElement, precCorr, tempCorr);
    thisElement = thisElement->GetNextElement();
  }
  //  thisElement=thisSubCatchment->GetLandScapeElement();
  //  while (thisElement) {
  //    fout << "  " << thisSubCatchment->GetIdentifier() << "  " << 
  //      thisElement->GetPrecipitationCorrection() << "  " << 
  //      thisElement->GetTemperatureCorrection() << endl;
  //    thisElement = thisElement->GetNextElement();
  //  }
}


void TraverseCorrectionLandScape(DistributedHbv * const thisElement, double precCorr, double tempCorr)
{
  thisElement->SetPrecipitationCorrection(precCorr);
  thisElement->SetTemperatureCorrection(tempCorr);
}


void TraverseSubCatchment(SubCatchment * const thisSubCatchment, int timeStep, ofstream &fout)
{
  int i;
  double accumulatedSum=0.0;
  double accumulatedSumLake=0.0;
  double accumulatedSumSnow=0.0;
  double accumulatedSumGlacier=0.0;
  double accumulatedSumHbv=0.0;
  double subCatchmentSum=0.0;
  double subCatchmentSumLake=0.0;
  double subCatchmentSumSnow=0.0;
  double subCatchmentSumGlacier=0.0;
  double subCatchmentSumHbv=0.0;
  double accumulatedDischarge=0.0;
  double accumulatedPrecipitation=0.0;
  double accumulatedTemperature=0.0;
  double accumulatedLakeStorage=0.0;
  double accumulatedSnowStore=0.0;
  double accumulatedEvapotranspiration=0.0;
  double accumulatedRunoff=0.0;
  double accumulatedHbvSoilMoisture=0.0;
  double accumulatedHbvSoilMoistureDeficit=0.0;
  double accumulatedHbvPercSoilUpper=0.0;
  double accumulatedHbvUpperZone=0.0;
  double accumulatedHbvLowerZone=0.0;
  double subCatchmentPrecipitation=0.0;
  double subCatchmentTemperature=0.0;
  double subCatchmentLakeStorage=0.0;
  double subCatchmentSnowStore=0.0;
  double accumulatedMeltWater=0.0;
  double accumulatedWaterOutput=0.0;
  double accumulatedSnowCoverFraction=0.0;
  double accumulatedGlacierMassBalance=0.0;
  double accumulatedGlacierIceMelt=0.0;
  double subCatchmentMeltWater=0.0;
  double subCatchmentWaterOutput=0.0;
  double subCatchmentSnowCoverFraction=0.0;
  double subCatchmentGlacierMassBalance=0.0;
  double subCatchmentGlacierIceMelt=0.0;
  double subCatchmentEvapotranspiration=0.0;
  double subCatchmentRunoff=0.0;
  double subCatchmentHbvSoilMoisture=0.0;
  double subCatchmentHbvSoilMoistureDeficit=0.0;
  double subCatchmentHbvPercSoilUpper=0.0;
  double subCatchmentHbvUpperZone=0.0;
  double subCatchmentHbvLowerZone=0.0;
  DistributedHbv * thisElement;

  for (i=0; i<thisSubCatchment->GetNumUpStream(); i++) {
    TraverseSubCatchment(thisSubCatchment->GetUpStream(i), timeStep, fout);
    // Fluxes accumulated
    accumulatedDischarge=accumulatedDischarge+thisSubCatchment->GetUpStream(i)->GetAccumulatedDischarge(timeStep);
    accumulatedPrecipitation=accumulatedPrecipitation+thisSubCatchment->GetUpStream(i)->GetAccumulatedPrecipitation(timeStep);
    accumulatedTemperature=accumulatedTemperature+thisSubCatchment->GetUpStream(i)->GetAccumulatedTemperature(timeStep);
    accumulatedEvapotranspiration=accumulatedEvapotranspiration+thisSubCatchment->GetUpStream(i)->GetAccumulatedEvapotranspiration(timeStep);
    accumulatedRunoff=accumulatedRunoff+thisSubCatchment->GetUpStream(i)->GetAccumulatedRunoff(timeStep);
    accumulatedSum=accumulatedSum+thisSubCatchment->GetUpStream(i)->GetAccumulatedSum(timeStep);
    // Lake state variables accumulated
    accumulatedLakeStorage=accumulatedLakeStorage+thisSubCatchment->GetUpStream(i)->GetAccumulatedLakeStorage(timeStep);
    accumulatedSumLake=accumulatedSumLake+thisSubCatchment->GetUpStream(i)->GetAccumulatedSumLake(timeStep);
    // Snow state variables accumulated
    // Snow state variables accumulated
    accumulatedSnowStore=accumulatedSnowStore+thisSubCatchment->GetUpStream(i)->GetAccumulatedSnowStore(timeStep);
    accumulatedMeltWater=accumulatedMeltWater+thisSubCatchment->GetUpStream(i)->GetAccumulatedMeltWater(timeStep);
    accumulatedWaterOutput=accumulatedWaterOutput+thisSubCatchment->GetUpStream(i)->GetAccumulatedWaterOutput(timeStep);
    accumulatedSnowCoverFraction=accumulatedSnowCoverFraction+thisSubCatchment->GetUpStream(i)->GetAccumulatedSnowCoverFraction(timeStep);
    accumulatedSumSnow=accumulatedSumSnow+thisSubCatchment->GetUpStream(i)->GetAccumulatedSumSnow(timeStep);
    // Glacier state variables accumulated
    accumulatedGlacierMassBalance=accumulatedGlacierMassBalance+thisSubCatchment->GetUpStream(i)->GetAccumulatedGlacierMassBalance(timeStep);
    accumulatedGlacierIceMelt=accumulatedGlacierIceMelt+thisSubCatchment->GetUpStream(i)->GetAccumulatedGlacierIceMelt(timeStep);
    accumulatedSumGlacier=accumulatedSumGlacier+thisSubCatchment->GetUpStream(i)->GetAccumulatedSumGlacier(timeStep);
    // Hbv state variables accumulated
    accumulatedHbvSoilMoisture=accumulatedHbvSoilMoisture+thisSubCatchment->GetUpStream(i)->GetAccumulatedHbvSoilMoisture(timeStep);
    accumulatedHbvSoilMoistureDeficit=accumulatedHbvSoilMoistureDeficit+thisSubCatchment->GetUpStream(i)->GetAccumulatedHbvSoilMoistureDeficit(timeStep);
    accumulatedHbvPercSoilUpper=accumulatedHbvPercSoilUpper+thisSubCatchment->GetUpStream(i)->GetAccumulatedHbvPercSoilUpper(timeStep);
    accumulatedHbvUpperZone=accumulatedHbvUpperZone+thisSubCatchment->GetUpStream(i)->GetAccumulatedHbvUpperZone(timeStep);
    accumulatedHbvLowerZone=accumulatedHbvLowerZone+thisSubCatchment->GetUpStream(i)->GetAccumulatedHbvLowerZone(timeStep);
    accumulatedSumHbv=accumulatedSumHbv+thisSubCatchment->GetUpStream(i)->GetAccumulatedSumHbv(timeStep);
  }
  thisElement=thisSubCatchment->GetLandScapeElement();
  while (thisElement) {
    TraverseLandScape(thisElement, timeStep, fout);
    // Fluxes accumulated
    accumulatedDischarge=accumulatedDischarge+thisElement->GetAccumulatedDischarge();
    accumulatedPrecipitation=accumulatedPrecipitation+thisElement->GetAccumulatedPrecipitation();
    accumulatedTemperature=accumulatedTemperature+thisElement->GetAccumulatedTemperature();
    accumulatedEvapotranspiration=accumulatedEvapotranspiration+thisElement->GetAccumulatedEvapotranspiration();
    accumulatedRunoff=accumulatedRunoff+thisElement->GetAccumulatedRunoff();
    accumulatedSum=accumulatedSum+thisElement->GetAccumulatedSum();
    subCatchmentPrecipitation=subCatchmentPrecipitation+thisElement->GetAccumulatedPrecipitation();
    subCatchmentTemperature=subCatchmentTemperature+thisElement->GetAccumulatedTemperature();
    subCatchmentEvapotranspiration=subCatchmentEvapotranspiration+thisElement->GetAccumulatedEvapotranspiration();
    subCatchmentRunoff=subCatchmentRunoff+thisElement->GetAccumulatedRunoff();
    subCatchmentSum=subCatchmentSum+thisElement->GetAccumulatedSum();
    // Lake state variables accumulated
    if (thisElement->GetAccumulatedSumLake() > 0) {
        accumulatedLakeStorage=accumulatedLakeStorage+thisElement->GetAccumulatedLakeStorage();
        accumulatedSumLake=accumulatedSumLake+thisElement->GetAccumulatedSumLake();
        subCatchmentLakeStorage=subCatchmentLakeStorage+thisElement->GetAccumulatedLakeStorage();
        subCatchmentSumLake=subCatchmentSumLake+thisElement->GetAccumulatedSumLake();
    }
    // Snow state variables accumulated
    if (thisElement->GetAccumulatedSumSnow() > 0) {
        accumulatedSnowStore=accumulatedSnowStore+thisElement->GetAccumulatedSnowStore();
        accumulatedMeltWater=accumulatedMeltWater+thisElement->GetAccumulatedMeltWater();
        accumulatedWaterOutput=accumulatedWaterOutput+thisElement->GetAccumulatedWaterOutput();
        accumulatedSnowCoverFraction=accumulatedSnowCoverFraction+thisElement->GetAccumulatedSnowCoverFraction();
        accumulatedSumSnow=accumulatedSumSnow+thisElement->GetAccumulatedSumSnow();
        subCatchmentSnowStore=subCatchmentSnowStore+thisElement->GetAccumulatedSnowStore();
        subCatchmentMeltWater=subCatchmentMeltWater+thisElement->GetAccumulatedMeltWater();
        subCatchmentWaterOutput=subCatchmentWaterOutput+thisElement->GetAccumulatedWaterOutput();
        subCatchmentSnowCoverFraction=subCatchmentSnowCoverFraction+thisElement->GetAccumulatedSnowCoverFraction();
        subCatchmentSumSnow=subCatchmentSumSnow+thisElement->GetAccumulatedSumSnow();
    }
    // Glacier state variables accumulated
    if (thisElement->GetAccumulatedSumGlacier() > 0) {
        accumulatedGlacierMassBalance=accumulatedGlacierMassBalance+thisElement->GetAccumulatedGlacierMassBalance();
        accumulatedGlacierIceMelt=accumulatedGlacierIceMelt+thisElement->GetAccumulatedGlacierIceMelt();
        accumulatedSumGlacier=accumulatedSumGlacier+thisElement->GetAccumulatedSumGlacier();
        subCatchmentGlacierMassBalance=subCatchmentGlacierMassBalance+thisElement->GetAccumulatedGlacierMassBalance();
        subCatchmentGlacierIceMelt=subCatchmentGlacierIceMelt+thisElement->GetAccumulatedGlacierIceMelt();
        subCatchmentSumGlacier=subCatchmentSumGlacier+thisElement->GetAccumulatedSumGlacier();
    }
    // Hbv state variables accumulated
    if (thisElement->GetAccumulatedSumHbv() > 0) {
        accumulatedHbvSoilMoisture=accumulatedHbvSoilMoisture+thisElement->GetAccumulatedHbvSoilMoisture();
        accumulatedHbvSoilMoistureDeficit=accumulatedHbvSoilMoistureDeficit+thisElement->GetAccumulatedHbvSoilMoistureDeficit();
        accumulatedHbvPercSoilUpper=accumulatedHbvPercSoilUpper+thisElement->GetAccumulatedHbvPercSoilUpper();
        accumulatedHbvUpperZone=accumulatedHbvUpperZone+thisElement->GetAccumulatedHbvUpperZone();
        accumulatedHbvLowerZone=accumulatedHbvLowerZone+thisElement->GetAccumulatedHbvLowerZone();
        accumulatedSumHbv=accumulatedSumHbv+thisElement->GetAccumulatedSumHbv();
        subCatchmentHbvSoilMoisture=subCatchmentHbvSoilMoisture+thisElement->GetAccumulatedHbvSoilMoisture();
        subCatchmentHbvSoilMoistureDeficit=subCatchmentHbvSoilMoistureDeficit+thisElement->GetAccumulatedHbvSoilMoistureDeficit();
        subCatchmentHbvPercSoilUpper=subCatchmentHbvPercSoilUpper+thisElement->GetAccumulatedHbvPercSoilUpper();
        subCatchmentHbvUpperZone=subCatchmentHbvUpperZone+thisElement->GetAccumulatedHbvUpperZone();
        subCatchmentHbvLowerZone=subCatchmentHbvLowerZone+thisElement->GetAccumulatedHbvLowerZone();
        subCatchmentSumHbv=subCatchmentSumHbv+thisElement->GetAccumulatedSumHbv();
    }
    thisElement = thisElement->GetNextElement();
  }
  // Fluxes accumulated
  thisSubCatchment->SetAccumulatedDischarge(timeStep, accumulatedDischarge);
  thisSubCatchment->SetAccumulatedPrecipitation(timeStep, accumulatedPrecipitation);
  thisSubCatchment->SetAccumulatedTemperature(timeStep, accumulatedTemperature);
  thisSubCatchment->SetAccumulatedEvapotranspiration(timeStep, accumulatedEvapotranspiration);
  thisSubCatchment->SetAccumulatedRunoff(timeStep, accumulatedRunoff);
  thisSubCatchment->SetAccumulatedSum(timeStep, accumulatedSum);
  thisSubCatchment->SetSubCatchmentPrecipitation(timeStep, subCatchmentPrecipitation/subCatchmentSum);
  thisSubCatchment->SetSubCatchmentTemperature(timeStep, subCatchmentTemperature/subCatchmentSum);
  thisSubCatchment->SetSubCatchmentEvapotranspiration(timeStep, subCatchmentEvapotranspiration/subCatchmentSum);
  thisSubCatchment->SetSubCatchmentRunoff(timeStep, subCatchmentRunoff/subCatchmentSum);
  // Lake state variables accumulated
  if (accumulatedSumLake > 0) {
    thisSubCatchment->SetAccumulatedLakeStorage(timeStep, accumulatedLakeStorage);
    thisSubCatchment->SetAccumulatedSumLake(timeStep, accumulatedSumLake);
  }
  else {
    thisSubCatchment->SetAccumulatedLakeStorage(timeStep, accumulatedLakeStorage);
    thisSubCatchment->SetAccumulatedSumLake(timeStep, accumulatedSumLake);
  }
  if (subCatchmentSumLake > 0) {
    thisSubCatchment->SetSubCatchmentLakeStorage(timeStep, subCatchmentLakeStorage/subCatchmentSumLake);
    //    thisSubCatchment->SetSubCatchmentLakeStorage(timeStep, subCatchmentLakeStorage/subCatchmentSum);
  }
  else {
    thisSubCatchment->SetSubCatchmentLakeStorage(timeStep, missingData);
  }
  // Snow state variables accumulated
  if (accumulatedSumSnow > 0) {
    thisSubCatchment->SetAccumulatedSnowStore(timeStep, accumulatedSnowStore);
    thisSubCatchment->SetAccumulatedMeltWater(timeStep, accumulatedMeltWater);
    thisSubCatchment->SetAccumulatedWaterOutput(timeStep, accumulatedWaterOutput);
    thisSubCatchment->SetAccumulatedSnowCoverFraction(timeStep, accumulatedSnowCoverFraction);
    thisSubCatchment->SetAccumulatedSumSnow(timeStep, accumulatedSumSnow);
  }
  else {
    thisSubCatchment->SetAccumulatedSnowStore(timeStep, accumulatedSnowStore);
    thisSubCatchment->SetAccumulatedMeltWater(timeStep, accumulatedMeltWater);
    thisSubCatchment->SetAccumulatedWaterOutput(timeStep, accumulatedWaterOutput);
    thisSubCatchment->SetAccumulatedSnowCoverFraction(timeStep, accumulatedSnowCoverFraction);
    thisSubCatchment->SetAccumulatedSumSnow(timeStep, accumulatedSumSnow);
  }
  if (subCatchmentSumSnow > 0) {
    thisSubCatchment->SetSubCatchmentSnowStore(timeStep, subCatchmentSnowStore/subCatchmentSumSnow);
    thisSubCatchment->SetSubCatchmentMeltWater(timeStep, subCatchmentMeltWater/subCatchmentSumSnow);
    thisSubCatchment->SetSubCatchmentWaterOutput(timeStep, subCatchmentWaterOutput/subCatchmentSumSnow);
    thisSubCatchment->SetSubCatchmentSnowCoverFraction(timeStep, subCatchmentSnowCoverFraction/subCatchmentSumSnow);
    //    thisSubCatchment->SetSubCatchmentSnowStore(timeStep, subCatchmentSnowStore/subCatchmentSum);
    //    thisSubCatchment->SetSubCatchmentMeltWater(timeStep, subCatchmentMeltWater/subCatchmentSum);
    //    thisSubCatchment->SetSubCatchmentWaterOutput(timeStep, subCatchmentWaterOutput/subCatchmentSum);
    //    thisSubCatchment->SetSubCatchmentSnowCoverFraction(timeStep, subCatchmentSnowCoverFraction/subCatchmentSum);
  }
  else {
    thisSubCatchment->SetSubCatchmentSnowStore(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentMeltWater(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentWaterOutput(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentSnowCoverFraction(timeStep, missingData);
  }
  // Glacier state variables accumulated
  if (accumulatedSumGlacier > 0) {
    thisSubCatchment->SetAccumulatedGlacierMassBalance(timeStep, accumulatedGlacierMassBalance);
    thisSubCatchment->SetAccumulatedGlacierIceMelt(timeStep, accumulatedGlacierIceMelt);
    thisSubCatchment->SetAccumulatedSumGlacier(timeStep, accumulatedSumGlacier);
  }
  else {
    thisSubCatchment->SetAccumulatedGlacierMassBalance(timeStep, accumulatedGlacierMassBalance);
    thisSubCatchment->SetAccumulatedGlacierIceMelt(timeStep, accumulatedGlacierIceMelt);
    thisSubCatchment->SetAccumulatedSumGlacier(timeStep, accumulatedSumGlacier);
  }
  if (subCatchmentSumGlacier > 0) {
    thisSubCatchment->SetSubCatchmentGlacierMassBalance(timeStep, subCatchmentGlacierMassBalance/subCatchmentSumGlacier);
    thisSubCatchment->SetSubCatchmentGlacierIceMelt(timeStep, subCatchmentGlacierIceMelt/subCatchmentSum);
    //    thisSubCatchment->SetSubCatchmentGlacierMassBalance(timeStep, subCatchmentGlacierMassBalance/subCatchmentSumGlacier);
    //    thisSubCatchment->SetSubCatchmentGlacierIceMelt(timeStep, subCatchmentGlacierIceMelt/subCatchmentSum);
  }
  else {
    thisSubCatchment->SetSubCatchmentGlacierMassBalance(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentGlacierIceMelt(timeStep, missingData);
  }
  // Hbv state variables accumulated
  if (accumulatedSumHbv > 0) {
    thisSubCatchment->SetAccumulatedHbvSoilMoisture(timeStep, accumulatedHbvSoilMoisture);
    thisSubCatchment->SetAccumulatedHbvSoilMoistureDeficit(timeStep, accumulatedHbvSoilMoistureDeficit);
    thisSubCatchment->SetAccumulatedHbvPercSoilUpper(timeStep, accumulatedHbvPercSoilUpper);
    thisSubCatchment->SetAccumulatedHbvUpperZone(timeStep, accumulatedHbvUpperZone);
    thisSubCatchment->SetAccumulatedHbvLowerZone(timeStep, accumulatedHbvLowerZone);
    thisSubCatchment->SetAccumulatedSumHbv(timeStep, accumulatedSumHbv);
  }
  else {
    thisSubCatchment->SetAccumulatedHbvSoilMoisture(timeStep, accumulatedHbvSoilMoisture);
    thisSubCatchment->SetAccumulatedHbvSoilMoistureDeficit(timeStep, accumulatedHbvSoilMoistureDeficit);
    thisSubCatchment->SetAccumulatedHbvPercSoilUpper(timeStep, accumulatedHbvPercSoilUpper);
    thisSubCatchment->SetAccumulatedHbvUpperZone(timeStep, accumulatedHbvUpperZone);
    thisSubCatchment->SetAccumulatedHbvLowerZone(timeStep, accumulatedHbvLowerZone);
    thisSubCatchment->SetAccumulatedSumHbv(timeStep, accumulatedSumHbv);
  }
  if (subCatchmentSumHbv > 0) {
    thisSubCatchment->SetSubCatchmentHbvSoilMoisture(timeStep, subCatchmentHbvSoilMoisture/subCatchmentSumHbv);
    thisSubCatchment->SetSubCatchmentHbvSoilMoistureDeficit(timeStep, subCatchmentHbvSoilMoistureDeficit/subCatchmentSumHbv);
    thisSubCatchment->SetSubCatchmentHbvPercSoilUpper(timeStep, subCatchmentHbvPercSoilUpper/subCatchmentSumHbv);
    thisSubCatchment->SetSubCatchmentHbvUpperZone(timeStep, subCatchmentHbvUpperZone/subCatchmentSumHbv);
    thisSubCatchment->SetSubCatchmentHbvLowerZone(timeStep, subCatchmentHbvLowerZone/subCatchmentSumHbv);
    //    thisSubCatchment->SetSubCatchmentHbvSoilMoisture(timeStep, subCatchmentHbvSoilMoisture/subCatchmentSum);
    //    thisSubCatchment->SetSubCatchmentHbvSoilMoistureDeficit(timeStep, subCatchmentHbvSoilMoistureDeficit/subCatchmentSum);
    //    thisSubCatchment->SetSubCatchmentHbvPercSoilUpper(timeStep, subCatchmentHbvPercSoilUpper/subCatchmentSum);
    //    thisSubCatchment->SetSubCatchmentHbvUpperZone(timeStep, subCatchmentHbvUpperZone/subCatchmentSum);
    //    thisSubCatchment->SetSubCatchmentHbvLowerZone(timeStep, subCatchmentHbvLowerZone/subCatchmentSum);
  }
  else {
    thisSubCatchment->SetSubCatchmentHbvSoilMoisture(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentHbvSoilMoistureDeficit(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentHbvPercSoilUpper(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentHbvUpperZone(timeStep, missingData);
    thisSubCatchment->SetSubCatchmentHbvLowerZone(timeStep, missingData);
  }
  //  cout << thisSubCatchment->GetIdentifier() << "  " << timeStep << "  " << "  Discharge: "  
  //       << accumulatedDischarge << endl;
  //       << accumulatedDischarge << "  " << thisSubCatchment->GetAccumulatedDischarge(timeStep) << endl;
}


void TraverseLandScape(DistributedHbv * const thisElement, int timeStep, ofstream &fout)
{
  int j;
  double accumulatedSum=0.0;
  double accumulatedSumLake=0.0;
  double accumulatedSumSnow=0.0;
  double accumulatedSumGlacier=0.0;
  double accumulatedSumHbv=0.0;
  double accumulatedDischarge=0.0;
  double accumulatedLowerDischarge=0.0;
  double accumulatedUpperDischarge=0.0;
  double accumulatedPrecipitation=0.0;
  double accumulatedTemperature=0.0;
  double accumulatedLakeStorage=0.0;
  double accumulatedSnowStore=0.0;
  double accumulatedMeltWater=0.0;
  double accumulatedWaterOutput=0.0;
  double accumulatedSnowCoverFraction=0.0;
  double accumulatedGlacierMassBalance=0.0;
  double accumulatedGlacierIceMelt=0.0;
  double accumulatedEvapotranspiration=0.0;
  double accumulatedRunoff=0.0;
  double accumulatedHbvSoilMoisture=0.0;
  double accumulatedHbvSoilMoistureDeficit=0.0;
  double accumulatedHbvPercSoilUpper=0.0;
  double accumulatedHbvUpperZone=0.0;
  double accumulatedHbvLowerZone=0.0;
  // Fluxes accumulated
  thisElement->SetAccumulatedDischarge(accumulatedDischarge+thisElement->GetDischarge());
  thisElement->SetAccumulatedPrecipitation(accumulatedPrecipitation+thisElement->GetPrecipitation()*thisElement->GetArea());
  thisElement->SetAccumulatedTemperature(accumulatedTemperature+thisElement->GetTemperature()*thisElement->GetArea());
  thisElement->SetAccumulatedEvapotranspiration(accumulatedEvapotranspiration+thisElement->GetEvapotranspiration()*thisElement->GetArea());
  thisElement->SetAccumulatedRunoff(accumulatedRunoff+thisElement->GetRunoff()*thisElement->GetArea());
  thisElement->SetAccumulatedSum(accumulatedSum+thisElement->GetArea());
  // Lake state variables accumulated
  if (thisElement->GetLakeStorage() != missingData) {
    thisElement->SetAccumulatedLakeStorage(accumulatedLakeStorage+thisElement->GetLakeStorage()*thisElement->GetLakeArea());
    thisElement->SetAccumulatedSumLake(accumulatedSumLake+thisElement->GetLakeArea());
  }
  else {
    thisElement->SetAccumulatedLakeStorage(accumulatedLakeStorage);
    thisElement->SetAccumulatedSumLake(accumulatedSumLake);
  }
  // Snow state variables accumulated
  if (thisElement->GetSnowStore() != missingData) {
    thisElement->SetAccumulatedSnowStore(accumulatedSnowStore+thisElement->GetSnowStore()*thisElement->GetLandArea());
    thisElement->SetAccumulatedMeltWater(accumulatedMeltWater+thisElement->GetMeltWater()*thisElement->GetLandArea());
    thisElement->SetAccumulatedWaterOutput(accumulatedWaterOutput+thisElement->GetWaterOutput()*thisElement->GetLandArea());
    thisElement->SetAccumulatedSnowCoverFraction(accumulatedSnowCoverFraction+thisElement->GetSnowCoverFraction()*thisElement->GetLandArea());
    thisElement->SetAccumulatedSumSnow(accumulatedSumSnow+thisElement->GetLandArea());
  }
  else {
    thisElement->SetAccumulatedSnowStore(accumulatedSnowStore);
    thisElement->SetAccumulatedMeltWater(accumulatedMeltWater);
    thisElement->SetAccumulatedWaterOutput(accumulatedWaterOutput);
    thisElement->SetAccumulatedSnowCoverFraction(accumulatedSnowCoverFraction);
    thisElement->SetAccumulatedSumSnow(accumulatedSumSnow);
  }
  // Glacier state variables accumulated
  if (thisElement->GetGlacierMassBalance() != missingData) {
    thisElement->SetAccumulatedGlacierMassBalance(accumulatedGlacierMassBalance+thisElement->GetGlacierMassBalance()*thisElement->GetGlacierArea());
    thisElement->SetAccumulatedGlacierIceMelt(accumulatedGlacierIceMelt+thisElement->GetGlacierIceMelt()*thisElement->GetGlacierArea());
    thisElement->SetAccumulatedSumGlacier(accumulatedSumGlacier+thisElement->GetGlacierArea());
  }
  else {
    thisElement->SetAccumulatedGlacierMassBalance(accumulatedGlacierMassBalance);
    thisElement->SetAccumulatedGlacierIceMelt(accumulatedGlacierIceMelt);
    thisElement->SetAccumulatedSumGlacier(accumulatedSumGlacier);
  }
   // Hbv state variables accumulated
  if (thisElement->GetHbvSoilMoisture() != missingData) {
    thisElement->SetAccumulatedHbvSoilMoisture(accumulatedHbvSoilMoisture+thisElement->GetHbvSoilMoisture()*thisElement->GetHbvArea());
    thisElement->SetAccumulatedSumHbv(accumulatedSumHbv+thisElement->GetHbvArea());
  }
  else {
    thisElement->SetAccumulatedHbvSoilMoisture(accumulatedHbvSoilMoisture);
    thisElement->SetAccumulatedSumHbv(accumulatedSumHbv);
  }
  if (thisElement->GetHbvSoilMoistureDeficit() != missingData) {
    thisElement->SetAccumulatedHbvSoilMoistureDeficit(accumulatedHbvSoilMoistureDeficit+thisElement->GetHbvSoilMoistureDeficit()*thisElement->GetHbvArea());
  }
  else {
    thisElement->SetAccumulatedHbvSoilMoistureDeficit(accumulatedHbvSoilMoistureDeficit);
  }
  if (thisElement->GetHbvPercSoilUpper() != missingData) {
    thisElement->SetAccumulatedHbvPercSoilUpper(accumulatedHbvPercSoilUpper+thisElement->GetHbvPercSoilUpper()*thisElement->GetHbvArea());
  }
  else {
    thisElement->SetAccumulatedHbvPercSoilUpper(accumulatedHbvPercSoilUpper);
  }
  if (thisElement->GetHbvUpperZone() != missingData) {
    thisElement->SetAccumulatedHbvUpperZone(accumulatedHbvUpperZone+thisElement->GetHbvUpperZone()*thisElement->GetHbvArea());
  }
  else {
    thisElement->SetAccumulatedHbvUpperZone(accumulatedHbvUpperZone);
  }
  if (thisElement->GetHbvLowerZone() != missingData) {
    thisElement->SetAccumulatedHbvLowerZone(accumulatedHbvLowerZone+thisElement->GetHbvLowerZone()*thisElement->GetHbvArea());
  }
  else {
    thisElement->SetAccumulatedHbvLowerZone(accumulatedHbvLowerZone);
  }
  // Time series for landscape elements
  for (j=0; j<thisElement->GetSelectedTimeSeriesElements()->GetNumberElements(); j++) {
    if (thisElement->GetLandIndex() == thisElement->GetSelectedTimeSeriesElements()->GetTimeSeriesElement(j)) {
      thisElement->SetDistributedHbvPrecipitation(timeStep,thisElement->GetPrecipitation());
      thisElement->SetDistributedHbvTemperature(timeStep,thisElement->GetTemperature());
      thisElement->SetDistributedHbvEvapotranspiration(timeStep,thisElement->GetEvapotranspiration());
      thisElement->SetDistributedHbvRunoff(timeStep,thisElement->GetRunoff());
      thisElement->SetDistributedHbvLakeStorage(timeStep,thisElement->GetLakeStorage());
      thisElement->SetDistributedHbvSnowStore(timeStep,thisElement->GetSnowStore());
      thisElement->SetDistributedHbvSnowCoverFraction(timeStep,thisElement->GetSnowCoverFraction());
      thisElement->SetDistributedHbvMeltWater(timeStep,thisElement->GetMeltWater());
      thisElement->SetDistributedHbvWaterOutput(timeStep,thisElement->GetWaterOutput());
      thisElement->SetDistributedHbvGlacierMassBalance(timeStep,thisElement->GetGlacierMassBalance());
      thisElement->SetDistributedHbvGlacierIceMelt(timeStep,thisElement->GetGlacierIceMelt());
      thisElement->SetDistributedHbvSoilMoisture(timeStep,thisElement->GetHbvSoilMoisture());
      thisElement->SetDistributedHbvSoilMoistureDeficit(timeStep,thisElement->GetHbvSoilMoistureDeficit());
      thisElement->SetDistributedHbvPercSoilUpper(timeStep,thisElement->GetHbvPercSoilUpper());
      thisElement->SetDistributedHbvUpperZone(timeStep,thisElement->GetHbvUpperZone());
      thisElement->SetDistributedHbvLowerZone(timeStep,thisElement->GetHbvLowerZone());
    }
  }
}


void TraverseMissingDataSubCatchment(SubCatchment * const thisSubCatchment, int timeStep, ofstream &fout)
{
  int i;
  DistributedHbv * thisElement;

  for (i=0; i<thisSubCatchment->GetNumUpStream(); i++) {
    TraverseMissingDataSubCatchment(thisSubCatchment->GetUpStream(i), timeStep, fout);
  }
  thisElement=thisSubCatchment->GetLandScapeElement();
  while (thisElement) {
    TraverseMissingDataLandScape(thisElement, timeStep, fout);
    thisElement = thisElement->GetNextElement();
  }
  thisSubCatchment->SetAccumulatedDischarge(timeStep, missingData);
  thisSubCatchment->SetAccumulatedPrecipitation(timeStep, missingData);
  thisSubCatchment->SetAccumulatedTemperature(timeStep, missingData);
  thisSubCatchment->SetAccumulatedLakeStorage(timeStep, missingData);
  thisSubCatchment->SetAccumulatedSnowStore(timeStep, missingData);
  thisSubCatchment->SetAccumulatedMeltWater(timeStep, missingData);
  thisSubCatchment->SetAccumulatedWaterOutput(timeStep, missingData);
  thisSubCatchment->SetAccumulatedSnowCoverFraction(timeStep, missingData);
  thisSubCatchment->SetAccumulatedGlacierMassBalance(timeStep, missingData);
  thisSubCatchment->SetAccumulatedGlacierIceMelt(timeStep, missingData);
  thisSubCatchment->SetAccumulatedEvapotranspiration(timeStep, missingData);
  thisSubCatchment->SetAccumulatedRunoff(timeStep, missingData);
  thisSubCatchment->SetAccumulatedHbvSoilMoisture(timeStep, missingData);
  thisSubCatchment->SetAccumulatedHbvSoilMoistureDeficit(timeStep, missingData);
  thisSubCatchment->SetAccumulatedHbvPercSoilUpper(timeStep, missingData);
  thisSubCatchment->SetAccumulatedHbvUpperZone(timeStep, missingData);
  thisSubCatchment->SetAccumulatedHbvLowerZone(timeStep, missingData);
  thisSubCatchment->SetAccumulatedSum(timeStep, missingData);
  thisSubCatchment->SetAccumulatedSumSnow(timeStep, missingData);
  thisSubCatchment->SetAccumulatedSumGlacier(timeStep, missingData);
  thisSubCatchment->SetAccumulatedSumHbv(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentPrecipitation(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentTemperature(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentSnowStore(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentMeltWater(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentWaterOutput(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentSnowCoverFraction(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentGlacierMassBalance(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentGlacierIceMelt(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentEvapotranspiration(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentRunoff(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentHbvSoilMoisture(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentHbvSoilMoistureDeficit(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentHbvPercSoilUpper(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentHbvUpperZone(timeStep, missingData);
  thisSubCatchment->SetSubCatchmentHbvLowerZone(timeStep, missingData);
}


void TraverseMissingDataLandScape(DistributedHbv * const thisElement, int timeStep, ofstream &fout)
{
  int j;
  thisElement->SetAccumulatedDischarge(missingData);
  thisElement->SetAccumulatedPrecipitation(missingData);
  thisElement->SetAccumulatedTemperature(missingData);
  thisElement->SetAccumulatedLakeStorage(missingData);
  thisElement->SetAccumulatedSnowStore(missingData);
  thisElement->SetAccumulatedMeltWater(missingData);
  thisElement->SetAccumulatedWaterOutput(missingData);
  thisElement->SetAccumulatedSnowCoverFraction(missingData);
  thisElement->SetAccumulatedGlacierMassBalance(missingData);
  thisElement->SetAccumulatedGlacierIceMelt(missingData);
  thisElement->SetAccumulatedEvapotranspiration(missingData);
  thisElement->SetAccumulatedRunoff(missingData);
  thisElement->SetAccumulatedHbvSoilMoisture(missingData);
  thisElement->SetAccumulatedHbvSoilMoistureDeficit(missingData);
  thisElement->SetAccumulatedHbvPercSoilUpper(missingData);
  thisElement->SetAccumulatedHbvUpperZone(missingData);
  thisElement->SetAccumulatedHbvLowerZone(missingData);
  thisElement->SetAccumulatedSum(missingData);
  thisElement->SetAccumulatedSumSnow(missingData);
  thisElement->SetAccumulatedSumGlacier(missingData);
  thisElement->SetAccumulatedSumHbv(missingData);
  for (j=0; j<thisElement->GetSelectedTimeSeriesElements()->GetNumberElements(); j++) {
    if (thisElement->GetLandIndex() == thisElement->GetSelectedTimeSeriesElements()->GetTimeSeriesElement(j)) {
      thisElement->SetDistributedHbvPrecipitation(timeStep,missingData);
      thisElement->SetDistributedHbvTemperature(timeStep,missingData);
      thisElement->SetDistributedHbvLakeStorage(timeStep,missingData);
      thisElement->SetDistributedHbvSnowStore(timeStep,missingData);
      thisElement->SetDistributedHbvSnowCoverFraction(timeStep,missingData);
      thisElement->SetDistributedHbvMeltWater(timeStep,missingData);
      thisElement->SetDistributedHbvWaterOutput(timeStep,missingData);
      thisElement->SetDistributedHbvGlacierMassBalance(timeStep,missingData);
      thisElement->SetDistributedHbvGlacierIceMelt(timeStep,missingData);
      thisElement->SetDistributedHbvEvapotranspiration(timeStep,missingData);
      thisElement->SetDistributedHbvRunoff(timeStep,missingData);
      thisElement->SetDistributedHbvSoilMoisture(timeStep,missingData);
      thisElement->SetDistributedHbvSoilMoistureDeficit(timeStep,missingData);
      thisElement->SetDistributedHbvPercSoilUpper(timeStep,missingData);
      thisElement->SetDistributedHbvUpperZone(timeStep,missingData);
      thisElement->SetDistributedHbvLowerZone(timeStep,missingData);
    }
  }
}


// For test purpose only
void WriteSubCatchmentIdentifier(SubCatchment * const CatchmentElement, int numWatc, ofstream &fout)
{
  int i;
  DistributedHbv * thisElement;

  // Sub-Catchment outlets
  ofstream subCatchmentOut("test_waterland.txt");  // Open for writing"
  for (i=0; i<numWatc; i++) {
    if (CatchmentElement[i].GetLandScapeElement()) {
      subCatchmentOut << "#  " << CatchmentElement[i].GetIdentifier() << "  #  " 
                   << CatchmentElement[i].GetNumLandScape() << endl;
      thisElement = CatchmentElement[i].GetLandScapeElement();
      while (thisElement) {
        subCatchmentOut << thisElement->GetLandIndex() << "  " << thisElement->GetGeoIndex() << endl;
        thisElement = thisElement->GetNextElement();
      }
    }
  }
  subCatchmentOut << endl;
  subCatchmentOut.close();
}


void WriteSubCatchmentDischarge(SubCatchment * CatchmentElement, int numWatc, DateTime startSimulationTime, 
                                DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                                int secondsPerTimeStep, bool modelCalibration, ofstream &fout)
{
  FILE *fpOut;
  char fileName[100];
  double sumValue;
  double nashSutcliffe;        /*  Nash-Sutcliffe efficiency criterion  */
  double rootMSError;          /*  Root mean square error criterion  */
  double biasVol;              /*  Bias (volume error) criterion  */
  double pearsCorr;            /*  Pearson's product-moment correlation coefficient  */
  int i,j,timeStep;
  DateTime datetime;
  double * observed = new double [numberTimeSteps];
  double * simulated = new double [numberTimeSteps];

  for (i=0; i<numWatc; i++) {
    sprintf(fileName,"hbv_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpOut = fopen(fileName, "w")) == NULL ) {
        printf("\n File %s not found!\n\n",fileName);
        exit(1);
    }
    sumValue = 0.0;
    j = 0;
    timeStep=initialTimeSteps;
    for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=secondsPerTimeStep) {
      if (CatchmentElement[i].GetAccumulatedDischarge(timeStep) > missingData) {
        fprintf(fpOut,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                CatchmentElement[i].GetAccumulatedDischarge(timeStep)*CatchmentElement[i].GetCorrection());
        if (modelCalibration && CatchmentElement[i].GetObsData(j) > missingData) 
          sumValue = sumValue+CatchmentElement[i].GetAccumulatedDischarge(timeStep);
      }
      else {
        fprintf(fpOut,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
                datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
      }
      observed[j] = CatchmentElement[i].GetObsData(j);
      simulated[j] = CatchmentElement[i].GetAccumulatedDischarge(timeStep)*CatchmentElement[i].GetCorrection();
      j++;
      timeStep++;
    }
    if (modelCalibration) fprintf(fpOut,"%29.6f\n", sumValue*CatchmentElement[i].GetCorrection());
    fclose(fpOut);
    if (timeStep != initialTimeSteps+numberTimeSteps) {
      cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
      exit(1);
    }
    ObjectiveCriteria(numberTimeSteps, observed, simulated, &nashSutcliffe, &rootMSError, &biasVol, &pearsCorr);
    cout.precision(2); cout.setf(ios::fixed); cout.setf(ios::showpoint); 
    cout << "\n    Sub-catchment :  " << CatchmentElement[i].GetIdentifier() << endl << endl;
    cout.width(10); 
    cout << "    Nash-Sutcliffe efficiency                           " << nashSutcliffe << endl;; 
    cout.width(10); 
    cout << "    Root mean square error                              " << rootMSError << endl;
    cout.width(10); 
    cout << "    Bias (volume error)                                 " << biasVol << endl;
    cout.width(10); 
    cout << "    Pearson's product-moment correlation coefficient    " << pearsCorr << endl;
    cout << endl;
    fout.precision(2); fout.setf(ios::fixed); fout.setf(ios::showpoint); 
    fout << "\n    Sub-catchment :  " << CatchmentElement[i].GetIdentifier() << endl << endl;
    fout.width(10); 
    fout << "    Nash-Sutcliffe efficiency                           " << nashSutcliffe << endl;; 
    fout.width(10); 
    fout << "    Root mean square error                              " << rootMSError << endl;
    fout.width(10); 
    fout << "    Bias (volume error)                                 " << biasVol << endl;
    fout.width(10); 
    fout << "    Pearson's product-moment correlation coefficient    " << pearsCorr << endl;
    fout << endl;
  }

  delete [] observed;
  delete [] simulated;
}


void WriteSubCatchmentWaterBalance(SubCatchment * CatchmentElement, int numWatc, DateTime startSimulationTime, 
                                   DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                                   int secondsPerTimeStep)
{
  FILE *fpPre,*fpTem,*fpLak,*fpSwe,*fpIns,*fpGmb,*fpGim,*fpEva,*fpRun,*fpHsd,*fpHsm,*fpHpe,*fpHuz,*fpHlz,*fpHgw;
  char fileName[100];
  int i,timeStep;
  DateTime datetime;
  for (i=0; i<numWatc; i++) {
    sprintf(fileName,"pre_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpPre = fopen(fileName, "w")) == NULL ) {
        printf("\n File %s not found!\n\n",fileName);
        exit(1);
    }
    sprintf(fileName,"tem_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpTem = fopen(fileName, "w")) == NULL ) {
        printf("\n File %s not found!\n\n",fileName);
        exit(1);
    }
    sprintf(fileName,"lak_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpLak = fopen(fileName, "w")) == NULL ) {
      printf("\n Filen %s ikke funnet!\n\n",fileName);
      exit(1);
    }
    sprintf(fileName,"swe_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpSwe = fopen(fileName, "w")) == NULL ) {
        printf("\n File %s not found!\n\n",fileName);
        exit(1);
    }
    sprintf(fileName,"ins_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpIns = fopen(fileName, "w")) == NULL ) {
        printf("\n Filen %s ikke funnet!\n\n",fileName);
        exit(1);
    }
    sprintf(fileName,"gmb_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpGmb = fopen(fileName, "w")) == NULL ) {
      printf("\n Filen %s ikke funnet!\n\n",fileName);
      exit(1);
    }
    sprintf(fileName,"gim_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpGim = fopen(fileName, "w")) == NULL ) {
      printf("\n Filen %s ikke funnet!\n\n",fileName);
      exit(1);
    }
    sprintf(fileName,"eva_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpEva = fopen(fileName, "w")) == NULL ) {
        printf("\n File %s not found!\n\n",fileName);
        exit(1);
    }
    sprintf(fileName,"run_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpRun = fopen(fileName, "w")) == NULL ) {
        printf("\n File %s not found!\n\n",fileName);
        exit(1);
    }
    sprintf(fileName,"hsd_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpHsd = fopen(fileName, "w")) == NULL ) {
      printf("\n Filen %s ikke funnet!\n\n",fileName);
      exit(1);
    }
    sprintf(fileName,"hsm_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpHsm = fopen(fileName, "w")) == NULL ) {
        printf("\n File %s not found!\n\n",fileName);
        exit(1);
    }
    sprintf(fileName,"hpe_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpHpe = fopen(fileName, "w")) == NULL ) {
        printf("\n Filen %s ikke funnet!\n\n",fileName);
        exit(1);
    }
    sprintf(fileName,"huz_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpHuz = fopen(fileName, "w")) == NULL ) {
        printf("\n File %s not found!\n\n",fileName);
        exit(1);
    }
    sprintf(fileName,"hlz_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpHlz = fopen(fileName, "w")) == NULL ) {
        printf("\n File %s not found!\n\n",fileName);
        exit(1);
    }
    sprintf(fileName,"hgw_%08d.var",CatchmentElement[i].GetIdentifier());
    if ((fpHgw = fopen(fileName, "w")) == NULL ) {
        printf("\n File %s not found!\n\n",fileName);
        exit(1);
    }
    timeStep=initialTimeSteps;
    for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=secondsPerTimeStep) {
      // Fluxes 
      if (CatchmentElement[i].GetSubCatchmentPrecipitation(timeStep) != missingData) {
        fprintf(fpPre,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                CatchmentElement[i].GetSubCatchmentPrecipitation(timeStep)*1000.0);
        fprintf(fpTem,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                CatchmentElement[i].GetSubCatchmentTemperature(timeStep));
        fprintf(fpEva,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                CatchmentElement[i].GetSubCatchmentEvapotranspiration(timeStep)*1000.0);
        fprintf(fpRun,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                CatchmentElement[i].GetSubCatchmentRunoff(timeStep)*1000.0);
      }
      else {
        fprintf(fpPre,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
        fprintf(fpTem,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
        fprintf(fpEva,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
        fprintf(fpRun,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      }
      // Lake state variables
      if (CatchmentElement[i].GetSubCatchmentLakeStorage(timeStep) != missingData) {
	fprintf(fpLak,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
		datetime.getDay(), datetime.getHour(), datetime.getMinute(),
		CatchmentElement[i].GetSubCatchmentLakeStorage(timeStep)*1000.0);
      }
      else {
	fprintf(fpLak,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
		datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      }
    // Snow state variables
      if (CatchmentElement[i].GetSubCatchmentSnowStore(timeStep) != missingData) {
        fprintf(fpSwe,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                (CatchmentElement[i].GetSubCatchmentSnowStore(timeStep)+CatchmentElement[i].GetSubCatchmentMeltWater(timeStep))*1000.0);
        fprintf(fpIns,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                CatchmentElement[i].GetSubCatchmentWaterOutput(timeStep)*1000.0);
      }
      else {
        fprintf(fpSwe,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
        fprintf(fpIns,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      }
      // Glacier mass balance and ice melt
      if (CatchmentElement[i].GetSubCatchmentGlacierMassBalance(timeStep) != missingData) {
	fprintf(fpGmb,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
		datetime.getDay(), datetime.getHour(), datetime.getMinute(),
		CatchmentElement[i].GetSubCatchmentGlacierMassBalance(timeStep)*1000.0);
	fprintf(fpGim,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
		datetime.getDay(), datetime.getHour(), datetime.getMinute(),
		CatchmentElement[i].GetSubCatchmentGlacierIceMelt(timeStep)*1000.0);
      }
      else {
	fprintf(fpGmb,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
		datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
	fprintf(fpGim,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
		datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      }
      // Hbv state variables
      if (CatchmentElement[i].GetSubCatchmentHbvSoilMoisture(timeStep) != missingData) { 
        fprintf(fpHsd,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                CatchmentElement[i].GetSubCatchmentHbvSoilMoistureDeficit(timeStep)*1000.0);
        fprintf(fpHsm,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                CatchmentElement[i].GetSubCatchmentHbvSoilMoisture(timeStep)*1000.0);
        fprintf(fpHpe,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                CatchmentElement[i].GetSubCatchmentHbvPercSoilUpper(timeStep)*1000.0);
        fprintf(fpHuz,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                CatchmentElement[i].GetSubCatchmentHbvUpperZone(timeStep)*1000.0);
        fprintf(fpHlz,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                CatchmentElement[i].GetSubCatchmentHbvLowerZone(timeStep)*1000.0);
        fprintf(fpHgw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                (CatchmentElement[i].GetSubCatchmentHbvUpperZone(timeStep)+CatchmentElement[i].GetSubCatchmentHbvLowerZone(timeStep))*1000.0);
      }
      else {
        fprintf(fpHsd,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
        fprintf(fpHsm,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
        fprintf(fpHpe,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
        fprintf(fpHuz,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
        fprintf(fpHlz,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
        fprintf(fpHgw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      }
      timeStep++;
    }
    fclose(fpPre);
    fclose(fpSwe);
    fclose(fpLak);
    fclose(fpIns);
    fclose(fpGmb);
    fclose(fpGim);
    fclose(fpEva);
    fclose(fpRun);
    fclose(fpHsd);
    fclose(fpHsm);
    fclose(fpHpe);
    fclose(fpHuz);
    fclose(fpHlz);
    if (timeStep != initialTimeSteps+numberTimeSteps) {
      cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
      exit(1);
    }
  }
}


void WriteDistributedHbvTimeSeries(DistributedHbv * const DistHbv, int numLand, 
                                   DateTime startSimulationTime, DateTime endSimulationTime, 
                                   int initialTimeSteps, int numberTimeSteps, int secondsPerTimeStep)
{
  FILE *fpHgw;
  char fileName[100];
  int i,j,timeStep;
  DateTime datetime;

  for (i=0; i<numLand; i++) {
    for (j=0; j<DistHbv[i].GetSelectedTimeSeriesElements()->GetNumberElements(); j++) {
      if (DistHbv[i].GetLandIndex() == DistHbv[i].GetSelectedTimeSeriesElements()->GetTimeSeriesElement(j)) {
        sprintf(fileName,"HBV_groundwater_%d.var",DistHbv[i].GetLandIndex());
        if ((fpHgw = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        timeStep=initialTimeSteps;
        for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=secondsPerTimeStep) {
          // Hbv state variables
          if (DistHbv[i].GetDistributedHbvUpperZone(timeStep) != missingData) { 
            fprintf(fpHgw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    (DistHbv[i].GetDistributedHbvUpperZone(timeStep)+DistHbv[i].GetDistributedHbvLowerZone(timeStep))*1000.0);
          }
          else {
            fprintf(fpHgw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
          }
          timeStep++;
        }
        fclose(fpHgw);
        if (timeStep != initialTimeSteps+numberTimeSteps) {
          cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
          exit(1);
        }
      }
    }
  }
}


void WriteBinaryGrid(DistributedHbv * const DistHbv, DateTime datetime, int numLand, int timeStep, int nRows,  
                     int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout, int writePET)
{
  int i,j,k;
  int writeRun,writeEva,writeTem,writePre,writeSwe,writeScf,writeGmb,writeLak;
  int writeWout,writeSmw,writeHsd,writeHsm,writeHgw,writeHlz,writeHuz;
  char fileName[100];
  char dirName[100];
  char timeName[100];
  ofstream fileTemp,filePre,fileEva,filePeva,fileSwe,fileScf,fileSmw,fileRun,fileHbv,fileHsd,fileHsm,fileHgw,fileHlz,fileHuz,fileLak,fileGmb,fileWou;
  
  writeRun=writeEva=writeSwe=1;
  writeWout=writeHsd=writeHsm=writeHgw=writeHuz=1;
  writeTem=writePre=0;
  writeGmb=writeLak=writeSmw=writeHlz=writeScf=0;
   
  sprintf(dirName,"./hbv_output/");
  sprintf(timeName,"_%04d_%02d_%02d.bil",datetime.getYear(),datetime.getMonth(),datetime.getDay());
  //cout << endl << " date  " << datetime.getYear() << datetime.getMonth() << datetime.getDay() << endl;
  // cout << endl << " numLand  " << numLand << endl;

  if (writeTem == 1) {   //  Temperature in landscape elements 
    short int * temp = new short int[numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"tem");
    strcat(fileName,timeName);
    fileTemp.open(fileName,ios::out | ios::binary);
    if (!fileTemp.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetTemperature() > missingData) 
	temp[k] = (short int)((DistHbv[k].GetTemperature()) * 100); //temperature C*100
      else
	temp[k] = (short int)noData;
      k++;
    }
    fileTemp.write(reinterpret_cast<char *>(temp), sizeof(short int)*numLand);
    fileTemp.close();
    delete [] temp;
  }
  
  if (writePre == 1) { //  Precipitation in landscape elements 
    float * pre = new float [numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"pre");
    strcat(fileName,timeName);
    filePre.open(fileName,ios::out | ios::binary);
    if (!filePre.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetPrecipitation() > missingData) 
	pre[k] = (float) ((DistHbv[k].GetPrecipitation())*1000.0);
      else
	pre[k] = (float) noData;
      k++;
    }
    filePre.write(reinterpret_cast<char *>(pre), sizeof(float)*numLand);
    filePre.close();
    delete [] pre;
  }
  
  if (writeEva == 1) {    //  Evapotranspiration from landscape elements
    float * eva = new float [numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"eva");
    strcat(fileName,timeName);
    fileEva.open(fileName,ios::out | ios::binary);
    if (!fileEva.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetEvapotranspiration() > missingData)
	eva[k] = (float) ((DistHbv[k].GetEvapotranspiration()*1000.));
      else
	eva[k] = (float) noData;
      k++;
    }
    fileEva.write(reinterpret_cast<char *>(eva), sizeof(float)*numLand);
    fileEva.close();
    delete [] eva;
  }
  
  if (writePET == 1) {  //  Potential Evapotranspiration from landscape elements
    unsigned short int   * peva = new unsigned short int [numLand];
    strcpy(fileName, dirName);
    strcat(fileName, "peva");
    strcat(fileName, timeName);
    filePeva.open(fileName, ios::out | ios::binary);
    if (!filePeva.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit(1);
    }
    k = 0;
    while (k<numLand) {
      if (DistHbv[k].GetPET() > missingData) 
	peva[k] = (unsigned short int)((DistHbv[k].GetPET())*1000.0*100.);  //mm/day *100 so that I can have two digital precision
      else
	peva[k] = (unsigned short int) noData;
      k++;
    }
    filePeva.write(reinterpret_cast<char *>(peva), sizeof(unsigned short int)*numLand);
    filePeva.close();
    delete[] peva;
  }
    
  if (writeSwe == 1) { //  Snow store in landscape elements
    unsigned short int * swe = new unsigned short int[numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"swe");
    strcat(fileName,timeName);
    fileSwe.open(fileName,ios::out | ios::binary);
    if (!fileSwe.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetSnowStore() > missingData) 
      swe[k] = (unsigned short int) ((DistHbv[k].GetSnowStore() + DistHbv[k].GetMeltWater()) * 1000); //mm (ingjerd: kan ikke ha *10, da kan du gå over unsigned short int grensa på 65535
    else
      swe[k] = (unsigned short int) noData;
    k++;
    }
    fileSwe.write(reinterpret_cast<char *>(swe), sizeof(unsigned short int)*numLand);
    fileSwe.close();
    delete [] swe;
  }

  if (writeWout == 1) {  //  water to soil (out of snowpack if snow, othwerwise prec
    unsigned short int * wou = new unsigned short int[numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"wou");
    strcat(fileName,timeName);
    fileWou.open(fileName,ios::out | ios::binary);
    if (!fileWou.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetWaterOutput() > missingData) 
	wou[k] = (unsigned short int) ((DistHbv[k].GetWaterOutput()) * 1000 * 10); //mm*10
      else
	wou[k] = (unsigned short int) noData;
      k++;
    }
    fileWou.write(reinterpret_cast<char *>(wou), sizeof(unsigned short int)*numLand);
    fileWou.close();
    delete [] wou;
  }

  if (writeScf == 1) {  //  Snowcover fraction in landscape elements
    unsigned short int * scf = new unsigned short int[numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"scf");
    strcat(fileName,timeName);
    fileScf.open(fileName,ios::out | ios::binary);
    if (!fileScf.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetSnowCoverFraction() > missingData) 
	scf[k] = (unsigned short int) ((DistHbv[k].GetSnowCoverFraction()) * 10000); //prosent*100
      else
	scf[k] = (unsigned short int) noData;
      k++;
    }
    fileScf.write(reinterpret_cast<char *>(scf), sizeof(unsigned short int)*numLand);
    fileScf.close();
    delete [] scf;
  }
   
  if (writeSmw == 1) { //  Snow meltwater in landscape elements
    unsigned short int * smw = new unsigned short int[numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"smw");
    strcat(fileName,timeName);
    fileSmw.open(fileName,ios::out | ios::binary);
    if (!fileSmw.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetMeltWater() > missingData) 
	smw[k] = (unsigned short int) ((DistHbv[k].GetMeltWater()) * 10000); //mm*10
      else
	smw[k] = (unsigned short int) noData;
      k++;
    }
    fileSmw.write(reinterpret_cast<char *>(smw), sizeof(unsigned short int)*numLand);
    fileSmw.close();
    delete [] smw;
  }
 
  if (writeRun == 1) { //  Runoff from landscape elements
    float * runoff = new float [numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"run");
    strcat(fileName,timeName);
    fileRun.open(fileName,ios::out | ios::binary);
    if (!fileRun.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetRunoff() > missingData) 
	runoff[k] = (float) ((DistHbv[k].GetRunoff())*1000.0);
      else
	runoff[k] = (float) noData;
      k++;
    }
    fileRun.write(reinterpret_cast<char *>(runoff), sizeof(float)*numLand);
    fileRun.close();
    delete [] runoff;
  }
  
  //  Discharge from landscape elements
  /* float * hbv = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"hbv");
  strcat(fileName,timeName);
  fileHbv.open(fileName,ios::out | ios::binary);
  if (!fileHbv.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /* while (k<numLand) {
  if (DistHbv[k].GetDischarge() > missingData) 
  hbv[k] = (float) (DistHbv[k].GetDischarge());
  else
  hbv[k] = (float) noData;
  k++;
  }
  /* fileHbv.write(reinterpret_cast<char *>(hbv), sizeof(float)*numLand);
  fileHbv.close();
  delete [] hbv; */
  
  if (writeHsd == 1) {   //  Soil moisture deficit in landscape elements
    unsigned short int * hsd = new unsigned short int [numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"hsd");
    strcat(fileName,timeName);
    fileHsd.open(fileName,ios::out | ios::binary);
    if (!fileHsd.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetHbvSoilMoistureDeficit() > missingData) {
	if (DistHbv[k].GetHbvSoilMoistureDeficit() >= 0.0)
	  hsd[k] = (unsigned short int) ((DistHbv[k].GetHbvSoilMoistureDeficit())*1000.0*10.);
	else
	  hsd[k] = 0;
      }
      else
	hsd[k] = (unsigned short int) noData;
      k++;
    }
    fileHsd.write(reinterpret_cast<char *>(hsd), sizeof(unsigned short int)*numLand);
    fileHsd.close();
    delete [] hsd;
  }
  
  if (writeHsm == 1) {  //  Soil moisture in landscape elements
    unsigned short int * hsm = new unsigned short int [numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"hsm");
    strcat(fileName,timeName);
    fileHsm.open(fileName,ios::out | ios::binary);
    if (!fileHsm.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetHbvSoilMoisture() > missingData) {
	if (DistHbv[k].GetHbvSoilMoisture() >= 0.0)
	  hsm[k] = (unsigned short int) ((DistHbv[k].GetHbvSoilMoisture())*1000.0*10.);
	else
	  hsm[k] = (unsigned short int) 0;
      }
      else
	hsm[k] = (unsigned short int) noData;
      k++;
    }
    fileHsm.write(reinterpret_cast<char *>(hsm), sizeof(unsigned short int)*numLand);
    fileHsm.close();
    delete [] hsm;
  }
  
  if (writeHgw == 1) {  //  Groundwater in landscape elements
    unsigned short int * hgw = new unsigned short int[numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"hgw");
    strcat(fileName,timeName);
    fileHgw.open(fileName,ios::out | ios::binary);
    if (!fileHgw.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetHbvUpperZone() > missingData && DistHbv[k].GetHbvLowerZone() > missingData) 
	hgw[k] = (unsigned short int) ((DistHbv[k].GetHbvUpperZone() + DistHbv[k].GetHbvLowerZone()) * 10000); //mm *10
      else
	hgw[k] = (unsigned short int) noData;
      k++;
    }
    fileHgw.write(reinterpret_cast<char *>(hgw), sizeof(unsigned short int)*numLand);
    fileHgw.close();
    delete [] hgw;
  }
  
  if (writeHlz == 1) {  //  Lower zone groundwater in landscape elements
    unsigned short int * hlz = new unsigned short int[numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"hlz");
    strcat(fileName,timeName);
    fileHlz.open(fileName,ios::out | ios::binary);
    if (!fileHlz.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetHbvLowerZone() > missingData) { 
	hlz[k] = (unsigned short int) ((DistHbv[k].GetHbvLowerZone()) * 10000); //mm *10
      }
      else {
	hlz[k] = (unsigned short int) noData;
      }
      k++;
    }
    fileHlz.write(reinterpret_cast<char *>(hlz), sizeof(unsigned short int)*numLand);
    fileHlz.close();
    delete [] hlz; 
  }

  if (writeHuz == 1) {  //  Upper zone groundwater in landscape elements
    unsigned short int * huz = new unsigned short int[numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"huz");
    strcat(fileName,timeName);
    fileHuz.open(fileName,ios::out | ios::binary);
    if (!fileHuz.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetHbvUpperZone() > missingData) { 
	huz[k] = (unsigned short int) ((DistHbv[k].GetHbvUpperZone()) * 10000); //mm *10
      }
      else {
	huz[k] = (unsigned short int) noData;
      }
      k++;
    }
    fileHuz.write(reinterpret_cast<char *>(huz), sizeof(unsigned short int)*numLand);
    fileHuz.close();
    delete [] huz; 
  }

  
  if (writeLak == 1) {  //Lake storage in landscape elements
    unsigned short int * lak = new unsigned short int[numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"lak");
    strcat(fileName,timeName);
    fileLak.open(fileName,ios::out | ios::binary);
    if (!fileLak.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetLakeStorage() > missingData)  
	lak[k] = (unsigned short int) ((DistHbv[k].GetLakeStorage()) * 10000); //mm *10
      else 
	lak[k] = (unsigned short int) noData;
      k++;
    }
    fileLak.write(reinterpret_cast<char *>(lak), sizeof(unsigned short int)*numLand);
    fileLak.close();
    delete [] lak; 
  }
  
   if (writeGmb == 1) {  //Glacier mass balance in landscape elements
    unsigned short int * gmb = new unsigned short int[numLand];
    strcpy(fileName,dirName);
    strcat(fileName,"gmb");
    strcat(fileName,timeName);
    fileGmb.open(fileName,ios::out | ios::binary);
    if (!fileGmb.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    k=0;
    while (k<numLand) {
      if (DistHbv[k].GetGlacierMassBalance() > missingData)  
	gmb[k] = (unsigned short int) ((DistHbv[k].GetGlacierMassBalance()) * 10000); //mm *10
      else 
	gmb[k] = (unsigned short int) noData;
      k++;
    }
    fileGmb.write(reinterpret_cast<char *>(gmb), sizeof(unsigned short int)*numLand);
    fileGmb.close();
    delete [] gmb; 
  }
    
}


void WriteAsciiGrid(DistributedHbv * const DistHbv, DateTime datetime, int numLand, int timeStep, int nRows,  
                    int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
  int i,j,k;
  char fileName[100];
  char timeName[100];

  sprintf(timeName,"_%04d_%02d_%02d.asc",datetime.getYear(),datetime.getMonth(),datetime.getDay());

  //  fout << "Temperature in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv_output/tem");
  strcat(fileName,timeName);
  ofstream fileTemp(fileName);
  if (!fileTemp.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
            fileTemp.width(15); fileTemp.precision(5); fileTemp.setf(ios::showpoint); fileTemp.setf(ios::fixed); 
            fileTemp <<  i << "\t" << j << "\t"  << k << "\t"  <<  (DistHbv[k].GetTemperature()) << endl;
          k++;
        }
      }
    }
  }
  fileTemp.close();

  fout << "Precipitation in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv_output/pre");
  strcat(fileName,timeName);
  ofstream filePrec(fileName);
  if (!filePrec.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  //  cout << endl << " file " << fileName << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
            filePrec.width(15); filePrec.precision(5); filePrec.setf(ios::showpoint); filePrec.setf(ios::fixed); 
            filePrec <<  i << "\t" << j << "\t"  << k << "\t"  <<  (DistHbv[k].GetPrecipitation())*1000.0 << endl;
          k++;
        }
      }
    }
  }
  filePrec.close();

  //  fout << "Evapotranspiration from landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv_output/eva");
  strcat(fileName,timeName);
  ofstream fileEvap(fileName);
  if (!fileEvap.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
            fileEvap.width(15); fileEvap.precision(5); fileEvap.setf(ios::showpoint); fileEvap.setf(ios::fixed); 
	    fileEvap <<  i << "\t" << j << "\t"  << k << "\t"  << (DistHbv[k].GetEvapotranspiration())*1000.0 << endl;
         k++;
        }
      }
    }
  }
  fileEvap.close();
  
  //  fout << "Snow store in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv_output/swe");
  strcat(fileName,timeName);
  ofstream fileSwe(fileName);
  if (!fileSwe.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
            fileSwe.width(15); fileSwe.precision(5); fileSwe.setf(ios::showpoint); fileSwe.setf(ios::fixed); 
            fileSwe <<  i << "\t" << j << "\t"  << k << "\t"  << (DistHbv[k].GetSnowStore()+DistHbv[k].GetMeltWater())*1000.0 << endl;
          k++;
        }
      }
    }
  }
  fileSwe.close();
  
  //  fout << "Snowcover fraction in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv_output/scf");
  strcat(fileName,timeName);
  ofstream fileScf(fileName);
  if (!fileScf.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
            fileScf.width(15); fileScf.precision(5); fileScf.setf(ios::showpoint); fileScf.setf(ios::fixed); 
            fileScf <<  i << "\t" << j << "\t"  << k << "\t"  << (DistHbv[k].GetSnowCoverFraction())*100.0 << endl;
          k++;
        }
      }
    }
  }
  fileScf.close();
  
  //  fout << "Snow meltwater in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv_output/smw");
  strcat(fileName,timeName);
  ofstream fileSmw(fileName);
  if (!fileSmw.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
            fileSmw.width(15); fileSmw.precision(5); fileSmw.setf(ios::showpoint); fileSmw.setf(ios::fixed); 
            fileSmw <<  i << "\t" << j << "\t"  << k << "\t"  << (DistHbv[k].GetMeltWater())*1000.0 << endl;
          k++;
        }
      }
    }
  }
  fileSmw.close();
  
  //  fout << "Runoff from landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv_output/run");
  strcat(fileName,timeName);
  ofstream fileRunoff(fileName);
  if (!fileRunoff.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
            fileRunoff.width(15); fileRunoff.precision(5); fileRunoff.setf(ios::showpoint); fileRunoff.setf(ios::fixed); 
            fileRunoff <<  i << "\t" << j << "\t"  << k << "\t"  << (DistHbv[k].GetRunoff())*1000.0 << endl;
          k++;
        }
      }
    }
  }
  fileRunoff.close();
  
  //  fout << "Discharge from landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv_output/hbv");
  strcat(fileName,timeName);
  ofstream fileDisch(fileName);
  if (!fileDisch.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
            fileDisch.width(15); fileDisch.precision(5); fileDisch.setf(ios::showpoint); fileDisch.setf(ios::fixed); 
            fileDisch <<  i << "\t" << j << "\t"  << k << "\t"  << DistHbv[k].GetDischarge() << endl;
          k++;
        }
      }
    }
  }
  fileDisch.close();
  
  //  fout << "Soil moisture deficit in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv_output/hsd");
  strcat(fileName,timeName);
  ofstream fileHsd(fileName);
  if (!fileHsd.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
            fileHsd.width(15); fileHsd.precision(5); fileHsd.setf(ios::showpoint); fileHsd.setf(ios::fixed);
            fileHsd <<  i << "\t" << j << "\t"  << k << "\t"  << (DistHbv[k].GetHbvSoilMoistureDeficit())*1000.0 << endl;
          k++;
        }
      }
    }
  }
  fileHsd.close();
  
  //  fout << "Soil moisture in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv_output/hsm");
  strcat(fileName,timeName);
  ofstream fileHsm(fileName);
  if (!fileHsm.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
            fileHsm.width(15); fileHsm.precision(5); fileHsm.setf(ios::showpoint); fileHsm.setf(ios::fixed);
            fileHsm <<  i << "\t" << j << "\t"  << k << "\t"  << (DistHbv[k].GetHbvSoilMoisture())*1000.0 << endl;
          k++;
        }
      }
    }
  }
  fileHsm.close();
  
  //  fout << "Upper zone in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv_output/huz");
  strcat(fileName,timeName);
  ofstream fileHuz(fileName);
  if (!fileHuz.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
            fileHuz.width(15); fileHuz.precision(5); fileHuz.setf(ios::showpoint); fileHuz.setf(ios::fixed);
            fileHuz <<  i << "\t" << j << "\t"  << k << "\t"  << (DistHbv[k].GetHbvUpperZone())*1000.0 << endl;
          k++;
        }
      }
    }
  }
  fileHuz << endl;
  fileHuz.close();
  
  //  fout << "Lower zone in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv_output/hlz");
  strcat(fileName,timeName);
  ofstream fileHlz(fileName);
  if (!fileHlz.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
            fileHlz.width(15); fileHlz.precision(5); fileHlz.setf(ios::showpoint); fileHlz.setf(ios::fixed);
            fileHlz <<  i << "\t" << j << "\t"  << k << "\t"  << (DistHbv[k].GetHbvLowerZone())*1000.0 << endl;
          k++;
        }
      }
    }
  }
  fileHlz.close();
}


void WriteAsciiGridWaterBalance(DistributedHbv * const DistHbv, DateTime startSimulationTime, DateTime endSimulationTime, int numLand, int nRows,  
				int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
    int i,j,k;
    double numberDaysPerYear=365.25;
    char fileName[100];
    char dirName[100];
    char timeName[100];
  
    sprintf(dirName,"./hbv_output/");
    sprintf(timeName,"_%04d_%02d_%02d_%04d_%02d_%02d.asc",startSimulationTime.getYear(),startSimulationTime.getMonth(),startSimulationTime.getDay(),
  	  endSimulationTime.getYear(),endSimulationTime.getMonth(),endSimulationTime.getDay());
  
    //  fout << "Mean annual precipitation from landscape elements :\n";
    strcpy(fileName,dirName);
    strcat(fileName,"meanAnnualPrecipitation");
    strcat(fileName,timeName);
    ofstream filePrec(fileName);
    if (!filePrec.is_open()) 
    {
        cout << endl << "Error opening file " << fileName << endl << endl;
        exit (1);
    }
    filePrec.precision(0); 
    filePrec.setf(ios::fixed); 
    filePrec << "ncols         " << nCols << endl;
    filePrec << "nrows         " << nRows << endl;
    filePrec << "xllcorner     " << xllCorner << endl;
    filePrec << "yllcorner     " << yllCorner << endl;
    filePrec << "cellsize      " << cellSize << endl;
    filePrec << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++) 
    {
        for (j = 0; j < nCols; j++) 
	{
            if (k < numLand) 
	    {
	        if (DistHbv[k].GetGeoIndex() == ELEMENT(i,j)) 
		{
		    if (DistHbv[k].GetSumPrecipitation() > missingData) 
		    {
		        filePrec.width(15); 
			filePrec.precision(5); 
			filePrec.setf(ios::showpoint); 
			filePrec.setf(ios::fixed); 
			filePrec << (DistHbv[k].GetSumPrecipitation()/DistHbv[k].GetNumberSum())*numberDaysPerYear*1000.0 << endl;
		    }
		    else 
		    {
  		        filePrec.width(15); filePrec << noData << endl;
		    }
		    k++;
		}
		else 
		{
		    filePrec.width(15); 
		    filePrec << noData << endl;
		}
	    }
	    else 
	    {
	        filePrec.width(15); 
		filePrec << noData << endl;
	    }
	}
    }
    filePrec.close();
  
    //  fout << "Mean annual evapotranspiration from landscape elements :\n";
    strcpy(fileName,dirName);
    strcat(fileName,"meanAnnualEvapotranspiration");
    strcat(fileName,timeName);
    ofstream fileEvap(fileName);
    if (!fileEvap.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    fileEvap.precision(0); 
    fileEvap.setf(ios::fixed); 
    fileEvap << "ncols         " << nCols << endl;
    fileEvap << "nrows         " << nRows << endl;
    fileEvap << "xllcorner     " << xllCorner << endl;
    fileEvap << "yllcorner     " << yllCorner << endl;
    fileEvap << "cellsize      " << cellSize << endl;
    fileEvap << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++) 
    {
        for (j = 0; j < nCols; j++) 
	{
            if (k < numLand) 
	    {
	        if (DistHbv[k].GetGeoIndex() == ELEMENT(i,j)) 
		{
		    if (DistHbv[k].GetSumEvapotranspiration() > missingData) 
		    {
		        fileEvap.width(15); 
			fileEvap.precision(5); 
			fileEvap.setf(ios::showpoint); 
			fileEvap.setf(ios::fixed); 
			fileEvap << (DistHbv[k].GetSumEvapotranspiration()/DistHbv[k].GetNumberSum())*numberDaysPerYear*1000.0 << endl;
		    }
		    else 
		    {
		        fileEvap.width(15); 
			fileEvap << noData << endl;
                     }
		    k++;
		}
		else 
		{
		    fileEvap.width(15); 
		    fileEvap << noData << endl;
		}
	    }
	    else 
	    {
	        fileEvap.width(15); 
		fileEvap << noData << endl;
	    }
	}
    }
    fileEvap.close();
  
    //  fout << "Mean annual runoff from landscape elements :\n";
    strcpy(fileName,dirName);
    strcat(fileName,"meanAnnualRunoff");
    strcat(fileName,timeName);
    ofstream fileRunoff(fileName);
    if (!fileRunoff.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    fileRunoff.precision(0); 
    fileRunoff.setf(ios::fixed); 
    fileRunoff << "ncols         " << nCols << endl;
    fileRunoff << "nrows         " << nRows << endl;
    fileRunoff << "xllcorner     " << xllCorner << endl;
    fileRunoff << "yllcorner     " << yllCorner << endl;
    fileRunoff << "cellsize      " << cellSize << endl;
    fileRunoff << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++) 
    {
        for (j = 0; j < nCols; j++) 
	{
	    if (k < numLand) 
	    {
	        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) 
		{
		    if (DistHbv[k].GetSumRunoff() > missingData) 
		    {
		        fileRunoff.width(15); 
			fileRunoff.precision(5); 
			fileRunoff.setf(ios::showpoint); 
			fileRunoff.setf(ios::fixed); 
			fileRunoff << (DistHbv[k].GetSumRunoff()/DistHbv[k].GetNumberSum())*numberDaysPerYear*1000.0 << endl;
		    }
		    else 
		    {
		        fileRunoff.width(15); 
			fileRunoff << noData << endl;
		    }
		    k++;
		}
		else 
		{
		    fileRunoff.width(15); 
		    fileRunoff << noData << endl;
		}
	    }
	    else 
	    {
	        fileRunoff.width(15); 
		fileRunoff << noData << endl;
	    }
	}
    }
    fileRunoff.close();
  
    //  fout << "Change in storage except glacier storage for landscape elements :\n";
    strcpy(fileName,dirName);
    strcat(fileName,"storageChange");
    strcat(fileName,timeName);
    ofstream fileChange(fileName);
    if (!fileChange.is_open()) {
      cout << endl << "Error opening file " << fileName << endl << endl;
      exit (1);
    }
    fileChange.precision(0); 
    fileChange.setf(ios::fixed); 
    fileChange << "ncols         " << nCols << endl;
    fileChange << "nrows         " << nRows << endl;
    fileChange << "xllcorner     " << xllCorner << endl;
    fileChange << "yllcorner     " << yllCorner << endl;
    fileChange << "cellsize      " << cellSize << endl;
    fileChange << "NODATA_value  " << noData << endl;
    k = 0;
    for (i = 0; i < nRows; i++) 
    {
        for (j = 0; j < nCols; j++) 
	{
            if (k < numLand) 
            {
                if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) 
                {
                    if (DistHbv[k].GetSumRunoff() > missingData) 
                    {
                        fileChange.width(15); 
                        fileChange.precision(5); 
                        fileChange.setf(ios::showpoint); 
                        fileChange.setf(ios::fixed); 
                        fileChange << (DistHbv[k].GetFinalStorage()-DistHbv[k].GetInitialStorage())*1000.0 << endl;
                    }
                    else 
                    {
                        fileChange.width(15); 
                        fileChange << noData << endl;
                    }
                    k++;
                }
                else 
                {
                    fileChange.width(15); fileChange << noData << endl;
                }
            }
            else 
            {
                fileChange.width(15); fileChange << noData << endl;
            }
	}
    }
    fileChange.close();
}


/*  Objective criteria of model performance  */
void ObjectiveCriteria(int numberTimeSteps, double *obs, double *sim, 
                        double *ns, double *rmse, double *bias, double *pears)
{
  int i, sum_count;
  double obs_mean, sim_mean;
  double sum_1;
  double sum_2;
  double sum_3;

  *ns = missingData;            /*  Nash-Sutcliffe efficiency criterion  */
  *rmse = missingData;          /*  Root mean square error criterion  */
  *bias = missingData;          /*  Bias (volume error) criterion  */
  *pears = missingData;         /*  Pearson's product-moment correlation coefficient  */
  /*    for (i = 0; i < numberTimeSteps; i++) sim[i] = obs[i];*/

  /* Mean values */
  obs_mean = 0; 
  sim_mean = 0;
  sum_count = 0;
  for (i = 0; i < numberTimeSteps; i++) {
    if (sim[i] > missingData && obs[i] > missingData) {
      obs_mean = obs_mean + obs[i];
      sim_mean = sim_mean + sim[i];
      sum_count++;
    }
  }
  if (sum_count > 0) {
    obs_mean = obs_mean / (double)sum_count;
    sim_mean = sim_mean / (double)sum_count;
    /* Nash-Sutcliffe efficiency criterion */
    sum_1 = 0;
    sum_2 = 0;
    for (i = 0; i < numberTimeSteps; i++) {
      if (sim[i] > missingData && obs[i] > missingData) {
        sum_1 = sum_1 + pow((obs[i] - sim[i]),2);
        sum_2 = sum_2 + pow((obs[i] - obs_mean),2);
      }
    }
    *ns = 1.0-(sum_1/sum_2);
    /* Root mean square error criterion */
    sum_1 = 0;
    for (i = 0; i < numberTimeSteps; i++) {
      if (sim[i] > missingData && obs[i] > missingData) {
        sum_1 = sum_1 + pow((obs[i] - sim[i]),2);
      }
    }
    *rmse = pow((sum_1/(double)i),0.5);
    /* Bias (volume error) criterion */
    sum_1 = 0;
    sum_2 = 0;
    for (i = 0; i < numberTimeSteps; i++) {
      if (sim[i] > missingData && obs[i] > missingData) {
        sum_1 = sum_1 + (sim[i] - obs[i]);
        sum_2 = sum_2 + obs[i];
      }
    }
    *bias = sum_1/sum_2;
    /* Pearson's product-moment correlation coefficient criterion */
    sum_1 = 0;
    sum_2 = 0;
    sum_3 = 0;
    for (i = 0; i < numberTimeSteps; i++) {
      if (sim[i] > missingData && obs[i] > missingData) {
        sum_1 = sum_1 + (obs[i] - obs_mean) * (sim[i] - sim_mean);
        sum_2 = sum_2 + pow((obs[i] - obs_mean),2);
        sum_3 = sum_3 + pow((sim[i] - sim_mean),2);
      }
    }
    *pears = pow((sum_1 / sqrt(sum_2 * sum_3)),2);
  }    
  return;
}


void WriteModelStateVariables(DistributedHbv * const DistHbv, int numLand, DateTime datetime, int timeStep, int secondsTimestep)
{
  FILE *fileStatesOutput;
  HbvAquifer *lastHbvAquifer;
  char fileName[100];
  int i;

  datetime=datetime-secondsTimestep;
  sprintf(fileName,"modelOutputStates.txt");
  if ((fileStatesOutput=fopen(fileName,"w"))==NULL ) {
    printf("\n File %s not found!\n\n",fileName);
    exit(1);
  }
  fprintf(fileStatesOutput," Time  %7d%02d%02d/%02d%02d\n", datetime.getYear(), datetime.getMonth(), datetime.getDay(), datetime.getHour(), datetime.getMinute());
  for (i=0; i<numLand; i++) {
    fprintf(fileStatesOutput," Element no. %10d\n", i);
    if (DistHbv[i].GetLake()) {
      fprintf(fileStatesOutput," lake  %21.10f", DistHbv[i].GetLake()->GetLakeWaterBalance()->GetWaterLevel());
      fprintf(fileStatesOutput," %21.10f\n", DistHbv[i].GetLake()->GetLakeWaterBalance()->GetTemperature());
    }
    if (DistHbv[i].GetGlacier()) {
      fprintf(fileStatesOutput," gmbal %21.10f\n", DistHbv[i].GetGlacier()->GetSurfaceMassBalance());
      fprintf(fileStatesOutput," gsnow ");
      if (DistHbv[i].GetGlacier()->GetSnow()->GetLandSurfacePar()->GetCV_SNOW() == 0.0) {
	fprintf(fileStatesOutput,"%21.10f", DistHbv[i].GetGlacier()->GetSnow()->GetSnowStore());
	fprintf(fileStatesOutput,"%21.10f", DistHbv[i].GetGlacier()->GetSnow()->GetMeltWater());
      } 
      else {
	for (int j=0; j<numberSnowClasses; j++) {
	  fprintf(fileStatesOutput,"%21.10f", DistHbv[i].GetGlacier()->GetSnow()->GetDistSnowStore(j));
	  fprintf(fileStatesOutput,"%21.10f", DistHbv[i].GetGlacier()->GetSnow()->GetDistMeltWater(j));
	} 
      } 
      fprintf(fileStatesOutput,"\n");
      fprintf(fileStatesOutput," gsoil %21.10f\n", DistHbv[i].GetGlacier()->GetHBV()->GetSoilMoisture());
      fprintf(fileStatesOutput," guzon %21.10f\n", DistHbv[i].GetGlacier()->GetHBV()->GetUpperZone());
      fprintf(fileStatesOutput," glzon %21.10f\n", DistHbv[i].GetGlacier()->GetHBV()->GetLowerZone());
    }
    if (DistHbv[i].GetHbvAquifer()) {
      lastHbvAquifer=DistHbv[i].GetHbvAquifer();
      while (lastHbvAquifer) {
	fprintf(fileStatesOutput," veg   %21.10f\n", lastHbvAquifer->GetVegetation()->GetInterceptionStore());
	fprintf(fileStatesOutput," snow  ");
	if (lastHbvAquifer->GetSnow()->GetLandSurfacePar()->GetCV_SNOW() == 0.0) {
	  fprintf(fileStatesOutput,"%21.10f", lastHbvAquifer->GetSnow()->GetSnowStore());
	  fprintf(fileStatesOutput,"%21.10f", lastHbvAquifer->GetSnow()->GetMeltWater());
	} 
	else {
	  for (int j=0; j<numberSnowClasses; j++) {
	    fprintf(fileStatesOutput,"%21.10f", lastHbvAquifer->GetSnow()->GetDistSnowStore(j));
	    fprintf(fileStatesOutput,"%21.10f", lastHbvAquifer->GetSnow()->GetDistMeltWater(j));
	  } 
	} 
	fprintf(fileStatesOutput,"\n");
	fprintf(fileStatesOutput," soil  %21.10f\n", lastHbvAquifer->GetHBV()->GetSoilMoisture());
	fprintf(fileStatesOutput," upzon %21.10f\n", lastHbvAquifer->GetHBV()->GetUpperZone());
	fprintf(fileStatesOutput," lozon %21.10f\n", lastHbvAquifer->GetHBV()->GetLowerZone());
	lastHbvAquifer=lastHbvAquifer->GetNextHbvAquifer();
      }
    }
  }

  fclose(fileStatesOutput);
}
  

void ReadModelStateVariables(DistributedHbv * const DistHbv, int numLand)
{
  char ch;
  char buffer[256];
  int i,j,landIndex;
  double lakeWaterLevel,lakeTemperature,surfaceMassBalance,snowStore,meltWater,soilWater,upperZone,lowerZone,interceptionStore,totalStorage;
  double distSnowStore[numberSnowClasses];   /*  Distributed snow store (m)  */
  double distMeltWater[numberSnowClasses];   /*  Distributed meltwater in snow (m)  */
  HbvAquifer *lastHbvAquifer;

  cout << "\n  Read model states\n\n";
  ifstream fileStatesInput("modelInputStates.txt");  // Open for reading
  if (!fileStatesInput.is_open()) {
    cout << endl << " Error opening file modelInputStates.txt" << endl << endl;
    exit(1);
  } 
  fileStatesInput.getline(buffer, 256);
  cout << buffer << endl;
  for (i=0; i<numLand; i++) {
    fileStatesInput.ignore(100,'.');
    fileStatesInput >> landIndex;
    fileStatesInput.ignore(256,'\n');
    cout << " Element no. " << landIndex << " " << i << endl;
    if (landIndex != i) {
      cout << endl << "Error reading file modelInputStates.txt " << "\telement no. " << i << "\tlandIndex" << landIndex << endl;
      exit (1);
    }
    if (DistHbv[i].GetLake()) {
      fileStatesInput >> buffer >> lakeWaterLevel >> lakeTemperature;
      DistHbv[i].GetLake()->GetLakeWaterBalance()->SetLakeValues(lakeTemperature, lakeWaterLevel);
      cout << " lakeWaterLevel " << lakeWaterLevel << endl;
      cout << " lakeTemperature " << lakeTemperature << endl;
    }
    if (DistHbv[i].GetGlacier()) {
      fileStatesInput >> buffer >> surfaceMassBalance;
      DistHbv[i].GetGlacier()->SetSurfaceMassBalance(surfaceMassBalance);
      cout << " surfaceMassBalance " << surfaceMassBalance << endl;
      if (DistHbv[i].GetGlacier()->GetSnow()->GetLandSurfacePar()->GetCV_SNOW() == 0.0) {
	fileStatesInput >> buffer >> snowStore >> meltWater;
	DistHbv[i].GetGlacier()->GetSnow()->SetSnowValues(snowStore, meltWater);
	cout << " snowStore " << snowStore << endl;
	cout << " meltWater " << meltWater << endl;
      } 
      else {
	fileStatesInput >> buffer;
	cout << " Distributed snow store " << buffer << " ";
	for (int j=0; j<numberSnowClasses; j++) {
	  fileStatesInput >> distSnowStore[j] >> distMeltWater[j];
	  DistHbv[i].GetGlacier()->GetSnow()->SeDistSnowStore(j, snowStore);
	  DistHbv[i].GetGlacier()->GetSnow()->SeDistMeltWater(j, meltWater);
	  cout << distSnowStore[j] << " " << distMeltWater[j] << " ";
	}
	cout << endl;
      } 
      fileStatesInput >> buffer >> soilWater;
      fileStatesInput >> buffer >> upperZone;
      fileStatesInput >> buffer >> lowerZone;
      DistHbv[i].GetGlacier()->GetHBV()->SetSubSurfaceHbvStore(soilWater, upperZone, lowerZone);
      cout << " soilWater " << soilWater << endl;
      cout << " upperZone " << upperZone << endl;
      cout << " lowerZone " << lowerZone << endl;
    }
    if (DistHbv[i].GetHbvAquifer()) {
      lastHbvAquifer=DistHbv[i].GetHbvAquifer();
      while (lastHbvAquifer) {
	fileStatesInput >> buffer >> interceptionStore;
	cout << " interceptionStore " << interceptionStore << endl;
	if (lastHbvAquifer->GetSnow()->GetLandSurfacePar()->GetCV_SNOW() == 0.0) {
	  fileStatesInput >> buffer >> snowStore >> meltWater;
	  lastHbvAquifer->GetSnow()->SetSnowValues(snowStore, meltWater);
	  cout << " snowStore " << snowStore << endl;
	  cout << " meltWater " << meltWater << endl;
	} 
	else {
	  fileStatesInput >> buffer;
	  cout << " Distributed snow store " << buffer << " ";
	  for (int j=0; j<numberSnowClasses; j++) {
	    fileStatesInput >> distSnowStore[j] >> distMeltWater[j];
	    lastHbvAquifer->GetSnow()->SeDistSnowStore(j, snowStore);
	    lastHbvAquifer->GetSnow()->SeDistMeltWater(j, meltWater);
	    cout << distSnowStore[j] << " " << distMeltWater[j] << " ";
	  }
	  cout << endl;
	} 
	fileStatesInput >> buffer >> soilWater;
	fileStatesInput >> buffer >> upperZone;
	fileStatesInput >> buffer >> lowerZone;
	lastHbvAquifer->GetHBV()->SetSubSurfaceHbvStore(soilWater, upperZone, lowerZone);
	cout << " soilWater " << soilWater << endl;
	cout << " upperZone " << upperZone << endl;
	cout << " lowerZone " << lowerZone << endl;
	lastHbvAquifer=lastHbvAquifer->GetNextHbvAquifer();
      }
    }
  }

  cout << endl;
  fileStatesInput.close();
}
