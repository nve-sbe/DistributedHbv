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
                      int numLand, int timeStep, int nRows, int nCols, DateTime datetime, char * metPath,
                      unsigned short int * precip10, unsigned short int * temp10K, unsigned short int * tmax10K, 
		      unsigned short int * tmin10K, unsigned short int * wind10, unsigned short int * solar10, 
	              unsigned short int * vp10,  bool * inputDataFound,
                      int * indexStore, int numberIndexStore, ofstream &fout);
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
  char fileNameInput[100];
  char fileObsStreamflow[100];
  char buffer[256];
  char modelRun = 'R';
  char modelStates = 'S';
  char inputFormat = 'F';
  char evaporationModellingControl = 'T';
  char ch;
  char *metPath;
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
  double seriesWeight;
  double sumWeight;
  double sumArea;
  double precStationsWeightedElevation, tempStationsWeightedElevation;
  double correction;
  double xllCorner, yllCorner, cellSize;
  double elementArea, elementLatitude, elementElevation, elementTimefalling, elementSlopeAngle, elementAspect, lakePercent, glacierPercent;
  double areaFraction[maximumNumberLandClasses];
  DateTime datetime;
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

  /*  while (modelRun != 'S' && modelRun != 'C' && modelRun != 's' && modelRun != 'c') {
      cout << " Type of model run, simulation(S) or calibration(C): ";
      cin >> modelRun;
      cout << endl;
      }*/   
  fileControl.ignore(100,':');
  fileControl >> modelRun;
  fileControl.ignore(256,'\n');
  if (modelRun != 'S' && modelRun != 'C' && modelRun != 's' && modelRun != 'c') {
    cout << "\n Type of model run, simulation(S) or calibration(C) \n\n";
    exit (1);
  }
  if (modelRun == 'C' || modelRun == 'c') modelCalibration = true;

  /*  while (modelStates != 'R' && modelStates != 'S' && modelStates != 'r' && modelStates != 's') {
      cout << " Model states, not in use(N), read(R), write(W) or both read and write(B): ";
      cin >> modelStates;
      cout << endl;
      }*/   
  /*  fileControl.ignore(100,':');
  fileControl >> modelStates;
  fileControl.ignore(256,'\n');
  if (modelStates != 'N' && modelStates != 'R' && modelStates != 'W' && modelStates != 'B' &&
      modelStates != 'n' && modelStates != 'r' && modelStates != 'w' && modelStates != 'b') {
    cout << "\n Model states, not in use(N), read(R), write(W) or both read and write(B) \n\n";
    exit (1);
  }
  if (modelStates == 'R' || modelStates == 'r' || modelStates == 'B' || modelStates == 'b') {
    readModelStates = true;
  }
  if (modelStates == 'W' || modelStates == 'w' || modelStates == 'B' || modelStates == 'b') { 
    writeModelStates = true;
    }*/

  /*  while (inputFormat != 'G' && inputFormat != 'T' && inputFormat != 'g' && inputFormat != 't') {
      cout << " Input data format, grid files(G) or time series file(T): ";
      cin >> inputFormat;
      cout << endl;
      }*/   
  fileControl.ignore(100,':');
  fileControl >> inputFormat;
  fileControl.ignore(256,'\n');
  if (inputFormat != 'G' && inputFormat != 'T' && inputFormat != 'g' && inputFormat != 't') {
    cout << "\n Input data format, grid files(G) or time series file(T) \n\n";
    exit (1);
  }

  /*  while (evaporationModellingControl != 'T' && evaporationModellingControl != 'M' && evaporationModellingControl != 't' && evaporationModellingControl != 't') {
      cout << " Potential evaporation, temperature index (T) or long-term mean monthly values (M): ";
      cin >> evaporationModellingControl;
      cout << endl;
      }*/
  /*  fileControl.ignore(256,':');
  fileControl >> evaporationModellingControl;
  fileControl.ignore(1024,'\n');
  if (evaporationModellingControl != 'T' && evaporationModellingControl != 'M' && evaporationModellingControl != 't' && evaporationModellingControl != 'm') {
    cout << "\n Potential evaporation, temperature index (T) or long-term mean monthly values (M) \n\n";
    exit (1);
    }*/
  // Object for storing evaporation modelling information
  EvaporationControl * EvaporationControlObject = new EvaporationControl;
  EvaporationControlObject->SetEvaporationModellingControl(evaporationModellingControl);

  /*  cout << " Output file: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ofstream fout(fileName);  // Open for writing
  
  fileControl.ignore(100, ':');
  fileControl >> writePET;
  fileControl.ignore(256, '\n');

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
  /*  cout << endl;
      cout << " Start simulation date and time (day, month, year, hour, minute): ";
      cin >> day >> mth >> year >> hour >> minute;*/
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
  //  initialTimeSteps = (int)((startSimulationTime-startModelTime)/(double)ParGeneralStore->GetSECONDS_TIMESTEP());
  //  numberTimeSteps = 1+(int)((endSimulationTime-startSimulationTime)/(double)ParGeneralStore->GetSECONDS_TIMESTEP());

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
    /*    finLandScape >> landIndex >> geoIndex >> elementArea >> elementElevation >> lakePercent >> glacierPercent;
    elementSlopeAngle = 0.0;
    elementAspect = 0.0;*/
    finLandScape >> landIndex >> geoIndex >> elementArea >> elementLatitude >> elementElevation >> elementTimefalling;
    finLandScape >> elementSlopeAngle >> elementAspect >> lakePercent >> glacierPercent;
    //  cout << landIndex << "  " << geoIndex << "  " << elementArea << "  " << elementElevation << "  " << elementTimefalling << "  ";
    //  cout << elementSlopeAngle << "  " << elementAspect << "  " << lakePercent << "  " << glacierPercent << "  ";
    sumArea = lakePercent + glacierPercent;
    for (j=0; j<maximumNumberLandClasses; j++) {
      finLandScape >> landSurf >> soil >> areaFraction[j]; 
      //      finLandScape >> landSurf >> areaFraction[j]; 
      landSurfType[j]=LANDSURFACE(landSurf);
      //      soil = landSurf;
      soilType[j]=SOIL(soil);
      sumArea = sumArea + areaFraction[j];
    }
    //    for (j=0; j<maximumNumberLandClasses; j++) {
    //      cout << landSurfType[j] << "  " << soilType[j] << "  " << areaFraction[j] << "  "; 
    //    }
    //    cout << endl;
    if (sumArea != 100.0) {
      lakePercent = lakePercent*100.0/sumArea;
      glacierPercent = glacierPercent*100.0/sumArea;
      for (j=0; j<maximumNumberLandClasses; j++) 
        areaFraction[j] = areaFraction[j]*100.0/sumArea;
    }
    //    cout << lakePercent << "  " << glacierPercent << "  ";
    //    for (j=0; j<maximumNumberLandClasses; j++) {
    //      cout << areaFraction[j] << "  "; 
    //    }
    //    cout << endl;
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
      //      cout << " tempStationsWeightedElevation = " << tempStationsWeightedElevation << endl;
        
      /*      for (j=0; j<ParGeneralStore->GetNUM_PREC_SERIES(); j++) {
              cout << "P " << DistHbv[landIndex].GetMetSeriesNumber(j) << "  " << DistHbv[landIndex].GetMetSeriesWeight(j) << endl;
              }
              for (j=0; j<ParGeneralStore->GetNUM_TEMP_SERIES(); j++) {
              cout << "T " << DistHbv[landIndex].GetMetSeriesNumber(ParGeneralStore->GetNUM_PREC_SERIES()+j) << 
              "  " << DistHbv[landIndex].GetMetSeriesWeight(ParGeneralStore->GetNUM_PREC_SERIES()+j) << endl;
              }*/
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
  //  printf("\n buffer = %s",buffer);
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
  //  printf("\n buffer = %s\n",buffer);
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

  // Time series format input data
  if (inputFormat == 'T' || inputFormat == 't') {
    // Object for storing input data for time series
    InputTimeSeries * InputTimeSeriesStore = new
      InputTimeSeries(initialTimeSteps+numberTimeSteps, MetStations->GetNumPrecStations()+MetStations->GetNumTempStations(),
                      startModelTime, endSimulationTime, ParGeneralStore->GetSECONDS_TIMESTEP());
    InputTimeSeriesStore->SetGeneralPar(ParGeneralStore);
    /*    cout << " File with meteorological input data: ";
          cin >> fileName;
          cout << endl;*/
    fileControl.ignore(100,':');
    fileControl >> fileNameInput;
    fileControl.ignore(256,'\n');
    //    cout << fileNameInput << endl;
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
	//WriteBinaryGrid(DistHbv, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout, writePET);
	//WriteAsciiGrid(DistHbv, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
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

  // Grid file format input data
  else { 
    // Path to grid files with meteorological input data
    metPath = getenv("METDATA");
    if (!metPath) {
      cout <<" Environment variable METDATA not defined " << endl << endl;
      exit(1);
    }
    unsigned short int * precip10 = new unsigned short int [nRows*nCols];
    unsigned short int * temp10K = new unsigned short int [nRows*nCols];
    unsigned short int * tmax10K = new unsigned short int[nRows*nCols];
    unsigned short int * tmin10K = new unsigned short int[nRows*nCols];
    unsigned short int * wind10 = new unsigned short int[nRows*nCols];
    unsigned short int * solar10 = new unsigned short int[nRows*nCols];
    unsigned short int * vp10 = new unsigned short int[nRows*nCols];

    // Water balance for all elements and time steps for for spin-up period and simulation period
    timeStep=0;
    for (datetime=startModelTime; datetime<=endSimulationTime; datetime+=ParGeneralStore->GetSECONDS_TIMESTEP()) {
      cout << " timeStep  " << timeStep << endl;
      //      cout << "  " << datetime.getYear() << "  " << datetime.getMonth() << "  " << datetime.getDay() << "  " 
      //	   << datetime.getHour() << "  " << datetime.getMinute() << "  " << datetime.getSecond() << endl;
      inputDataFound=true;
      WaterBalanceGrid(DistHbv, ParGeneralStore, InputElementStore, initialTimeSteps, numberTimeSteps, 
                       numLand, timeStep, nRows, nCols, 
                       datetime, metPath, precip10, temp10K, tmax10K, tmin10K, wind10, solar10, vp10, 
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
      //      if (timeStep == (int)(initialTimeSteps+numberTimeSteps/2.0) || timeStep == initialTimeSteps+numberTimeSteps-1) {
      WriteBinaryGrid(DistHbv, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout, writePET);
      //        WriteAsciiGrid(DistHbv, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
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
    delete [] precip10;
    delete [] temp10K;
    delete [] tmax10K;
    delete [] tmin10K;
    delete [] wind10;
    delete [] solar10;
    delete [] vp10;
  }

  // Write model state variables
  if (writeModelStates) { 
    WriteModelStateVariables(DistHbv, numLand, datetime, timeStep-1, ParGeneralStore->GetSECONDS_TIMESTEP());
  }

  // Write water balance grid
  //  WriteAsciiGridWaterBalance(DistHbv, startSimulationTime, endSimulationTime, numLand, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
 
  // Write state variable time series for landscape elements selected for output
  WriteDistributedHbvTimeSeries(DistHbv, numLand, startSimulationTime, endSimulationTime, 
                                initialTimeSteps, numberTimeSteps, ParGeneralStore->GetSECONDS_TIMESTEP());


  /*  k=0;
      fout << "\nPrecipitation correction grid:\n";
      for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {
      if (k<numLand) {
      if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
      fout.width(10); 
      fout << DistHbv[k].GetPrecipitationCorrection() << "  ";
      k++;
      }
      else {
      fout.width(10); fout << noData << "  ";
      }
      }
      else {
      fout.width(10); fout << noData << "  ";
      }
      }
      fout << endl;
      }
      fout << endl;
      k=0;
      fout << "\nTemperature correction grid:\n";
      for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {
      if (k<numLand) {
      if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
      fout.width(10); 
      fout << DistHbv[k].GetTemperatureCorrection() << "  ";
      k++;
      }
      else {
      fout.width(10); fout << noData << "  ";
      }
      }
      else {
      fout.width(10); fout << noData << "  ";
      }
      }
      fout << endl;
      }*/

  
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
    //    cout << indexPrec << "  "  << indexTemp << "  " << metStationElevation << "  " << weight << "  " << elementTemperature << endl;
    // Snow store is removed at day no. DAY_SNOW_ZERO
    /*if (ParGeneralStore->GetDAY_SNOW_ZERO() > 0 &&
        dayNumber(InputTimeSeriesStore->GetDateTime(timeStep).getYear(),
                  InputTimeSeriesStore->GetDateTime(timeStep).getMonth(),
                  InputTimeSeriesStore->GetDateTime(timeStep).getDay()) == 
		  ParGeneralStore->GetDAY_SNOW_ZERO()+leapYear(InputTimeSeriesStore->GetDateTime(timeStep).getYear())) {*/
      /*        cout << InputTimeSeriesStore->GetDateTime(timeStep).getYear() << " " 
             << InputTimeSeriesStore->GetDateTime(timeStep).getMonth() << " "
             << InputTimeSeriesStore->GetDateTime(timeStep).getDay() << " "
             << InputTimeSeriesStore->GetDateTime(timeStep).getHour() << " "
             << InputTimeSeriesStore->GetDateTime(timeStep).getMinute() << " "
             << "    " << DistHbv[i].GetSnowStore() << endl;*/
      //        DistHbv[i].SetSnowStore(0.0);
        //      DistHbv[i].SetSubSurfaceHbvStore(0.2,0.0,0.05);
    //    }
    if (elementPrecipitation > missingData && elementTemperature > missingData) {
//      *inputDataFound=true;
      InputElementStore->SetInput(0,elementPrecipitation*DistHbv[i].GetPrecipitationCorrection());
      InputElementStore->SetInput(1,elementTemperature+DistHbv[i].GetTemperatureCorrection());
      //      cout << elementPrecipitation*1000 << "  " << elementTemperature << endl;
      /*      if (dayNumber(InputTimeSeriesStore->GetYear(timeStep),InputTimeSeriesStore->GetMth(timeStep),
              InputTimeSeriesStore->GetDay(timeStep)) == ParGeneralStore->GetDAY_SNOW_ZERO()+
              leapYear(InputTimeSeriesStore->GetYear(timeStep))) {
              InputElementStore->SetInput(0,0.0);
              }*/
      dayofyear_PM2 = dayNumber(InputTimeSeriesStore->GetDateTime(timeStep).getYear(),
				InputTimeSeriesStore->GetDateTime(timeStep).getMonth(),
				InputTimeSeriesStore->GetDateTime(timeStep).getDay());
      DistHbv[i].WaterBalance(timeStep,datetime,initialTimeSteps,numberTimeSteps,dayofyear_PM2);
      /*      DistHbv[i].SetSumWaterBalance();
      if (timeStep == initialTimeSteps) DistHbv[i].SetInitialStorage();
      if (timeStep == numberTimeSteps) DistHbv[i].SetFinalStorage();*/
    }
    else {
      *inputDataFound=false;
    }
  }
}


void WaterBalanceGrid(DistributedHbv * DistHbv,  ParametersGeneral * ParGeneralStore, InputElement * InputElementStore,
                      int initialTimeSteps, int numberTimeSteps, 
                      int numLand, int timeStep, int nRows, int nCols, DateTime datetime, char * metPath,
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

  strcpy(precFileName,metPath);
  strcpy(tempFileName,metPath);
  strcpy(TmaxFileName, metPath);
  strcpy(TminFileName, metPath);
  strcpy(windFileName, metPath);
  strcpy(solarFileName, metPath);
  strcpy(VPFileName, metPath);

  strcat(precFileName,"/rr/");
  strcat(tempFileName,"/tm/");
  strcat(TmaxFileName, "/tmax/");
  strcat(TminFileName, "/tmin/");
  strcat(windFileName, "/wind/");
  strcat(solarFileName, "/srad/");
  strcat(VPFileName, "/vp/");
  
  //  if (datetime.getMonth() < 9)
    sprintf(hydYear,"%04d",datetime.getYear());
    //  else
    //    sprintf(hydYear,"%04d",datetime.getYear()+1);
  strcat(precFileName,hydYear);
  strcat(tempFileName,hydYear);
  strcat(TmaxFileName, hydYear);
  strcat(TminFileName, hydYear);
  strcat(windFileName, hydYear);
  strcat(solarFileName, hydYear);
  strcat(VPFileName, hydYear);

  sprintf(fileName,"/tm_%04d_%02d_%02d.bil",datetime.getYear(),datetime.getMonth(),datetime.getDay());
  strcat(tempFileName,fileName);
  //  if (datetime.getMonth() != 8) {

  sprintf(fileName,"/rr_%04d_%02d_%02d.bil",datetime.getYear(),datetime.getMonth(),datetime.getDay());
  /*  }
      else {
      sprintf(fileName,"/rr_%04d_%02d_%02d.bil",datetime.getYear(),7,datetime.getDay());
      }*/
  strcat(precFileName,fileName);
  //  cout << precFileName << "  " << tempFileName << endl;

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

  //  filePrec.read((unsigned short int*) precip10, sizeof(unsigned short int)*nRows*nCols);
  //  fileTemp.read((unsigned short int*) temp10K, sizeof(unsigned short int)*nRows*nCols);
  streamoff newPosition;
  k=0;
  //  for (i=0; i<nRows; i++) {
  //    for (j=0; j<nCols; j++) {
  for (i = 0; i<numberIndexStore; i++) {
     if (k<numLand) {
       //  if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
        while (DistHbv[k].GetGeoIndex()<indexStore[i]) {
	    cout << " element not found " << i << " " << indexStore[i] << "    " << k << " " << DistHbv[k].GetGeoIndex() << endl;
	    k++;
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

			  preci_corr = precipitation/((0.82-(0.81*exp((temperature-0.66)/1.07))/(1+exp((temperature-0.66)/1.07)))*exp(-pow((wind/4.24),1.81))
                                  +(0.81*exp((temperature-0.66)/1.07))/(1+exp((temperature-0.66)/1.07))+0.18);

			  if(temperature < 0.5) {
			    preci_corr = ParGeneralStore->GetPREC_CORR_SNOW()*preci_corr ;
			  } else {
			    preci_corr = ParGeneralStore->GetPREC_CORR_RAIN()*preci_corr ;
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
          // Snow store is removed at day no. DAY_SNOW_ZERO
          /*if (ParGeneralStore->GetDAY_SNOW_ZERO() > 0 &&
              dayNumber(datetime.getYear(),datetime.getMonth(),datetime.getDay()) == 
              ParGeneralStore->GetDAY_SNOW_ZERO()+leapYear(datetime.getYear())) {
              //      InputElementStore->SetInput(0,0.0);
              DistHbv[k].SetSnowStore(0.0);
              //      DistHbv[k].SetSubSurfaceHBVStore(0.2,0.0,0.05);
          }*/
		  dayofyear_PM2 = dayNumber(datetime.getYear(), datetime.getMonth(), datetime.getDay());



		  DistHbv[k].WaterBalance(timeStep,datetime,initialTimeSteps,numberTimeSteps,dayofyear_PM2);
	  // Water balance grid values
	  /*	  DistHbv[i].SetSumWaterBalance();
	  if (timeStep == initialTimeSteps) DistHbv[i].SetInitialStorage();
	  if (timeStep == numberTimeSteps) DistHbv[i].SetFinalStorage();*/
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
  char fileName[100];
  char dirName[100];
  char timeName[100];
  //  unsigned short int noData16=65535;
  ofstream fileTemp, filePre, fileEva, filePeva, fileSwe, fileScf, fileSmw, fileRun, fileHbv, fileHsd, fileHsm, fileHgw, fileHlz;
  sprintf(dirName,"./%04d/",datetime.getYear());
  sprintf(timeName,"_%04d_%02d_%02d.bil",datetime.getYear(),datetime.getMonth(),datetime.getDay());
   //cout << endl << " date  " << datetime.getYear() << datetime.getMonth() << datetime.getDay() << endl;
   // cout << endl << " numLand  " << numLand << endl;

    if (writePET > 0) {
  //  Temperature in landscape elements 
  // float * temp = new float [numLand];
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
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//    if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
//      if (DistHbv[k].GetAccumulatedTemperature() > missingData) 
      if (DistHbv[k].GetTemperature() > missingData) 
       // temp[k] = (float) ((DistHbv[k].GetTemperature()));
		  temp[k] = (short int)((DistHbv[k].GetTemperature()) * 100); //temperature C*100
      else
        // temp[k] = (float) noData;
		  temp[k] = (short int)noData;
      k++;
  }
/*        else {
          temp[ELEMENT(i,j)] = noData;
        }
      }
      else {
        temp[ELEMENT(i,j)] = noData;
      }
      }
      }*/
  // fileTemp.write(reinterpret_cast<char *>(temp), sizeof(float)*numLand);
     fileTemp.write(reinterpret_cast<char *>(temp), sizeof(short int)*numLand);
  fileTemp.close();
  delete [] temp;
	}
  
    if (writePET > 0) {
  //  Precipitation in landscape elements 
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
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//    if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
//      if (DistHbv[k].GetAccumulatedPrecipitation() > missingData) 
      if (DistHbv[k].GetPrecipitation() > missingData) 
        pre[k] = (float) ((DistHbv[k].GetPrecipitation())*1000.0);
      else
        pre[k] = (float) noData;
      k++;
  }
/*        else {
          pre[ELEMENT(i,j)] = noData;
          }
          }
          else {
          pre[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  filePre.write(reinterpret_cast<char *>(pre), sizeof(float)*numLand);
  filePre.close();
  delete [] pre;
	}
  
    if (writePET > 0) {
  //  Evapotranspiration from landscape elements
  unsigned short int * eva = new unsigned short int [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"eva");
  strcat(fileName,timeName);
  fileEva.open(fileName,ios::out | ios::binary);
  if (!fileEva.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//      if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
    if (DistHbv[k].GetEvapotranspiration() > missingData)
    //DistHbv[k].GetInterceptionLoss() > missingData &&
        //DistHbv[k].GetTranspSoilEvap() > missingData &&
        //DistHbv[k].GetLakeEvap() > missingData)
      //      eva[k] = (float) ((DistHbv[k].GetInterceptionLoss() + DistHbv[k].GetTranspSoilEvap() + 
      //                         DistHbv[k].GetLakeEvap())*1000.0);
      eva[k] = (unsigned short int) ((DistHbv[k].GetEvapotranspiration())*1000.0*100.);
     // cout << endl << " EP =  " << eva[k] ;
      
    else
       eva[k] = (unsigned short int) noData;
       // cout << endl << " EP =  " << eva[k] ;
        
      

    k++;
  }
/*        else {
          eva[ELEMENT(i,j)] = noData;
          }
          }
          else {
          eva[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  fileEva.write(reinterpret_cast<char *>(eva), sizeof(unsigned short int)*numLand);
  fileEva.close();
  delete [] eva;
	}

   // if (writePET < 1) {
  //  Potential Evapotranspiration from landscape elements
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
  /*  for (i=0; i<nRows; i++) {
  for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
	  //      if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
	  if (DistHbv[k].GetPET() > missingData) 
		  //&& DistHbv[k].GetLakeEvap() > missingData)
		  //      eva[k] = (float) ((DistHbv[k].GetInterceptionLoss() + DistHbv[k].GetTranspSoilEvap() + 
		  //                         DistHbv[k].GetLakeEvap())*1000.0);
		   peva[k] = (unsigned short int)((DistHbv[k].GetPET())*1000.0*100.);              //mm/day *100 so that I can have two digital precision
        // if(k==4244)      cout <<  " ETP =  " << peva[k] << endl;
                
	  else
		   peva[k] = (unsigned short int)(0);
     // cout << endl << " ETP =  " << peva[k] ;
    
	  k++;

  }

  /*        else {
  eva[ELEMENT(i,j)] = noData;
  }
  }
  else {
  eva[ELEMENT(i,j)] = noData;
  }
  }
  }*/
  filePeva.write(reinterpret_cast<char *>(peva), sizeof(unsigned short int)*numLand);
  filePeva.close();
  delete[] peva;
	//}
  
    if (writePET > 0) {
  //  Snow store in landscape elements
  // float * swe = new float [numLand];
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
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (DistHbv[k].GetSnowStore() > missingData) 
          // swe[k] = (float) ((DistHbv[k].GetSnowStore()+DistHbv[k].GetMeltWater())*1000.0);
		  swe[k] = (unsigned short int) ((DistHbv[k].GetSnowStore() + DistHbv[k].GetMeltWater()) * 1000); //mm (ingjerd: kan ikke ha *10, da kan du gå over unsigned short int grensa på 65535
      else
          // swe[k] = (float) noData;
		  swe[k] = (unsigned short int) (0);
      k++;
  }
/*        else {
          swe[ELEMENT(i,j)] = noData;
          }
          }
          else {
          swe[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  // fileSwe.write(reinterpret_cast<char *>(swe), sizeof(float)*numLand);
  fileSwe.write(reinterpret_cast<char *>(swe), sizeof(unsigned short int)*numLand);
  fileSwe.close();
  delete [] swe;
	}
  
    if (writePET > 0) {
  //  Snowcover fraction in landscape elements
  // float * scf = new float [numLand];
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
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (DistHbv[k].GetSnowCoverFraction() > missingData) 
          // scf[k] = (float) ((DistHbv[k].GetSnowCoverFraction())*100.0);
	      scf[k] = (unsigned short int) ((DistHbv[k].GetSnowCoverFraction()) * 10000); //prosent*100
      else
          // scf[k] = (float) noData;
		  scf[k] = (unsigned short int) noData;
      k++;
  }
/*        else {
          scf[ELEMENT(i,j)] = noData;
          }
          }
          else {
          scf[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  // fileScf.write(reinterpret_cast<char *>(scf), sizeof(float)*numLand);
  fileScf.write(reinterpret_cast<char *>(scf), sizeof(unsigned short int)*numLand);
  fileScf.close();
  delete [] scf;
	}
   
    if (writePET > 0) {
  //  Snow meltwater in landscape elements
  // float * smw = new float [numLand];
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
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (DistHbv[k].GetMeltWater() > missingData) 
          // smw[k] = (float) ((DistHbv[k].GetMeltWater())*1000.0);
	      smw[k] = (unsigned short int) ((DistHbv[k].GetMeltWater()) * 10000); //mm*10
      else
          // smw[k] = (float) noData;
	      smw[k] = (unsigned short int) noData;
      k++;
  }
/*        else {
          smw[ELEMENT(i,j)] = noData;
          }
          }
          else {
          smw[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  // fileSmw.write(reinterpret_cast<char *>(smw), sizeof(float)*numLand);
  fileSmw.write(reinterpret_cast<char *>(smw), sizeof(unsigned short int)*numLand);
  fileSmw.close();
  delete [] smw;
	}
 
    if (writePET > 0) {
  //  Runoff from landscape elements
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
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (DistHbv[k].GetRunoff() > missingData) 
          runoff[k] = (float) ((DistHbv[k].GetRunoff())*1000.0);
      else
          runoff[k] = (float) noData;
      k++;
  }
/*        else {
          runoff[ELEMENT(i,j)] = noData;
          }
          }
          else {
          runoff[ELEMENT(i,j)] = noData;
          }
          }
          }*/
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
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  /* while (k<numLand) {
//        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (DistHbv[k].GetDischarge() > missingData) 
          hbv[k] = (float) (DistHbv[k].GetDischarge());
      else
          hbv[k] = (float) noData;
      k++;
  }
/*        else {
          hbv[ELEMENT(i,j)] = noData;
          }
          }
          else {
          hbv[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  /* fileHbv.write(reinterpret_cast<char *>(hbv), sizeof(float)*numLand);
  fileHbv.close();
  delete [] hbv; */
  
    if (writePET > 0) {
  //  Soil moisture deficit in landscape elements
  float * hsd = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"hsd");
  strcat(fileName,timeName);
  fileHsd.open(fileName,ios::out | ios::binary);
  if (!fileHsd.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
    if (DistHbv[k].GetHbvSoilMoistureDeficit() > missingData) {
      if (DistHbv[k].GetHbvSoilMoistureDeficit() >= 0.0)
        hsd[k] = (float) ((DistHbv[k].GetHbvSoilMoistureDeficit())*1000.0);
      else
        hsd[k] = 0.0;
    }
    else
      hsd[k] = (float) noData;
    k++;
  }
/*        else {
          hsd[ELEMENT(i,j)] = noData;
          }
          }
          else {
          hsd[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  fileHsd.write(reinterpret_cast<char *>(hsd), sizeof(float)*numLand);
  fileHsd.close();
  delete [] hsd;
	}
  
    if (writePET > 0) {
  //  Soil moisture in landscape elements
  float * hsm = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"hsm");
  strcat(fileName,timeName);
  fileHsm.open(fileName,ios::out | ios::binary);
  if (!fileHsm.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
    if (DistHbv[k].GetHbvSoilMoistureDeficit() > missingData) {
      if (DistHbv[k].GetHbvSoilMoistureDeficit() >= 0.0)
        hsm[k] = (float) ((DistHbv[k].GetHbvSoilMoisture())*1000.0);
      else
        hsm[k] = 0.0;
    }
    else
      hsm[k] = (float) noData;
    k++;
  }
/*        else {
          hsm[ELEMENT(i,j)] = noData;
          }
          }
          else {
          hsm[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  fileHsm.write(reinterpret_cast<char *>(hsm), sizeof(float)*numLand);
  fileHsm.close();
  delete [] hsm;
	}
  
    if (writePET > 0) {
  //  Groundwater in landscape elements
  // float * hgw = new float [numLand];
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
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
    //        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (DistHbv[k].GetHbvUpperZone() > missingData && DistHbv[k].GetHbvLowerZone() > missingData) 
          // hgw[k] = (float) ((DistHbv[k].GetHbvUpperZone() + DistHbv[k].GetHbvLowerZone())*1000.0);
	      hgw[k] = (unsigned short int) ((DistHbv[k].GetHbvUpperZone() + DistHbv[k].GetHbvLowerZone()) * 10000); //mm *10
      else
          // hgw[k] = (float) noData;
	      hgw[k] = (unsigned short int) noData;
      k++;
  }
/*        else {
          hgw[ELEMENT(i,j)] = noData;
          }
          }
          else {
          hgw[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  // fileHgw.write(reinterpret_cast<char *>(hgw), sizeof(float)*numLand);
  fileHgw.write(reinterpret_cast<char *>(hgw), sizeof(unsigned short int)*numLand);
  fileHgw.close();
  delete [] hgw;
	}
  
  //  Lower zone groundwater in landscape elements
  /*float * hlz = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"hlz");
  strcat(fileName,timeName);
  fileHlz.open(fileName,ios::out | ios::binary);
  if (!fileHlz.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  /* while (k<numLand) {
//        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (DistHbv[k].GetHbvLowerZone() > missingData) 
          hlz[k] = (float) ((DistHbv[k].GetHbvLowerZone())*1000.0);
      else
          hlz[k] = (float) noData;
      k++;
  }
/*        else {
          hlz[ELEMENT(i,j)] = noData;
          }
          }
          else {
          hlz[ELEMENT(i,j)] = noData;
          }
          }
  }*/
  /* fileHlz.write(reinterpret_cast<char *>(hlz), sizeof(float)*numLand);
  fileHlz.close();
  delete [] hlz; */
}


void WriteAsciiGrid(DistributedHbv * const DistHbv, DateTime datetime, int numLand, int timeStep, int nRows,  
                    int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
  int i,j,k;
  char fileName[100];
  char timeName[100];

  sprintf(timeName,"_%04d_%02d_%02d.asc",datetime.getYear(),datetime.getMonth(),datetime.getDay());

  //  fout << "Temperature in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"tem");
  strcat(fileName,timeName);
  ofstream fileTemp(fileName);
  if (!fileTemp.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileTemp << "ncols         " << nCols << endl;
  fileTemp << "nrows         " << nRows << endl;
  fileTemp << "xllcorner     " << xllCorner << endl;
  fileTemp << "yllcorner     " << yllCorner << endl;
  fileTemp << "cellsize      " << cellSize << endl;
  fileTemp << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (DistHbv[k].GetTemperature() > missingData) {
            fileTemp.width(15); fileTemp.precision(5); fileTemp.setf(ios::showpoint); fileTemp.setf(ios::fixed); 
            fileTemp << (DistHbv[k].GetTemperature()) << endl;
          }
          else {
            fileTemp.width(15); fileTemp << noData << endl;
          }
          k++;
        }
        else {
          fileTemp.width(15); fileTemp << noData << endl;
        }
      }
      else {
        fileTemp.width(15); fileTemp << noData << endl;
      }
    }
  }
  fileTemp.close();

  //  fout << "Precipitation in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"pre");
  strcat(fileName,timeName);
  ofstream filePrec(fileName);
  if (!filePrec.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  filePrec << "ncols         " << nCols << endl;
  filePrec << "nrows         " << nRows << endl;
  filePrec << "xllcorner     " << xllCorner << endl;
  filePrec << "yllcorner     " << yllCorner << endl;
  filePrec << "cellsize      " << cellSize << endl;
  filePrec << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (DistHbv[k].GetPrecipitation() > missingData) {
            filePrec.width(15); filePrec.precision(5); filePrec.setf(ios::showpoint); filePrec.setf(ios::fixed); 
            filePrec << (DistHbv[k].GetPrecipitation())*1000.0 << endl;
          }
          else {
            filePrec.width(15); filePrec << noData << endl;
          }
          k++;
        }
        else {
          filePrec.width(15); filePrec << noData << endl;
        }
      }
      else {
        filePrec.width(15); filePrec << noData << endl;
      }
    }
  }
  filePrec.close();

  //  fout << "Evapotranspiration from landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"eva");
  strcat(fileName,timeName);
  ofstream fileEvap(fileName);
  if (!fileEvap.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileEvap << "ncols         " << nCols << endl;
  fileEvap << "nrows         " << nRows << endl;
  fileEvap << "xllcorner     " << xllCorner << endl;
  fileEvap << "yllcorner     " << yllCorner << endl;
  fileEvap << "cellsize      " << cellSize << endl;
  fileEvap << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (DistHbv[k].GetInterceptionLoss() > missingData &&
              DistHbv[k].GetTranspSoilEvap() > missingData &&
              DistHbv[k].GetLakeEvap() > missingData) {
            fileEvap.width(15); fileEvap.precision(5); fileEvap.setf(ios::showpoint); fileEvap.setf(ios::fixed); 
	    //            fileEvap << (DistHbv[k].GetInterceptionLoss() + DistHbv[k].GetTranspSoilEvap() + 
	    //                         DistHbv[k].GetLakeEvap())*1000.0 << endl;
	    fileEvap << (DistHbv[k].GetEvapotranspiration())*1000.0 << endl;
          }
          else {
            fileEvap.width(15); fileEvap << noData << endl;
          }
         k++;
        }
        else {
          fileEvap.width(15); fileEvap << noData << endl;
        }
      }
      else {
        fileEvap.width(15); fileEvap << noData << endl;
      }
    }
  }
  fileEvap.close();
  
  //  fout << "Snow store in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"swe");
  strcat(fileName,timeName);
  ofstream fileSwe(fileName);
  if (!fileSwe.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileSwe << "ncols         " << nCols << endl;
  fileSwe << "nrows         " << nRows << endl;
  fileSwe << "xllcorner     " << xllCorner << endl;
  fileSwe << "yllcorner     " << yllCorner << endl;
  fileSwe << "cellsize      " << cellSize << endl;
  fileSwe << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (DistHbv[k].GetSnowStore() > missingData) {
            fileSwe.width(15); fileSwe.precision(5); fileSwe.setf(ios::showpoint); fileSwe.setf(ios::fixed); 
            fileSwe << (DistHbv[k].GetSnowStore()+DistHbv[k].GetMeltWater())*1000.0 << endl;
          }
          else {
            fileSwe.width(15); fileSwe << noData << endl;
          }
          k++;
        }
        else {
          fileSwe.width(15); fileSwe << noData << endl;
        }
      }
      else {
        fileSwe.width(15); fileSwe << noData << endl;
      }
    }
  }
  fileSwe.close();
  
  //  fout << "Snowcover fraction in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"scf");
  strcat(fileName,timeName);
  ofstream fileScf(fileName);
  if (!fileScf.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileScf << "ncols         " << nCols << endl;
  fileScf << "nrows         " << nRows << endl;
  fileScf << "xllcorner     " << xllCorner << endl;
  fileScf << "yllcorner     " << yllCorner << endl;
  fileScf << "cellsize      " << cellSize << endl;
  fileScf << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (DistHbv[k].GetSnowCoverFraction() > missingData) {
            fileScf.width(15); fileScf.precision(5); fileScf.setf(ios::showpoint); fileScf.setf(ios::fixed); 
            fileScf << (DistHbv[k].GetSnowCoverFraction())*100.0 << endl;
          }
          else {
            fileScf.width(15); fileScf << noData << endl;
          }
          k++;
        }
        else {
          fileScf.width(15); fileScf << noData << endl;
        }
      }
      else {
        fileScf.width(15); fileScf << noData << endl;
      }
    }
  }
  fileScf.close();
  
  //  fout << "Snow meltwater in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"smw");
  strcat(fileName,timeName);
  ofstream fileSmw(fileName);
  if (!fileSmw.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileSmw << "ncols         " << nCols << endl;
  fileSmw << "nrows         " << nRows << endl;
  fileSmw << "xllcorner     " << xllCorner << endl;
  fileSmw << "yllcorner     " << yllCorner << endl;
  fileSmw << "cellsize      " << cellSize << endl;
  fileSmw << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (DistHbv[k].GetMeltWater() > missingData) {
            fileSmw.width(15); fileSmw.precision(5); fileSmw.setf(ios::showpoint); fileSmw.setf(ios::fixed); 
            fileSmw << (DistHbv[k].GetMeltWater())*1000.0 << endl;
          }
          else {
            fileSmw.width(15); fileSmw << noData << endl;
          }
          k++;
        }
        else {
          fileSmw.width(15); fileSmw << noData << endl;
        }
      }
      else {
        fileSmw.width(15); fileSmw << noData << endl;
      }
    }
  }
  fileSmw.close();
  
  //  fout << "Runoff from landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"run");
  strcat(fileName,timeName);
  ofstream fileRunoff(fileName);
  if (!fileRunoff.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileRunoff << "ncols         " << nCols << endl;
  fileRunoff << "nrows         " << nRows << endl;
  fileRunoff << "xllcorner     " << xllCorner << endl;
  fileRunoff << "yllcorner     " << yllCorner << endl;
  fileRunoff << "cellsize      " << cellSize << endl;
  fileRunoff << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (DistHbv[k].GetRunoff() > missingData) {
            fileRunoff.width(15); fileRunoff.precision(5); fileRunoff.setf(ios::showpoint); fileRunoff.setf(ios::fixed); 
            fileRunoff << (DistHbv[k].GetRunoff())*1000.0 << endl;
          }
          else {
            fileRunoff.width(15); fileRunoff << noData << endl;
          }
          k++;
        }
        else {
          fileRunoff.width(15); fileRunoff << noData << endl;
        }
      }
      else {
        fileRunoff.width(15); fileRunoff << noData << endl;
      }
    }
  }
  fileRunoff.close();
  
  //  fout << "Discharge from landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hbv");
  strcat(fileName,timeName);
  ofstream fileDisch(fileName);
  if (!fileDisch.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileDisch << "ncols         " << nCols << endl;
  fileDisch << "nrows         " << nRows << endl;
  fileDisch << "xllcorner     " << xllCorner << endl;
  fileDisch << "yllcorner     " << yllCorner << endl;
  fileDisch << "cellsize      " << cellSize << endl;
  fileDisch << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (DistHbv[k].GetDischarge() > missingData) {
            fileDisch.width(15); fileDisch.precision(5); fileDisch.setf(ios::showpoint); fileDisch.setf(ios::fixed); 
            fileDisch << DistHbv[k].GetDischarge() << endl;
          }
          else {
            fileDisch.width(15); fileDisch << noData << endl;
          }
          k++;
        }
        else {
          fileDisch.width(15); fileDisch << noData << endl;
        }
      }
      else {
        fileDisch.width(15); fileDisch << noData << endl;
      }
    }
  }
  fileDisch.close();
  
  //  fout << "Soil moisture deficit in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hsd");
  strcat(fileName,timeName);
  ofstream fileHsd(fileName);
  if (!fileHsd.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileHsd << "ncols         " << nCols << endl;
  fileHsd << "nrows         " << nRows << endl;
  fileHsd << "xllcorner     " << xllCorner << endl;
  fileHsd << "yllcorner     " << yllCorner << endl;
  fileHsd << "cellsize      " << cellSize << endl;
  fileHsd << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (DistHbv[k].GetHbvSoilMoisture() > missingData) {
            fileHsd.width(15); fileHsd.precision(5); fileHsd.setf(ios::showpoint); fileHsd.setf(ios::fixed);
            fileHsd << (DistHbv[k].GetHbvSoilMoistureDeficit())*1000.0 << endl;
          }
          else {
            fileHsd.width(15); fileHsd << noData << endl;
          }
          k++;
        }
        else {
          fileHsd.width(15); fileHsd << noData << endl;
        }
      }
      else {
        fileHsd.width(15); fileHsd << noData << endl;
      }
    }
  }
  fileHsd.close();
  
  //  fout << "Soil moisture in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hsm");
  strcat(fileName,timeName);
  ofstream fileHsm(fileName);
  if (!fileHsm.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileHsm << "ncols         " << nCols << endl;
  fileHsm << "nrows         " << nRows << endl;
  fileHsm << "xllcorner     " << xllCorner << endl;
  fileHsm << "yllcorner     " << yllCorner << endl;
  fileHsm << "cellsize      " << cellSize << endl;
  fileHsm << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (DistHbv[k].GetHbvSoilMoisture() > missingData) {
            fileHsm.width(15); fileHsm.precision(5); fileHsm.setf(ios::showpoint); fileHsm.setf(ios::fixed);
            fileHsm << (DistHbv[k].GetHbvSoilMoisture())*1000.0 << endl;
          }
          else {
            fileHsm.width(15); fileHsm << noData << endl;
          }
          k++;
        }
        else {
          fileHsm.width(15); fileHsm << noData << endl;
        }
      }
      else {
        fileHsm.width(15); fileHsm << noData << endl;
      }
    }
  }
  fileHsm.close();
  
  //  fout << "Upper zone in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"huz");
  strcat(fileName,timeName);
  ofstream fileHuz(fileName);
  if (!fileHuz.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileHuz << "ncols         " << nCols << endl;
  fileHuz << "nrows         " << nRows << endl;
  fileHuz << "xllcorner     " << xllCorner << endl;
  fileHuz << "yllcorner     " << yllCorner << endl;
  fileHuz << "cellsize      " << cellSize << endl;
  fileHuz << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (DistHbv[k].GetHbvUpperZone() > missingData) {
            fileHuz.width(15); fileHuz.precision(5); fileHuz.setf(ios::showpoint); fileHuz.setf(ios::fixed);
            fileHuz << (DistHbv[k].GetHbvUpperZone())*1000.0 << endl;
          }
          else {
            fileHuz.width(15); fileHuz << noData << endl;
          }
          k++;
        }
        else {
          fileHuz.width(15); fileHuz << noData << endl;
        }
      }
      else {
        fileHuz.width(15); fileHuz << noData << endl;
      }
    }
    fileHuz << endl;
  }
  fileHuz << endl;
  fileHuz.close();
  
  //  fout << "Lower zone in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,"hlz");
  strcat(fileName,timeName);
  ofstream fileHlz(fileName);
  if (!fileHlz.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileHlz << "ncols         " << nCols << endl;
  fileHlz << "nrows         " << nRows << endl;
  fileHlz << "xllcorner     " << xllCorner << endl;
  fileHlz << "yllcorner     " << yllCorner << endl;
  fileHlz << "cellsize      " << cellSize << endl;
  fileHlz << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (DistHbv[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (DistHbv[k].GetHbvLowerZone() > missingData) {
            fileHlz.width(15); fileHlz.precision(5); fileHlz.setf(ios::showpoint); fileHlz.setf(ios::fixed);
            fileHlz << (DistHbv[k].GetHbvLowerZone())*1000.0 << endl;
          }
          else {
            fileHlz.width(15); fileHlz << noData << endl;
          }
          k++;
        }
        else {
          fileHlz.width(15); fileHlz << noData << endl;
        }
      }
      else {
        fileHlz.width(15); fileHlz << noData << endl;
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
  
    sprintf(dirName,"./");
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
