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

void ReadLandUse(DistributedHbv * const DistHbv, int nRows, int nCols, int
                 noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout);
void ReadLandUseGeneral(DistributedHbv * const DistHbv, int nRows, int nCols, int
                 noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout);
void WriteLandScapeElements(DistributedHbv * const DistHbv, ParametersGeneral * const ParGeneralStore,
                            MeteorologicalStations * const MetStations, int landIndex, int nRows, int nCols, 
                            int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void DistanceSort(DistributedHbv  * const distElement, ParametersGeneral * const ParGeneralStore, 
                  MeteorologicalStations * const MetStations,
                  int * precStations, int * tempStations, double * precWeights, double * tempWeights);
void WriteCatchmentIdentifier(DistributedHbv * const DistHbv, SubCatchment * const CatchmentElement, int numWatc, ofstream &fout);
void FindCatchmentIdentifier(DistributedHbv * const DistHbv, SubCatchment * const CatchmentElement, int numWatc, 
                             int nRows, int nCols, int noData, double xllCorner, double yllCorner, 
                             double cellSize, ifstream &fileControl, ofstream &fout);
void GetRowCol(int elementNo, int nCols, int &row, int &col); 
void SetGeneralParameters(ParametersGeneral * const ParGeneralStore, ifstream &fileControl, ofstream &fout);

int main(int argc, char *argv[])
{
  std::cout << "\n\n Preprocessing for distributed HBV water balance model \n\n";

  char fileName[80];
  char buffer[256];
  char ch;
  int i,j,k;
  int landIndex;
  int nRows, nCols, noData;
  int numWatc, numWatcUp, numWatcOut;
  int subCatchmentId;
  double value;
  double correction;
  double xllCorner, yllCorner, cellSize;

  if (argc != 2) {
    cout << " " << argv[0] << "  <control file name>\n\n";
    exit(1);
  }
  ifstream fileControl(argv[1]);
  if (!fileControl.is_open()) {
    cout << " Error opening file " << argv[1] << endl << endl;
    exit (1);
  }

  /*  cout << " Output file: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ofstream fout(fileName);  // Open for writing

  // Object for storing meteorological station information
  MeteorologicalStations * MetStations = new MeteorologicalStations;
  MetStations->SetMeteorologicalStations(fileControl, fout);

  // Read common parameters file and set parameter values
  ParametersGeneral * ParGeneralStore = new ParametersGeneral;
  SetGeneralParameters(ParGeneralStore, fileControl, fout);
  // End read common parameters

  // Read landscape element file and generate landscape element objects  
  /*  cout << " File with geographical analysis area: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finGeo(fileName);  // Open for reading
  if (!finGeo.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  } 
  finGeo >> buffer >> nCols;
  finGeo >> buffer >> nRows;
  finGeo >> buffer >> xllCorner;
  finGeo >> buffer >> yllCorner;
  finGeo >> buffer >> cellSize;
  finGeo >> buffer >> noData;
  cout << endl << nCols << endl;
  cout << nRows << endl;
  cout << xllCorner << endl;
  cout << yllCorner << endl;
  cout << cellSize << endl;
  cout << buffer << "    " << noData << endl;
  landIndex = 0;
  DistributedHbv * DistHbv = new DistributedHbv [nRows*nCols];
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      DistHbv[ELEMENT(i,j)].SetGeoIndex(ELEMENT(i,j));
      finGeo >> value;
      if (value!=noData) {
        DistHbv[ELEMENT(i,j)].SetLandIndex(landIndex++);
      } else {
        DistHbv[ELEMENT(i,j)].SetLandIndex(noData);
      }
      DistHbv[ELEMENT(i,j)].SetXCoord(xllCorner + (j+0.5)*cellSize);
      DistHbv[ELEMENT(i,j)].SetYCoord(yllCorner + (nRows-i-0.5)*cellSize);
    }
  }
  finGeo.close();
  // Read end

  // Read land use for landscape elements
  //  ReadLandUse(DistHbv, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fileControl, fout);
  ReadLandUseGeneral(DistHbv, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fileControl, fout);

  // Open file with sub-catchment hierarchy
  /*  cout << " File with sub-catchment hierarchy: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream fileWCo(fileName);
  if (!fileWCo.is_open()) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  
  // Read sub-catchment information and generate sub-catchment objects  
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
    cout << " Sub-catchment index " << i << "  " << "Sub-catchment identifier " << CatchmentElement[i].GetIdentifier() << endl;
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
  // Read end

  // Write landscape element information to input file to dew
  WriteLandScapeElements(DistHbv, ParGeneralStore, MetStations, landIndex, nRows, nCols, noData, 
                         xllCorner, yllCorner, cellSize, fout);
  
  // Find sub-catchment identifiers for landscape elements and 
  // connect landscape elements to sub-catchment outlets
  FindCatchmentIdentifier(DistHbv, CatchmentElement, numWatc, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fileControl, fout);
  
  // Write information about sub-catchment outlets and landscape elements to input file to dew.cpp
  WriteCatchmentIdentifier(DistHbv, CatchmentElement, numWatc, fout);
  
  delete [] Outlet;
  delete [] CatchmentElement;
  
  fout.close();
  fileControl.close();
  delete [] DistHbv;
  return 0;
}


// Land surface classes based on potential tree level
void ReadLandUse(DistributedHbv * const DistHbv, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[256];
  int i,j,k;
  int nRo, nCo, noDa;
  double value, notAssigned;
  double xllC, yllC, cellS;

  // Area of grid cells 
  for (i=0; i<nRows*nCols; i++) {
    DistHbv[i].SetArea(cellSize*cellSize);
  }

  // Read elevation of grid cells
  /*  cout << " File with grid cell elevations: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finElevation(fileName);  // Open for reading
  if (!finElevation .is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finElevation >> buffer >> nCo;
  finElevation >> buffer >> nRo;
  finElevation >> buffer >> xllC;
  finElevation >> buffer >> yllC;
  finElevation >> buffer >> cellS;
  finElevation >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for elevation!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finElevation >> value;
    if (value<0.0) value=0.0;
    DistHbv[i].SetElevation(value);
  }
  finElevation.close();
  // Read end

  // Read slope angle of grid cells
  /*  cout << " File with slope angles: ";
      cin >> fileName;
      cout << endl;*/
  /*  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finSlope(fileName);  // Open for reading
  if (!finSlope.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finSlope >> buffer >> nCo;
  finSlope >> buffer >> nRo;
  finSlope >> buffer >> xllC;
  finSlope >> buffer >> yllC;
  finSlope >> buffer >> cellS;
  finSlope >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for slope!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finSlope >> value;
    if (value<0) value=0;
    DistHbv[i].SetSlopeAngle(value);
  }
  finSlope.close();*/
  // Read end

  // Read slope aspect of grid cells
  /*  cout << " File with slope aspects: ";
      cin >> fileName;
      cout << endl;*/
  /*  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finAspect(fileName);  // Open for reading
  if (!finAspect.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finAspect >> buffer >> nCo;
  finAspect >> buffer >> nRo;
  finAspect >> buffer >> xllC;
  finAspect >> buffer >> yllC;
  finAspect >> buffer >> cellS;
  finAspect >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for aspect!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finAspect >> value;
    if (value<0) value=0;
    DistHbv[i].SetAspect(value);
  }
  finAspect.close();*/
  // Read end

  // Read percentage of grid cells covered by lakes
  /*  cout << " File with lake percentage: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finLakePercent(fileName);  // Open for reading
  if (!finLakePercent.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finLakePercent >> buffer >> nCo;
  finLakePercent >> buffer >> nRo;
  finLakePercent >> buffer >> xllC;
  finLakePercent >> buffer >> yllC;
  finLakePercent >> buffer >> cellS;
  finLakePercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for lake percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finLakePercent >> value;
    if (value<0.0) value=0.0;
    DistHbv[i].SetLakePercent(value);
  }
  finLakePercent.close();
  // Read end

  // Read percentage of grid cells covered by forest
  /*  cout << " File with forest percentage: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finForestPercent(fileName);  // Open for reading
  if (!finForestPercent.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finForestPercent >> buffer >> nCo;
  finForestPercent >> buffer >> nRo;
  finForestPercent >> buffer >> xllC;
  finForestPercent >> buffer >> yllC;
  finForestPercent >> buffer >> cellS;
  finForestPercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for forest percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finForestPercent >> value;
    if (value<0.0) value=0.0;
    DistHbv[i].SetForestPercent(value);
  }
  finForestPercent.close();
  // Read end

  // Read percentage of grid cells covered by bogs
  /*  cout << " File with bog percentage: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finBogPercent(fileName);  // Open for reading
  if (!finBogPercent.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finBogPercent >> buffer >> nCo;
  finBogPercent >> buffer >> nRo;
  finBogPercent >> buffer >> xllC;
  finBogPercent >> buffer >> yllC;
  finBogPercent >> buffer >> cellS;
  finBogPercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for bog percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finBogPercent >> value;
    if (value<0.0) value=0.0;
    DistHbv[i].SetBogPercent(value);
  }
  finBogPercent.close();
  // Read end

  // Read percentage of grid cells covered by glaciers
  /*  cout << " File with glacier percentage: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finGlacierPercent(fileName);  // Open for reading
  if (!finGlacierPercent.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finGlacierPercent >> buffer >> nCo;
  finGlacierPercent >> buffer >> nRo;
  finGlacierPercent >> buffer >> xllC;
  finGlacierPercent >> buffer >> yllC;
  finGlacierPercent >> buffer >> cellS;
  finGlacierPercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for glacier percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finGlacierPercent >> value;
    if (value<0.0) value=0.0;
    DistHbv[i].SetGlacierPercent(value);
  }
  finGlacierPercent.close();
  // Read end

  // Read tree level for grid cells
  /*  cout << " File with tree levels: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finTreeLevel(fileName);  // Open for reading
  if (!finTreeLevel.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finTreeLevel >> buffer >> nCo;
  finTreeLevel >> buffer >> nRo;
  finTreeLevel >> buffer >> xllC;
  finTreeLevel >> buffer >> yllC;
  finTreeLevel >> buffer >> cellS;
  finTreeLevel >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for glacier percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finTreeLevel >> value;
    if (value<0.0) value=0.0;
    DistHbv[i].SetTreeLevel(value);
  }
  finTreeLevel.close();
  // Read end

  // Assign land use classes based on potential tree level
  for (i=0; i<nRows*nCols; i++) {
    if (DistHbv[i].GetLandIndex()!=noData) { 
      notAssigned = 100.0 - DistHbv[i].GetLakePercent() - DistHbv[i].GetGlacierPercent() -
          DistHbv[i].GetForestPercent() - DistHbv[i].GetBogPercent();
      if (notAssigned < 0.0) notAssigned = 0.0;
      if (DistHbv[i].GetElevation() > DistHbv[i].GetTreeLevel()) {
        if (DistHbv[i].GetForestPercent() > 0.0) {
          if (notAssigned > 0.0) {
            DistHbv[i].SetAlpinePercent(DistHbv[i].GetForestPercent() + notAssigned*0.2);
            DistHbv[i].SetHeatherPercent(notAssigned*0.8);
            DistHbv[i].SetForestPercent(0.0);
          } 
          else {
            DistHbv[i].SetAlpinePercent(DistHbv[i].GetForestPercent());
            DistHbv[i].SetForestPercent(0.0);
          }
        }
        else if (DistHbv[i].GetElevation() <= DistHbv[i].GetTreeLevel()+100.0) {
          DistHbv[i].SetAlpinePercent(notAssigned*0.2);
          DistHbv[i].SetHeatherPercent(notAssigned*0.8);
        }
        else if (DistHbv[i].GetElevation() <= DistHbv[i].GetTreeLevel()+200.0) {
          DistHbv[i].SetHeatherPercent(notAssigned*0.5);
          DistHbv[i].SetBedrockPercent(notAssigned*0.5);
        }
        else { 
          DistHbv[i].SetBedrockPercent(notAssigned);
        }
      }
      else if (DistHbv[i].GetElevation() > DistHbv[i].GetTreeLevel()*0.9) {
        if (notAssigned > 0.0) {
          DistHbv[i].SetAlpinePercent(DistHbv[i].GetForestPercent() + notAssigned*0.2);
          DistHbv[i].SetHeatherPercent(notAssigned*0.8);
          DistHbv[i].SetForestPercent(0.0);
        }
        else {
          DistHbv[i].SetAlpinePercent(DistHbv[i].GetForestPercent()*0.5);
          DistHbv[i].SetForestPercent(DistHbv[i].GetForestPercent()*0.5);
        }
      }
      else if (DistHbv[i].GetElevation() > DistHbv[i].GetTreeLevel()*0.8) {
        if (notAssigned > 0.0) {
          DistHbv[i].SetOpenLandPercent(notAssigned);
        }
        else {
          DistHbv[i].SetAlpinePercent(DistHbv[i].GetForestPercent()*0.2);
          DistHbv[i].SetForestPercent(DistHbv[i].GetForestPercent()*0.8);
        }
      }
      else {
        if (notAssigned > 0.0) {
          DistHbv[i].SetOpenLandPercent(notAssigned);
        }
      }
    }
  }

  // Control and correct area fractions of land use classes
  for (i=0; i<nRows*nCols; i++) {
    if (DistHbv[i].GetLandIndex()!=noData) { 
      value=DistHbv[i].GetForestPercent()+DistHbv[i].GetAlpinePercent()+
        DistHbv[i].GetHeatherPercent()+DistHbv[i].GetBedrockPercent()+
        DistHbv[i].GetLakePercent()+DistHbv[i].GetGlacierPercent()+DistHbv[i].GetBogPercent();

      if (value > 100.0) {
        DistHbv[i].SetForestPercent(DistHbv[i].GetForestPercent()*100.0/value);
        DistHbv[i].SetAlpinePercent(DistHbv[i].GetAlpinePercent()*100.0/value);
        DistHbv[i].SetHeatherPercent(DistHbv[i].GetHeatherPercent()*100.0/value);
        DistHbv[i].SetBedrockPercent(DistHbv[i].GetBedrockPercent()*100.0/value);
        DistHbv[i].SetLakePercent(DistHbv[i].GetLakePercent()*100.0/value);
        DistHbv[i].SetGlacierPercent(DistHbv[i].GetGlacierPercent()*100.0/value);
        DistHbv[i].SetBogPercent(DistHbv[i].GetBogPercent()*100.0/value);

        value=DistHbv[i].GetForestPercent()+DistHbv[i].GetAlpinePercent()+
        DistHbv[i].GetHeatherPercent()+DistHbv[i].GetBedrockPercent()+
        DistHbv[i].GetLakePercent()+DistHbv[i].GetGlacierPercent()+DistHbv[i].GetBogPercent();
      }

      if (value > 100.0) {
        GetRowCol(i, nCols, j, k);
        cout << "\n Sum of land use percentages = " << value << "    Element : " << j << "," << k << "\n";
        DistHbv[i].SetOpenLandPercent(0.0);
      } else {
        DistHbv[i].SetOpenLandPercent(100.0-value);
      }
    }
  }
}


// Land surface classes general
void ReadLandUseGeneral(DistributedHbv * const DistHbv, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[1024];
  int i,j,k;
  int nRo, nCo, noDa;
  double value, notAssigned, previousValue;
  double xllC, yllC, cellS;

  // Area of grid cells 
  for (i=0; i<nRows*nCols; i++) {
    DistHbv[i].SetArea(cellSize*cellSize);
  }

  // Read elevation of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with grid cell elevations: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finElevation(fileName);  // Open for reading
  if (!finElevation.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finElevation >> buffer >> nCo;
  finElevation >> buffer >> nRo;
  finElevation >> buffer >> xllC;
  finElevation >> buffer >> yllC;
  finElevation >> buffer >> cellS;
  finElevation >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for elevation!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finElevation >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    DistHbv[i].SetElevation(value);
    //    previousValue=value;
  }
  finElevation.close();
  // Read end

  // Read slope angle of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with slope angles: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finSlopeAngle(fileName);  // Open for reading
  if (!finSlopeAngle.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finSlopeAngle >> buffer >> nCo;
  finSlopeAngle >> buffer >> nRo;
  finSlopeAngle >> buffer >> xllC;
  finSlopeAngle >> buffer >> yllC;
  finSlopeAngle >> buffer >> cellS;
  finSlopeAngle >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for slope!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finSlopeAngle >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    DistHbv[i].SetSlopeAngle(value);
    //    previousValue=value;
  }
  finSlopeAngle.close();
  // Read end

  // Read slope aspect of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with slope aspects: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finAspect(fileName);  // Open for reading
  if (!finAspect.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finAspect >> buffer >> nCo;
  finAspect >> buffer >> nRo;
  finAspect >> buffer >> xllC;
  finAspect >> buffer >> yllC;
  finAspect >> buffer >> cellS;
  finAspect >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for aspect!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finAspect >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    DistHbv[i].SetAspect(value);
    //    previousValue=value;
  }
  finAspect.close();
  // Read end
  
  // Read latitude of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with latitude of grid cells: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finLatitude(fileName);  // Open for reading
  if (!finLatitude.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finLatitude >> buffer >> nCo;
  finLatitude >> buffer >> nRo;
  finLatitude >> buffer >> xllC;
  finLatitude >> buffer >> yllC;
  finLatitude >> buffer >> cellS;
  finLatitude >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for latitude!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finLatitude >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    DistHbv[i].SetLatitude(value);
    //    previousValue=value;
  }
  finLatitude.close();
  // Read end
  
  // Read time falling values of grid cells
  //  previousValue=-9999.0;
  /*  cout << " File with time falling values: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finT1(fileName);  // Open for reading
  if (!finT1.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finT1 >> buffer >> nCo;
  finT1 >> buffer >> nRo;
  finT1 >> buffer >> xllC;
  finT1 >> buffer >> yllC;
  finT1 >> buffer >> cellS;
  finT1 >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for T1!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finT1 >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    DistHbv[i].SetT1(value);
    //    previousValue=value;
  }
  finT1.close();
  // Read end

  // Read percentage of grid cells covered by lakes
  //  previousValue=-9999.0;
  /*cout << " File with lake percentage: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finLakePercent(fileName);  // Open for reading
  if (!finLakePercent.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finLakePercent >> buffer >> nCo;
  finLakePercent >> buffer >> nRo;
  finLakePercent >> buffer >> xllC;
  finLakePercent >> buffer >> yllC;
  finLakePercent >> buffer >> cellS;
  finLakePercent >> buffer >> noDa;
  //  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize) {
    cout << "\n\n Error in grid header information for lake percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finLakePercent >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    DistHbv[i].SetLakePercent(value);
    //    previousValue=value;
  }
  finLakePercent.close();
  // Read end

  // Read percentage of grid cells covered by all land surface classes, number = numberLandSurfaceClasses-1
  for (j=0; j<numberLandSurfaceClasses-1; j++) {
    //    previousValue=-9999.0;
    cout << " File with land surface class percentage: " << j << endl;
    /*    cin >> fileName;
      cout << endl; */
    fileControl.ignore(100,':');
    fileControl >> fileName;
    fileControl.ignore(1024,'\n');
    ifstream finLandSurfacePercent(fileName);  // Open for reading
    if (!finLandSurfacePercent.is_open()) {
      cout << endl << " Error opening file " << fileName << endl << endl;
      exit (1);
    }
    cout << fileName << endl;
    finLandSurfacePercent >> buffer >> nCo;
    finLandSurfacePercent >> buffer >> nRo;
    finLandSurfacePercent >> buffer >> xllC;
    finLandSurfacePercent >> buffer >> yllC;
    finLandSurfacePercent >> buffer >> cellS;
    finLandSurfacePercent >> buffer >> noDa;
    //    if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize) {
      cout << "\n\n Error in grid header information for land surface class percent! " << j << "\n\n";
      cout << " nCols        " << nCo << "\t  " << nCols << endl;
      cout << " nRows        " << nRo << "\t  " << nRows << endl;
      cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
      cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
      cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
      cout << " noData       " << noDa << "\t  " << noData << endl << endl;
      exit (1);
    }
    for (i=0; i<nRows*nCols; i++) {
      finLandSurfacePercent >> value;
      /*      if (j==0) {
	if (i==48921 || i==49341 || i==49760 || i==49761 || i==50180 || i==50181) cout << "land " << i << "  " << value << "\n";
	}*/
      if (value<0.0) value=0.0;
      //    if (value<0.0 && previousValue>=0.0) value=previousValue;
      //      else if (value<0.0) value=0.0;
      DistHbv[i].SetLandSurfacePercent(j,value);
      //      previousValue=value;
    }
    finLandSurfacePercent.close();
  }
  // Read end

  // Read percentage of grid cells covered by glaciers
  //  previousValue=-9999.0;
  /*cout << " File with glacier percentage: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finGlacierPercent(fileName);  // Open for reading
  if (!finGlacierPercent.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finGlacierPercent >> buffer >> nCo;
  finGlacierPercent >> buffer >> nRo;
  finGlacierPercent >> buffer >> xllC;
  finGlacierPercent >> buffer >> yllC;
  finGlacierPercent >> buffer >> cellS;
  finGlacierPercent >> buffer >> noDa;
  //  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize) {
    cout << "\n\n Error in grid header information for glacier percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finGlacierPercent >> value;
    if (value<0.0) value=0.0;
    //    if (value<0.0 && previousValue>=0.0) value=previousValue;
    //    else if (value<0.0) value=0.0;
    DistHbv[i].SetGlacierPercent(value);
    //    previousValue=value;
  }
  finGlacierPercent.close();
  // Read end
  
  // Read soil type for landscape elements
  fileControl.ignore(100, ':');
  fileControl >> fileName;
  fileControl.ignore(1024, '\n');
  ifstream finSoil(fileName);  // Open for reading
  if (!finSoil.is_open()) {
	  cout << endl << " Error opening file " << fileName << endl << endl;
	  exit(1);
  }
  finSoil >> buffer >> nCo;
  finSoil >> buffer >> nRo;
  finSoil >> buffer >> xllC;
  finSoil >> buffer >> yllC;
  finSoil >> buffer >> cellS;
  finSoil >> buffer >> noDa;
  if (nCo != nCols || nRo != nRows || xllC != xllCorner || yllC != yllCorner || cellS != cellSize || noDa != noData) {
	  cout << "\n\n Error in grid header information for soil percent!\n\n";
	  cout << " nCols        " << nCo << "\t  " << nCols << endl;
	  cout << " nRows        " << nRo << "\t  " << nRows << endl;
	  cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
	  cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
	  cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
	  cout << " noData       " << noDa << "\t  " << noData << endl << endl;
	  exit(1);
  }
  for (i = 0; i<nRows*nCols; i++) {
	  finSoil >> value;
	  if (value<0.0) value = 0.0;
	  //    if (value<0.0 && previousValue>=0.0) value=previousValue;
	  //    else if (value<0.0) value=0.0;
	  DistHbv[i].SetSoil(value);
	  //    previousValue=value;
  }
  finSoil.close();
  // Read end


  // Control and correct area fractions of land use classes
  for (i=0; i<nRows*nCols; i++) {
    if (DistHbv[i].GetLandIndex()!=noData) { 
      //      if (DistHbv[i].GetLakePercent() > 0) cout << "\n  cell no. " << i << " glac " << DistHbv[i].GetLakePercent() << "\n";
      //      if (DistHbv[i].GetGlacierPercent() > 0) cout << "\n cell no. " << i << " glac " << DistHbv[i].GetGlacierPercent() << "\n";
      value=DistHbv[i].GetLakePercent()+DistHbv[i].GetGlacierPercent();
      for (j=0; j<numberLandSurfaceClasses-1; j++) {
	//	if (value > 0 && DistHbv[i].GetLandSurfacePercent(j) > 0) cout << " cell no. " << i << " valu " << value << " land " << DistHbv[i].GetLandSurfacePercent(j) << "\n";
	value=value+DistHbv[i].GetLandSurfacePercent(j);
      }
      //      if (value > 100) cout << " cell no. " << i << "      " << value << "\n";
      if (value < 1.0) {
        DistHbv[i].SetLandSurfacePercent(0,100.0);
        value = 100.0;
      }
      /*      else if (value < 100.0) {
              DistHbv[i].SetLandSurfacePercent(2,100.0-value);
              value=value+DistHbv[i].GetLandSurfacePercent(2);
              }
	      if (value != 100.0) {*/
      else if (value != 100.0) {
        DistHbv[i].SetLakePercent(DistHbv[i].GetLakePercent()*100.0/value);
        DistHbv[i].SetGlacierPercent(DistHbv[i].GetGlacierPercent()*100.0/value);
        for (j=0; j<numberLandSurfaceClasses-1; j++) DistHbv[i].SetLandSurfacePercent(j,DistHbv[i].GetLandSurfacePercent(j)*100.0/value);

        value=DistHbv[i].GetLakePercent()+DistHbv[i].GetGlacierPercent();
        for (j=0; j<numberLandSurfaceClasses-1; j++) value=value+DistHbv[i].GetLandSurfacePercent(j);
      }

      if (value != 100.0) {
        GetRowCol(i, nCols, j, k);
        cout << "\n Sum of land use percentages = " << value << "    Element : " << j << "," << k << "\n";
      }
    }
  }
}


void WriteLandScapeElements(DistributedHbv * const DistHbv, ParametersGeneral * const ParGeneralStore,
                            MeteorologicalStations * const MetStations, int landIndex, int nRows, int nCols, 
                            int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
  int i, k, n, maxIndex;
  double lakePer, glacPer, areaCorr[maximumNumberLandClasses], totalArea, tempArea, areaFraction[numberLandSurfaceClasses-1];
  LANDSURFACE tempLandSurf, landSurfType[numberLandSurfaceClasses-1];
  SOIL soilType[numberSoilClasses-1];
  int * precipitationStation = new int [ParGeneralStore->GetNUM_PREC_SERIES()];
  int * temperatureStation = new int [ParGeneralStore->GetNUM_TEMP_SERIES()];
  double * precipitationWeight = new double [ParGeneralStore->GetNUM_PREC_SERIES()];
  double * temperatureWeight = new double [ParGeneralStore->GetNUM_TEMP_SERIES()];

  // Landscape elements  
  ofstream landScapeOut("hbv_landscape.txt");  // Open for writing
  ofstream geoIndexOut("hbv_grid_index.txt");  // Open for writing
  landScapeOut << "ncols         " << nCols << endl;
  landScapeOut << "nrows         " << nRows << endl;
  landScapeOut << "xllcorner     " << xllCorner << endl;
  landScapeOut << "yllcorner     " << yllCorner << endl;
  landScapeOut << "cellsize      " << cellSize << endl;
  landScapeOut << "NODATA_value  " << noData << endl;
  landScapeOut << "# Number of landscape elements :  " << landIndex << endl;
  geoIndexOut << " Index of grid cells relative to upper left corner of rectangle with no. of rows = " << nRows << " and no. of columns = " << nCols << endl;
  for (i=0; i<nRows*nCols; i++) {
    if (DistHbv[i].GetLandIndex()!=noData) {
      // Sort land surface types based on area
      for (k=0; k<numberLandSurfaceClasses-1; k++) landSurfType[k]=LANDSURFACE(k);
      for (k=0; k<numberLandSurfaceClasses-1; k++) areaFraction[k] = 0.0;
      // Area fraction for land surface classes based on potential tree level
      /*      areaFraction[0]=DistHbv[i].GetOpenLandPercent();      // Open land (meadows, agriculture)
      areaFraction[1]=DistHbv[i].GetBogPercent();           // Bogs
      areaFraction[2]=DistHbv[i].GetForestPercent();        // Forest
      areaFraction[3]=DistHbv[i].GetAlpinePercent();        // Alpine forest
      areaFraction[4]=DistHbv[i].GetHeatherPercent();       // Low mountain with vegetation
      areaFraction[5]=DistHbv[i].GetBedrockPercent();       // Exposed bedrock (high mountain)*/
      // Area fraction for land surface classes general
      for (k=0; k<numberLandSurfaceClasses-1; k++) areaFraction[k] = DistHbv[i].GetLandSurfacePercent(k);
      // End area fraction for land surface classes
      for (k=0; k<numberLandSurfaceClasses-1; k++) {
        maxIndex=k;
        for (n=k+1; n<numberLandSurfaceClasses-1; n++) {
          if (areaFraction[n] > areaFraction[maxIndex]) {
            tempLandSurf=landSurfType[maxIndex];
            landSurfType[maxIndex]=landSurfType[n];
            landSurfType[n]=tempLandSurf;
            tempArea=areaFraction[maxIndex];
            areaFraction[maxIndex]=areaFraction[n];
            areaFraction[n]=tempArea; 
          }
        }
      }
      // Read explicit soil types for Penman-Monteith version of code
      //      for (k=0; k<numberLandSurfaceClasses-1; k++) soilType[k]=SOIL(landSurfType[k]);
      lakePer=DistHbv[i].GetLakePercent();
      glacPer=DistHbv[i].GetGlacierPercent();
      areaCorr[0]=areaFraction[0];
      areaCorr[1]=areaFraction[1];
      areaCorr[2] = areaFraction[2];
      //      cout << lakePer << "    " << glacPer << "    "  << areaCorr[0] << "    "
      //           << areaCorr[1] << endl;
      totalArea = lakePer + glacPer + areaCorr[0] + areaCorr[1] + areaCorr[2];
      if (totalArea!=100.0) {
        //      if (areaFraction[0]==0.0 && areaFraction[1]==0.0) {
        //        lakePer = lakePer*100/totalArea;
        //        glacPer = glacPer*100/totalArea;
        //        } else {
        areaCorr[0]=areaFraction[0]+areaFraction[0]*(100.0-totalArea)/(areaFraction[0]+areaFraction[1]+areaFraction[2]);
        areaCorr[1]=areaFraction[1]+areaFraction[1]*(100.0-totalArea)/(areaFraction[0]+areaFraction[1]+areaFraction[2]);
	areaCorr[2] = areaFraction[2] + areaFraction[2] * (100.0 - totalArea) / (areaFraction[0] + areaFraction[1] + areaFraction[2]);
        //       }
        totalArea = lakePer + glacPer + areaCorr[0] + areaCorr[1] + areaCorr[2];
        if (totalArea!=100.0) cout << endl << "totalArea    " << totalArea << endl;
      }

      /* Find meteorological stations */
      DistanceSort(&DistHbv[i], ParGeneralStore, MetStations, 
                   precipitationStation, temperatureStation, precipitationWeight, temperatureWeight);

      landScapeOut.precision(4); landScapeOut.setf(ios::fixed); landScapeOut.setf(ios::showpoint);
      landScapeOut.width(10); landScapeOut << DistHbv[i].GetLandIndex() << "  ";
      landScapeOut.width(10); landScapeOut << DistHbv[i].GetGeoIndex() << "  ";
      landScapeOut.width(15); landScapeOut << DistHbv[i].GetArea() << "  ";
      landScapeOut.precision(6); 
      landScapeOut.width(15); landScapeOut << DistHbv[i].GetLatitude() << "  ";
      landScapeOut.precision(4); 
      landScapeOut.width(10); landScapeOut << DistHbv[i].GetElevation() << "  ";
      landScapeOut.width(10); landScapeOut << DistHbv[i].GetT1() << "  ";
      landScapeOut.width(10); landScapeOut << DistHbv[i].GetSlopeAngle() << "  ";
      landScapeOut.width(10); landScapeOut << DistHbv[i].GetAspect() << "  ";
      landScapeOut.width(10); landScapeOut << lakePer << "  ";
      landScapeOut.width(10); landScapeOut << glacPer << "  ";
      for (k=0; k<maximumNumberLandClasses; k++) {
        landScapeOut.width(10); landScapeOut << landSurfType[k] << "  ";
        landScapeOut.width(10); landScapeOut << DistHbv[i].GetSoil()  << "  ";
        landScapeOut.width(10); landScapeOut << areaCorr[k] << "  ";
      }
      for (k=0; k<ParGeneralStore->GetNUM_PREC_SERIES(); k++) {
        landScapeOut.width(10); landScapeOut << precipitationStation[k] << "  ";
        landScapeOut.width(10); landScapeOut << precipitationWeight[k] << "  ";
      }
      for (k=0; k<ParGeneralStore->GetNUM_TEMP_SERIES(); k++) {
        landScapeOut.width(10); landScapeOut << temperatureStation[k] << "  ";
        landScapeOut.width(10); landScapeOut << temperatureWeight[k] << "  ";
      }
      landScapeOut<< endl;
      geoIndexOut.width(10); geoIndexOut << "  " << DistHbv[i].GetGeoIndex() << "  " << endl;
    }
  }
  landScapeOut << endl;
  landScapeOut.close();
  geoIndexOut.close();

  delete [] precipitationStation;
  delete [] temperatureStation;
  delete [] precipitationWeight;
  delete [] temperatureWeight;
}


void DistanceSort(DistributedHbv * const distElement, ParametersGeneral * const ParGeneralStore, 
                  MeteorologicalStations * const MetStations,
                  int * precipitationStation, int * temperatureStation, double * precipitationWeight, double * temperatureWeight)
{
  int i, j;
  int numP, minIndex;
  int temporaryStation;
  double temporaryDistance;
  double sumWeight;
  int * precipitationSort = new int [MetStations->GetNumPrecStations()];
  int * temperatureSort = new int [MetStations->GetNumTempStations()];
  double * distancePrec = new double [MetStations->GetNumPrecStations()];
  double * distanceTemp = new double [MetStations->GetNumTempStations()];

  for (i=0; i<MetStations->GetNumPrecStations(); i++) {
    precipitationSort[i] = i;
    distancePrec[i] = sqrt(pow((distElement->GetXCoord() - MetStations->GetStationCoordX(i)),2.0) 
                           + pow((distElement->GetYCoord() - MetStations->GetStationCoordY(i)),2.0));
    //    cout << "Prec " << precipitationSort[i] << "  " << distancePrec[i] << endl;
  }
  for (i=0; i<MetStations->GetNumPrecStations(); i++) {
    minIndex = i;
    for (j=i+1; j<MetStations->GetNumPrecStations(); j++) {
      if (distancePrec[j] < distancePrec[minIndex]) {
        temporaryDistance = distancePrec[minIndex];
        distancePrec[minIndex] = distancePrec[j];
        distancePrec[j] = temporaryDistance;
        temporaryStation = precipitationSort[minIndex];
        precipitationSort[minIndex] = precipitationSort[j];
        precipitationSort[j] = temporaryStation;
      }
    }
  }
  //  cout << "PREC  " << precipitationSort[0] << "  " << distancePrec[0] << endl;

  numP = MetStations->GetNumPrecStations();
  for (i=0; i<MetStations->GetNumTempStations(); i++) {
    temperatureSort[i] = i;
    distanceTemp[i] = sqrt(pow((distElement->GetXCoord() - MetStations->GetStationCoordX(numP+i)),2.0) 
                    + pow((distElement->GetYCoord() - MetStations->GetStationCoordY(numP+i)),2.0));
    //    cout << "Temp " << temperatureSort[i] << "  " << distanceTemp[i] << endl;
  }
  for (i=0; i<MetStations->GetNumTempStations(); i++) {
    minIndex = i;
    for (j=i+1; j<MetStations->GetNumTempStations(); j++) {
      if (distanceTemp[j] < distanceTemp[minIndex]) {
        temporaryDistance = distanceTemp[minIndex];
        distanceTemp[minIndex] = distanceTemp[j];
        distanceTemp[j] = temporaryDistance;
        temporaryStation = temperatureSort[minIndex];
        temperatureSort[minIndex] = temperatureSort[j];
        temperatureSort[j] = temporaryStation;
      }
    }
  }
  //  cout << "TEMP  "  << temperatureSort[0] << "  " << distanceTemp[0] << endl;

  sumWeight = 0.0;
  for (i=0; i<ParGeneralStore->GetNUM_PREC_SERIES(); i++) {
    precipitationStation[i] = precipitationSort[i];
    precipitationWeight[i] = 1.0/distancePrec[i];
    sumWeight = sumWeight + precipitationWeight[i];
  }
  if (sumWeight != 1.0) {
    for (i=0; i<ParGeneralStore->GetNUM_PREC_SERIES(); i++) {
      precipitationWeight[i] = precipitationWeight[i]/sumWeight;
    }
  }

  sumWeight = 0.0;
  for (i=0; i<ParGeneralStore->GetNUM_TEMP_SERIES(); i++) {
    temperatureStation[i] = temperatureSort[i];
    temperatureWeight[i] = 1.0/distanceTemp[i];
    sumWeight = sumWeight + temperatureWeight[i];
  }
  if (sumWeight != 1.0) {
    for (i=0; i<ParGeneralStore->GetNUM_TEMP_SERIES(); i++) {
      temperatureWeight[i] = temperatureWeight[i]/sumWeight;
    }
  }

  //  cout << distElement->GetLandIndex() << "  " << distElement->GetGeoIndex() << "  "
  //       << distElement->GetXCoord() << "  " << distElement->GetYCoord() << endl;

  delete [] precipitationSort;
  delete [] temperatureSort;
  delete [] distancePrec;
  delete [] distanceTemp;

}


void WriteCatchmentIdentifier(DistributedHbv * const DistHbv, SubCatchment * const CatchmentElement, int numWatc, ofstream &fout)
{
  int i;
  DistributedHbv * thisElement;

  // Sub-catchment outlets
  ofstream subCatchmentOut("hbv_waterland.txt");  // Open for writing
  for (i=0; i<numWatc; i++) {
    if (CatchmentElement[i].GetLandScapeElement()) {
      subCatchmentOut << "#  " << CatchmentElement[i].GetIdentifier() << "  #  " << CatchmentElement[i].GetNumLandScape() << endl;
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


void FindCatchmentIdentifier(DistributedHbv * const DistHbv, SubCatchment * const CatchmentElement, int numWatc, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[256];
  int i,j;
  int nRo, nCo, noDa;
  int subCatchmentId;
  double xllC, yllC, cellS;
  bool subCatchmentFound;
  bool noElement=false;
  DistributedHbv * lastElement;
  
  // Read sub-catchment identifiers for landscape elements
  /*  cout << "\n File with sub-catchment identifiers: ";
      cin >> file;
      cout << endl; */
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finId(fileName);  // Open for reading
  if (!finId.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finId >> buffer >> nCo;
  finId >> buffer >> nRo;
  finId >> buffer >> xllC;
  finId >> buffer >> yllC;
  finId >> buffer >> cellS;
  finId >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for sub-catchment identifier one!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  // Read end

  // Connect landscape elements to sub-catchment outlets
  for (i=0; i<nRows*nCols; i++) {
    finId >> subCatchmentId;
    DistHbv[i].SetSubCatchmentValue(subCatchmentId);
    j=0;
    subCatchmentFound = false;
    while (j<numWatc) {
      if (subCatchmentId==CatchmentElement[j].GetIdentifier()) {
        if (!subCatchmentFound) {
          subCatchmentFound=true;
          if (CatchmentElement[j].GetLandScapeElement()) {
            lastElement = CatchmentElement[j].GetLandScapeElement();
            while (lastElement->GetNextElement()) lastElement = lastElement->GetNextElement();
            lastElement->SetNextElement(&DistHbv[i]);
          } else {
            CatchmentElement[j].SetLandScapeElement(&DistHbv[i]);
          }
          CatchmentElement[j].SetNumLandScape(CatchmentElement[j].GetNumLandScape()+1);
        } else {
          cout << "\n\n Error in sub-catchment identifiers!   ";
          cout << CatchmentElement[j].GetIdentifier() << endl << endl;
          exit (1);
        }
      } 
      j++;
    }
  }

  for (j=0; j<numWatc; j++) {
    if (!(CatchmentElement[j].GetLandScapeElement())) {
      cout << "\n No landscape elements for sub-catchment index " << j << "  "  << "Sub-catchment identifier " << CatchmentElement[j].GetIdentifier() << endl;
      noElement=true;
    }
  }
  if (noElement) {
    cout << "\n Program is terminated! " << endl;
    exit(1);
  }
    
  finId.close();
}


void GetRowCol(int elementNo, int nCols, int &row, int &col) 
{
  row = elementNo/nCols; 
  col = elementNo%nCols; 
}


void SetGeneralParameters(ParametersGeneral * const ParGeneralStore, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  int numPrec, numTemp;

  /*  cout << " File with common parameters: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finGeneralPar(fileName);  // Open for reading
  if (!finGeneralPar.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finGeneralPar.ignore(256,'\n');
  finGeneralPar.ignore(100,':'); finGeneralPar >> numPrec; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> numTemp; 
  ParGeneralStore->SetNUM_PREC_SERIES(numPrec);
  ParGeneralStore->SetNUM_TEMP_SERIES(numTemp);
  finGeneralPar.close();
  fout << "Common parameters: \n";
  fout << ParGeneralStore->GetNUM_PREC_SERIES() << endl;
  fout << ParGeneralStore->GetNUM_TEMP_SERIES() << endl;
  fout << endl << endl;
}
