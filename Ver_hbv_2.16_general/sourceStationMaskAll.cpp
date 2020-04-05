/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments,          *
 *  river networks, lakes and landscape elements.                                       *
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

const int maximumStations=1000;

using namespace std;

int main(int argc, char *argv[])
{
  std::cout << "\n\n Preprocessing for Distributed Element Water Model \n\n";

  char fileName[80];
  char buffer[1024];
  char ch;
  int i,j;
  int nRo, nCo, noDa;
  int numCatch;
  int value, valueArea;
  double xllC, yllC, cellS;
  bool catchmentFound;
  bool noElement=false;
  
  int * stationId = new int[maximumStations];

  if (argc != 2) {
    cout << " " << argv[0] << "  <control file name>\n\n";
    exit(1);
  }
  ifstream fileControl(argv[1]);
  if (!fileControl.is_open()) {
    cout << " Error opening file " << argv[1] << endl << endl;
    exit (1);
  }

  /*  cout << "\n File with sub-catchment identifiers: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream fileWCo(fileName);
  if (!fileWCo.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }  fileWCo.ignore(100,':');
  fileWCo >> numCatch;
  cout << " # Number of sub-catchments: " << numCatch << endl;
  for (i=0; i<numCatch; i++) {
    fileWCo >> j >> ch >> stationId[i];
    if (j != i) {
      cout << endl << "Error reading file " << fileName << "\t" << i << "\t" << j << endl;
      exit (1);
    }
    fileWCo.ignore(1024,'\n');
  }

  /*  cout << "\n File with grid cell sub-catchment identifiers: ";
      cin >> fileName;
      cout << endl; */
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finId(fileName); 
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

  /*  cout << "\n File with model domain area: ";
      cin >> fileName;
      cout << endl; */
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finArea(fileName); 
  if (!finArea.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finArea >> buffer >> nCo;
  finArea >> buffer >> nRo;
  finArea >> buffer >> xllC;
  finArea >> buffer >> yllC;
  finArea >> buffer >> cellS;
  finArea >> buffer >> noDa;

  /*  cout << " Output file: ";
  ci|n >> fileName;
  cou|t << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ofstream fout(fileName);  
  fout.precision(0); fout.setf(ios::fixed); 
  fout << "ncols         " << nCo << endl;
  fout << "nrows         " << nRo << endl;
  fout << "xllcorner     " << xllC << endl;
  fout << "yllcorner     " << yllC << endl;
  fout << "cellsize      " << cellS << endl;
  fout << "NODATA_value  " << noDa << endl;

  for (i=0; i<nRo*nCo; i++) {
    finId >> value;
    finArea >> valueArea;
    if (value > noDa) { 
      j=0;
      catchmentFound = false;
      while (j<numCatch && !catchmentFound) {
        //      cout << value << "  " << j << "  " << stationId[j] << endl;
        if (stationId[j]==value) catchmentFound=true;
        j++;
      } 
      if (catchmentFound) {
        fout << value << endl;
      } 
      else if (valueArea > noDa) {
	fout << 0 << endl;
      }
      else {
        fout << noDa << endl;
      }
    }
    else if (valueArea > noDa) {
      fout << 0 << endl;
    }
    else {
      fout << noDa << endl;
    }
  }

  fileWCo.close();
  finId.close();
  fout.close();
  fileControl.close();
  delete [] stationId;

  return 0;
}

