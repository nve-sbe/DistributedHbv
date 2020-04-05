#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int maximumStations=1000;

using namespace std;

int main()
{
  char fileName[80];
  char buffer[256];
  char ch;
  int i,j;
  int nRo, nCo, noDa;
  int numCatch;
  int value;
  int domain;
  double xllC, yllC, cellS;
  bool catchmentFound;
  bool noElement=false;
  
  int * stationId = new int[maximumStations];

  cout << "\n File with sub-catchment hierarchy: ";
  cin >> fileName;
  cout << endl;
  ifstream fileWCo(fileName);
  if (fileWCo == NULL) {
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
    fileWCo.ignore(256,'\n');
  }

  cout << "\n File with model domain: ";
  cin >> fileName;
  cout << endl; 
  ifstream finDom(fileName); 
  if (finDom == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finDom >> buffer >> nCo;
  finDom >> buffer >> nRo;
  finDom >> buffer >> xllC;
  finDom >> buffer >> yllC;
  finDom >> buffer >> cellS;
  finDom >> buffer >> noDa;

  cout << "\n File with sub-catchment identifiers: ";
  cin >> fileName;
  cout << endl; 
  ifstream finId(fileName); 
  if (finId == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finId >> buffer >> nCo;
  finId >> buffer >> nRo;
  finId >> buffer >> xllC;
  finId >> buffer >> yllC;
  finId >> buffer >> cellS;
  finId >> buffer >> noDa;

  cout << " Output file: ";
  cin >> fileName;
  cout << endl;
  ofstream fout(fileName);  
  fout << "ncols         " << nCo << endl;
  fout << "nrows         " << nRo << endl;
  fout << "xllcorner     " << xllC << endl;
  fout << "yllcorner     " << yllC << endl;
  fout << "cellsize      " << cellS << endl;
  fout << "NODATA_value  " << noDa << endl;

  for (i=0; i<nRo*nCo; i++) {
    finDom >> domain;
    finId >> value;
    if (domain > noDa && value > noDa) { 
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
      else {
        fout << noDa << endl;
      }
    }
    else {
      fout << noDa << endl;
    }
  }

  fileWCo.close();
  finDom.close();
  finId.close();
  fout.close();
  delete [] stationId;

  return 0;
}

