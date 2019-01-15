#pragma once

#include "classDistHbv.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

void SetLandSurfaceParameters(ParametersLandSurface * const ParLandSurfaceStore, ifstream &fileControl, ofstream &fout);
void SetSnowDistribution(ParametersLandSurface * const thisParLandSurface, double cvSnow);
void SetSubSurfaceHbvParameters(ParametersSubSurfaceHbv * const ParSubSurfaceHbvStore, ifstream &fileControl, ofstream &fout);
void SetGeneralParameters(ParametersGeneral * const ParGeneralStore, ifstream &fileControl, ofstream &fout);
