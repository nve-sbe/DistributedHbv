#pragma once

#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

double power(double base, double exponent);
int leapYear(int year);
int dayNumber(int year, int month, int day);
void dayNo2Date(int dayNo, int year, int * month, int * day);
double potentialEvapTemperatureIndex(double temp, double epotPar);
double HBVTranspSoilEvapTemperatureIndex(double soilMoist, double temp, double epotPar, double fieldCapacity, double fcDel);
double potentialEvapLongTermMean(double temp, double potentialEvaporation);
double HBVTranspSoilEvapLongTermMean(double soilMoist, double temp, double epotPar, double fieldCapacity, double fcDel, double potentialEvaporation);
double MaxBasWeight(double, double, double);
double MaxBasFunction(double, double);
double absValue(double);
