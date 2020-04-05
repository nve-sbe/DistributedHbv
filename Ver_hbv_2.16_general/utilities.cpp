#include "classDistHbv.h"
#include "utilities.h"


// Power function 
double power(double base, double exponent)
{
  if (base > 0.0) 
    return pow(base,exponent);
  else if (base == 0.0) 
    return 0.0;
  else if (base >= -epsilon)
    return 0.0;
  else {
    printf("\n\n    *****     function power     base = %f     exponent = %f     *****\n\n",base,exponent);
    return 0.0;
  }
}


// Leap year
int leapYear(int year)
{
  if (year%400==0) 
    return 1;
  else if (year%100==0) 
    return 0;
  else if (year%4==0) 
    return 1;
  else
    return 0;
}


// Find day number
int dayNumber(int year, int month, int day)
{
  int accDays[] = {0,31,59,90,120,151,181,212,243,273,304,334};
  if (month<3) 
    return accDays[month-1]+day;
  else
    return accDays[month-1]+day+leapYear(year);
}


// Day number to date
void dayNo2Date(int dayNo, int year, int * month, int * day)
{      
  int i;
  int accDays[] = {0,31,59,90,120,151,181,212,243,273,304,334};
  if (dayNo>365+leapYear(year)) dayNo=dayNo-365-leapYear(year);
  if (leapYear(year)) 
    for (i=2; i<12; i++) accDays[i]++;
  i=12;
  while (dayNo<=accDays[i-1]) i--;
  *month=i;
  *day=dayNo-accDays[i-1];
}


// Potential evaporation based on temperature index
double potentialEvapTemperatureIndex(double temp, double epotPar)
{
  if (temp > 0.0)
    return epotPar * temp;
  else
    return 0.0;
}


// Actual evapotranspiration HBV based on temperature index potential evapotranspiration
double HBVTranspSoilEvapTemperatureIndex(double soilMoist, double temp, double epotPar, double fieldCapacity, double fcDel)
{
  double epot;
  epot = potentialEvapTemperatureIndex(temp,epotPar);            
  if (soilMoist > fcDel * fieldCapacity) {
    return epot;
  }
  else if (soilMoist > 0.0) {
    return epot * soilMoist / (fieldCapacity * fcDel);
  }
  else {
    return 0.0;
  }
}


// Potential evaporation based on long-term mean potential evapotranspiration
double potentialEvapLongTermMean(double temp, double potentialEvaporation)
{
  if (temp > 0.0)
    return potentialEvaporation;
  else
    return 0.0;
}


// Actual evapotranspiration HBV based on long-term mean potential evapotranspiration
double HBVTranspSoilEvapLongTermMean(double soilMoist, double temp, double epotPar, double fieldCapacity, double fcDel, double potentialEvaporation)
{
  double epot;
  //  epot = potentialEvapLongTermMean(temp, potentialEvaporation);
  epot = potentialEvaporation;
  if (soilMoist > fcDel * fieldCapacity) {
    return epot;
  }
  else if (soilMoist > 0.0) {
    return epot * soilMoist / (fieldCapacity * fcDel);
  }
  else {
    return 0.0;
  }
}


// MaxBas triangular weight
double MaxBasWeight(double lower, double upper, double maxBas)
{
  return (MaxBasFunction(lower,maxBas) + MaxBasFunction(upper,maxBas)) * (upper-lower)/2.0;
}


// MaxBas triangular function
double MaxBasFunction(double argument, double maxBas)
{
  return (2.0/maxBas) - absValue(argument-(maxBas/2.0)) * (4.0/(maxBas*maxBas));
}


// Absolute value
double absValue(double value)
{
  if (value >= 0) return value;
  else return -value;
}
