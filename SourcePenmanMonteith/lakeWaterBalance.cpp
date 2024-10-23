#include "classDistHbv.h"
#include "utilities.h"

void LakeWaterBalance::WaterBalance(int timeStep, DateTime datetime, int dayofyear_PM2)
{
  double tempMemory;
  double stage, newStage;
  double kLake;
  double deltaLevel;
  double nLake;
  double waterInput;
  int i;
  precipitation=GetInputElement()->GetInput(0);
  temp = GetInputElement()->GetInput(1);
  tempMax = GetInputElement()->GetInput(2);
  tempMin = GetInputElement()->GetInput(3);
  wind1 = GetInputElement()->GetInput(4);
  radiationS = GetInputElement()->GetInput(5);
  vp = GetInputElement()->GetInput(6);
  /*  cout << "    timeStep " << timeStep << "          Lake precipitation " << precipitation;
      cout << "    Temperature " << temp << endl;*/

  waterInput = precipitation;

  /*  Lake temperature and evaporation  */
  tempMemory = commonPar->GetDAY_TEMP_MEMORY() * 86400.0/commonPar->GetSECONDS_TIMESTEP();
  lakeTemp = lakeTemp*(1.0-(1.0/tempMemory))+temp/tempMemory;
  if (lakeTemp > 0.0) {
    if(GetLandScapeElement()->GetEvaporationControlObj()->GetEvaporationModellingControl() == 'P')
      lakeEvaporation = potentialEvapLake(temp, tempMax, tempMin, radiationS, vp, 
					  wind1, dayofyear_PM2,GetLandScapeElement()->GetLatitude(),
					  GetLandScapeElement()->GetElevation(),
					  commonPar->GetHeight_WIND_INSTRUMENT() ) / 1000.; // convert from mm/day to m/day
     else
       lakeEvaporation =  potentialEvapTemperatureIndex(lakeTemp, commonPar->GetLAKE_EPOT_PAR()); //hvis ikke P, så er det T. ikke dele på 1000 her! (?)
  }
  else
    lakeEvaporation = 0.0;
  if (waterLevel < (-1)*deltaLevel) {
    lakeEvaporation = 0.0;
  }

  /*  Lake rating curve parameters  */
  deltaLevel = commonPar->GetDELTA_LEVEL();
  nLake = commonPar->GetNLAKE();
  kLake = commonPar->GetKLAKE();

  /*  Initial lake water level  */
  stage = waterLevel + waterInput - lakeEvaporation;

  /*  Lake outflow  */
  if (stage < (-1)*deltaLevel) {
    runoff = 0.0;
  }
  else {
    runoff = kLake * pow((stage+deltaLevel),nLake);
  }

  /*  Final runoff  */
  newStage = stage - runoff;
  //if (runoff > 0.0) {                                  // New test added in order to allow lake evaporation to draw water below -deltaLevel
  //if (newStage + deltaLevel < 0.0) {
  //  runoff = runoff + newStage + deltaLevel;
  //  if (runoff < 0.0) runoff = 0.0;
  //  newStage = (-1)*deltaLevel;
  //}
  //}

  /*  Final lake water level and lake water level change  */
  waterLevelChange = newStage - waterLevel;
  waterLevel = newStage;
  discharge = runoff*(GetLandScapeElement()->GetArea()*GetLandScapeElement()->GetLake()->GetAreaFraction()/100.0)/commonPar->GetSECONDS_TIMESTEP();

}

