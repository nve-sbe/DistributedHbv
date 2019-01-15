#include "classDistHbv.h"
#include "utilities.h"

void Vegetation::WaterBalance(int timeStep, DateTime datetime)
{
  double timeResolution=1.0;
  int i;
  precipitation=GetInputElement()->GetInput(0);
  temp = GetInputElement()->GetInput(1);
  //  cout << "    timeStep " << timeStep << "    Vegetation precipitation " << precipitation;
  //  cout << "    Temperature " << temp << endl;

  /* Potential evapotranspiration */
  wetPeriod = 0.0;
  dryPeriod = timeResolution;
  throughFall = 0.0;
  interceptionLoss = 0.0;
  if (GetLandScapeElement()->GetEvaporationControlObj()->GetEvaporationModellingControl() == 'M' || 
      GetLandScapeElement()->GetEvaporationControlObj()->GetEvaporationModellingControl() == 'm') {
    i = datetime.getMonth()-1;
    potev = potentialEvapLongTermMean(temp, GetLandScapeElement()->GetEvaporationControlObj()->GetEvaporationArray(i)/1000.0);
    /*for (i=0; i<numberPotentialEvaporationValuesPerYear; i++) {
      cout << " " << i << "    " << GetLandScapeElement()->GetEvaporationControlObj()->GetEvaporationArray(i) << "\n";
      }*/
  }
  else {
    potev = potentialEvapTemperatureIndex(temp, landSurfacePar->GetEPOT_PAR());
  }

  /* Water input from precipitation and/or snowmelt > 0 or interception remaining from previous time step */
  if (precipitation > 0.0 || prevInterception > 0.0) {
    interceptionStore = prevInterception + precipitation;
    if (potev > 0.0) {
      wetPeriod = interceptionStore / potev;
    }
    else {
      wetPeriod = 0.0;
    }
    //    cout << interceptionStore << " " << potev  << " " << timeResolution  << " " << wetPeriod << " " << dryPeriod << endl;
    if (wetPeriod < timeResolution) {
      interceptionLoss = potev * wetPeriod;
      dryPeriod = timeResolution - (wetPeriod*landSurfacePar->GetWET_PER_CORR());
    }
    else {
      interceptionLoss = potev * timeResolution;
      dryPeriod = 0.0;
    }
    interceptionStore = interceptionStore - interceptionLoss;

    /* If interception store > INTER_MAX, surplus water is infiltrated through the soil surface
       Soil moistured deficit and percolation is calculated */
    //    cout << interceptionStore << " " << potev  << " " << timeResolution  << " " << wetPeriod << " " << dryPeriod << endl;
    if (interceptionStore > landSurfacePar->GetINTER_MAX()) {
      throughFall = interceptionStore - landSurfacePar->GetINTER_MAX();
      interceptionStore = landSurfacePar->GetINTER_MAX();
    }

    /* If interception store < 0 following evaporation at potential rate then dry period > 0, 
       transpiration and soil evaporation will be calculated for dryPeriod */
    if (interceptionStore < 0.0-epsilon*epsilon) {
      printf("    timeStep = %d      interceptionStore = %f\n",timeStep,interceptionStore);
      /*        exit(1);*/
    }
    if (wetPeriod < 0.0) {
      printf("    timeStep = %d      wetPeriod = %f\n",timeStep,wetPeriod);
      /*        exit(1);*/
    }
    if (dryPeriod < 0.0 || dryPeriod > timeResolution) {
      printf("    timeStep = %d      dryPeriod = %f\n",timeStep,dryPeriod);
      /*        exit(1);*/
    }

    /* If 0 <= interception store <= INTER_MAX, no action is necessary
       Infiltration through soil surface = 0, soil moisture deficit and depth of saturated zone is unchanged */
  }
  prevInterception = interceptionStore;
  //  cout << "throughFall " << throughFall << endl;
}

