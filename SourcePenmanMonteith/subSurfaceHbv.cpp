#include "classDistHbv.h"
#include "utilities.h"

void HBV::WaterBalance(int timeStep, DateTime datetime, double waterInput, double snowCoverFraction, double dryPeriod, double potev)
{
  int i, hoursPerTimeStep;
  double inSoilHour, inSoil, outSoil;
  double lowerRunoff, upperRunoff, lowerPercolation, upperPercolation, upperTemporary;
  double lowerZoneMax, drawUp;
  temp = GetInputElement()->GetInput(1);

  /* Maximum lower zone storage */
  lowerZoneMax = subSurfacePar->GetPERC()/subSurfacePar->GetKLZ();
  /*if(timeStep==0) {
    cout << "timeStep " << timeStep << " sM " << soilMoisture << " uZ " << upperZone << " lZ " << lowerZone << " " << endl;
    cout << "subSurfacePar->GetFC() " << subSurfacePar->GetFC() << " subSurfacePar->GetKUZ()  " << subSurfacePar->GetKUZ() << endl;
    cout << "subSurfacePar->GetPerc() " << subSurfacePar->GetPERC() << " subSurfacePar->GetKLZ() " << subSurfacePar->GetKLZ() << endl;
     cout << "timeStep " << timeStep << "  lowerzonemax " << lowerZoneMax << endl;
     }*/
  
  /* Capillary rise from lower zone to soil moisture zone */
  drawUp = (2.0*subSurfacePar->GetDRAW())*(lowerZone/lowerZoneMax)*(subSurfacePar->GetFC()-soilMoisture)/(subSurfacePar->GetFC());    
  lowerZone = lowerZone - drawUp;
  if (lowerZone < 0.0) {
    drawUp = drawUp + lowerZone;
    lowerZone = 0.0;
  }
  soilMoisture = soilMoisture + drawUp;
    
  /* Water exceeding lower zone maximum to upper zone */
  if (lowerZone > lowerZoneMax) {
    upperZone = upperZone + lowerZone - lowerZoneMax;
    lowerZone = lowerZoneMax;
  }

  /* Soil moisture deficit in use for glacier-free areas */
  if (dryPeriod > missingData) {
    /* Infiltration to soil moisture zone */
    if (waterInput > subSurfacePar->GetINFMAX()) {
      upperPercolation = waterInput - subSurfacePar->GetINFMAX();
      inSoil = subSurfacePar->GetINFMAX();
    }
    else {
      upperPercolation = 0.0;
      inSoil = waterInput;
    }

//    cout << " subSurfaceHbv   snowCoverFraction   " << snowCoverFraction << endl;
    /* Transpiration and soil evaporation for snow free areas */
    if (dryPeriod > 0.0) {
//    if (dryPeriod > 0.0 && GetLandScapeElement()->GetSnowStore() <= 0.0) {
//    if (dryPeriod > 0.0 && GetLandScapeElement()->GetSnowStore() > 0.0) {
      transpSoilEvap = dryPeriod * (1.0 - snowCoverFraction) * 
                       HBVTranspSoilEvap (soilMoisture,  potev, 
                                          subSurfacePar->GetFC(), subSurfacePar->GetFCDEL());
      /*      if (GetLandScapeElement()->GetEvaporationControlObj()->GetEvaporationModellingControl() == 'M' || 
	  GetLandScapeElement()->GetEvaporationControlObj()->GetEvaporationModellingControl() == 'm') {
	i = datetime.getMonth()-1;
	transpSoilEvap = dryPeriod * (1.0 - snowCoverFraction) * 
	  HBVTranspSoilEvapLongTermMean(soilMoisture, temp, landSurfacePar->GetEPOT_PAR(), 
					    subSurfacePar->GetFC(), subSurfacePar->GetFCDEL(),
					    GetLandScapeElement()->GetEvaporationControlObj()->GetEvaporationArray(i)/1000.0);*/
	/*for (i=0; i<numberPotentialEvaporationValuesPerYear; i++) {
	  cout << " " << i << "    " << GetLandScapeElement()->GetEvaporationControlObj()->GetEvaporationArray(i) << "\n";
	  }*/
      /*      }
      else {
	transpSoilEvap = dryPeriod * (1.0 - snowCoverFraction) * 
	  HBVTranspSoilEvapTemperatureIndex(soilMoisture, temp, landSurfacePar->GetEPOT_PAR(), 
					    subSurfacePar->GetFC(), subSurfacePar->GetFCDEL());
					    }*/
    }
    else {
      transpSoilEvap = 0.0;
    }
    soilMoisture = soilMoisture - transpSoilEvap;
    if (soilMoisture < 0.0) {
      transpSoilEvap = transpSoilEvap + soilMoisture;
      soilMoisture = 0.0;
    }
    
    /* Water balance for soil moisture zone per hour, percolation to upper groundwater zone */
    if (inSoil > 0.0) {
      hoursPerTimeStep = commonPar->GetSECONDS_TIMESTEP() / minimumTimeStep;
      inSoilHour = inSoil / hoursPerTimeStep;
      //      cout << " hoursPerTimeStep : " << hoursPerTimeStep << " * " << minimumTimeStep << " = " << commonPar->GetSECONDS_TIMESTEP() << endl;
      //      cout << " inSoilHour : " << inSoilHour << " * " << hoursPerTimeStep << " = " << inSoil << endl;
      for (i=0; i<hoursPerTimeStep; i++) {
        if (soilMoisture < subSurfacePar->GetFC())
          outSoil = inSoilHour*pow(soilMoisture/subSurfacePar->GetFC(),subSurfacePar->GetBETA());
        else
          outSoil = inSoilHour;
        soilMoisture = soilMoisture + inSoilHour - outSoil;
        if (soilMoisture > subSurfacePar->GetFC()) {
          outSoil = outSoil + soilMoisture - subSurfacePar->GetFC();
          soilMoisture = subSurfacePar->GetFC();
        }
        upperPercolation = upperPercolation + outSoil;
      }
    }
    else {
      outSoil = 0.0;
    }
    //    cout << "upperPercolation not Glacier  " << upperPercolation << endl;
  }
  
  /* Soil moisture deficit not in use */
  else {
    /*    transpSoilEvap = 0.0;
          upperPercolation = waterInput;
          soilMoisture = subSurfacePar->GetFC();*/
    transpSoilEvap = 0.0;
    upperPercolation = 0.0;
    inSoil = waterInput;

    /* Water balance for soil moisture zone, percolation to upper groundwater zone */
    if (inSoil > 0.0) {
      if (soilMoisture < subSurfacePar->GetFC())
        outSoil = inSoil*power(soilMoisture/subSurfacePar->GetFC(),subSurfacePar->GetBETA());
      else
        outSoil = inSoil;
      soilMoisture = soilMoisture + inSoil - outSoil;
      if (soilMoisture > subSurfacePar->GetFC()) {
        outSoil = outSoil + soilMoisture - subSurfacePar->GetFC();
        soilMoisture = subSurfacePar->GetFC();
      }
      upperPercolation = upperPercolation + outSoil;
    }
    else {
      outSoil = 0.0;
    }
    //    cout << "Glacier:   upperPercolation " << upperPercolation << "  ";
  }

  /* Water balance for upper groundwater zone, percolation to lower groundwater zone */
  lowerPercolation = subSurfacePar->GetPERC();
  //  if (dryPeriod == missingData) cout << " subSurfacePar->GetPERC() " << subSurfacePar->GetPERC() << "  ";
  //  if (dryPeriod == missingData) cout << " lowerPercolation " << lowerPercolation << "  ";
  upperTemporary = upperZone + 0.5*(upperPercolation-lowerPercolation);
  //  if (dryPeriod == missingData) cout << " upperTemporary " << upperTemporary << "  ";
  if (lowerPercolation > upperZone+upperPercolation) {
    lowerPercolation = upperZone + upperPercolation;
    upperRunoff = 0.0;
  }
  else {
    upperRunoff = subSurfacePar->GetKUZ() * pow(upperTemporary,subSurfacePar->GetALFA());
  }
  upperZone = upperZone + upperPercolation - lowerPercolation - upperRunoff;
  if (upperZone < 0.0) {
    upperRunoff = upperRunoff + upperZone;
    if (upperRunoff < 0.0) {
      lowerPercolation = lowerPercolation + upperRunoff;
      upperRunoff = 0.0;
    }
    if (lowerPercolation < 0.0) lowerPercolation = 0.0;
    upperZone = 0.0;
  }
  //  if (dryPeriod == missingData) cout << " lowerPercolation " << lowerPercolation << "  " << endl;

  /* Water balance for lower groundwater zone */
  lowerZone = lowerZone + lowerPercolation;
  if (lowerZone < 0.0) lowerZone = 0.0;
  lowerRunoff = lowerZone * subSurfacePar->GetKLZ();
  lowerZone = lowerZone - lowerRunoff;

  /* Add direct runoff to upper runoff */
  //  upperRunoff = upperRunoff + directRunoff;

  /* Runoff from both zones */
  runoff = lowerRunoff + upperRunoff;

  /* Percolation from soil moisture zone to upper zone */
  percSoilUpper = upperPercolation;
}

