#include "classDistHbv.h"
#include "utilities.h"

void Vegetation::WaterBalance(int timeStep, DateTime datetime, int dayofyear_PM2, double snowCoverFraction)
{
  double timeResolution=1.0;
  double delg = 0.;
  double par = 0.;
  double ddm = 0.;
  double f = 0.;
  double ff= 0.;
  double deltalai = 0.;
  double yy = 0.;
  double zz = 0.;
  double int_max = 0.;

  int t1, t2;
  precipitation=GetInputElement()->GetInput(0);
  temp = GetInputElement()->GetInput(1);
  tempMax = GetInputElement()->GetInput(2);
  tempMin = GetInputElement()->GetInput(3);
  wind1 = GetInputElement()->GetInput(4);
  radiationS = GetInputElement()->GetInput(5);
  vp = GetInputElement()->GetInput(6);


  //  cout << "    timeStep " << timeStep << "    Vegetation precipitation " << precipitation;
  //  cout << "    Temperature " << temp << endl;

  /* Potential evapotranspiration */
  wetPeriod = 0.0;
  dryPeriod = timeResolution;
  throughFall = 0.0;
  interceptionLoss = 0.0;
  
  // LAI for deciduous forest
  		  t1 = GetLandScapeElement()->GetTimeFalling();
		    t2 = t1 + 13;
    //      cout << "    t1, t2 " << t1 << " " << t2 << endl;
      
 if (landSurfacePar->GetDECIDUOUS_SHARE() == 1) {      //deciduous forest
	  // calculate the budburst time
	  if (dayofyear_PM2 == 1) {
		  sft = 0.;
		  budburst = 0;

	  }

  if(budburst==0) calc_lai = 0;

	  if (dayofyear_PM2 >= 65 && temp > 0 && sft <= 198.) {
		  sft = sft + temp;
	  }

	  if ((sft > 198. && budburst < 15) | dayofyear_PM2 >= 181) {
		  budburst = budburst + 1;
      calc_lai = double(budburst) / 15. * landSurfacePar->GetTREE_LAI();
	  }
      //  cout << "    budburst " << budburst << endl;
        

	  if (budburst >= 15) {
		  calc_lai = landSurfacePar->GetTREE_LAI();
	  }

	  // calculate the autumn leaf fall
	  if (dayofyear_PM2 > 182) {
		  calc_lai = landSurfacePar->GetTREE_LAI() *(1 - 1 / (1 + exp(2.2 / (t2 - t1)*(t2 - dayofyear_PM2))));
	  }

        calc_lai = calc_lai+ 1.5*exp(-0.5*calc_lai);

  }

    if(landSurfacePar->GetDECIDUOUS_SHARE() == 2) {        // crops
	  if(dayofyear_PM2 >120 && dayofyear_PM2 < 300) {
// COMPUTE DAILY INCREASE IN HEAT UNITS delg
      delg = (temp-0.)/1900.;        // 0. is the base temperature for spring barley to grow, and 1900 is the potential heat units for the maturity of the spring barley
      if (delg < 0.) delg = 0.;
      g = g + delg ;           
      if (g > 1.) g = 1. ;  
      
// GROWTH SEASON
      if (g <= 1.) {
      
//   CALC daily biomass increase: ddm
        par = .02092 * (radiationS/0.04184)  * (1.-exp(-.65*(calc_lai+.05)));            // 0.04184 conversion factor from MJ/m2 to Ly
        ddm = 30. * par ;       // 30 is the biomass-energy ratio for spring barley
        if (ddm < 0.) ddm = 0. ;
        
//   CALC BIOMASS dm(), ROOT WEIGHT rwt(9 
      
        dm = dm + ddm ;
        
//  CALC f, ff, huharv()
        f = g / (g + exp(5.41325-18.10167*g));    // for spring barley only
        ff = f - huharv;
        huharv = f ;

        
        
// CALC LAI and adjust for lower limit of LAI for forest alnm()
        if (g <= 0.9) {
          if (calc_lai > 6.) calc_lai = 6. ;
          deltalai = ff * 6.0 * (1.-exp(5.*(calc_lai - 6.))) ;
          calc_lai = calc_lai + deltalai   ;
          if (calc_lai > 6.)  calc_lai = 6.  ;
          olai = calc_lai ;
          }
        else   {
          yy = sqrt(1.- g) ;
          zz = 1./sqrt(1.-0.9)  ;
          calc_lai = zz * olai * yy ;
       }
           
      }
      }   else {       // dormant season
      g = 0. ;
      calc_lai = 0. ;
      dm  = 0. ;
       huharv = 0. ;
       olai  = 0. ;
      }
  }   
  
  
    if(landSurfacePar->GetDECIDUOUS_SHARE() == 0) {      // other vegetation
	  calc_lai = landSurfacePar->GetTREE_LAI();
  }

    if(landSurfacePar->GetDECIDUOUS_SHARE() == 3) {      // evergreen forests
	  calc_lai = landSurfacePar->GetTREE_LAI()+ 1.5*exp(-0.5*landSurfacePar->GetTREE_LAI());
  }



     if(calc_lai < 0.001) {
     calc_lai = 0.;
     }

    // Maximum intercetion capacity as a function of leaf area index
     if(landSurfacePar->GetTREE_LAI()>0.01) {
     int_max = landSurfacePar->GetINTER_MAX() * calc_lai/landSurfacePar->GetTREE_LAI();
     } 
       else {
     int_max = landSurfacePar->GetINTER_MAX();
     }
     
  potev = potentialEvap (temp, tempMax, tempMin, radiationS, vp, 
	                     wind1, dayofyear_PM2, GetLandScapeElement()->GetLatitude(), GetLandScapeElement()->GetElevation(),
	  commonPar->GetHeight_WIND_INSTRUMENT(), commonPar->GetHeight_HUMI_INSTRUMENT(),
	  landSurfacePar->GetTREE_HEIGHT(), calc_lai, landSurfacePar->GetTREE_LAI_CORR(), 
	  landSurfacePar->GetALBEDO(), landSurfacePar->GetBULK_RESISTANCE(),landSurfacePar->GetALBEDO_SNOW(),
	  landSurfacePar->GetDECIDUOUS_SHARE(), snowCoverFraction,landSurfacePar->GetWind_H(),landSurfacePar->GetTopen_min(), landSurfacePar->GetTclose_min(), landSurfacePar->GetVPDclose(), landSurfacePar->GetVPDopen(),
			 landSurfacePar->GetGH(), landSurfacePar->GetCL(), landSurfacePar->GetZ0G(), landSurfacePar->GetGSMAX(), landSurfacePar->GetCR(), landSurfacePar->GetD50(),   landSurfacePar->GetQ50()) / 1000.;  // convert from mm/day to m/day 

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
    // if (interceptionStore > landSurfacePar->GetINTER_MAX()) {
    //  throughFall = interceptionStore - landSurfacePar->GetINTER_MAX();
    //  interceptionStore = landSurfacePar->GetINTER_MAX();
    // }

    if (interceptionStore > int_max) {
      throughFall = interceptionStore - int_max;
      interceptionStore = int_max;
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

