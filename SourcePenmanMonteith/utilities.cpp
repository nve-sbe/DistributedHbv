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


// Actual evapotranspiration HBV
double HBVTranspSoilEvap(double soilMoist, double epot, double fieldCapacity, double fcDel)
{
  // double epot;
  // epot = potentialEvap(temp,epotPar);            
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

// Potential evapotranspiration (Penmann Monteith)
double potentialEvap(double temp, double tempMax, double tempMin, double radiationS, double vp, 
		     double wind1, int dayofyear, double elemLatitute, double elemElevation, double windHeight, double humidityHeight,
		     double treeHeight, double treeLai, double treeLaiCorr, double treeAlbedo, double bulkResist,
		     double albedo_snow,double deciduous, double snowCoverFraction, double tree_wind_H, double Topen_min, double Tclose_min, double VPDclose,
		     double VPDopen, double gh, double cl, double z0g,
		     double gsmax, double CR, double D50, double Q50)
{
  // variables
  float Rns;
  double delta, vs_Tmax, vs_Tmin, vas, vabar, lat_rad;
  double d_r2, delta2, w_s, R_a, R_so, R_nl, Rn,sun_hour;
  double lv, pz, gamma, airDensity, albedo;
  double aeroResist, surResist, ETP, Windspeed_tree;
  double mTmin, mVPD, mRad, sradWm2, Rcorr, Gs, z0c, cd, d, z0;
  double Rns_snow, ETP_snow, Rn_snow, ETP_tot;
  
  
  //if (temp > 0.0) {    // need to consider this precondition later, but there are many errors for high latitudes without it
  /*####### slope of the saturation vapour pressure calculated by VIC. Unit: 1000Pa
    
    svp <-A_svp*exp((B_svp * Tair) / (C_svp + Tair))
    if (Tair < 0) {
    svp <-svp * (1.0 + 0.00972 * Tair + 0.000042 * Tair * Tair)   ## need to check
    }
    
    delta <-(B_svp*C_svp) / ((C_svp + Tair) * (C_svp + Tair))* svp */
  
  //###### slope of the slope of the saturation vapour pressure calculated using S2.4. Unit: 1000Pa
  
  delta = 4098 * (0.6108 * exp((17.27 * temp) / (temp + 237.3))) / pow((temp + 237.3), 2);
  
  // std::cout << std::endl << " delta =  " << delta << std::endl << std::endl;
  
  //###### The results from these two methods are quite similar when Tair >0 ########
  //###### need to decide why VIC has concern for Tair < 0 ##########################################
  
  
  //###### Actual vapour pressure is estimated by Mtclim in VIC   ####
  //###### It can also be estimated using the maximum and minimum humidity  using S2.7 and FAO method, unit:1000Pa
  
  vs_Tmax = 0.6108 * exp(17.27 * tempMax / (tempMax + 237.3));
  vs_Tmin = 0.6108 * exp(17.27 * tempMin / (tempMin + 237.3));
  vas = (vs_Tmax + vs_Tmin) / 2;        // saturation vapour pressure(kPa)
  vabar = vp;
  // vabar = (vs_Tmin * humidityMax / 100 + vs_Tmax * humidityMin / 100) / 2;
  
  // std::cout << std::endl << " vas =  " << vas << std::endl << std::endl;
  //std::cout << std::endl << " vabar =  " << vabar << std::endl << std::endl;
  //## if RHmin has large errors, then only use RHmax, RHmean can also be used only when RHmax is not available
  // vabar <-vs_Tmin * RHmax / 100
  
  //###### if humidity is not availabe, vabar can also be estimated by dew point temperature
  
  // vabar <-0.6108 * exp(17.27 * Td / (Td + 237.3))
  
  //###### Albedo estimation based on snow water equalvalent
  
  albedo = treeAlbedo;
  //if(treeLai < 4. & treeHeight > 0.1) { 
  //  albedo = 0.1 + 0.25 * treeLai * (treeAlbedo - 0.1);  // equation from AMOR model (D1), 0.1 is bare soil albedo
  //}	
  
  //######  net radiation using S3.2 - S3.8. Unit: MJ/m2/day
  Rns = (1. - albedo) * radiationS;
  Rns_snow = (1. - albedo_snow) * radiationS;
  
  lat_rad = elemLatitute * M_PI / 180.;   // lat unit : radian
  d_r2 = 1. + 0.033 * cos(2 * M_PI / 365 * dayofyear);  // inverse relative distance Earth - Sun
  delta2 = 0.409 * sin(2 * M_PI / 365 * dayofyear - 1.39);   // solar declination(rad)
  
  sun_hour = -tan(lat_rad) * tan(delta2);
  if (sun_hour < -1.) {
    sun_hour = -1.;
  }
  if (sun_hour > 1.) {
    sun_hour = 1.;
  }
  w_s = acos(sun_hour);   // sunset hour angle(rad)
  // N = 24 / pi * w_s
  R_a = (1440. / M_PI) * d_r2 * 0.082 * (w_s * sin(lat_rad) * sin(delta2) + cos(lat_rad) * cos(delta2) *sin(w_s));
  R_so = (0.75 + (2e-5) * elemElevation) * R_a;
  
  if (R_so < radiationS) {
    R_so = radiationS;
  }
  
  if (R_so < 0.01) {     //radiationS can be 0 in high latitude in winter
    R_so = 0.01;
  }
  
  R_nl = 4.903e-09 * (0.34 - 0.14 * sqrt(vabar)) * (pow((tempMax + 273.2), 4) + pow((tempMin + 273.2), 4)) / 2 * (1.35 * radiationS / R_so - 0.35);
  //Rs / R_so <= 1
  Rn = Rns - R_nl;
  Rn_snow = Rns_snow - R_nl;
  
  //std::cout << std::endl << " R_a =  " << R_a << std::endl << std::endl;
  
  //std::cout << std::endl << " Rn =  " << Rn << std::endl << std::endl;
  //###### latent heat of vaporization estimated by VIC, unit: MJ/kg; can be a constant 2.45 MJ/kg
  
  lv = 2.501000 - 0.002361 * temp;
  
  // std::cout << std::endl << " lv =  " << lv << std::endl << std::endl;
  
  /*###### Atmospheric pressure estimated by VIC, unit: 1000Pa
    ###########
    
    scaleH <-287 / 9.81*((Tair + 273.15) + 0.5 * Elev*LAPSE_PM)   # Elev should be in float
    pz <-101.3 * exp(-Elev / scaleH) */
  
  
  //########## S2.10 estimates pz as below, unit : 1000Pa
  pz = 101.3 * pow(((293. - 0.0065*elemElevation) / 293.), 5.26);
  
  // std::cout << std::endl << " pz =  " << pz << std::endl << std::endl;
  
  //################# The results are also similar, S2.10 is simpler
  
  //###### gamma, using the lv from VIC and pz from S2.10
  gamma = 0.00163 * pz / lv;
  
  // std::cout << std::endl << " gamma =  " << gamma << std::endl << std::endl;
  
  //###### air density from VIC, unit: kg/m3
  airDensity = 3.486 * pz / (275. + temp);
  
  //###### aeroResist using dynamic resistance from Stephanie
  
  Windspeed_tree = (log(tree_wind_H / 0.03) / log(windHeight / 0.03))*wind1;
  
  if(treeHeight <= 1.) {
    z0c = 0.13 * treeHeight;
  } else if(treeHeight >= 10.) {
    z0c = 0.05 * treeHeight; 
  } else  {
    z0c = 0.089 + 0.041 * treeHeight;
  }
  
  if(treeHeight < 0.001)  {
    cd = 0.0014;
  } else {
    cd = pow((-1. + exp(0.909-3.03*z0c/treeHeight)),4)/4;
  }
  
  if(treeLai < 4.) {
    d = 1.1 * treeHeight * log(1 + pow((cd*treeLai),0.25));
  } else {
    d = treeHeight - z0c/0.3;
  }
  
  z0 = min(0.3*(treeHeight - d), z0g + 0.3* treeHeight*pow(cd*treeLai,0.5)) ;
  
  aeroResist = log((tree_wind_H - d) / z0)*log((tree_wind_H - d) / (0.1*z0)) / 0.41 / 0.41 / Windspeed_tree ;
  
  //##### surResist based on equation 12 in Leuning et al.
  
  Rcorr = 1 / (101300 / pz / 1000 * pow(((temp + 273.15) / 293.15), 1.75));
  
  if (treeLai > 0.01) {
    
    sradWm2 = radiationS / 86400 *1.0e6; // ### convert radiationS to unit W m-2
    
    mRad = 1 / CR * log((0.8*sradWm2+Q50) / (0.8*sradWm2*exp(-1*CR*treeLai)+Q50));
    
    mVPD = 1 / (1+(vas-vabar)/D50);
    
    if (temp >= Topen_min) mTmin = 1;
    if (temp > Tclose_min && temp < Topen_min) mTmin = (temp - Tclose_min) / (Topen_min - Tclose_min);
    if (temp <= Tclose_min) mTmin = 0;
    
    Gs = gsmax*mRad*mVPD*mTmin*Rcorr;
    
    if (Gs > 0) {
      surResist = 1 / Gs;
    }
    else {
      surResist = 10000;
    }
    
  }
  else {
    surResist = Rcorr*107;   // bare soil resistance from Mu et al.
  }
  
  // std::cout << std::endl << " surResist =  " << surResist << std::endl << std::endl;
  //##### Penmann Monteith equation (mm/day)
  
  //Cp <-gamma*0.622*lv / pz
  
  ETP =  (delta*Rn + numberSecondsDay*airDensity*0.001013*(vas - vabar) / aeroResist) / lv / (delta + gamma*(1 + surResist / aeroResist));
  ETP_snow =  (delta*Rn_snow + numberSecondsDay*airDensity*0.001013*(vas - vabar) / aeroResist) / lv / (delta + gamma*(1 + surResist / aeroResist));
  
  if(snowCoverFraction > 0.) 
    ETP_tot = ETP * (1.0 - snowCoverFraction) + ETP_snow * snowCoverFraction;
  else 
    ETP_tot = ETP;

  if ((vas - vabar) >= 0.0 && ETP < 0.0) 
    ETP = 0.0;
  
  if(ETP < 0.00001) 
    ETP = 0.0;
  
  return ETP;	
}

// Potential evapotranspiration (Penmann Monteith) for lakes
double potentialEvapLake(double temp, double tempMax, double tempMin, double radiationS, double vp, 
	double wind1, int dayofyear, double elemLatitute, double elemElevation, double windHeight)
{
	// variables
	float Rns;
	double delta, vs_Tmax, vs_Tmin, vas, vabar, lat_rad;
	double d_r2, delta2, w_s, R_a, R_so, R_nl, Rn, sun_hour;
	double lv, pz, gamma, airDensity;
	double aeroResist, surResist, ETP;

	//if (temp > 0.0) {    // need to consider this precondition later, but there are many errors for high latitudes without it
						 /*####### slope of the saturation vapour pressure calculated by VIC. Unit: 1000Pa

						 svp <-A_svp*exp((B_svp * Tair) / (C_svp + Tair))
						 if (Tair < 0) {
						 svp <-svp * (1.0 + 0.00972 * Tair + 0.000042 * Tair * Tair)   ## need to check
						 }

						 delta <-(B_svp*C_svp) / ((C_svp + Tair) * (C_svp + Tair))* svp */

						 //###### slope of the slope of the saturation vapour pressure calculated using S2.4. Unit: 1000Pa

		delta = 4098 * (0.6108 * exp((17.27 * temp) / (temp + 237.3))) / pow((temp + 237.3), 2);

		// std::cout << std::endl << " delta =  " << delta << std::endl << std::endl;

		//###### The results from these two methods are quite similar when Tair >0 ########
		//###### need to decide why VIC has concern for Tair < 0 ##########################################


		//###### Actual vapour pressure is estimated by Mtclim in VIC   ####
		//###### It can also be estimated using the maximum and minimum humidity  using S2.7 and FAO method, unit:1000Pa

		vs_Tmax = 0.6108 * exp(17.27 * tempMax / (tempMax + 237.3));
		vs_Tmin = 0.6108 * exp(17.27 * tempMin / (tempMin + 237.3));
		vas = (vs_Tmax + vs_Tmin) / 2;        // saturation vapour pressure(kPa)
		vabar = vp;
		//vabar = (vs_Tmin * humidityMax / 100 + vs_Tmax * humidityMin / 100) / 2;

		// std::cout << std::endl << " vas =  " << vas << std::endl << std::endl;
		//std::cout << std::endl << " vabar =  " << vabar << std::endl << std::endl;
		//## if RHmin has large errors, then only use RHmax, RHmean can also be used only when RHmax is not available
		// vabar <-vs_Tmin * RHmax / 100

		//###### if humidity is not availabe, vabar can also be estimated by dew point temperature

		// vabar <-0.6108 * exp(17.27 * Td / (Td + 237.3))

		//######  net radiation using S3.2 - S3.8. Unit: MJ/m2/day
		Rns = (1. - 0.22) * radiationS;     // 0.22 is from http://agsys.cra-cin.it/tools/SolarRadiation/help/Albedo.html, for water at 10 degree. The snow albedo should be added later

		lat_rad = elemLatitute * M_PI / 180.;   // lat unit : radian
		d_r2 = 1. + 0.033 * cos(2 * M_PI / 365 * dayofyear);  // inverse relative distance Earth - Sun
		delta2 = 0.409 * sin(2 * M_PI / 365 * dayofyear - 1.39);   // solar declination(rad)

		sun_hour = -tan(lat_rad) * tan(delta2);
		if (sun_hour < -1.) {
			sun_hour = -1.;
		}
		if (sun_hour > 1.) {
			sun_hour = 1.;
		}
		w_s = acos(sun_hour);   // sunset hour angle(rad)
								// N = 24 / pi * w_s
		R_a = (1440. / M_PI) * d_r2 * 0.082 * (w_s * sin(lat_rad) * sin(delta2) + cos(lat_rad) * cos(delta2) *sin(w_s));
		R_so = (0.75 + (2e-5) * elemElevation) * R_a;

		if (R_so < radiationS) {
			R_so = radiationS;
		}

		if (R_so < 0.01) {     //radiationS can be 0 in high latitude in winter
			R_so = 0.01;
		}

		R_nl = 4.903e-09 * (0.34 - 0.14 * sqrt(vabar)) * (pow((tempMax + 273.2), 4) + pow((tempMin + 273.2), 4)) / 2 * (1.35 * radiationS / R_so - 0.35);
		//Rs / R_so <= 1
		Rn = Rns - R_nl;

		//std::cout << std::endl << " R_a =  " << R_a << std::endl << std::endl;

		//std::cout << std::endl << " Rn =  " << Rn << std::endl << std::endl;
		//###### latent heat of vaporization estimated by VIC, unit: MJ/kg; can be a constant 2.45 MJ/kg

		lv = 2.501000 - 0.002361 * temp;

		// std::cout << std::endl << " lv =  " << lv << std::endl << std::endl;

		/*###### Atmospheric pressure estimated by VIC, unit: 1000Pa
		###########

		scaleH <-287 / 9.81*((Tair + 273.15) + 0.5 * Elev*LAPSE_PM)   # Elev should be in float
		pz <-101.3 * exp(-Elev / scaleH) */


		//########## S2.10 estimates pz as below, unit : 1000Pa
		pz = 101.3 * pow(((293. - 0.0065*elemElevation) / 293.), 5.26);

		// std::cout << std::endl << " pz =  " << pz << std::endl << std::endl;

		//################# The results are also similar, S2.10 is simpler

		//###### gamma, using the lv from VIC and pz from S2.10
		gamma = 0.00163 * pz / lv;

		// std::cout << std::endl << " gamma =  " << gamma << std::endl << std::endl;

		//###### air density from VIC, unit: kg/m3
		airDensity = 3.486 * pz / (275. + temp);

		//###### aeroResist using equation S5.3 and S5.7, in HBV the lake and glaicer was calculated separately
		//	if (waterID > 0) {
				aeroResist = 4.72*pow(pow(log(windHeight / 1.37e-3), 2), 2) / (1. + 0.54*wind1);
		//	}
		//	else {
		// aeroResist = log((windHeight - 0.67*treeHeight) / 0.123*treeHeight)*log((humidityHeight - 0.67*treeHeight) / 0.0123*treeHeight) / 0.41 / 0.41 / wind1;
		//	}

		// std::cout << std::endl << " aeroResist =  " << aeroResist << std::endl << std::endl;

		//##### surResist using equation S5.9
		// if (treeLai < 0.00001) {
		//	surResist = 1.e20;
		// }
		// else {
		//	surResist = bulkResist / (treeLaiCorr*treeLai);
		// }

		//	if (waterID > 0) {
				surResist = 0.;
		//	}    //surResist = 0 for water surface

		// std::cout << std::endl << " surResist =  " << surResist << std::endl << std::endl;
		//##### Penmann Monteith equation (mm/day)

		//Cp <-gamma*0.622*lv / pz

		ETP = (delta*Rn + numberSecondsDay*airDensity*0.001013*(vas - vabar) / aeroResist) / lv / (delta + gamma*(1 + surResist / aeroResist));
	//}
	//else {

		//ETP = 0.0;
	//}


	if ((vas - vabar) >= 0.0 && ETP < 0.0) {
		ETP = 0.0;
	}
	
		if(ETP < 0.00001) {
		ETP = 0.0;
	}

	/* cout << endl << " Tair =  " << temp << endl;
	cout << endl << " Tmax =  " << tempMax << endl;
	cout << endl << " Tmin =  " << tempMin << endl;
	cout << endl << " Rs =  " << radiationS << endl;
	cout << endl << " Rn_l =  " << R_nl << endl;
	cout << endl << " vp =  " << vp << endl;
	cout << endl << " vas =  " << vas << endl;
	cout << endl << " Wind =  " << wind1 << endl;
	cout << endl << " DoY =  " << dayofyear << endl;
	cout << endl << " LAT =  " << elemLatitute << endl;
	cout << endl << " Elev =  " << elemElevation << endl;
	cout << endl << " ETP =  " << ETP << endl << endl;    */

	return ETP;

}
