#include "parameters.h"
#include "classDistHbv.h"


void SetLandSurfaceParameters(ParametersLandSurface * const ParLandSurfaceStore, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[512];
  int i,j,k;
  double interMax, epotPar, wetPerCorr;
  double accTemp, meltTemp, snowMeltRate, iceMeltRate, freezeEff, maxRel, albedo, cvSnow;
  double tree_height, tree_lai, tree_lai_corr, bulk_resistance,albedo_snow,deciduous_share;
  double wind_h,Topen_min, Tclose_min, VPDclose, VPDopen, gh, cl, z0g;
  double gsmax, CR, D50, Q50;

  /*  cout << " File with landsurface parameters: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(512,'\n');
  ifstream finLandSurfacePar(fileName);  // Open for reading
  if (!finLandSurfacePar.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finLandSurfacePar.getline(buffer, 512);
  for (i=0; i<numberLandSurfaceClasses; i++) {
    finLandSurfacePar >> buffer >> j >> interMax >> epotPar >> wetPerCorr >> accTemp >> meltTemp ;
    finLandSurfacePar >> snowMeltRate >> iceMeltRate >> freezeEff >> maxRel >> albedo >> cvSnow;
    finLandSurfacePar >> tree_height >> tree_lai >> tree_lai_corr >> bulk_resistance >> albedo_snow >> deciduous_share;
    finLandSurfacePar >> wind_h >> Topen_min >> Tclose_min >> VPDclose >> VPDopen >> gh >> cl >> z0g;
    finLandSurfacePar >> gsmax >> CR >> D50 >> Q50;
    if (i!=j) {
      cout << endl << " Error in land surface parameter file, parameter no. " 
           << i << endl << endl;
      exit (1);
    } 
    ParLandSurfaceStore[i].SetINTER_MAX(interMax);
    ParLandSurfaceStore[i].SetEPOT_PAR(epotPar);
    ParLandSurfaceStore[i].SetWET_PER_CORR(wetPerCorr);
    ParLandSurfaceStore[i].SetACC_TEMP(accTemp);
    ParLandSurfaceStore[i].SetMELT_TEMP(meltTemp);
    ParLandSurfaceStore[i].SetSNOW_MELT_RATE(snowMeltRate);
    ParLandSurfaceStore[i].SetICE_MELT_RATE(iceMeltRate);
    ParLandSurfaceStore[i].SetFREEZE_EFF(freezeEff);
    ParLandSurfaceStore[i].SetMAX_REL(maxRel);
    ParLandSurfaceStore[i].SetALBEDO(albedo);
    ParLandSurfaceStore[i].SetCV_SNOW(cvSnow);
    ParLandSurfaceStore[i].SetTREE_HEIGHT(tree_height);
    ParLandSurfaceStore[i].SetTREE_LAI(tree_lai);
    ParLandSurfaceStore[i].SetTREE_LAI_CORR(tree_lai_corr);
    ParLandSurfaceStore[i].SetBULK_RESISTANCE(bulk_resistance);
    ParLandSurfaceStore[i].SetALBEDO_SNOW(albedo_snow);
    ParLandSurfaceStore[i].SetDECIDUOUS_SHARE(deciduous_share);
    ParLandSurfaceStore[i].SetWind_H(wind_h);
    ParLandSurfaceStore[i].SetTopen_min(Topen_min);
    ParLandSurfaceStore[i].SetTclose_min(Tclose_min);
    ParLandSurfaceStore[i].SetVPDclose(VPDclose);
    ParLandSurfaceStore[i].SetVPDopen(VPDopen);
    ParLandSurfaceStore[i].SetGH(gh);
    ParLandSurfaceStore[i].SetCL(cl);
    ParLandSurfaceStore[i].SetZ0G(z0g);
    ParLandSurfaceStore[i].SetGSMAX(gsmax);
    ParLandSurfaceStore[i].SetCR(CR);
    ParLandSurfaceStore[i].SetD50(D50);
    ParLandSurfaceStore[i].SetQ50(Q50);
    // Set snow distribution parameters
    SetSnowDistribution(&ParLandSurfaceStore[i], cvSnow);
  }
  finLandSurfacePar.close();

  fout << "Land surface parameters: \n";
  for (i=0; i<numberLandSurfaceClasses; i++) {
    fout << ParLandSurfaceStore[i].GetINTER_MAX() << "    ";
    fout << ParLandSurfaceStore[i].GetEPOT_PAR() << "    ";
    fout << ParLandSurfaceStore[i].GetWET_PER_CORR() << "    ";
    fout << ParLandSurfaceStore[i].GetACC_TEMP() << "    ";
    fout << ParLandSurfaceStore[i].GetMELT_TEMP() << "    ";
    fout << ParLandSurfaceStore[i].GetSNOW_MELT_RATE() << "    ";
    fout << ParLandSurfaceStore[i].GetICE_MELT_RATE() << "    ";
    fout << ParLandSurfaceStore[i].GetFREEZE_EFF() << "    ";
    fout << ParLandSurfaceStore[i].GetMAX_REL() << "    ";
    fout << ParLandSurfaceStore[i].GetALBEDO() << "    ";
    fout << ParLandSurfaceStore[i].GetCV_SNOW() << endl;
    fout << ParLandSurfaceStore[i].GetTREE_HEIGHT() << endl;
    fout << ParLandSurfaceStore[i].GetTREE_LAI() << endl;
    fout << ParLandSurfaceStore[i].GetTREE_LAI_CORR() << endl;
    fout << ParLandSurfaceStore[i].GetBULK_RESISTANCE() << endl;
    fout << ParLandSurfaceStore[i].GetALBEDO_SNOW() << endl;
    fout << ParLandSurfaceStore[i].GetDECIDUOUS_SHARE() << endl;
    fout << ParLandSurfaceStore[i].GetWind_H() << endl;
    fout << ParLandSurfaceStore[i].GetTopen_min() << endl;
    fout << ParLandSurfaceStore[i].GetTclose_min() << endl;
    fout << ParLandSurfaceStore[i].GetVPDclose() << endl;
    fout << ParLandSurfaceStore[i].GetVPDopen() << endl;
    fout << ParLandSurfaceStore[i].GetGH() << endl;
    fout << ParLandSurfaceStore[i].GetCL() << endl;
    fout << ParLandSurfaceStore[i].GetZ0G() << endl;
    fout << ParLandSurfaceStore[i].GetGSMAX() << endl;
    fout << ParLandSurfaceStore[i].GetCR() << endl;
    fout << ParLandSurfaceStore[i].GetD50() << endl;
    fout << ParLandSurfaceStore[i].GetQ50() << endl;
    for (k=0; k<numberSnowClasses; k++) {
      fout << ParLandSurfaceStore[i].GetSNOW_WEIGHT(k) << "  ";
    } 
    fout << endl;
  }
  fout << endl;
}


void SetSnowDistribution(ParametersLandSurface * thisParLandSurface, double cvSnow)
{
  int k;
  double stdDevNorm, meanNorm, sumNorm, sumNorm2;
  double stdNormVar[numberSnowClasses], logNormWeight[numberSnowClasses];
  stdNormVar[0]=-2.326347;
  stdNormVar[1]=-1.644853476;
  stdNormVar[2]=-1.036433474;
  stdNormVar[3]=-0.385320604;
  stdNormVar[4]=0.385320604;
  stdNormVar[5]=1.036433474;
  stdNormVar[6]=1.644853476;
  stdNormVar[7]=2.326347;
  stdNormVar[8]=3.71909027;
  //  stdNormVar[8]=4.265043367;
  meanNorm = 0.5*log(1.0/(1.0+cvSnow*cvSnow));
  stdDevNorm = sqrt(log(1.0+cvSnow*cvSnow));
  sumNorm = 0.0;
  for (k=0; k<numberSnowClasses; k++) {
    logNormWeight[k] = exp(stdNormVar[k]*stdDevNorm + meanNorm);
    sumNorm = sumNorm + logNormWeight[k]*probNorm[k];
  } 
  sumNorm2 = 0.0;
  for (k=0; k<numberSnowClasses; k++) {
    logNormWeight[k] = logNormWeight[k]/sumNorm;
    sumNorm2 = sumNorm2 + logNormWeight[k]*probNorm[k];
    thisParLandSurface->SetSNOW_WEIGHT(k,logNormWeight[k]);
    //    cout << k << "  " << logNormWeight[k] << endl;
  } 
  if (sumNorm2 < 1.0-epsilon || sumNorm2 > 1.0+epsilon) {
    cout << endl << " Sum of snow distribution weights = " << sumNorm2 << endl << endl;
    exit(1);
  }
}


void SetSubSurfaceHbvParameters(ParametersSubSurfaceHbv * const ParSubSurfaceHbvStore, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[256];
  int i, j;
  double fc, fcdel, beta, infmax; 
  double kuz, alfa, perc, klz, draw;

  /*  cout << " File with HBV subsurface parameters: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finSubSurface(fileName);  // Open for reading
  if (!finSubSurface.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finSubSurface.getline(buffer, 256);
  for (i=0; i<numberSoilClasses; i++) {
    finSubSurface >> buffer >> j >> fc >> fcdel >> beta >> infmax;
    finSubSurface >> kuz >> alfa >> perc >> klz >> draw;
    if (i!=j) {
      cout << endl << " Error in subsurface parameter file, parameter no. " 
           << i << endl << endl;
      exit (1);
    } 
    ParSubSurfaceHbvStore[i].SetFC(fc);
    ParSubSurfaceHbvStore[i].SetFCDEL(fcdel);
    ParSubSurfaceHbvStore[i].SetBETA(beta);
    ParSubSurfaceHbvStore[i].SetINFMAX(infmax);
    ParSubSurfaceHbvStore[i].SetKUZ(kuz);
    ParSubSurfaceHbvStore[i].SetALFA(alfa);
    ParSubSurfaceHbvStore[i].SetPERC(perc);
    ParSubSurfaceHbvStore[i].SetKLZ(klz);
    ParSubSurfaceHbvStore[i].SetDRAW(draw);
  }
  finSubSurface.close();
  fout << "Subsurface parameters: \n";
  for (i=0; i<numberSoilClasses; i++) {
    fout << ParSubSurfaceHbvStore[i].GetFC() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetFCDEL() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetBETA() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetINFMAX() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetKUZ() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetALFA() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetPERC() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetKLZ() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetDRAW() << "\n";
  }
  fout << endl;
}


void SetGeneralParameters(ParametersGeneral * const ParGeneralStore, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  double precGradLow, precGradHigh, gradChangeAltitude, lapseDry, lapseWet; 
  double precCorrRain, precCorrSnow;
  double dayTempMem, lakeEpotPar, kLake, deltaLevel, nLake;
  double initialSoilMoisture, initialUpperZone, initialLowerZone;
  double initialLakeTemp, initialLakeLevel, initialSnow;
  int secondsTimestep, numPrec, numTemp, daySnowZero;
  float height_wind_instrument, height_humi_instrument;

  /*  cout << " File with common parameters: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finGeneralPar(fileName);  // Open for reading
  if (!finGeneralPar.is_open()) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finGeneralPar.ignore(100,':'); finGeneralPar >> secondsTimestep; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> numPrec; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> numTemp; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> precGradLow;
  finGeneralPar.ignore(100,':'); finGeneralPar >> precGradHigh;
  finGeneralPar.ignore(100,':'); finGeneralPar >> gradChangeAltitude;
  finGeneralPar.ignore(100,':'); finGeneralPar >> precCorrRain;
  finGeneralPar.ignore(100,':'); finGeneralPar >> precCorrSnow;
  finGeneralPar.ignore(100,':'); finGeneralPar >> lapseDry;
  finGeneralPar.ignore(100,':'); finGeneralPar >> lapseWet; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> dayTempMem;
  finGeneralPar.ignore(100,':'); finGeneralPar >> lakeEpotPar;
  finGeneralPar.ignore(100,':'); finGeneralPar >> kLake; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> deltaLevel; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> nLake; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> initialSoilMoisture; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> initialUpperZone; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> initialLowerZone; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> initialLakeTemp; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> initialLakeLevel; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> initialSnow; 
  finGeneralPar.ignore(100,':'); finGeneralPar >> daySnowZero; 
  finGeneralPar.ignore(100, ':'); finGeneralPar >> height_wind_instrument;
  finGeneralPar.ignore(100, ':'); finGeneralPar >> height_humi_instrument;
  ParGeneralStore->SetSECONDS_TIMESTEP(secondsTimestep);
  ParGeneralStore->SetNUM_PREC_SERIES(numPrec);
  ParGeneralStore->SetNUM_TEMP_SERIES(numTemp);
  ParGeneralStore->SetPREC_GRAD_LOW(precGradLow);
  ParGeneralStore->SetPREC_GRAD_HIGH(precGradHigh);
  ParGeneralStore->SetGRAD_CHANGE_ALT(gradChangeAltitude);
  ParGeneralStore->SetPREC_CORR_RAIN(precCorrRain);
  ParGeneralStore->SetPREC_CORR_SNOW(precCorrSnow);
  ParGeneralStore->SetLAPSE_DRY(lapseDry);
  ParGeneralStore->SetLAPSE_WET(lapseWet);
  ParGeneralStore->SetDAY_TEMP_MEMORY(dayTempMem);
  ParGeneralStore->SetLAKE_EPOT_PAR(lakeEpotPar);
  ParGeneralStore->SetKLAKE(kLake);
  ParGeneralStore->SetDELTA_LEVEL(deltaLevel);
  ParGeneralStore->SetNLAKE(nLake);
  ParGeneralStore->SetINITIAL_SOIL_MOISTURE(initialSoilMoisture);
  ParGeneralStore->SetINITIAL_UPPER_ZONE(initialUpperZone); 
  ParGeneralStore->SetINITIAL_LOWER_ZONE(initialLowerZone);
  ParGeneralStore->SetINITIAL_LAKE_TEMP(initialLakeTemp);
  ParGeneralStore->SetINITIAL_LAKE_LEVEL(initialLakeLevel);
  ParGeneralStore->SetINITIAL_SNOW(initialSnow);
  ParGeneralStore->SetDAY_SNOW_ZERO(daySnowZero);
  ParGeneralStore->SetHeight_WIND_INSTRUMENT(height_wind_instrument);
  ParGeneralStore->SetHeight_HUMI_INSTRUMENT(height_humi_instrument);
  //  ParGeneralStore->SetNumStations(finGeneralPar, numPrec, numTemp); 
  finGeneralPar.close();
  fout << "Common parameters: \n";
  fout << ParGeneralStore->GetSECONDS_TIMESTEP() << endl;
  fout << ParGeneralStore->GetNUM_PREC_SERIES() << endl;
  fout << ParGeneralStore->GetNUM_TEMP_SERIES() << endl;
  fout << ParGeneralStore->GetPREC_GRAD_LOW() << endl;
  fout << ParGeneralStore->GetPREC_GRAD_HIGH() << endl;
  fout << ParGeneralStore->GetGRAD_CHANGE_ALT() << endl;
  fout << ParGeneralStore->GetPREC_CORR_RAIN() << endl;
  fout << ParGeneralStore->GetPREC_CORR_SNOW() << endl;
  fout << ParGeneralStore->GetLAPSE_DRY() << endl;
  fout << ParGeneralStore->GetLAPSE_WET() << endl;
  fout << ParGeneralStore->GetDAY_TEMP_MEMORY() << endl;
  fout << ParGeneralStore->GetLAKE_EPOT_PAR() << endl;
  fout << ParGeneralStore->GetKLAKE() << endl;
  fout << ParGeneralStore->GetDELTA_LEVEL() << endl;
  fout << ParGeneralStore->GetNLAKE() << endl;
  fout << ParGeneralStore->GetINITIAL_SOIL_MOISTURE() << endl;
  fout << ParGeneralStore->GetINITIAL_UPPER_ZONE() << endl;
  fout << ParGeneralStore->GetINITIAL_LOWER_ZONE() << endl;
  fout << ParGeneralStore->GetINITIAL_LAKE_TEMP() << endl;
  fout << ParGeneralStore->GetINITIAL_LAKE_LEVEL() << endl;
  fout << ParGeneralStore->GetINITIAL_SNOW() << endl;
  fout << ParGeneralStore->GetDAY_SNOW_ZERO() << endl;
  fout << ParGeneralStore->GetHeight_WIND_INSTRUMENT() << endl;
  fout << ParGeneralStore->GetHeight_HUMI_INSTRUMENT() << endl;
  /*  fout << " Prec.    "; 
      for (i=0; i<ParGeneralStore->GetNUM_PREC_SERIES(); i++) {
      fout << ParGeneralStore->GetSTATION_ALTITUDE(i) << "  ";
      fout << ParGeneralStore->GetSTATION_WEIGHT(i) << "  ";
      }
      fout << endl<< " Temp.    ";
      for (i=0; i<ParGeneralStore->GetNUM_TEMP_SERIES(); i++) {
      fout << ParGeneralStore->GetSTATION_ALTITUDE(ParGeneralStore->GetNUM_PREC_SERIES()+i) << "  ";
      fout << ParGeneralStore->GetSTATION_WEIGHT(ParGeneralStore->GetNUM_PREC_SERIES()+i) << "  ";
      }*/
  fout << endl;
}



