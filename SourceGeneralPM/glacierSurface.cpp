#include "classDistHbv.h"

void GlacierSurface::WaterBalance(int timeStep, DateTime datetime)
{
  //  cout << "start GlacierSurfaceWaterBalance " << endl;
  precipitation = GetInputElement()->GetInput(0);
  temp = GetInputElement()->GetInput(1);
  /*  cout << "    timeStep " << timeStep << "       Glacier precipitation " << precipitation;
      cout << "    Temperature " << temp << endl;*/

  /*  Glacier ice melt  */
//  if (temp > landSurfacePar->GetMELT_TEMP() && GetLandScapeElement()->GetGlacier()->GetSnowStore() == 0.0)
  if (temp > landSurfacePar->GetMELT_TEMP())
    iceMelt = landSurfacePar->GetICE_MELT_RATE() * landSurfacePar->GetSNOW_MELT_RATE() * (temp - landSurfacePar->GetMELT_TEMP());
  else
    iceMelt = 0.0;
  //  cout << temp << "  " << landSurfacePar->GetMELT_TEMP() << "  " << landSurfacePar->GetICE_MELT_RATE() << "  " << GetLandScapeElement()->GetSnowStore() << endl;
  //  cout << "  end GlacierSurfaceWaterBalance " << endl;
}

