OBJSD = sourceDistHbv.o classDistHbv.o DateTime.o Date.o CTime.o utilities.o parameters.o vegetation.o snow.o glacierSurface.o lakeWaterBalance.o subSurfaceHbv.o
OBJSP = sourcePreHbv.o classPreHbv.o
C++C=g++

hbv: 	$(OBJSD)
	$(C++C) $(C++FLAGS) $(OBJSD) -o hbv -I/usr/local/include -lnetcdf

prehbv:	$(OBJSP)
	$(C++C) $(C++FLAGS) $(OBJSP) -o prehbv 

stationMask:	sourceStationMask.cpp 
		$(C++C) $(C++FLAGS) sourceStationMask.cpp -o stationMask

mask:	sourceMask.cpp 
	$(C++C) $(C++FLAGS) sourceMask.cpp -o mask

clean:
	rm -f $(OBJSD) $(OBJSP) *.o core
