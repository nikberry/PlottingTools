LIBS=`root-config --libs`
CFLAGS=`root-config --cflags`
CC=g++

#set compliler options
#	-g = debugging
#	-O# = optimisation
COPT=-g

default:
	$(CC) $(COPT) src/TTbarXSectionPlottingTools.cpp src/GlobalVariables.cpp src/Samples/AllSamples.cpp src/Samples/Sample.cpp src/Variables/Variable.cpp src/Objects/CutFlow.cpp src/Objects/multiFit.cpp src/Objects/Jets.cpp  src/Objects/MET.cpp src/Objects/Muon.cpp src/Objects/Object.cpp  src/Objects/TdrStyle.cpp -o TTbarXSectionPlottingTools $(LIBS) $(CFLAGS)
	
clean:
	rm -rf *o TTbarXSectionPlottingTools
