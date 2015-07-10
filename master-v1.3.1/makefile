DEBUG ?= no

SOURCE = CenResGrp.cpp CenterResidue.cpp Common.cpp Match.cpp ProteinStruct.cpp QueryStruct.cpp SearchOptions.cpp SearchResults.cpp TargetStruct.cpp
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCE))

MASTER = master
CREATE = createPDS
EXTRACT = extractPDB

MSLLIB = mslib
MSLINCL = mslib/src
CC = g++
FLAGS = -O3
LINKFLAGS = -static

ifeq ($(strip $(DEBUG)), yes)
	FLAGS = -O0 -g
endif

all: $(MASTER) $(CREATE) $(EXTRACT)

$(OBJECTS): %.o : %.cpp %.h
	$(CC) $(FLAGS) -c $< -o $@
$(MASTER): % : %.cpp $(OBJECTS)
	$(CC) $(FLAGS) $(OBJECTS) $< $(LINKFLAGS) -o $@
$(CREATE): % : %.cpp 
	$(CC) $(FLAGS) $< -I$(MSLINCL) -L$(MSLLIB) -lmsl $(LINKFLAGS) -o $@
$(EXTRACT): % : %.cpp
	$(CC) $(FLAGS) $< $(LINKFLAGS) -o $@

clean :
	-rm -f $(OBJECTS) $(MASTER) $(CREATE) $(EXTRACT)
