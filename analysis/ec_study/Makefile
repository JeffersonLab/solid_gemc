# (modified version of $(ROOTSYS)/test/Makefile)
#

#   Tue Jun 28 16:14:11 EDT 2016      imported from Remoll simulation analysis scripts to use in GEMC based analysis scripts : Rakitha S Beminiwattha, rakithab@jlab.org

ObjSuf        = o
SrcSuf        = cc
ExeSuf        =
DllSuf        = so
EVENTLIB      = $(EVENTO)
OutPutOpt     = -o  

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)   -lNew  -lMinuit
ROOTGLIBS     = $(shell root-config --glibs)


MYOS := $(subst -,,$(shell uname))

CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

EXE := ecal_trigger_analysis2_withBkg_PVDISeAll_bkgcheck_2023_Eend 
# Linux with egcs	
#ifeq ($(MYOS),Linux)
	CXX           = g++ 
	CXXFLAGS      = -O2 -Wall -fPIC 
#-I$(G0ANALYSIS)/Analysis/include
	LD            = g++
	LDFLAGS       = -O
	SOFLAGS       = -shared

	CXXFLAGS     += $(ROOTCFLAGS)
	CXXFLAGS     += -I./include
	LIBS          = $(ROOTLIBS)
#-L$(G0ANALYSIS)/lib -lG0 -lmyevio -L$(CODA_LIB) -let
	GLIBS         = $(ROOTGLIBS)
#endif

all:   $(EXE) 
	@echo "---------------------> $(EXE) done <---------------------- "

ecal_trigger_analysis2_withBkg_PVDISeAll_bkgcheck_2023_Eend: ecal_trigger_analysis2_withBkg_PVDISeAll_bkgcheck_2023_Eend.o
	$(CXX) $(CXXFLAGS)  $(LIBS) $(GLIBS) -o $@ $<
clean:
	@rm -rf *.o
	@rm -rf *~
	@rm -rf  core

cleanall:
	@rm -f $(OBJS) $(EXE) *~  core

.PHONY : clean cleanall all
