#compiler
ifeq ($(origin comp), undefined) 
comp       = pgfortran
endif

ifeq ($(origin oflags), undefined)
oflags := -fast -Mvect=sse,simd -Minline -Mconcur -Minform=warn -Minfo=accel,inline,intensity,loop,mp,opt,par,vect -Mcuda=fastmath,cuda8.0,cc2x 
endif

# include directeries
LIBDIRD    := $(WORKSPACE)/LIB/$(ConfigName)/
incdir     := $(WORKSPACE)/LIB/$(ConfigName)/

##########################################################
#sorce directories
ifeq ($(origin MDLIBDIRS), undefined) 
MDLIBDIRS := ./sor
endif

#objective name
objname := CalForceTest

#sorce directories
sor  := $(MDLIBDIRS)/

#target directories

ifeq ($(ConfigName), Debug)
tgt  := $(WORKSPACE)/applications/$(ConfigName)/
else
tgt  := $(WORKSPACE)/applications/
endif

#executable name
exename  := $(tgt)$(objname).exe
##########################################################
mlist    :=  calforcetest
nlist    :=  CalForceTest      
objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
CUFfiles := $(foreach n, $(nlist), $(sor)$(n).CUF)

libnames    := MiniUtilities  RandGenerators
mdlibnames  := Appshell CommonGPU  Common
potlibnames := EM_TB_WangJun_W-HE_2010\
               EM_TB_Ackland_Vitek_PRB41_10324\
               EM_TB_AWB_PHIL_MAG_A71_1995_553\
               EM_TB_Cleri_Rosato_PRB48_22 

libs     := $(foreach n, $(potlibnames), $(LIBDIRD)libPOT_$(n).$(LIB_EXT))
libs     += $(foreach n, $(mdlibnames),  $(LIBDIRD)libMD_$(n).$(LIB_EXT))
libs     += $(foreach n, $(libnames),    $(LIBDIRD)lib$(n).$(LIB_EXT))

liblist  := $(foreach n, $(potlibnames), -L$(LIBDIRD) -lPOT_$(n))
liblist  += $(foreach n, $(mdlibnames),  -L$(LIBDIRD)  -lMD_$(n))
liblist  += $(foreach n, $(libnames),    -L$(LIBDIRD)  -l$(n))
#######################################################
$(exename) : $(tgt) CalForceTest.mk $(tgt)CalForceTest.o $(libs) 
	$(comp)  $(oflags) $(tgt)CalForceTest.o $(liblist) -o $(exename)
	rm $(tgt)CalForceTest.o

$(tgt) : 
	-mkdir $(WORKSPACE)/applications/
	-mkdir $(tgt)

$(tgt)CalForceTest.o : $(sor)CalForceTest.f90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) $< -o $@ 


clean:
	-rm $(tgt)CalForceTest.o $(exename) 
