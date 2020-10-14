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
objname := NeighboresListTest

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
mlist    :=  neighboreslisttest 
nlist    :=  NeighboresListTest      
objects  :=   $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
CUFfiles := $(foreach n, $(nlist), $(sor)$(n).CUF)

libs     := $(foreach n, $(potlibnames), $(LIBDIRD)libPOT_$(n).$(LIB_EXT))
libs     += $(foreach n, $(mdlibnames),  $(LIBDIRD)libMD_$(n).$(LIB_EXT))
libs     += $(foreach n, $(libnames),    $(LIBDIRD)lib$(n).$(LIB_EXT))

liblist  := $(foreach n, $(potlibnames), -L$(LIBDIRD) -lPOT_$(n))
liblist  += $(foreach n, $(mdlibnames),  -L$(LIBDIRD)  -lMD_$(n))
liblist  += $(foreach n, $(libnames),    -L$(LIBDIRD)  -l$(n))
#######################################################
$(exename) : $(tgt) NeighboresListTest.mk $(tgt)NeighboresListTest.o $(libs) 
	$(comp)  $(oflags) $(tgt)NeighboresListTest.o $(liblist) -o $(exename)
	rm $(tgt)NeighboresListTest.o

$(tgt) : 
	-mkdir $(WORKSPACE)/applications/
	-mkdir $(tgt)

$(tgt)NeighboresListTest.o : $(sor)NeighboresListTest.f90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) $< -o $@


clean:
	-rm $(tgt)NeighboresListTest.o $(exename) 
