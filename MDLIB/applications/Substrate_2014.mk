#compiler
ifeq ($(origin comp), undefined) 
comp       = pgfortran
endif

ifeq ($(origin oflags), undefined)
oflags := -fast -Mvect=sse -Minline -Mconcur -Minform=warn -Minfo=accel,inline,intensity,loop,mp,opt,par,vect 
endif
oflags += -Mcuda=fastmath

# include directeries
LIBDIRD    := $(WORKSPACE)/LIB/$(ConfigName)/
incdir     := $(WORKSPACE)/LIB/$(ConfigName)/

##########################################################
#sorce directories
ifeq ($(origin MDLIBDIRS), undefined) 
MDLIBDIRS := ./sor
endif

#objective name
objname := Substrate_2014

#sorce directories
sor  := $(MDLIBDIRS)/$(objname)/

#target directories
ifeq ($(ConfigName), Debug)
tgt  := $(WORKSPACE)/applications/$(ConfigName)/
else
tgt  := $(WORKSPACE)/applications/
endif

#executable name
exename  := $(tgt)$(objname).exe
##########################################################
mlist    :=  substrate_2014
nlist    :=  SUBSTRATE_2014_Main_GPU      
objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
CUFfiles := $(foreach n, $(nlist), $(sor)$(n).CUF)

libs     := $(foreach n, $(potlibnames), $(LIBDIRD)libPOT_$(n).$(LIB_EXT))
libs     += $(foreach n, $(mdlibnames),  $(LIBDIRD)libMD_$(n).$(LIB_EXT))
libs     += $(foreach n, $(libnames),    $(LIBDIRD)lib$(n).$(LIB_EXT))

liblist  := $(foreach n, $(utilibnames), -L$(LIBDIRD) -lMD_$(n))
liblist  += $(foreach n, $(potlibnames), -L$(LIBDIRD) -lPOT_$(n))
liblist  += $(foreach n, $(mdlibnames),  -L$(LIBDIRD)  -lMD_$(n))
liblist  += $(foreach n, $(libnames),    -L$(LIBDIRD)  -l$(n))

#######################################################
$(exename) : $(tgt) $(tgt)substrate_2014.o Substrate_2014.mk $(libs)
	$(comp)  $(oflags) $(tgt)substrate_2014.o  $(liblist) -o $(exename)
	rm $(tgt)substrate_2014.o

$(tgt) : 
	-mkdir $(WORKSPACE)/applications/
	-mkdir $(tgt)

$(tgt)substrate_2014.o : $(sor)SUBSTRATE_2014_Main_GPU.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@ 


clean:
	-rm $(tgt)substrate_2014.o $(tgt)substrate_2014.mod $(exename) 
	     
