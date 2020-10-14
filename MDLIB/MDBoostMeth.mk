#compiler
ifeq ($(origin comp), undefined) 
comp       = pgfortran
endif

ifeq ($(origin oflags), undefined)
oflags := -fast -Mvect=sse,simd -Minline -Mconcur -Minform=warn -Minfo=accel,inline,intensity,loop,mp,opt,par,vect -Mcuda=fastmath,cuda8.0,cc2x
endif

#include dir
incdir := $(LIBDIRD)

##########################################################
#sorce dir name
objname := BoostMeths

#sorce directories
ifeq ($(origin MDLIBDIRS), undefined) 
MDLIBDIRS := MDLIB/sor/
endif
sor     := $(MDLIBDIRS)$(objname)/

#target directories
ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

ifeq ($(origin LIBDIRD), undefined) 
LIBDIRD := $(WORKSPACE)/LIB/$(ConfigName)
endif
tgt  := $(LIBDIRD)

#target lib name
libname  := libMD_$(objname).$(LIB_EXT)

##########################################################
mlist    :=  md_springboost_gpu  \
             md_boostmeth_gpu
 

nlist    :=  MD_SpringBoost_GPU  \
             MD_BoostMeth_GPU
                
objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
CUFfiles := $(foreach n, $(nlist), $(sor)$(n).CUF)

#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
######### for EPC
$(tgt)MD_BoostMeth_GPU.o : $(sor)MD_BoostMeth_GPU.F90
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_SpringBoost_GPU.o : $(sor)MD_SpringBoost_GPU.F90
	$(comp) -c $(oflags) -I$(incdir) -I$(MSMLIBDIRS) -module $(tgt) $< -o $@
	

clean:
	-rm $(objects) $(libname) $(modules) 
