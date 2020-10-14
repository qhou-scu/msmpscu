#compiler
ifeq ($(origin comp), undefined) 
comp       = pgfortran
endif

ifeq ($(origin oflags), undefined)
oflags := -fast
endif

#include dir
incdir := $(LIBDIRD)

##########################################################
#sorce dir name
objname := EAM_WHeH_Bonny_JPCM26_2014

#sorce directories
ifeq ($(origin POTLIBDIRS), undefined) 
POTLIBDIRS := MDLIB/sor/Potentials/
endif
sor  := $(POTLIBDIRS)$(objname)/

#target directories
ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

ifeq ($(origin LIBDIRD), undefined) 
LIBDIRD := $(WORKSPACE)/LIB/$(ConfigName)
endif
tgt  := $(LIBDIRD)

#target lib name
libname  := libPOT_$(objname).$(LIB_EXT)

#######################################################
mlist    := eam1_wheh_bonny_jpcm26_2014    \
            eam_wheh_forcetable_bonny_jpcm26_2014

nlist    := EAM1_WHeH_Bonny_JPCM26_2014 \
            EAM_ForceTable_Bonny_JPCM26_2014
          

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)EAM_ForceTable_Bonny_JPCM26_2014.o : $(sor)EAM_ForceTable_Bonny_JPCM26_2014.F90  \
                                              $(tgt)EAM1_WHeH_Bonny_JPCM26_2014.o           
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)EAM1_WHeH_Bonny_JPCM26_2014.o : $(sor)EAM1_WHeH_Bonny_JPCM26_2014.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
