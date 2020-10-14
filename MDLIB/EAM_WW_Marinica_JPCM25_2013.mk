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
objname := EAM_WW_Marinica_JPCM25_2013

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
mlist    := eam2_ww_marinica_jpcm25_2013    \
            eam3_ww_marinica_jpcm25_2013    \
            eam4_ww_marinica_jpcm25_2013    \
            eam_ww_forcetable_marinica_jpcm25_2013

nlist    := EAM2_WW_Marinica_JPCM25_2013 \
            EAM3_WW_Marinica_JPCM25_2013 \
            EAM4_WW_Marinica_JPCM25_2013 \
            EAM_ForceTable_Marinica_JPCM25_2013
          

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)EAM_ForceTable_Marinica_JPCM25_2013.o : $(sor)EAM_ForceTable_Marinica_JPCM25_2013.F90  \
                                              $(tgt)EAM2_WW_Marinica_JPCM25_2013.o           \
                                              $(tgt)EAM3_WW_Marinica_JPCM25_2013.o           \
                                              $(tgt)EAM4_WW_Marinica_JPCM25_2013.o           
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)EAM2_WW_Marinica_JPCM25_2013.o : $(sor)EAM2_WW_Marinica_JPCM25_2013.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)EAM3_WW_Marinica_JPCM25_2013.o : $(sor)EAM3_WW_Marinica_JPCM25_2013.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)EAM4_WW_Marinica_JPCM25_2013.o : $(sor)EAM4_WW_Marinica_JPCM25_2013.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
