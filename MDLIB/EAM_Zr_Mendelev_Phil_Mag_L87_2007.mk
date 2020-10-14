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
objname := EAM_Zr_Mendelev_Phil_Mag_L87_2007

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
mlist    := eam_zrzr_mendelev_phil_mag_l87_2007    \
            eam_zr_forcetable_mendelev_phil_mag_l87_2007

nlist    := EAM_ZrZr_Mendelev_Phil_Mag_L87_2007 \
            EAM_Zr_ForceTable_Mendelev_Phil_Mag_L87_2007
          

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)EAM_Zr_ForceTable_Mendelev_Phil_Mag_L87_2007.o : \
            $(sor)EAM_Zr_ForceTable_Mendelev_Phil_Mag_L87_2007.F90  \
            $(tgt)EAM_ZrZr_Mendelev_Phil_Mag_L87_2007.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)EAM_ZrZr_Mendelev_Phil_Mag_L87_2007.o : $(sor)EAM_ZrZr_Mendelev_Phil_Mag_L87_2007.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
