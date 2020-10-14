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
objname := EM_TB_WangJun_Ti-HE_2007

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
mlist    :=  em_tb_forcetable_wangjun_ti_he_2007  \
             pot_lj_he_ti

nlist    :=  EM_TB_ForcetTable_WangJun_Ti_He_2007 \
             LJ_He_Ti

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)EM_TB_ForcetTable_WangJun_Ti_He_2007.o : $(sor)EM_TB_ForcetTable_WangJun_Ti_He_2007.F90  \
             $(tgt)TB_CR_TiTi.o                                          \
             $(tgt)LJ_He_Ti.o                                          


	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)LJ_He_Ti.o : $(sor)LJ_He_Ti.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@


######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
