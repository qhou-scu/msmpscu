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
objname := EM_TB_Cleri_Rosato_PRB48_22

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
mlist    :=  em_tb_forcetable_cr_2014             \
             pot_tb_cr_alal                       \
             pot_tb_cr_auau                       \
             pot_tb_cr_titi

nlist    :=  EM_TB_ForcetTable_CR_2014           \
             TB_CR_AlAl                          \
             TB_CR_AuAu                          \
             TB_CR_TiTi

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)EM_TB_ForcetTable_CR_2014.o : $(sor)EM_TB_ForcetTable_CR_2014.F90  \
             $(tgt)TB_CR_AlAl.o                                          \
             $(tgt)TB_CR_AuAu.o                                          \
             $(tgt)TB_CR_TiTi.o                                          

	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)TB_CR_AlAl.o : $(sor)TB_CR_AlAl.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)TB_CR_AuAu.o : $(sor)TB_CR_AuAu.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)TB_CR_TiTi.o : $(sor)TB_CR_TiTi.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@


######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
