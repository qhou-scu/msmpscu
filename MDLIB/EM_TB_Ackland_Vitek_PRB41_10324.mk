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
objname := EM_TB_Ackland_Vitek_PRB41_10324

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
libname  :=libPOT_$(objname).$(LIB_EXT)

#######################################################
mlist    :=  em_tb_auau_ackland_vitek_prb41_10324  \
             em_tb_cuau_ackland_vitek_prb41_10324  \
             em_tb_cucu_ackland_vitek_prb41_10324  \
             em_tb_forcetable_ackland_vitek_prb41_10324

nlist    := EM_TB_AuAu_Ackland_Vitek_PRB41_10324 \
            EM_TB_CuAu_Ackland_Vitek_PRB41_10324 \
            EM_TB_CuCu_Ackland_Vitek_PRB41_10324 \
            EM_TB_ForceTable_Ackland_Vitek_PRB41_10324            

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)EM_TB_ForceTable_Ackland_Vitek_PRB41_10324.o : $(sor)EM_TB_ForceTable_Ackland_Vitek_PRB41_10324.F90  \
             $(tgt)EM_TB_AuAu_Ackland_Vitek_PRB41_10324.o \
             $(tgt)EM_TB_CuAu_Ackland_Vitek_PRB41_10324.o \
             $(tgt)EM_TB_CuCu_Ackland_Vitek_PRB41_10324.o           
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)EM_TB_AuAu_Ackland_Vitek_PRB41_10324.o : $(sor)EM_TB_AuAu_Ackland_Vitek_PRB41_10324.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
$(tgt)EM_TB_CuAu_Ackland_Vitek_PRB41_10324.o : $(sor)EM_TB_CuAu_Ackland_Vitek_PRB41_10324.F90 \
             $(tgt)EM_TB_AuAu_Ackland_Vitek_PRB41_10324.o \
             $(tgt)EM_TB_CuCu_Ackland_Vitek_PRB41_10324.o           
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
$(tgt)EM_TB_CuCu_Ackland_Vitek_PRB41_10324.o : $(sor)EM_TB_CuCu_Ackland_Vitek_PRB41_10324.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
