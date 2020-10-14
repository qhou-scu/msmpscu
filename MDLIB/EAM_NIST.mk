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
objname := EAM_NIST

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
mlist    :=  nist_forcetable                    \
             filedatas_func_lspt                \
             filedatas_func_moldy               \
             filedatas_func_setfl               \
             my_miniUtilities_qmj               

nlist    :=  NIST_ForceTable           \
             Filedatas_Func_Lspt       \
             Filedatas_Func_Moldy      \
             Filedatas_Func_Setfl      \
             My_MiniUtilities_QMJ  

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)NIST_ForceTable.o : $(sor)NIST_ForceTable.F90  \
             $(tgt)Filedatas_Func_Lspt.o              \
             $(tgt)Filedatas_Func_Moldy.o             \
             $(tgt)Filedatas_Func_Setfl.o             \
             $(tgt)My_MiniUtilities_QMJ.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)Filedatas_Func_Lspt.o : $(sor)Filedatas_Func_Lspt.F90  $(tgt)My_MiniUtilities_QMJ.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)Filedatas_Func_Moldy.o : $(sor)Filedatas_Func_Moldy.F90  $(tgt)My_MiniUtilities_QMJ.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)Filedatas_Func_Setfl.o : $(sor)Filedatas_Func_Setfl.F90  $(tgt)My_MiniUtilities_QMJ.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)My_MiniUtilities_QMJ.o : $(sor)My_MiniUtilities_QMJ.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@


######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
