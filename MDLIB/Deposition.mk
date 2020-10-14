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
objname := Deposition

#sorce directories
ifeq ($(origin MDLIBDIRS), undefined) 
MDLIBDIRS := MDLIB/sor/
endif
sor  := $(MDLIBDIRS)$(objname)/

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

#######################################################
mlist    :=  deposition_typedef_ctrlparam_2010 \
             deposition_common_2010            \
             md_typedef_projectile
             
nlist    :=  DEPOSITION_TypeDef_CtrlParam_10  \
             DEPOSITION_COMMON_10             \
             MD_TypeDef_Projectile  

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
	
$(tgt)DEPOSITION_TypeDef_CtrlParam_10.o : $(sor)DEPOSITION_TypeDef_CtrlParam_10.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@


$(tgt)DEPOSITION_COMMON_10.o : $(sor)DEPOSITION_COMMON_10.F90 \
                               $(tgt)DEPOSITION_TypeDef_CtrlParam_10.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

$(tgt)MD_TypeDef_Projectile.o : $(sor)MD_TypeDef_Projectile.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@


clean:
	-rm $(objects) $(libname) $(modules) 
