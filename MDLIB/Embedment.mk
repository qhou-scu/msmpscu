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
objname := Embedment

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
mlist    :=  embedment_common_2010   \
             embedment_typedef_ctrlparam_2010 \
             embedment_Sandwich
             
nlist    :=  EMBEDMENT_COMMON_2010              \
             EMBEDMENT_TypeDef_CtrlParam_10     \
             EMBEDMENT_Sandwich

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
	
$(tgt)EMBEDMENT_TypeDef_CtrlParam_10.o : $(sor)EMBEDMENT_TypeDef_CtrlParam_10.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@


$(tgt)EMBEDMENT_COMMON_2010.o : $(sor)EMBEDMENT_COMMON_2010.F90 \
                                $(tgt)EMBEDMENT_TypeDef_CtrlParam_10.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

$(tgt)EMBEDMENT_Sandwich.o : $(sor)EMBEDMENT_Sandwich.F90 \
                             $(tgt)EMBEDMENT_TypeDef_CtrlParam_10.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

clean:
	-rm $(objects) $(libname) $(modules) 
