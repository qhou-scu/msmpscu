#compiler
ifeq ($(origin comp), undefined) 
comp       = pgfortran
endif

ifeq ($(origin oflags), undefined)
oflags := -fast
endif

ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

#include dir
ifeq ($(origin LIBDIRD), undefined)
LIBDIRD := $(WORKSPACE)/LIB/$(ConfigName)
endif
incdir := $(LIBDIRD)

##########################################################
#sorce dir name
objname := Common

#sorce directories
ifeq ($(origin MSMLIBDIRS), undefined) 
MSMLIBDIRS := $(MSMPSCUSOR)/MSMLIB/sor/
endif
sor  := $(MSMLIBDIRS)$(objname)/

#target directories
ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

ifeq ($(origin LIBDIRD), undefined) 
LIBDIRD := $(WORKSPACE)/LIB/$(ConfigName)
endif
tgt  := $(LIBDIRD)

#target lib name
libname  := libMSM_$(objname).$(LIB_EXT)

#######################################################
mlist    :=  msm_constants                   \
             msm_memallocator                \
             msm_typedf_recordstamp          \
             msm_typedef_datapad             \
             msm_typedef_InputPaser
           
             
nlist    :=  MSM_Const                       \
             MSM_MemAllocator                \
             MSM_TypeDef_RecordStamp         \
             MSM_TypeDef_DataPad             \
             MSM_TypeDef_InputPaser
             
inclist  :=  MSM_MemAllocate_Temp.F90       

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)

#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects)
	mv $(libname) $(tgt)
	
$(tgt)MSM_Const.o : $(sor)MSM_Const.F90 
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MSM_MemAllocator.o : $(sor)MSM_MemAllocator.F90              \
                                $(tgt)MSM_Const.o                  \
                                $(MSMLIBDIRS)MSM_MemAllocate_Temp.F90               
	$(comp) -c $(oflags) -I$(incdir) -I$(MSMLIBDIRS) -module $(tgt) $< -o $@
	
$(tgt)MSM_TypeDef_RecordStamp.o : $(sor)MSM_TypeDef_RecordStamp.F90    \
                                  $(tgt)MSM_Const.o                  
	$(comp) -c $(oflags) -I$(incdir) -I$(MSMLIBDIRS) -module $(tgt) $< -o $@
	
	
$(tgt)MSM_TypeDef_DataPad.o : $(sor)MSM_TypeDef_DataPad.F90        \
                                $(tgt)MSM_Const.o             
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MSM_TypeDef_InputPaser.o : $(sor)MSM_TypeDef_InputPaser.F90   \
                                   $(tgt)MSM_Const.o                  
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@



######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
