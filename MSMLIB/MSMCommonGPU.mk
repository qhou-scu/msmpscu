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
objname := CommonGPU

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
mlist    :=  msm_memallocator_gpu \
             msm_multigpu_basic   \
             msm_multigpu_datlist
           
             
nlist    :=  MSM_MemAllocator_GPU  \
             MSM_MultiGPU_Basic    \
             MSM_MultiGPU_DatList  \
			 MSM_MultiGPU_DatList_DFMat \
			 MSM_MultiGPU_DatList_DFVec \
			 MSM_MultiGPU_DatList_IVec  \
			 MSM_MultiGPU_DatList_IMat

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)MSM_MemAllocator_GPU.o : $(sor)MSM_MemAllocator_GPU.F90 \
                               $(MSMLIBDIRS)MSM_MemAllocate_Temp.F90  
	$(comp) -c $(oflags) -I$(incdir) -I$(MSMLIBDIRS) -module $(tgt) $< -o $@

$(tgt)MSM_MultiGPU_Basic.o : $(sor)MSM_MultiGPU_Basic.F90 
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MSM_MultiGPU_DatList_DFMat.o : $(sor)MSM_MultiGPU_DatList_DFMat.F90 \
                                     $(MSMLIBDIRS)MSM_MultiGPU_DatList_Header.F90 \
							                  		 $(MSMLIBDIRS)MSM_MultiGPU_DatList_Body.F90 \
									                   $(tgt)MSM_MultiGPU_Basic.o
	$(comp) -c $(oflags) -I$(incdir) -I$(MSMLIBDIRS) -module $(tgt) $< -o $@

$(tgt)MSM_MultiGPU_DatList_IMat.o :  $(sor)MSM_MultiGPU_DatList_IMat.F90 \
                                     $(MSMLIBDIRS)MSM_MultiGPU_DatList_Header.F90 \
									                   $(MSMLIBDIRS)MSM_MultiGPU_DatList_Body.F90 \
									                   $(tgt)MSM_MultiGPU_Basic.o
	$(comp) -c $(oflags) -I$(incdir) -I$(MSMLIBDIRS) -module $(tgt) $< -o $@

$(tgt)MSM_MultiGPU_DatList_DFVec.o : $(sor)MSM_MultiGPU_DatList_DFVec.F90 \
                                     $(MSMLIBDIRS)MSM_MultiGPU_DatList_Header.F90 \
									                   $(MSMLIBDIRS)MSM_MultiGPU_DatList_Body.F90 \
								                  	 $(tgt)MSM_MultiGPU_Basic.o
	$(comp) -c $(oflags) -I$(incdir) -I$(MSMLIBDIRS) -module $(tgt) $< -o $@

$(tgt)MSM_MultiGPU_DatList_IVec.o : $(sor)MSM_MultiGPU_DatList_IVec.F90 \
                                    $(MSMLIBDIRS)MSM_MultiGPU_DatList_Header.F90 \
							                   		$(MSMLIBDIRS)MSM_MultiGPU_DatList_Body.F90 \
									                  $(tgt)MSM_MultiGPU_Basic.o
	$(comp) -c $(oflags) -I$(incdir) -I$(MSMLIBDIRS) -module $(tgt) $< -o $@

$(tgt)MSM_MultiGPU_DatList.o : $(sor)MSM_MultiGPU_DatList.F90 \
                               $(tgt)MSM_MultiGPU_Basic.o     \
							                 $(tgt)MSM_MultiGPU_DatList_DFMat.o \
						               	   $(tgt)MSM_MultiGPU_DatList_IMat.o  \
							                 $(tgt)MSM_MultiGPU_DatList_DFVec.o \
							                 $(tgt)MSM_MultiGPU_DatList_IVec.o
	$(comp) -c $(oflags) -I$(incdir) -I$(MSMLIBDIRS) -module $(tgt) $< -o $@


######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
