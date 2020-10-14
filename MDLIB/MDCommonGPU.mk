#compiler
ifeq ($(origin comp), undefined) 
comp       = pgfortran
endif

ifeq ($(origin oflags), undefined)
oflags := -fast -Mvect=sse,simd -Minline -Mconcur -Minform=warn -Minfo=accel,inline,intensity,loop,mp,opt,par,vect -Mcuda=fastmath,cuda5.5,cc35
endif

#include dir
incdir := $(LIBDIRD)

##########################################################
#sorce dir name
objname := CommonGPU

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

##########################################################
mlist    :=  md_activeregion_gpu             \
             md_bp_coupling_gpu              \
             md_cgscheme_gpu                 \
             md_diffscheme_gpu               \
             md_eam_force_table_gpu          \
             md_fs_force_table_gpu           \
             md_forceclass_registergpu       \
             md_globle_variables_gpu         \
             md_lbfgsscheme_gpu              \
             md_multigpu_basic               \
             md_neighborslist_gpu            \
             md_simboxarray_gpu              \
             md_steepestscheme_gpu           \
             md_utilities_gpu

nlist    :=  MD_ActiveRegion_GPU              \
             MD_BP_Coupling_GPU               \
             MD_CGScheme_GPU                  \
             MD_DiffScheme_GPU                \
             MD_EAM_ForceTable_GPU            \
             MD_FS_ForceTable_GPU             \
             MD_ForceClass_Register_GPU       \
             MD_Globle_Variables_GPU          \
             MD_LBFGSScheme_GPU               \
             MD_MultiGPU_Basic                \
             MD_NeighborsList_GPU             \
             MD_SimBoxArray_GPU               \
             MD_SteepestScheme_GPU            \
             MD_Utilities_GPU

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
	
$(tgt)MD_ActiveRegion_GPU.o : $(sor)MD_ActiveRegion_GPU.F90 \
                                            $(tgt)MD_Globle_Variables_GPU.o \
                                            $(tgt)MD_NeighborsList_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_BP_Coupling_GPU.o : $(sor)MD_BP_Coupling_GPU.F90    \
                                $(tgt)MD_Globle_Variables_GPU.o 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_CGScheme_GPU.o : $(sor)MD_CGScheme_GPU.F90 \
                                $(tgt)MD_Globle_Variables_GPU.o \
                                $(tgt)MD_ForceClass_Register_GPU.o \
                                $(tgt)MD_MultiGPU_Basic.o   
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_DiffScheme_GPU.o : $(sor)MD_DiffScheme_GPU.F90 \
                                $(tgt)MD_Globle_Variables_GPU.o    
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_EAM_ForceTable_GPU.o : $(sor)MD_EAM_ForceTable_GPU.F90 \
                                            $(tgt)MD_Globle_Variables_GPU.o \
                                            $(tgt)MD_NeighborsList_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_FS_ForceTable_GPU.o : $(sor)MD_FS_ForceTable_GPU.F90 \
                                            $(tgt)MD_Globle_Variables_GPU.o \
                                            $(tgt)MD_NeighborsList_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_ForceClass_Register_GPU.o : $(sor)MD_ForceClass_Register_GPU.F90 \
                                        $(tgt)MD_EAM_ForceTable_GPU.o     \
                                        $(tgt)MD_FS_ForceTable_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_Globle_Variables_GPU.o : $(sor)MD_Globle_Variables_GPU.F90 \
                                  $(tgt)MD_MultiGPU_Basic.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_LBFGSScheme_GPU.o : $(sor)MD_LBFGSScheme_GPU.F90 \
                                       $(tgt)MD_Globle_Variables_GPU.o  \
                                       $(tgt)MD_ForceClass_Register_GPU.o \
                                       $(tgt)MD_MultiGPU_Basic.o     
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_SteepestScheme_GPU.o : $(sor)MD_SteepestScheme_GPU.F90 \
                                       $(tgt)MD_Globle_Variables_GPU.o  \
                                       $(tgt)MD_ForceClass_Register_GPU.o \
                                       $(tgt)MD_MultiGPU_Basic.o     
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
	
$(tgt)MD_MultiGPU_Basic.o : $(sor)MD_MultiGPU_Basic.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_NeighborsList_GPU.o : $(sor)MD_NeighborsList_GPU.F90 \
                                    $(tgt)MD_Globle_Variables_GPU.o    
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_SimBoxArray_GPU.o : $(sor)MD_SimBoxArray_GPU.F90 \
                                       $(tgt)MD_Globle_Variables_GPU.o   
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_Utilities_GPU.o : $(sor)MD_Utilities_GPU.F90 \
                              $(tgt)MD_Globle_Variables_GPU.o \
                              $(tgt)MD_NeighborsList_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

clean:
	-rm $(objects) $(libname) $(modules) 
