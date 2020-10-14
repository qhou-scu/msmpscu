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
objname := Appshell

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
mlist    :=  md_forcelib_factory_cpu             \
             md_forcelib_factory_gpu             \
             md_forcelib_register                \
             md_methodclass_factory_gpu          \
             md_method_art_gpu                   \
             md_method_bst_gpu                   \
             md_method_genericmd_gpu             \
             md_method_neb_gpu                   \
             md_method_parrep_gpu                \
             md_method_tad_gpu                   \
             md_simboxarray_appshell_16_gpu      \
             md_simboxarray_toolshell_14_cpu     \
             md_simboxarray_toolshell_14_gpu     
             
nlist    :=  MD_ForceLib_Factory_CPU           \
             MD_ForceLib_Factory_GPU           \
             MD_ForceLib_Register              \
             MD_MethodClass_Factory_GPU        \
             MD_Method_ART_GPU                 \
             MD_Method_BST_GPU                 \
             MD_Method_GenericMD_GPU           \
             MD_Method_NEB_GPU                 \
             MD_Method_ParRep_GPU              \
             MD_Method_TAD_GPU                 \
             MD_SimBoxArray_AppShell_16_GPU    \
             MD_SimBoxArray_ToolShell_14_CPU   \
             MD_SimBoxArray_ToolShell_14_GPU   

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)MD_ForceLib_Factory_CPU.o : $(sor)MD_ForceLib_Factory_CPU.F90 \
                                  $(tgt)MD_ForceLib_Register.o   
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

$(tgt)MD_ForceLib_Factory_GPU.o : $(sor)MD_ForceLib_Factory_GPU.F90 \
                                  $(tgt)MD_ForceLib_Register.o   
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

$(tgt)MD_ForceLib_Register.o : $(sor)MD_ForceLib_Register.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

$(tgt)MD_MethodClass_Factory_GPU.o : $(sor)MD_MethodClass_Factory_GPU.F90    \
                                     $(tgt)MD_Method_ART_GPU.o               \
                                     $(tgt)MD_Method_BST_GPU.o               \
                                     $(tgt)MD_Method_GenericMD_GPU.o         \
                                     $(tgt)MD_Method_NEB_GPU.o               \
                                     $(tgt)MD_Method_ParRep_GPU.o            \
                                     $(tgt)MD_Method_TAD_GPU.o          
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

$(tgt)MD_Method_ART_GPU.o : $(sor)MD_Method_ART_GPU.F90                     \
                            $(tgt)MD_ForceLib_Factory_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@
	
$(tgt)MD_Method_BST_GPU.o : $(sor)MD_Method_BST_GPU.F90                     \
                                    $(tgt)MD_ForceLib_Factory_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@
	

$(tgt)MD_Method_GenericMD_GPU.o : $(sor)MD_Method_GenericMD_GPU.F90         \
                                  $(tgt)MD_ForceLib_Factory_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

$(tgt)MD_Method_NEB_GPU.o : $(sor)MD_Method_NEB_GPU.F90                     \
                            $(tgt)MD_ForceLib_Factory_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

$(tgt)MD_Method_ParRep_GPU.o : $(sor)MD_Method_ParRep_GPU.F90                \
                               $(tgt)MD_ForceLib_Factory_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

$(tgt)MD_Method_TAD_GPU.o : $(sor)MD_Method_TAD_GPU.F90                      \
                            $(tgt)MD_ForceLib_Factory_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

$(tgt)MD_SimBoxArray_AppShell_16_GPU.o : $(sor)MD_SimBoxArray_AppShell_16_GPU.F90     \
                                           $(tgt)MD_ForceLib_Factory_GPU.o            \
                                           $(tgt)MD_MethodClass_Factory_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

$(tgt)MD_SimBoxArray_ToolShell_14_CPU.o : $(sor)MD_SimBoxArray_ToolShell_14_CPU.F90   \
                                           $(tgt)MD_ForceLib_Factory_CPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@

$(tgt)MD_SimBoxArray_ToolShell_14_GPU.o : $(sor)MD_SimBoxArray_ToolShell_14_GPU.F90   \
                                           $(tgt)MD_ForceLib_Factory_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(deplibs)) -module $(tgt) $< -o $@



######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
