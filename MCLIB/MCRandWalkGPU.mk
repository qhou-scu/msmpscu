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
objname := RandWalk

#sorce directories
ifeq ($(origin MCLIBDIRS), undefined) 
MCLIBDIRS := MCLIB/sor/
endif
sor  := $(MCLIBDIRS)$(objname)/

#target directories
ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

ifeq ($(origin LIBDIRD), undefined) 
LIBDIRD := $(WORKSPACE)/LIB/$(ConfigName)
endif
tgt  := $(LIBDIRD)

#target lib name
libname  := libMC_$(objname).$(LIB_EXT)

##########################################################
mlist    :=  mc_randwalk_const               \
             mc_randwalk_exp_gpu             \
             mc_randwalk_gv                  \
             mc_randwalk_pow_gpu             \
             mc_randWalk_register            \
             mc_typedef_randwalkctrl         \
             mc_typedef_randwalker_base      \
             mc_typedef_randWalkerstate      \
             mc_typedef_randWalkerstatNumb   \
             mc_randwalk_evolve

nlist    :=  MC_RandWalk_Const               \
             MC_RandWalk-Exp-GPU             \
             MC_RandWalk_GV_GPU              \
             MC_RandWalk-Pow-GPU             \
             MC_RandWalk_Register            \
             MC_TypeDef_RandWalkCtrl         \
             MC_TypeDef_RandWalker_Base      \
             MC_TypeDef_RandWalkerState      \
             MC_TypeDef_RandWalkerStatNumb   \
             MC_RandWalk_Evolve


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

$(tgt)MC_RandWalk_Const.o : $(sor)MC_RandWalk_Const.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
	
$(tgt)MC_RandWalk_GV_GPU.o : $(sor)MC_RandWalk_GV_GPU.F90 $(tgt)MC_TypeDef_RandWalker_Base.o \
                                                          $(tgt)MC_TypeDef_RandWalkCtrl.o    \
                                                          $(tgt)MC_TypeDef_RandWalkerState.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
	
	
$(tgt)MC_RandWalk-Exp-GPU.o : $(sor)MC_RandWalk-Exp-GPU.F90 $(tgt)MC_RandWalk_GV_GPU.o 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
	

$(tgt)MC_RandWalk-Pow-GPU.o : $(sor)MC_RandWalk-Pow-GPU.F90 $(tgt)MC_RandWalk_GV_GPU.o 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
	
	
$(tgt)MC_RandWalk_Register.o : $(sor)MC_RandWalk_Register.F90 $(tgt)MC_RandWalk_GV_GPU.o     \
                                                              $(tgt)MC_RandWalk-Exp-GPU.o    \
                                                              $(tgt)MC_RandWalk-Pow-GPU.o


	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MC_RandWalk_Evolve.o : $(sor)MC_RandWalk_Evolve.F90 $(tgt)MC_TypeDef_RandWalker_Base.o \
                                                          $(tgt)MC_RandWalk_GV_GPU.o         \
                                                          $(tgt)MC_RandWalk_Register.o       \
                                                          $(tgt)MC_TypeDef_RandWalkCtrl.o    \
                                                          $(tgt)MC_TypeDef_RandWalkerState.o \
                                                          $(tgt)MC_TypeDef_RandWalkerStatNumb.o
                                                          
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MC_TypeDef_RandWalkCtrl.o : $(sor)MC_TypeDef_RandWalkCtrl.F90 $(tgt)MC_RandWalk_Const.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MC_TypeDef_RandWalker_Base.o : $(sor)MC_TypeDef_RandWalker_Base.F90 $(tgt)MC_RandWalk_Const.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MC_TypeDef_RandWalkerState.o : $(sor)MC_TypeDef_RandWalkerState.F90 $(tgt)MC_TypeDef_RandWalker_Base.o \
                                                                          $(tgt)MC_TypeDef_RandWalkCtrl.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
                                                                          

$(tgt)MC_TypeDef_RandWalkerStatNumb.o : $(sor)MC_TypeDef_RandWalkerStatNumb.F90 $(tgt)MC_TypeDef_RandWalker_Base.o \
                                                                                $(tgt)MC_TypeDef_RandWalkCtrl.o    \
                                                                                $(tgt)MC_TypeDef_RandWalkerState.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

clean:
	-rm $(objects) $(libname) $(modules) 
