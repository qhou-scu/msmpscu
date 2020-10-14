#compiler
ifeq ($(origin comp), undefined) 
comp       = pgfortran
endif

ifeq ($(origin oflags), undefined)
oflags := -fast -Mvect=sse,simd -Minline -Mconcur -Minform=warn -Minfo=accel,inline,intensity,loop,mp,opt,par,vect -Mcuda=fastmath,cuda8.0,cc2x
endif

#include dir
incdir := $(LIBDIRD)

##########################################################
#sorce dir name
objname := LocalTempCtrlMeths

#sorce directories
ifeq ($(origin MDLIBDIRS), undefined) 
MDLIBDIRS := MDLIB/sor/
endif
sor     := $(MDLIBDIRS)$(objname)/
sorepc  := $(sor)EPC/
sorst   := $(sor)Stopping/

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
mlist    :=  md_ep_coupling_cpu              \
             md_ep_coupling_gpu              \
             md_localtempmethod_gpu          \
             stop_b_module                   \
             stop_srim_module                \
             stop_z85_module                 \
             stop_z95_module 

nlist    :=  MD_EP_Coupling_CPU               \
             MD_EP_Coupling_GPU               \
             MD_LocalTempMethod_GPU           \
             MD_ST_Coupling_GPU               \
             MD_STPLib_Register               \
             Stop_B                           \
             Stop_Srim                        \
             Stop_Z85                         \
             Stop_Z95   
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
	
######### for EPC
$(tgt)MD_EP_Coupling_CPU.o : $(sorepc)MD_EP_Coupling_CPU.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_EP_Coupling_GPU.o : $(sorepc)MD_EP_Coupling_GPU.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

######### for Stopping
$(tgt)Stop_B.o : $(sorst)Stop_B.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)Stop_Srim.o : $(sorst)Stop_Srim.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)Stop_Z85.o : $(sorst)Stop_Z85.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)Stop_Z95.o : $(sorst)Stop_Z95.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_STPLib_Register.o : $(sorst)MD_STPLib_Register.F90 \
                                   $(tgt)Stop_B.o           \
                                   $(tgt)Stop_Srim.o        \
                                   $(tgt)Stop_Z85.o         \
                                   $(tgt)Stop_Z95.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)MD_ST_Coupling_GPU.o : $(sorst)MD_ST_Coupling_GPU.F90 \
                                   $(tgt)MD_STPLib_Register.o           
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

######### for LocalTempCtrl
$(tgt)MD_LocalTempMethod_GPU.o : $(sor)MD_LocalTempMethod_GPU.F90 \
                                   $(tgt)MD_EP_Coupling_GPU.o      \
                                   $(tgt)MD_ST_Coupling_GPU.o

	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

clean:
	-rm $(objects) $(libname) $(modules) 
