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
objname := Analysis

#sorce directories
ifeq ($(origin ANALYLIBDIRS), undefined) 
ANALYLIBDIRS := $(MDPSCUSOR)/MDLIB/sor/Analysis/
endif
sor  := $(ANALYLIBDIRS)

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
mlist    :=  \
             atomicvirialstress_gpu     \
             atomselectorcommon         \
             atomtrackcommon            \
             boxselector                \
             centrosymparam_for_bcc_gpu \
             centrosymparam_for_fcc_gpu \
             cfg_track_cluster_gpu      \
             coordnum_gpu               \
             createoctalattice_bcc      \
             createoctalattice_bcc      \
             defectcluster_analysis_bcc \
             defectcluster_analysis_hcp \
             embedment_diffusion1       \
             embedment_diffusion2       \
             latticeclustertrackcommon  \
             propcluster_atomType       \
             propclustercommon_14_gpu   \
             refenergy_gpu              \
             refstrainfield_gpu         \
             refvoronoivacancy_2014_gpu \
             voronoitessellation_13_cpu \
             voronoitessellation_13_gpu   

nlist    :=  \
             AtomicVirialStress_GPU     \
             AtomSelectorCommon         \
             AtomTrackCommon            \
             BoxSelector                \
             CentroSymParam_for_BCC_GPU \
             CentroSymParam_for_FCC_GPU \
             Cfg_Track_Cluster_GPU      \
             CoordNum_GPU               \
             CreateOctaLattice_BCC      \
             CreateTetraLattice_BCC     \
             DefectCluster_Analysis_bcc \
             DefectCluster_Analysis_hcp \
             EMBEDMENT_Diffusion1       \
             EMBEDMENT_Diffusion2       \
             LatticeClusterTrackCommon  \
             PropCluster_AtomType       \
             PropClusterCommon_GPU      \
             RefEnergy_GPU              \
             RefStrainField_GPU         \
             RefVoronoiVA_GPU           \
             VoronoiTessellationM_CPU   \
             VoronoiTessellationM_GPU                       

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)AtomicVirialStress_GPU.o : $(sor)AtomicVirialStress_GPU.F90 $(tgt)VoronoiTessellationM_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)AtomSelectorCommon.o : $(sor)AtomSelectorCommon.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)AtomTrackCommon.o : $(sor)AtomTrackCommon.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)BoxSelector.o : $(sor)BoxSelector.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)CentroSymParam_for_BCC_GPU.o : $(sor)CentroSymParam_for_BCC_GPU.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)CentroSymParam_for_FCC_GPU.o : $(sor)CentroSymParam_for_FCC_GPU.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)Cfg_Track_Cluster_GPU.o : $(sor)Cfg_Track_Cluster_GPU.F90 $(tgt)RefVoronoiVA_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)CoordNum_GPU.o : $(sor)CoordNum_GPU.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)CreateOctaLattice_BCC.o : $(sor)CreateOctaLattice_BCC.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)CreateTetraLattice_BCC.o : $(sor)CreateTetraLattice_BCC.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)DefectCluster_Analysis_bcc.o : $(sor)DefectCluster_Analysis_bcc.F90 $(tgt)LatticeClusterTrackCommon.o  
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)DefectCluster_Analysis_hcp.o : $(sor)DefectCluster_Analysis_hcp.F90 $(tgt)LatticeClusterTrackCommon.o  
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)EMBEDMENT_Diffusion1.o : $(sor)EMBEDMENT_Diffusion1.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)EMBEDMENT_Diffusion2.o : $(sor)EMBEDMENT_Diffusion2.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)LatticeClusterTrackCommon.o : $(sor)LatticeClusterTrackCommon.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)PropCluster_AtomType.o : $(sor)PropCluster_AtomType.F90 $(tgt)PropClusterCommon_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)PropClusterCommon_GPU.o : $(sor)PropClusterCommon_GPU.F90 $(tgt)VoronoiTessellationM_GPU.o
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)RefEnergy_GPU.o : $(sor)RefEnergy_GPU.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)RefStrainField_GPU.o : $(sor)RefStrainField_GPU.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)RefVoronoiVA_GPU.o : $(sor)RefVoronoiVA_GPU.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)VoronoiTessellationM_CPU.o : $(sor)VoronoiTessellationM_CPU.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)VoronoiTessellationM_GPU.o : $(sor)VoronoiTessellationM_GPU.F90
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
