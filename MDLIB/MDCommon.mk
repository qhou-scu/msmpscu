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
ifeq ($(origin MDLIBDIRS), undefined) 
MDLIBDIRS := $(MDPSCUSOR)/MDLIB/sor/
endif
sor  := $(MDLIBDIRS)$(objname)/

#target lib name
tgt  := $(LIBDIRD)
libname  := libMD_$(objname).$(LIB_EXT)

#######################################################
mlist    :=  md_analysistools               \
             md_bp_coupling                 \
             md_constants                   \
             md_corepotent                  \
             md_forceclass_register         \
             md_fs_force_table              \
             md_globle_variables            \
             md_neighborslist               \
             md_pot_eam_utilities           \
             md_simulationbox_array         \
             md_simctrlparam_art            \
             md_simctrlparam_bst            \
             md_simctrlparam_gmd            \
             md_simctrlparam_neb            \
             md_simulationctrlparam_parrep  \
             md_simctrlparam_tad            \
             md_swope_scheme                \
             md_typedef_clusterlist         \
             md_typedef_datapad             \
             md_typedef_epcctrl             \
             md_typedef_forcetable          \
             md_typedef_InputPaser          \
             md_typdef_printlist            \
             md_typedef_recordlist          \
             md_typedef_recordstamp         \
             md_typedef_simulationbox       \
             md_typedef_simulationctrlparam
           
             
nlist    :=  MD_ANALYSISTOOLS               \
             MD_BP_Coupling                 \
             MD_CorePotent                  \
             MD_Const                       \
             MD_ForceClass_Register         \
             MD_FS_ForceTable               \
             MD_Gvar                        \
             MD_NeighborsList               \
             MD_Pot_EAM_Utilities           \
             MD_SimBoxArray                 \
             MD_SimCtrlParam_ART            \
             MD_SimCtrlParam_BST            \
             MD_SimCtrlParam_GMD            \
             MD_SimCtrlParam_NEB            \
             MD_SimCtrlParam_PARREP         \
             MD_SimCtrlParam_TAD            \
             MD_SwopeScheme                 \
             MD_TypeDef_ClusterList         \
             MD_TypeDef_DataPad             \
             MD_TypeDef_EPCCtrl             \
             MD_TypeDef_ForceTable          \
             MD_TypeDef_InputPaser          \
             MD_TypeDef_PrintList           \
             MD_TypeDef_RecordList          \
             MD_TypeDef_RecordStamp         \
             MD_TypeDef_SimBox              \
             MD_TypeDef_SimCtrlParam        \
             MD_TypeDef_StpRangTable  

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)MD_Const.o : $(sor)MD_Const.F90 
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_ANALYSISTOOLS.o : $(sor)MD_ANALYSISTOOLS.F90 $(tgt)MD_Const.o
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_BP_Coupling.o : $(sor)MD_BP_Coupling.F90  \
                            $(tgt)MD_Const.o          \
                            $(tgt)MD_TypeDef_SimBox.o \
                            $(tgt)MD_TypeDef_SimCtrlParam.o
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_CorePotent.o : $(sor)MD_CorePotent.F90  \
                                 $(tgt)MD_Const.o          
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@


$(tgt)MD_ForceClass_Register.o : $(sor)MD_ForceClass_Register.F90         \
                                        $(tgt)MD_Const.o                  \
                                        $(tgt)MD_TypeDef_SimBox.o         \
                                        $(tgt)MD_FS_ForceTable.o          \
                                        $(tgt)MD_TypeDef_ForceTable.o
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@


$(tgt)MD_FS_ForceTable.o : $(sor)MD_FS_ForceTable.F90                  \
                                        $(tgt)MD_Const.o                  \
                                        $(tgt)MD_TypeDef_SimBox.o         \
                                        $(tgt)MD_TypeDef_SimCtrlParam.o   \
                                        $(tgt)MD_TypeDef_ForceTable.o     \
                                        $(tgt)MD_NeighborsList.o
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@


$(tgt)MD_Gvar.o : $(sor)MD_Gvar.F90                     \
                     $(tgt)MD_Const.o                   \
                     $(tgt)MD_TypeDef_SimBox.o          \
                     $(tgt)MD_TypeDef_SimCtrlParam.o    \
                     $(tgt)MD_SimCtrlParam_ART.o        \
                     $(tgt)MD_SimCtrlParam_BST.o        \
                     $(tgt)MD_SimCtrlParam_GMD.o        \
                     $(tgt)MD_SimCtrlParam_NEB.o        \
                     $(tgt)MD_SimCtrlParam_PARREP.o     \
                     $(tgt)MD_SimCtrlParam_TAD.o        \
                     $(tgt)MD_TypeDef_PrintList.o
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_NeighborsList.o : $(sor)MD_NeighborsList.F90      \
                              $(tgt)MD_Const.o                \
                              $(tgt)MD_TypeDef_SimBox.o       \
                              $(tgt)MD_TypeDef_SimCtrlParam.o   
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_Pot_EAM_Utilities.o : $(sor)MD_Pot_EAM_Utilities.F90      \
                              $(tgt)MD_Const.o      
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_SimBoxArray.o : $(sor)MD_SimBoxArray.F90           \
                            $(tgt)MD_Const.o                \
                            $(tgt)MD_TypeDef_SimBox.o       \
                            $(tgt)MD_TypeDef_SimCtrlParam.o   
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_SimCtrlParam_ART.o : $(sor)MD_SimCtrlParam_ART.F90     \
                            $(tgt)MD_Const.o                    \
                            $(tgt)MD_TypeDef_SimBox.o           \
                            $(tgt)MD_TypeDef_SimCtrlParam.o   
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@
	
$(tgt)MD_SimCtrlParam_BST.o : $(sor)MD_SimCtrlParam_BST.F90     \
                            $(tgt)MD_Const.o                    \
                            $(tgt)MD_TypeDef_SimBox.o           \
                            $(tgt)MD_TypeDef_SimCtrlParam.o   
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@
	

$(tgt)MD_SimCtrlParam_GMD.o : $(sor)MD_SimCtrlParam_GMD.F90     \
                            $(tgt)MD_Const.o                    \
                            $(tgt)MD_TypeDef_SimBox.o           \
                            $(tgt)MD_TypeDef_SimCtrlParam.o   
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_SimCtrlParam_NEB.o : $(sor)MD_SimCtrlParam_NEB.F90     \
                            $(tgt)MD_Const.o                    \
                            $(tgt)MD_TypeDef_SimBox.o           \
                            $(tgt)MD_TypeDef_SimCtrlParam.o   
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_SimCtrlParam_PARREP.o : $(sor)MD_SimCtrlParam_PARREP.F90   \
                            $(tgt)MD_Const.o                        \
                            $(tgt)MD_TypeDef_SimBox.o               \
                            $(tgt)MD_TypeDef_SimCtrlParam.o   
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_SimCtrlParam_TAD.o : $(sor)MD_SimCtrlParam_TAD.F90   \
                            $(tgt)MD_Const.o                  \
                            $(tgt)MD_TypeDef_SimBox.o         \
                            $(tgt)MD_TypeDef_SimCtrlParam.o   
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@


$(tgt)MD_SwopeScheme.o : $(sor)MD_SwopeScheme.F90           \
                            $(tgt)MD_Const.o                \
                            $(tgt)MD_TypeDef_SimBox.o       \
                            $(tgt)MD_TypeDef_SimCtrlParam.o   
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_TypeDef_ClusterList.o : $(sor)MD_TypeDef_ClusterList.F90   \
                                   $(tgt)MD_Const.o                  
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_TypeDef_DataPad.o : $(sor)MD_TypeDef_DataPad.F90    \
                                $(tgt)MD_Const.o             
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_TypeDef_EPCCtrl.o : $(sor)MD_TypeDef_EPCCtrl.F90    \
                                $(tgt)MD_Const.o             \
                                $(tgt)MD_TypeDef_InputPaser.o
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_TypeDef_ForceTable.o : $(sor)MD_TypeDef_ForceTable.F90   \
                                   $(tgt)MD_Const.o                  \
                                   $(tgt)MD_TypeDef_SimBox.o         \
                                   $(tgt)MD_TypeDef_SimCtrlParam.o   
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_TypeDef_InputPaser.o : $(sor)MD_TypeDef_InputPaser.F90   \
                                   $(tgt)MD_Const.o                  
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_TypeDef_PrintList.o : $(sor)MD_TypeDef_PrintList.F90   \
                                   $(tgt)MD_Const.o                  
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@


$(tgt)MD_TypeDef_RecordList.o : $(sor)MD_TypeDef_RecordList.F90      \
                                   $(tgt)MD_Const.o                  \
                                   $(tgt)MD_TypeDef_SimBox.o         \
                                   $(tgt)MD_TypeDef_SimCtrlParam.o   
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_TypeDef_RecordStamp.o : $(sor)MD_TypeDef_RecordStamp.F90    \
                                   $(tgt)MD_Const.o                  
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_TypeDef_SimBox.o : $(sor)MD_TypeDef_SimBox.F90       \
                               $(tgt)MD_Const.o               \
                               $(tgt)MD_TypeDef_DataPad.o     \
                               $(tgt)MD_TypeDef_InputPaser.o  \
                               $(tgt)MD_TypeDef_RecordStamp.o          
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@


$(tgt)MD_TypeDef_SimCtrlParam.o : $(sor)MD_TypeDef_SimCtrlParam.F90  \
                                     $(tgt)MD_Const.o                \
                                     $(tgt)MD_TypeDef_InputPaser.o   \
                                     $(tgt)MD_TypeDef_SimBox.o       \
                                     $(tgt)MD_TypeDef_EPCCtrl.o    
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)MD_TypeDef_StpRangTable.o : $(sor)MD_TypeDef_StpRangTable.F90  \
                                     $(tgt)MD_Const.o                 
      
	$(comp) -c $(oflags) -I$(incdir) -module $(tgt) $< -o $@


######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
