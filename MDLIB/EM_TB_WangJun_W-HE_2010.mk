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
objname := EM_TB_WangJun_W-HE_2010

#sorce directories
ifeq ($(origin POTLIBDIRS), undefined) 
POTLIBDIRS := MDLIB/sor/Potentials/
endif
sor  := $(POTLIBDIRS)$(objname)/

#target directories
ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

ifeq ($(origin LIBDIRD), undefined) 
LIBDIRD := $(WORKSPACE)/LIB/$(ConfigName)
endif
tgt  := $(LIBDIRD)

#target lib name
libname  := libPOT_$(objname).$(LIB_EXT)

#######################################################
mlist    :=  em_tb_forcetable_wangjun_w_he_2010  \
             pot_fs_ack_ww                       \
             pot_fs_pm_ww                        \
             pot_lj_hehe                         \
             pot_zbl_abintio_he_w_513            \
             pot_zbl_abintio_he_w_614            \
             pot_zbl_abintio_he_w_818            \
             pot_zbl_abintio_he_w_930            \
             pot_zbl_abintio_he_w_930            \
             pot_zbl_abintio_he_w                \
             pot_zbl_exp6_hehe                   \
             pot_zbl_he_w                

nlist    :=  EM_TB_ForceTable_WangJun_W_HE_2010  \
             FS_Ackland_WW                       \
             FS_PM_WW                            \
             LJ_HeHe                             \
             ZBL_ABINTIO_He_W_513                \
             ZBL_ABINTIO_He_W_614                \
             ZBL_ABINTIO_He_W_818                \
             ZBL_ABINTIO_He_W_930                \
             ZBL_ABINTIO_He_W                    \
             ZBL_EXP6_HeHe                       \
             ZBL_He_W                                         

objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)

$(tgt)EM_TB_ForceTable_WangJun_W_HE_2010.o : $(sor)EM_TB_ForceTable_WangJun_W_HE_2010.F90  \
             $(tgt)FS_Ackland_WW.o                       \
             $(tgt)FS_PM_WW.o                            \
             $(tgt)LJ_HeHe.o                             \
             $(tgt)ZBL_ABINTIO_He_W_513.o                \
             $(tgt)ZBL_ABINTIO_He_W_614.o                \
             $(tgt)ZBL_ABINTIO_He_W_818.o                \
             $(tgt)ZBL_ABINTIO_He_W_930.o                \
             $(tgt)ZBL_ABINTIO_He_W.o                    \
             $(tgt)ZBL_EXP6_HeHe.o                       \
             $(tgt)ZBL_He_W.o              

	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@

$(tgt)FS_Ackland_WW.o : $(sor)FS_Ackland_WW.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
$(tgt)FS_PM_WW.o : $(sor)FS_PM_WW.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
$(tgt)LJ_HeHe.o : $(sor)LJ_HeHe.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
$(tgt)ZBL_ABINTIO_He_W_513.o : $(sor)ZBL_ABINTIO_He_W_513.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
$(tgt)ZBL_ABINTIO_He_W_614.o : $(sor)ZBL_ABINTIO_He_W_614.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
$(tgt)ZBL_ABINTIO_He_W_818.o : $(sor)ZBL_ABINTIO_He_W_818.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
$(tgt)ZBL_ABINTIO_He_W_930.o : $(sor)ZBL_ABINTIO_He_W_930.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
$(tgt)ZBL_ABINTIO_He_W.o : $(sor)ZBL_ABINTIO_He_W.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
$(tgt)ZBL_EXP6_HeHe.o : $(sor)ZBL_EXP6_HeHe.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@
$(tgt)ZBL_He_W.o : $(sor)ZBL_He_W.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $< -o $@


######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
