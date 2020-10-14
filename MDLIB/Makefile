#choice compiler and compiling flags
export comp       := pgfortran

ifeq ($(origin oflags), undefined)
export oflags := = -fast -tp sandybridge-64  -Mvect=sse,simd -Minline -Mconcur -Minform=warn -Minfo=accel,inline,intensity,loop,mp,opt,par,vect \
                     -Mcuda=fastmath,cuda8.0,cc2x 
endif

ifeq ($(origin ConfigName), undefined) 
export ConfigName := Release
endif

#####################################
export LIBDIRS    := $(MSMPSCUSOR)/LIB/sor/f/
libnames          := MiniUtilities RandGenerators CudaRanGeneratorC2F LBFGS MATH_LBFGSB MATH90A

export MSMLIBDIRS := $(MSMPSCUSOR)/MSMLIB/sor/
msmlibnames       := Common CommonGPU


export MDLIBDIRS  := $(MSMPSCUSOR)/MDLIB/sor/
mdlibnames        := Common CommonGPU Appshell

export POTLIBDIRS := $(MSMPSCUSOR)/MDLIB/sor/Potentials/
potlibnames       := EAM_NIST\
                     EM_TB_WangJun_W-HE_2010\
                     EM_TB_Ackland_Vitek_PRB41_10324   \
                     EM_TB_AWB_Phil_Mag_A71_1995_553   \
                     EM_TB_Cleri_Rosato_PRB48_22       \
                     EAM_WW_Marinica_JPCM25_2013       \
                     EAM_WHeH_Bonny_JPCM26_2014        \
                     EAM_Zr_Mendelev_Phil_Mag_L87_2007 \

export APPLIBDIRS := $(MSMPSCUSOR)/MDLIB/sor/
applibnames       := Embedment \
                     Deposition

export ANALYLIBDIRS := $(MSMPSCUSOR)/MDLIB/sor/Analysis/
analysislibnames    := Analysis

export LIBDIR0    := $(WORKSPACE)/LIB/
export LIBDIRD    := $(LIBDIR0)$(ConfigName)/
             
all:  
#the math libs:
	$(MAKE) -f $(MSMPSCUSOR)/LIB/Makefile

#the common MSM libs:
	$(MAKE) -f $(MSMPSCUSOR)/MSMLIB/Makefile

#the common MD libs:
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/MDCommon.mk
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/MDCommonGPU.mk

#user supplied potntials
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/EAM_NIST.mk
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/EM_TB_WangJun_W-HE_2010.mk
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/EM_TB_Ackland_Vitek_PRB41_10324.mk
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/EM_TB_AWB_Phil_Mag_A71_1995_553.mk
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/EM_TB_Cleri_Rosato_PRB48_22.mk
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/EAM_WW_Marinica_JPCM25_2013.mk
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/EAM_WHeH_Bonny_JPCM26_2014.mk
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/EAM_Zr_Mendelev_Phil_Mag_L87_2007.mk
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/EM_TB_WangJun_Ti-HE_2007.mk

#the LocalTempCtrl libs:
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/MDLocalTempCtrl.mk
	

#the BoostMethCtrl libs:
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/MDBoostMeth.mk	

#the app shells
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/MDAppshell.mk

#some utility libs
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/Embedment.mk
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/Deposition.mk

#some analysis libs
	$(MAKE) -f $(MSMPSCUSOR)/MDLIB/Analysis.mk

############################################
.PHONY: clean
clean:
	$(MAKE) clean -f  $(MSMPSCUSOR)/LIB/MiniUtilities.mk
	$(MAKE) clean -f  $(MSMPSCUSOR)/LIB/RandGenerators.mk
	$(MAKE) clean -f  $(MSMPSCUSOR)/MDLIB/MDCommon.mk	
	$(MAKE) clean -f  $(MSMPSCUSOR)/MDLIB/MDCommonGPU.mk	
	$(MAKE) clean -f  $(MSMPSCUSOR)/MDLIB/MDAppshell.mk
	$(MAKE) clean -f  $(MSMPSCUSOR)/MDLIB/MDAppshellGPU.mk	
	$(MAKE) clean -f  $(MSMPSCUSOR)/MDLIB/EM_TB_WangJun_W-HE_2010.mk
	$(MAKE) clean -f  $(MSMPSCUSOR)/MDLIB/EM_TB_Ackland_Vitek_PRB41_10324.mk
	$(MAKE) clean -f  $(MSMPSCUSOR)/MDLIB/EM_TB_AWB_PHIL_MAG_A71_1995_553.mk
	$(MAKE) clean -f  $(MSMPSCUSOR)/MDLIB/EM_TB_Cleri_Rosato_PRB48_22.mk
	$(MAKE) clean -f  $(MSMPSCUSOR)/MDLIB/Embedment.mk
	$(MAKE) clean -f  $(MSMPSCUSOR)/MDLIB/Deposition.mk
	
