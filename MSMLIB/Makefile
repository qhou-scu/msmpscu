#choice compiler and compiling flags
export comp       := pgfortran

ifeq ($(origin oflags), undefined)
export oflags := = -fast -tp sandybridge-64  -Mvect=sse,simd -Minline -Mconcur -Minform=warn -Minfo=accel,inline,intensity,loop,mp,opt,par,vect \
                     -Mcuda=fastmath,cuda8.0,cc2x 
endif

ifeq ($(origin ConfigName), undefined) 
export ConfigName := Release
endif

export LIBDIRS    := $(MSMPSCUSOR)/LIB/sor/f/
libnames          := MiniUtilities RandGenerators CudaRanGeneratorC2F LBFGS MATH_LBFGSB MATH90A

export MSMLIBDIRS := $(MSMPSCUSOR)/MSMLIB/sor/
msmlibnames       := Common CommonGPU

export LIBDIR0    := $(WORKSPACE)/LIB/
export LIBDIRD    := $(LIBDIR0)$(ConfigName)/
             
all:  
#the math libs:
	$(MAKE) -f $(MSMPSCUSOR)/LIB/Makefile

#the common MSM libs:
	$(MAKE) -f $(MSMPSCUSOR)/MSMLIB/MSMCommon.mk
	$(MAKE) -f $(MSMPSCUSOR)/MSMLIB/MSMCommonGPU.mk

############################################
.PHONY: clean
clean:
	$(MAKE) clean -f $(MSMPSCUSOR)/LIB/Makefile
	$(MAKE) clean -f $(MSMPSCUSOR)/MSMLIB/MSMCommon.mk
	$(MAKE) clean -f $(MSMPSCUSOR)/MSMLIB/MSMCommonGPU.mk

	
