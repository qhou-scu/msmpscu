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


export LIBDIR0    := $(WORKSPACE)/LIB/
export LIBDIRD    := $(LIBDIR0)$(ConfigName)/
             
all:  
#the math libs:
	$(MAKE) -f $(MSMPSCUSOR)/LIB/LBFGS.mk
	$(MAKE) -f $(MSMPSCUSOR)/LIB/MATH_LBFGSB.mk
	$(MAKE) -f $(MSMPSCUSOR)/LIB/MATH90A.mk
	$(MAKE) -f $(MSMPSCUSOR)/LIB/MiniUtilities.mk
	$(MAKE) -f $(MSMPSCUSOR)/LIB/RandGenerators.mk
	$(MAKE) -f $(MSMPSCUSOR)/LIB/CudaRanGeneratorC2F.mk


############################################
.PHONY: clean
clean:
	$(MAKE) clean -f $(MSMPSCUSOR)/LIB/LBFGS.mk
	$(MAKE) clean -f $(MSMPSCUSOR)/LIB/MATH_LBFGSB.mk
	$(MAKE) clean -f $(MSMPSCUSOR)/LIB/MATH90A.mk
	$(MAKE) clean -f $(MSMPSCUSOR)/LIB/MiniUtilities.mk
	$(MAKE) clean -f $(MSMPSCUSOR)/LIB/RandGenerators.mk
	$(MAKE) clean -f $(MSMPSCUSOR)/LIB/CudaRanGeneratorC2F.mk
	
