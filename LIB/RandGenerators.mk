#compiler
ifeq ($(origin comp), undefined) 
comp       = pgfortran
endif

ifeq ($(origin oflags), undefined)
oflags := -fast
endif

##########################################################
#sorce directories
ifeq ($(origin LIBDIRS), undefined) 
LIBDIRS := $(MDPSCUSOR)/LIB/sor/f/
endif
sor  := $(LIBDIRS)/RandGenerators/

#target directories
ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

ifeq ($(origin LIBDIRD), undefined) 
LIBDIRD := $(WORKSPACE)/LIB/$(ConfigName)
endif
tgt  := $(LIBDIRD)

nlist    :=  DRAND32 DRAND32SEEDLIB
mlist    :=  rand32_module rand32seedlib_module
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)

libname  := libRandGenerators.$(LIB_EXT)
##########################################################
$(libname) : $(objects) 
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)DRAND32.o : $(sor)DRAND32.F90 
	$(comp) -c $(oflags)  -module $(tgt) $(sor)DRAND32.F90 -o $(tgt)DRAND32.o 

$(tgt)DRAND32SEEDLIB.o : $(sor)DRAND32SEEDLIB.F90 
	$(comp) -c $(oflags)  -module $(tgt) $(sor)DRAND32SEEDLIB.F90 -o $(tgt)DRAND32SEEDLIB.o

clean:
	-rm $(objects) $(libname) $(modules) 
