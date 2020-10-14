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
sor  := $(LIBDIRS)/LBFGS/

#target directories
ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

ifeq ($(origin LIBDIRD), undefined) 
LIBDIRD := $(WORKSPACE)/LIB/$(ConfigName)
endif
tgt  := $(LIBDIRD)

nlist    :=  LBFGS
mlist    :=  lbfgs_module
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)

libname  := libLBFGS.$(LIB_EXT)
##########################################################
$(libname) : $(objects) 
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)

$(tgt)LBFGS.o : $(sor)LBFGS.F90 
	$(comp) -c $(oflags)  -module $(tgt) $(sor)LBFGS.F90 -o $(tgt)LBFGS.o 
	

	

clean:
	-rm $(objects) $(libname) $(modules) 
