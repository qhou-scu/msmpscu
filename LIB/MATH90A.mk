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
sor  := $(LIBDIRS)/MATH90A/

#target directories
ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

ifeq ($(origin LIBDIRD), undefined) 
LIBDIRD := $(WORKSPACE)/LIB/$(ConfigName)
endif
tgt  := $(LIBDIRD)

nlist    :=  MATH90A
mlist    :=  math90a_module
modules  := $(foreach n, $(mlist), $(tgt)$(n).mod)
objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)

libname  := libMATH90A.$(LIB_EXT)
##########################################################
$(libname) : $(objects) 
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)MATH90A.o : $(sor)MATH90A.F 
	$(comp) -c $(oflags)  -module $(tgt) $(sor)MATH90A.F -o $(tgt)MATH90A.o 

clean:
	-rm $(objects) $(libname) $(modules) 
