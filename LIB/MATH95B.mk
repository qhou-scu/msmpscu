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
LIBDIRS := LIB/sor/f/
endif
sor  := $(LIBDIRS)/MATH95B/

#target directories
ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

ifeq ($(origin LIBDIRD), undefined) 
LIBDIRD := $(WORKSPACE)/LIB/$(ConfigName)
endif
tgt  := $(LIBDIRD)

nlist  :=   ADLINT  DAVINT DCHFIE DGAUS8 DPCHIA        \
            DPCHID  DQAGE DQAG DQAGSE DQAGS            \
            DQELG  DQK21 DQK31 DQK41 DQK51             \
            DQK61 DQNC79 DQNG DQPSRT DSIMP             \
            DSIMPF 
modules  := $(foreach n, $(nlist), $(tgt)$(n).mod)
objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)

libname  := libMATH95B.$(LIB_EXT)
##########################################################
$(libname) : $(objects)
	ar -rcs $(libname) $(objects) 
	mv $(libname) $(tgt)
	
$(tgt)%.o : $(sor)%.F
	$(comp) -c $(oflags)  $< -o $@
clean:
	-rm $(objects) $(libname)
