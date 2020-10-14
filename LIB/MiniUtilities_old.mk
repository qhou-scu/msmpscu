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
sor  := $(LIBDIRS)/MiniUtilities/

#target directories
ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

ifeq ($(origin LIBDIRD), undefined) 
LIBDIRD := $(WORKSPACE)/LIB/$(ConfigName)
endif
tgt  := $(LIBDIRD)

nlist  :=   ARRTOSTR  AvailableIOUnit CATSTR CLOCK1 COMBSTR  \
            CSTRTOF DINTF2 DRSTR DTOS  EXNUMB                \
            EXOPTSTR EXSUBSTR1 EXSUBSTR FINDMVAL FSTRROC     \
            GETKEYWORD GETFNAME GETINPUTSTRLINE GETPATH ISTR \
            LOWCASE NBIT1 NSTOKEN NTOKENS ONWARNING          \
            RSTR STOD STRCATI TIMER TOFNAME UPCASE
modules  := $(foreach n, $(nlist), $(tgt)$(n).mod)
objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)

libname  := $(tgt)libOthers.$(LIB_EXT)
##########################################################
$(libname) : $(objects)
	ar -rcs $(libname) $(objects) 

$(tgt)%.o : $(sor)%.F
	$(comp) -c $(oflags)  $< -o $@
clean:
	-rm $(objects) $(libname)
