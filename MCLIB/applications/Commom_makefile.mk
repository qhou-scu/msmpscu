
# include directeries
LIBDIRD    := $(WORKSPACE)/LIB/$(ConfigName)/
incdir     := $(WORKSPACE)/LIB/$(ConfigName)/

#  external dependence lib directories
DEPENDLIBD := $(MSMPSCUSOR)/Dependence/


##########################################################

#objective name
objname := $(OBJNAME)

#sorce name
# NOTE: in CYGWIN, we have use the full name for the source
#       otherwise, the compiler cannot find the source
sormain := $(MSMPSCUSOR)/MCLIB/sor/Applications/$(objname).F90

#target directories
ifeq ($(ConfigName), Debug)
tgt  := $(WORKSPACE)/applications/$(ConfigName)/
else
tgt  := $(WORKSPACE)/applications/
endif

#executable name
exename   := $(tgt)$(objname).exe

##########################################################
libs     := $(foreach n, $(msmcomlibnames),  $(LIBDIRD)libMSM_$(n).$(LIB_EXT))
libs     += $(foreach n, $(libnames),        $(LIBDIRD)lib$(n).$(LIB_EXT))
libs     += $(foreach n, $(dependlibname),   $(DEPENDLIBD)lib$(n).$(LIB_EXT))


liblist  := $(foreach n, $(mclibnames),      -L$(LIBDIRD) -lMC_$(n))
liblist  += $(foreach n, $(msmcomlibnames),  -L$(LIBDIRD) -lMSM_$(n))
liblist  += $(foreach n, $(libnames),        -L$(LIBDIRD) -l$(n))
liblist  += $(foreach n, $(dependlibname),   -L$(DEPENDLIBD) -l$(n))


#######################################################
#######################################################
# Note 1: in Linux $(sormain) must intermediately follows  $(oflags)
#         otherwise, there is a link error
# 
$(exename) : $(objmain) $(libs) $(sormain)
	$(comp)  $(oflags) $(sormain) $(patsubst %, -I%, $(incdir)) $(obj) $(liblist) -o $(exename)

#$(objmain) : $(sormain) $(obj) 
#	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $(sormain) -o $(objmain) 

clean:
	-rm  $(obj) $(objmain) $(tgt)$(objname).mod $(exename) 
