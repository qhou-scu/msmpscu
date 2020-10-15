
# include directeries
LIBDIRD    := $(WORKSPACE)/LIB/$(ConfigName)/
incdir     := $(WORKSPACE)/LIB/$(ConfigName)/

#  external dependence lib directories
DEPENDLIBD := $(MSMPSCUSOR)/Dependence/


##########################################################

#objective name
objname := $(OBJNAME)
exsor   := $(EXOBJS)

#sorce name
# NOTE: in CYGWIN, we have use the full name for the source
#       otherwise, the compiler cannot find the source
sormain := $(APPDIR)/$(objname).F90

#target directories
ifeq ($(ConfigName), Debug)
tgt  := $(WORKSPACE)/applications/$(ConfigName)/
else
tgt  := $(WORKSPACE)/applications/
endif

#executable name
exename   := $(tgt)$(objname).exe
##########################################################
sorsF90  := $(foreach n, $(exsor),   $(APPDIR)/$(n).F90)
objs     := $(foreach n, $(exsor),   $(tgt)$(n).o)
linkobjs := $(foreach n, $(objs),    -L$(n))


################### MD libs #######################################
mdlibs     := $(foreach n, $(utilibnames),        $(LIBDIRD)libMD_$(n).$(LIB_EXT))
mdlibs     += $(foreach n, $(boostmethslibnames), $(LIBDIRD)libMD_$(n).$(LIB_EXT))
mdlibs     += $(foreach n, $(ltempctlibnames),    $(LIBDIRD)libMD_$(n).$(LIB_EXT))
mdlibs     += $(foreach n, $(potlibnames),        $(LIBDIRD)libPOT_$(n).$(LIB_EXT))
mdlibs     += $(foreach n, $(potlibnames),        $(LIBDIRD)libPOT_$(n).$(LIB_EXT))
mdlibs     += $(foreach n, $(mdlibnames),         $(LIBDIRD)libMD_$(n).$(LIB_EXT))
mdlibs     += $(foreach n, $(msmcomlibnames),     $(LIBDIRD)libMSM_$(n).$(LIB_EXT))
mdlibs     += $(foreach n, $(libnames),           $(LIBDIRD)lib$(n).$(LIB_EXT))
mdlibs     += $(foreach n, $(dependlibname),      $(DEPENDLIBD)lib$(n).$(LIB_EXT))


mdliblist  := $(foreach n, $(utilibnames),        -L$(LIBDIRD) -lMD_$(n))
mdliblist  += $(foreach n, $(boostmethslibnames), -L$(LIBDIRD) -lMD_$(n))
mdliblist  += $(foreach n, $(ltempctlibnames),    -L$(LIBDIRD) -lMD_$(n))
mdliblist  += $(foreach n, $(potlibnames),        -L$(LIBDIRD) -lPOT_$(n))
mdliblist  += $(foreach n, $(mdlibnames),         -L$(LIBDIRD) -lMD_$(n))
mdliblist  += $(foreach n, $(msmcomlibnames),     -L$(LIBDIRD) -lMSM_$(n))
mdliblist  += $(foreach n, $(libnames),           -L$(LIBDIRD) -l$(n))
mdliblist  += $(foreach n, $(dependlibname),      -L$(DEPENDLIBD) -l$(n))

#######################################################

################## MC libs ########################################
mclibs     := $(foreach n, $(mclibnames),          $(LIBDIRD)libMC_$(n).$(LIB_EXT))
mclibs     += $(foreach n, $(msmcomlibnames),      $(LIBDIRD)libMSM_$(n).$(LIB_EXT))
mclibs     += $(foreach n, $(libnames),            $(LIBDIRD)lib$(n).$(LIB_EXT))
mclibs     += $(foreach n, $(dependlibname),       $(DEPENDLIBD)lib$(n).$(LIB_EXT))


mcliblist  := $(foreach n, $(randwalklibnames),   -L$(LIBDIRD) -lMC_$(n))
mcliblist  += $(foreach n, $(msmcomlibnames),     -L$(LIBDIRD) -lMSM_$(n))
mcliblist  += $(foreach n, $(libnames),           -L$(LIBDIRD) -l$(n))
mcliblist  += $(foreach n, $(dependlibname),      -L$(DEPENDLIBD) -l$(n))

#######################################################

#######################################################
# Note 1: in Linux $(sormain) must intermediately follows  $(oflags)
#         otherwise, there is a link error
# 
$(exename) : $(objmain) $(objs) $(libs) $(sormain)
	$(comp)  $(oflags) $(sormain) $(patsubst %, -I%, $(incdir)) $(mdliblist) $(mcliblist) $(linkobjs) -o $(exename)

$(objs) : $(sorsF90)
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) $< -o $@

#$(objmain) : $(sormain) $(obj) 
#	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) -module $(tgt) $(sormain) -o $(objmain) 

clean:
	-rm  $(obj) $(objmain) $(tgt)$(objname).mod $(exename) 
