#compiler
ifeq ($(origin comp), undefined) 
comp       = pgfortran
endif

ifeq ($(origin oflags), undefined)
oflags := -fast -Mvect=sse -Minline -Mconcur -Minform=warn -Minfo=accel,inline,intensity,loop,mp,opt,par,vect 
endif
oflags += -Mcuda=fastmath

# include directeries
LIBDIRD    := $(WORKSPACE)/LIB/$(ConfigName)/
incdir     := $(WORKSPACE)/LIB/$(ConfigName)/

##########################################################
#sorce directories
ifeq ($(origin MDLIBDIRS), undefined) 
MDLIBDIRS := ./sor
endif

#objective name
objname := Make_Sandwich

#sorce directories
sor  := $(MDLIBDIRS)/

#target directories
ifeq ($(ConfigName), Debug)
tgt  := $(WORKSPACE)/applications/$(ConfigName)/
else
tgt  := $(WORKSPACE)/applications/
endif

#executable name
exename  := $(tgt)$(objname).exe
##########################################################
nlist    := Make_Sandwich_Main   

libs     := $(foreach n, $(utilibnames), $(LIBDIRD)libMD_$(n).$(LIB_EXT))
libs     += $(foreach n, $(potlibnames), $(LIBDIRD)libPOT_$(n).$(LIB_EXT))
libs     += $(foreach n, $(mdlibnames),  $(LIBDIRD)libMD_$(n).$(LIB_EXT))
libs     += $(foreach n, $(libnames),    $(LIBDIRD)lib$(n).$(LIB_EXT))

liblist  := $(foreach n, $(utilibnames), -L$(LIBDIRD) -lMD_$(n))
liblist  += $(foreach n, $(potlibnames), -L$(LIBDIRD) -lPOT_$(n))
liblist  += $(foreach n, $(mdlibnames),  -L$(LIBDIRD)  -lMD_$(n))
liblist  += $(foreach n, $(libnames),    -L$(LIBDIRD)  -l$(n))

#######################################################
$(exename) : $(tgt) $(tgt)Make_Sandwich_Main.o $(libs)
	$(comp)  $(oflags) $(tgt)Make_Sandwich_Main.o  $(liblist) -o $(exename)
	rm $(tgt)Make_Sandwich_Main.o

$(tgt) : 
	-mkdir $(WORKSPACE)/applications/
	-mkdir $(tgt)

$(tgt)Make_Sandwich_Main.o : $(sor)Make_Sandwich_Main.F90 
	$(comp) -c $(oflags) $(patsubst %, -I%, $(incdir)) $< -o $@ 


clean:
	-rm $(tgt)Make_Sandwich_Main.o Make_Sandwich_Main.mod $(exename) 
