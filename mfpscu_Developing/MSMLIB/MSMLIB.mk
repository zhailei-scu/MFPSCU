#compiler
ifeq ($(origin comp), undefined) 
comp       := pgfortran
endif

ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

oflags_this := $(oflags)

###########################################################
#sorce dir name
objname := MSMLIB

#sorce directories
ifeq ($(origin MSMLIBDIRS), undefined) 
MSMLIBDIRS := $(mcpscusor)$(Segment)MSMLIB$(Segment)sor
endif
sor  := $(MSMLIBDIRS)

#include dir
ifeq ($(origin LIBDIRD), undefined)
LIBDIRD := $(mcworkspace)$(Segment)LIB$(Segment)$(ConfigName)
endif
incdir := $(LIBDIRD)

#target directories
tgt  := $(LIBDIRD)

#target lib name
libname  := lib_$(objname).$(LIB_EXT)

#######################################################          
nlist :=  MSM_Const		  	  \
	        MSM_TYPEDEF_DataPad	  \
	        MSM_TYPEDEF_InputPaser	  \
	        MSM_MultiGPU_Basic
             
objects  := $(foreach n, $(nlist), $(tgt)$(Segment)$(n).o)
modules  := $(foreach n, $(nlist), $(tgt)$(Segment)$(n).mod)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects)
	mv $(libname) $(tgt)


$(tgt)$(Segment)MSM_Const.o : $(sor)$(Segment)Common$(Segment)MSM_Const.F90
	$(comp) -c $(oflags_this) -I$(incdir) -module $(tgt) $< -o $@


$(tgt)$(Segment)MSM_TYPEDEF_DataPad.o : $(sor)$(Segment)Common$(Segment)MSM_TypeDef_DataPad.F90	\
			                                  $(tgt)$(Segment)MSM_Const.o
	$(comp) -c $(oflags_this) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)$(Segment)MSM_TYPEDEF_InputPaser.o : $(sor)$(Segment)Common$(Segment)MSM_TypeDef_InputPaser.F90	\
			      	                             $(tgt)$(Segment)MSM_Const.o
	$(comp) -c $(oflags_this) -I$(incdir) -module $(tgt) $< -o $@

$(tgt)$(Segment)MSM_MultiGPU_Basic.o : $(sor)$(Segment)CommonGPU$(Segment)MSM_MultiGPU_Basic.F90	\
			                                 $(tgt)$(Segment)MSM_Const.o
	$(comp) -c $(oflags_this) -I$(incdir) -module $(tgt) $< -o $@

######################################################################
clean:
	-rm $(objects) $(libname) $(modules) 
