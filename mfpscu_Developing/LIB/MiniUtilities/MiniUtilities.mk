#*********************************************************************************
#--- Description: 
#--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
#--- Email : zhaileiytp@163.com
#--- data: From 2017-09 to 2021-12
#--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
#--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
#---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
#*********************************************************************************

#compiler
ifeq ($(origin comp), undefined)
comp       := pgfortran
endif

ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

oflags_this := $(oflags)

##########################################################
#sorce dir name
objname := MiniUtilities

#sorce directories
ifeq ($(origin LIBDIRS), undefined) 
LIBDIRS := $(mcpscusor)$(Segment)LIB$(Segment)sor$(Segment)f
endif
sor  := $(LIBDIRS)$(Segment)$(objname)

ifeq ($(origin LIBDIRD), undefined) 
LIBDIRD := $(mcworkspace)$(Segment)LIB$(Segment)$(ConfigName)
endif

#target directories
tgt  := $(LIBDIRD)

nlist    :=  MiniUtilities
modules  := $(foreach n, $(nlist), $(tgt)$(Segment)$(n).mod)
objects  := $(foreach n, $(nlist), $(tgt)$(Segment)$(n).o)
ffiles   := $(foreach n, $(nlist), $(sor)$(Segment)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(Segment)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(Segment)$(n).F90)

libname  := lib_$(objname).$(LIB_EXT)
##########################################################
$(libname) : $(objects) 
	ar -rcs $(libname) $(objects)
	mv $(libname) $(tgt) 

$(tgt)$(Segment)MiniUtilities.o : $(sor)$(Segment)MiniUtilities.F 
	$(comp) -c $(oflags_this)  -module $(tgt) $< -o $@ 

clean:
	-rm $(objects) $(libname) $(modules) 
