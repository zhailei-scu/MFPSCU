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
#sorce directories

ifeq ($(origin Ccomp), undefined) 
  Ccomp   := g++
  ifeq ($(operateSystem), CYGWIN)
    Ccomp   := cl
  endif
  
  ifeq ($(operateSystem), LINUX)
    Ccomp   := g++
  endif
endif

ifeq ($(origin ConfigName), undefined) 
  ConfigName := Release
endif

CP := -c
OutF := -o
IncF := -I
ifeq ($(origin Coflags), undefined)
  ifeq ($(operateSystem), CYGWIN)
      ifeq ($(ConfigName), Release)
        Coflags := /Wall /O2
      endif

      ifeq ($(ConfigName), Debug)
        Coflags := /Wall /Zi
      endif

      CP := /c
      OutF := /Fo:
      IncF := /I	
  endif

  ifeq ($(operateSystem), LINUX)
      ifeq ($(ConfigName), Release)
        Coflags := -Wall -Mconcur -O2 -Minform=warn -Minfo=all
      endif

      ifeq ($(ConfigName), Debug)
        Coflags := -Wall -g -O0 -Minform=warn -Minfo=all
      endif

      CP := -c
      OutF := -o
      IncF := -I
  endif

endif

###########################################################
#sorce dir name
objname := TCCLIB

ifeq ($(origin TCCINTERDIRS), undefined) 
  TCCINTERDIRS := $(mcpscusor)$(Segment)TCCLIB$(Segment)sor
endif
sor  := $(TCCINTERDIRS)

#include dir
ifeq ($(origin LIBDIRD), undefined)
LIBDIRD := $(mcworkspace)$(Segment)LIB$(Segment)$(ConfigName)
endif
incdir := $(LIBDIRD)

#target directories
tgt  := $(LIBDIRD)

#target lib name
libname  := lib_$(objname).$(LIB_EXT)

#the Macro

ifeq ($(operateSystem), CYGWIN)
  Macro += /DLIBTCCPATH=\"$(tccpath)$(Segment)win32\" /DLIBTCCHEADERF=\"$(tccpath)$(Segment)win32$(Segment)libtcc$(Segment)libtcc.h\" /DDIFFUSORDEFPATH=\"$(TCCINTERDIRS)\" /D$(operateSystem)
endif

ifeq ($(operateSystem), LINUX)
  Macro += -DLIBTCCPATH=\"$(tccpath)\" -DLIBTCCHEADERF=\"$(tccpath)/libtcc.h\" -DDIFFUSORDEFPATH=\"$(TCCINTERDIRS)\" -D$(operateSystem)
endif

#######################################################          
nlist    :=  DiffusorList	\
             ReactionsList	\
	     TccInterp
             
objects  := $(foreach n, $(nlist), $(tgt)$(Segment)$(n).o)
cfiles   := $(foreach n, $(nlist), $(sor)$(Segment)$(n).cpp)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects)
	mv $(libname) $(tgt)

$(tgt)$(Segment)DiffusorList.o : $(sor)$(Segment)DiffusorList.cpp
	$(Ccomp) $(CP) $< $(Coflags) $(IncF)$(incdir) $(OutF) $@

$(tgt)$(Segment)ReactionsList.o : $(sor)$(Segment)ReactionsList.cpp
	$(Ccomp) $(CP) $< $(Coflags) $(IncF)$(incdir) $(OutF) $@

$(tgt)$(Segment)TccInterp.o : $(sor)$(Segment)TccInterp.cpp
	$(Ccomp) $(CP) $< $(Coflags) $(Macro) $(IncF)$(incdir) $(OutF) $@

######################################################################
clean:
	-rm $(objects) $(libname) 
