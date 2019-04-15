#compiler
ifeq ($(origin comp), undefined)
comp       := pgfortran
endif

ifeq ($(origin oflags), undefined)
oflags := -fast -tp sandybridge-64 -Minline -Mconcur -Minform=warn -Minfo=all
endif

##########################################################
#sorce dir name
objname := MiniUtilities

#sorce directories
ifeq ($(origin LIBDIRS), undefined) 
LIBDIRS := $(MCPSCUSOR)/LIB/sor/f/
endif
sor  := $(LIBDIRS)$(objname)/

#target directories
ifeq ($(origin ConfigName), undefined) 
ConfigName := Release
endif

ifeq ($(origin LIBDIRD), undefined) 
LIBDIRD := $(MCWORKSPACE)/LIB/$(ConfigName)/
endif
tgt  := $(LIBDIRD)

nlist    :=  MiniUtilities
modules  := $(foreach n, $(nlist), $(tgt)$(n).mod)
objects  := $(foreach n, $(nlist), $(tgt)$(n).o)
ffiles   := $(foreach n, $(nlist), $(sor)$(n).f)
Ffiles   := $(foreach n, $(nlist), $(sor)$(n).F)
F90files := $(foreach n, $(nlist), $(sor)$(n).F90)

libname  := $(tgt)lib_$(objname).a
##########################################################
$(libname) : $(objects) 
	ar -rcs $(libname) $(objects) 

$(tgt)MiniUtilities.o : $(sor)MiniUtilities.F 
	$(comp) -c $(oflags)  -module $(tgt) $< -o $@ 

clean:
	-rm $(objects) $(libname) $(modules) 
