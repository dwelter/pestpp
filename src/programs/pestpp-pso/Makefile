# This file is part of PEST++
top_builddir = ../..
include $(top_builddir)/global.mak

EXE := psopp$(EXE_EXT)
OBJECTS := \
    bestsubs \
    estimation \
    objsubs \
    pareto \
    paretosubs \
    parsubs \
    pertsubs \
    psomain \
    readsubs \
    rmsubs \
    utilsubs \
    writesubs
OBJECTS	:= $(addsuffix $(OBJ_EXT),$(OBJECTS))


all: $(EXE)

# Compile readsubs.f90 first, to get psodat.mod used by everything else
$(filter-out readsubs$(OBJ_EXT),$(OBJECTS)): readsubs$(OBJ_EXT)

$(EXE): $(OBJECTS)
	$(LD) $(LDFLAGS) $^ $(PESTPP_LIBS) -o $@

install: $(EXE)
	$(CP) $< $(bindir)

clean:
	$(RM) $(OBJECTS) $(EXE) *.mod

.PHONY: all install clean
