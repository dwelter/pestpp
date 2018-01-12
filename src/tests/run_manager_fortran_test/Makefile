# This file is part of PEST++
top_builddir = ../..
include $(top_builddir)/global.mak

EXE := fortran_test$(EXE_EXT)
OBJECTS := run_manager_fortran_test$(OBJ_EXT)


all: $(EXE)

$(EXE): $(OBJECTS)
	$(LD) $(LDFLAGS) $^ $(PESTPP_LIBS) -o $@

clean:
	$(RM) $(OBJECTS) $(EXE) *.mod

.PHONY: all test clean
