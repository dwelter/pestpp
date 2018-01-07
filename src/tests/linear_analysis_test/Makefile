# This file is part of PEST++
include ../../global.mak

EXE := linear_analysis_test$(EXE_EXT)
OBJECTS := \
    linear_analysis_test \
    stdafx
OBJECTS := $(addsuffix $(OBJ_EXT),$(OBJECTS))


all: $(EXE)

$(EXE): $(OBJECTS)
	$(LD) $(LDFLAGS) $^ $(PESTPP_LIBS) -o $@

clean:
	$(RM) $(OBJECTS) $(EXE)

.PHONY: all test clean
