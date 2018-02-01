# This file is part of PEST++
include ../../../global.mak

EXE := concurrency_tester$(EXE_EXT)
OBJECTS := \
    concurrency_tester \
    thread_struct
OBJECTS := $(addsuffix $(OBJ_EXT),$(OBJECTS))


all: $(EXE)

$(EXE): $(OBJECTS)
	$(LD) $(LDFLAGS) $^ $(PESTPP_LIBS) -o $@

clean:
	$(RM) $(OBJECTS) $(EXE)

.PHONY: all test clean
