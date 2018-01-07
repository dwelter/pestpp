# This file is part of PEST++
include ../../../global.mak

LIB := $(LIB_PRE)rm_wrappers$(LIB_EXT)
OBJECTS := \
    RunManagerCWrapper \
    RunManagerFortranWrapper
OBJECTS := $(addsuffix $(OBJ_EXT),$(OBJECTS))


all: $(LIB)

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

clean:
	$(RM) $(OBJECTS) $(LIB)

.PHONY: all clean
