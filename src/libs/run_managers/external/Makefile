# This file is part of PEST++
include ../../../global.mak

LIB := $(LIB_PRE)rm_external$(LIB_EXT)
OBJECTS := RunManagerExternal$(OBJ_EXT)


all: $(LIB)

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

clean:
	$(RM) $(OBJECTS) $(LIB)

.PHONY: all clean
