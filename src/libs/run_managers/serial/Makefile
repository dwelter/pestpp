# This file is part of PEST++
include ../../../global.mak

LIB := $(LIB_PRE)rm_serial$(LIB_EXT)
OBJECTS := RunManagerSerial$(OBJ_EXT)


all: $(LIB)

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

clean:
	$(RM) $(OBJECTS) $(LIB)

.PHONY: all clean
