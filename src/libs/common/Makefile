# This file is part of PEST++
include ../../global.mak

LIB := $(LIB_PRE)common$(LIB_EXT)
OBJECTS := \
    fortran_wrappers \
    network_package \
    network_wrapper \
    pest_error \
    system_variables \
    Transformable \
    utilities
OBJECTS := $(addsuffix $(OBJ_EXT),$(OBJECTS))


all: $(LIB)

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

clean:
	$(RM) $(OBJECTS) $(LIB)

.PHONY: all clean
