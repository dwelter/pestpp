# This file is part of PEST++
include ../../global.mak

EXE := gsa$(EXE_EXT)
OBJECTS := \
    main \
    GsaAbstractBase \
    MorrisMethod \
    PooledVariance \
    sobol \
    Stats \
    TornadoPlot
OBJECTS := $(addsuffix $(OBJ_EXT),$(OBJECTS))


all: $(EXE)

$(EXE): $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $^ $(PESTPP_LIBS)

install: $(EXE)
	$(MKDIR) $(bindir)
	$(CP) $< $(bindir)

clean:
	$(RM) $(OBJECTS) $(EXE)

.PHONY: all install clean
