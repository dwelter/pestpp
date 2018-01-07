# This file is part of PEST++
include ../global.mak

SUBDIRS := \
    common \
    linear_analysis \
    mio \
    pestpp_common \
    pest_routines \
    propack \
    run_managers \
    opt


all:	$(foreach d,$(SUBDIRS),$(d)-target)

clean:	$(foreach d,$(SUBDIRS),$(d)-clean)

.PHONY: all clean
