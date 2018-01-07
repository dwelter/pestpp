# This file is part of PEST++
include ../global.mak

SUBDIRS := \
    gsa \
    pest++ \
    pestpp-pso

# Does not seem to compile
# SUBDIRS += pestpp-ies

all:	$(foreach d,$(SUBDIRS),$(d)-target)

install:	$(foreach d,$(SUBDIRS),$(d)-install)

clean:	$(foreach d,$(SUBDIRS),$(d)-clean)

.PHONY: all install clean
