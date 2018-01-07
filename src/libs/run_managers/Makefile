# This file is part of PEST++
include ../../global.mak

SUBDIRS := \
    abstract_base \
    external \
    genie_wrapper \
    serial \
    wrappers \
    yamr

ifeq ($(SYSTEM),win)
SUBDIRS += genie
endif


all:	$(foreach d,$(SUBDIRS),$(d)-target)

clean:	$(foreach d,$(SUBDIRS),$(d)-clean)

.PHONY: all clean
