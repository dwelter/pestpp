# This file is part of PEST++
top_builddir = ..
include $(top_builddir)/global.mak

SUBDIRS := run_manager_fortran_test

ifeq ($(SYSTEM),win)
SUBDIRS += linear_analysis_test
endif

all:	$(foreach d,$(SUBDIRS),$(d)-target)

clean:	$(foreach d,$(SUBDIRS),$(d)-clean)

.PHONY: all install clean
