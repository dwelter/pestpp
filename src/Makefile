# Top-level GNU Make file for PEST++

# Customize options in global.mak
include ./global.mak


all: libs-target programs-target utilities-target

libs-target:
	(cd libs; $(MAKE))

programs-target:
	(cd programs; $(MAKE))

utilities-target:
	(cd utilities; $(MAKE))

install:
	(cd programs; $(MAKE) install)
	(cd utilities; $(MAKE) install)

clean:
	(cd programs; $(MAKE) clean)
	(cd libs; $(MAKE) clean)
	(cd utilities; $(MAKE) clean)

.PHONY: all libs-target programs-target utilities-target install clean
