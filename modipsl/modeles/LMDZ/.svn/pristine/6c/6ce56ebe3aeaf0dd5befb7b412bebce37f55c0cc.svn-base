# This is the makefile for "Max_diff_nc" with its libraries. For GNU
# make.

SUBDIRS = Max_diff_nc NetCDF95 NR_util Jumble

.PHONY: all clean $(SUBDIRS)

all clean: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ ${MAKECMDGOALS}

Max_diff_nc: NetCDF95 Jumble
Jumble: NR_util
