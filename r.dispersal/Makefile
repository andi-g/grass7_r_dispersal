#MODULE_TOPDIR = /home/wegmann/software/grass/grass6/

#MODULE_TOPDIR = ../grass7_trunk
MODULE_TOPDIR = ../..

PGM = r.dispersal

#CXX="clang++ -std=c++11 -stdlib=libc++"
#LIBES = $(SITESLIB) $(VECTLIB) $(GISLIB) $(SEGMENTLIB) $(DATETIMELIB)
LIBES = $(RASTERLIB) $(GISLIB) $(MATHLIB)
#DEPENDENCIES = $(VECTDEP) $(GISDEP) $(SEGMENTDEP) $(DATETIMEDEP)
DEPENDENCIES = $(RASTERDEP) $(GISDEP)  
EXTRA_INC = $(VECT_INC)
EXTRA_CFLAGS = $(VECT_CFLAGS) -H
EXTRA_LDFLAGS = -L/usr/lib -Wl -t 

include $(MODULE_TOPDIR)/include/Make/Module.make

default: cmd
