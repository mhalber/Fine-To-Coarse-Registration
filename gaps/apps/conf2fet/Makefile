#
# Application name and list of source files.
#

NAME=conf2fet
CCSRCS=$(NAME).cpp 



#
# With Ceres
# 

#USER_CFLAGS=-I/usr/local/suitesparse -I/usr/local/include/eigen3 -DRN_USE_CERES -DRN_USE_SPLM -DRN_USE_CSPARSE
#USER_LIBS=-L/usr/local/lib -lceres -lsplm -lCSparse -lminpack -lccolamd -lcamd -lcolamd -lamd -Wl,-Bdynamic -lcxsparse -lcholmod -lccolamd -lcamd -lcolamd -lamd -lsuitesparseconfig -llapack -lblas -Wl,-Bstatic -lcxsparse -Wl,-Bdynamic -lgomp -lglog -lgflags 


#
# Without Ceres
# 

USER_CFLAGS=-DRN_USE_SPLM -DRN_USE_CSPARSE
USER_LIBS=-lsplm -lCSparse



#
# Pkg libraries
#

PKG_LIBS=-lFET -lR3Surfels -lR3Shapes -lR2Shapes -lRNMath -lRNBasics -ljpeg -lpng


#
# R3 application makefile
#

include ../../makefiles/Makefile.apps


