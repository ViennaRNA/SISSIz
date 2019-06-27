AM_CPPFLAGS = \
    -I$(top_srcdir)/librna \
    -I$(top_srcdir)/sissiz \
    -I$(top_srcdir)/levmar \
    -I$(top_srcdir)

bin_PROGRAMS = SISSIz

SISSIz_SOURCES = \
    SISSIz.c \
    rnaz_utils.c \
    expfit.c \
    bionj_custom.c \
    bionj_custom.h \
    sissiz.h \
    rnaz_utils.h \
    expfit.h \
    cmdline_sissiz.c \
    cmdline_sissiz.h \
    treeML.h \
    treeML.c

LDADD = \
    -L$(top_builddir)/librna -lRNA \
    -lm \
    -L$(top_builddir)/levmar -llevmar \
    -L$(top_builddir)/sissiz -lSISSI \
    -L$(top_builddir)/phyml -lPHYML

