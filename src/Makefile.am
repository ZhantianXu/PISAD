bin_PROGRAMS = pisadCount pisadEval

# Global compilation and linking flags
AM_CXXFLAGS = -std=c++17 -fopenmp
AM_LDFLAGS = -fopenmp  # Ensure OpenMP is enabled during linking

# pisadCount compilation rules
pisadCount_CXXFLAGS = $(AM_CXXFLAGS)  # Inherit global AM_CXXFLAGS
pisadCount_CPPFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/vendor \
	-I$(top_srcdir)
pisadCount_LDFLAGS = $(AM_LDFLAGS)  # Inherit global AM_LDFLAGS
pisadCount_LDADD = -lz  # Add zlib library for compression support
pisadCount_SOURCES = ntSeqMatchCount.cpp \
	FingerPrint.hpp \
	Options.cpp Options.h \
	Util.h \
	FastxParser.cpp

# pisadEval compilation rules
pisadEval_CXXFLAGS = $(AM_CXXFLAGS)  # Inherit global AM_CXXFLAGS
pisadEval_CPPFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/vendor \
	-I$(top_srcdir)
pisadEval_LDFLAGS = $(AM_LDFLAGS)  # Inherit global AM_LDFLAGS
pisadEval_SOURCES = ntSeqMatchEval.cpp \
	CompareCounts.hpp \
	Options.cpp Options.h \
	KDTreeUtil.h \
	Util.h
