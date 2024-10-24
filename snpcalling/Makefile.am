# Makefile.am for PISAD/snpcalling

# Global compilation and linking flags
AM_CXXFLAGS = -g -gdwarf-4 -O3 -I$(CONDA_PREFIX)/include
AM_LDFLAGS = -L $(CONDA_PREFIX)/lib
LIBS = -lboost_system -lboost_filesystem -lboost_graph -lhdf5_cpp -lhdf5 -ltcmalloc

# Define programs and scripts
bin_PROGRAMS = kmer2snp snp_vs_vcf create
bin_SCRIPTS = run.sh 

# kmer2snp compilation rules
kmer2snp_SOURCES = kmer2snp.cpp
kmer2snp_CXXFLAGS = $(AM_CXXFLAGS) -std=c++17 -march=native -fopenmp  # Inherit global AM_CXXFLAGS
kmer2snp_CPPFLAGS = -I$(top_srcdir)/vendor  # Add specific include paths
kmer2snp_LDFLAGS = $(AM_LDFLAGS) -Wl,-rpath=$(CONDA_PREFIX)/lib # Inherit global AM_LDFLAGS
kmer2snp_LDADD = $(LIBS)

# snp_vs_vcf compilation rules
snp_vs_vcf_SOURCES = snp_vs_vcf.cpp
snp_vs_vcf_CXXFLAGS = $(AM_CXXFLAGS) -std=c++11  # Inherit global AM_CXXFLAGS
snp_vs_vcf_CPPFLAGS = -I$(top_srcdir)/vendor 
snp_vs_vcf_LDFLAGS = $(AM_LDFLAGS)  # Inherit global AM_LDFLAGS
snp_vs_vcf_LDADD = 

# create compilation rules
create_SOURCES = create.cpp
create_CXXFLAGS = $(AM_CXXFLAGS) -std=c++17 -march=native -fopenmp  # Inherit global AM_CXXFLAGS
create_CPPFLAGS = -I$(top_srcdir)/vendor 
create_LDFLAGS = $(AM_LDFLAGS) -Wl,-rpath=$(CONDA_PREFIX)/lib  # Inherit global AM_LDFLAGS
create_LDADD = $(LIBS)
