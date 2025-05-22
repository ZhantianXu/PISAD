# Makefile.am for PISAD/snpcalling

# Global compilation and linking flags
# Use Conda-provided environment variables for portability
AM_CXXFLAGS = -I${PREFIX}/include
AM_LDFLAGS = -L${PREFIX}/lib
LIBS = -lboost_system -lboost_filesystem -lboost_graph -lhdf5_cpp -lhdf5 -ltcmalloc

# Define programs and scripts
bin_PROGRAMS = kmer2snp snp_vs_vcf create
bin_SCRIPTS = run.sh 

# kmer2snp compilation rules
kmer2snp_SOURCES = kmer2snp.cpp
kmer2snp_CXXFLAGS = $(AM_CXXFLAGS) -std=c++17 -fopenmp  # Remove -march=native for portability
kmer2snp_CPPFLAGS = -I$(top_srcdir)/vendor
# Inherit AM_LDFLAGS and append Conda LDFLAGS via environment variable
kmer2snp_LDFLAGS = $(AM_LDFLAGS)
kmer2snp_LDADD = $(LIBS)

# snp_vs_vcf compilation rules
snp_vs_vcf_SOURCES = snp_vs_vcf.cpp
snp_vs_vcf_CXXFLAGS = $(AM_CXXFLAGS) -std=c++17
snp_vs_vcf_CPPFLAGS = -I$(top_srcdir)/vendor
snp_vs_vcf_LDFLAGS = $(AM_LDFLAGS)
# Ensure snp_vs_vcf links required libraries
snp_vs_vcf_LDADD = $(LIBS)

# create compilation rules
create_SOURCES = create.cpp
create_CXXFLAGS = $(AM_CXXFLAGS) -std=c++17 -fopenmp  # Remove -march=native for portability
create_CPPFLAGS = -I$(top_srcdir)/vendor
create_LDFLAGS = $(AM_LDFLAGS)
create_LDADD = $(LIBS)