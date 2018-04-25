.PHONY: clean all

all: simple_sampling

CCXX := g++
CXXFLAGS := -std=c++11 -O2 -Wall -Wextra -g

ar_he_pes.o: ./ar_he_pes.cpp
	$(CCXX) $(CXXFLAGS) -c $< -o $@

ar_he_dip_buryak_fit.o: ./ar_he_dip_buryak_fit.cpp
	$(CCXX) $(CXXFLAGS) -c $< -o $@

simple_sampling.o: ./simple_sampling.cpp
	$(CCXX) $(CXXFLAGS) -c $< -o $@

simple_sampling: simple_sampling.o ar_he_pes.o ar_he_dip_buryak_fit.o
	$(CCXX) $(CXXFLAGS) $^ -o $@


