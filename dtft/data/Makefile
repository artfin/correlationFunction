.PHONY: all clean cleandata

CCXX = g++

CXXFLAGS = -std=c++11 -O2 -Wall -Wextra

dtft: DTFT.cpp
	$(CCXX) $(CXXFLAGS) $< -o $@

symm: symm.cpp
	$(CCXX) $(CXXFLAGS) $< -o $@

all: dtft symm

clean:
	@echo "Cleaning directory from binaries..."
	@rm -f dtft
	@rm -f symm

cleandata:
	rm -f input.txt
	rm -f symm.txt
	rm -f spectrum.txt
	rm -f specfunc.txt
