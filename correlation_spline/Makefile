.PHONY: all clean cleandata

all: spline

CCXX := g++
CXXFLAGS := -std=c++11 -O2 -Wall -Wextra -g

GSL := -lgsl -lgslcblas 

spline: spline.cpp
	$(CCXX) $(CXXFLAGS) $(GSL) $< -o $@

clean: 
	rm spline

cleandata:
	rm -f input.txt
	rm -f output.txt
