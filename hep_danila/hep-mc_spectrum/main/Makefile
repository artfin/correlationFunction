TARGET := main 

CCXX := g++
MPICCXX := mpic++

all: $(TARGET) 
	@echo Compilation finished

CXXFLAGS := -O2 -std=gnu++11 
FORFLAGS := -O1

INCLUDEHEP := -I"/home/artfin/Downloads/hep-mc-0.6/include"
#LIBS := -lgfortran -lmsmpi -lgsl -lgslcblas
LIBS := -lgfortran

OBJ = monte_carlo.o ar_he_pes.o ar_he_dip_buryak_fit.o basis_r.o awp.o vmblock.o fgauss.o     gear.o  ar_he_dip_buryak_fit.o \
 			   diatomic_equations_2D.o ar_he_pes_der.o  

$(TARGET): $(OBJ) 
	$(MPICCXX) $(CXXFLAGS) $^ -o $@  $(LIBMPI) $(LIBS)
 
monte_carlo.o: monte_carlo.cpp
	$(MPICCXX) $(CXXFLAGS) -c $^ $(INCLUDEMPI) $(INCLUDEHEP)

ar_he_dip_buryak_fit.o: ar_he_dip_buryak_fit.cpp
	$(CCXX) $(CXXFLAGS) -c $^

ar_he_pes.o: ar_he_pes.cpp
	$(CCXX) $(CXXFLAGS) -c $^

awp.o: 
	$(CCXX) $(CXXFLAGS) -c ../jean-pierre/awp.cpp

basis_r.o: 
	$(CCXX) $(CXXFLAGS) -c  ../jean-pierre/basis_r.cpp 

gear.o: 
	$(CCXX) $(CXXFLAGS) -c  ../jean-pierre/gear.cpp

t_dgls.o: 
	$(CCXX) $(CXXFLAGS) -c  ../jean-pierre/t_dgls.cpp

vmblock.o: 
	$(CCXX) $(CXXFLAGS) -c  ../jean-pierre/vmblock.cpp

fgauss.o: 
	$(CCXX) $(CXXFLAGS) -c  ../jean-pierre/fgauss.cpp

ar_he_pes_der.o: 
	$(CCXX) $(CXXFLAGS) -c  ar_he_pes_der.cpp

diatomic_equations_2D.o:  
	$(CCXX) $(CXXFLAGS) -c diatomic_equations_2D.cpp

clean:
	rm *.o
