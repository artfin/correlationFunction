TARGET := monte_carlo.exe

all: $(TARGET) clean
	@echo Compilation finished

CXXFLAGS := -O1 -std=gnu++11 
FORFLAGS := -O1

INCLUDEHEP := -I"hep-mc-0.5/include"
INCLUDEMPI := -I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include" 
LIBMPI := -L"C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86"
LIBS := -lgfortran -lmsmpi -lgsl -lgslcblas

#SRC_C = main_mudt_transformed.cpp co2co2_hamiltonian_deriv.cpp co2co2_dipole_lr.cpp constant.cpp n2h2_dip_derivatives.cpp   
#SRC_F = co2-h2-pes.f
OBJ = monte_carlo.o ar_he_pes.o ar_he_dip_buryak_fit.o basis_r.o awp.o vmblock.o fgauss.o     gear.o  ar_he_dip_buryak_fit.o \
 			   diatomic_equations_2D.o ar_he_pes_der.o  

$(TARGET): $(OBJ) 
	g++ $(CXXFLAGS) $^ -o $@  $(LIBMPI) $(LIBS)

 
monte_carlo.o: monte_carlo.cpp
	g++ $(CXXFLAGS) -c $^ $(INCLUDEMPI) $(INCLUDEHEP)

ar_he_dip_buryak_fit.o: ar_he_dip_buryak_fit.cpp
	g++ $(CXXFLAGS) -c $^

ar_he_pes.o: ar_he_pes.cpp
	g++ $(CXXFLAGS) -c $^

awp.o: 
	g++ $(CXXFLAGS) -c ../jean-pierre/awp.cpp

basis_r.o: 
	g++ $(CXXFLAGS) -c  ../jean-pierre/basis_r.cpp 

gear.o: 
	g++ $(CXXFLAGS) -c  ../jean-pierre/gear.cpp

t_dgls.o: 
	g++ $(CXXFLAGS) -c  ../jean-pierre/t_dgls.cpp

vmblock.o: 
	g++ $(CXXFLAGS) -c  ../jean-pierre/vmblock.cpp

fgauss.o: 
	g++ $(CXXFLAGS) -c  ../jean-pierre/fgauss.cpp


ar_he_pes_der.o: 
	g++ $(CXXFLAGS) -c  ar_he_pes_der.cpp

diatomic_equations_2D.o:  
	g++ $(CXXFLAGS) -c diatomic_equations_2D.cpp

clean:
	rm *.o
