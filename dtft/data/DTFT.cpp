#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*This code works with already symmetrized correlation
function. That means, that first 
(N-1)/2 times are negative, (N-1)/2+1 is zero, and all (N-1)/2 times left are positive.
 Integrate calculates real cosine DTFT*/

const double LSP = 3E10;
const double NL = 2.6867805E+25;
const double hbar = 1.0545718E-34;

int reader(FILE * fp1, double * time, double * corr_val, int N)
{
	for (int i = 0; i < (N-1)/2; ++i)
	{
		fscanf(fp1,"%lg %lg", &(time[i]), &(corr_val[i]));
		//time[i] *= -1;
		//fscanf(fp1, "%lf", &f);
		//fprintf(stderr, "%s\n", );
	} 
	for (int i = (N-1)/2; i < N; ++i)
	{
		fscanf(fp1,"%lg %lg", &(time[i]), &(corr_val[i]));
		//time[i] *= -1;
		//fscanf(fp1, "%lf", &f);
		//fprintf(stderr, "%s\n", );
	} 
	// for (int i = 0; i < N; ++i)
	// {
	// 	fscanf(fp1,"%lg %lg", &(time[i]), &(corr_val[i]));
	// 	//time[i] *= -1;
	// 	//fscanf(fp1, "%lf", &f);
	// 	//fprintf(stderr, "%s\n", );
	// } 
}


double integrand(double time, double func, double omega)
{
	return func*cos(time*LSP*omega*2.0*M_PI);
}

double integrate(double omega, double * time, double * corr_val, int N)
{
	double int_val = 0.0;
	double step = (time[N-1] - time[0])/N;
	std::cout << "step: " << step << std::endl;

	for (int i = 1; i < N-1; i+=2)
	{
		int_val += integrand(time[i-1], corr_val[i-1], omega) + 
					4*integrand(time[i], corr_val[i], omega) +
					integrand(time[i+1], corr_val[i+1], omega);
	}
	std::cout << "integral: " << int_val << std::endl;
	int_val *= step/3.0;
	return int_val/2.0/M_PI;
}



int main(int argc, char const *argv[])
{
	int Npoints = 180*2-1;
	double * time = new double[Npoints];
	double * correlation = new double[Npoints];

	FILE * fp = fopen("danila_sym.txt","r");
	FILE * out = fopen("danila_specfunc.txt","w");
	FILE * spectrum = fopen("danila_spectrum.txt","w");
	reader(fp, time, correlation, Npoints);

	int Nomega = 200;
	double omega0  = 0.0;
	double omegastep = 800.0/Nomega;

	//double specfun[Nomega];
	double spfunval;
	double specval;
	for (int i = 0; i < Nomega; ++i)
	{
		spfunval = integrate(omega0, time, correlation, Npoints);
		std::cout << "omega: " << omega0 << "; spfunval: " << spfunval << std::endl;
		fprintf(out, "%f  %e\n", omega0, spfunval);
		//specfun[i] = spfunval;

		specval = pow(2*M_PI,3.0)*NL*NL/3.0/hbar*omega0*(1-exp(-2*0.7193875*omega0/295.0))*spfunval;
		fprintf(spectrum, "%f  %e\n", omega0, specval);
		omega0 += omegastep;
	}

	fclose(out);
	fclose(fp);
	fclose(spectrum);

	delete [] time;
	delete [] correlation;
	system("pause");
	return 0;
}
