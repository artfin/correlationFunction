#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <gsl/gsl_spline.h>

// atomic time unit to seconds
const double ATU = 2.418884326505 * std::pow( 10, -17 );

void readFile( std::string const & filename,
			   std::vector<double> & time,
			   std::vector<double> & correlation )
{
	std::ifstream in(filename);
	if ( !in )
		throw std::invalid_argument("Can't open input file");

	const int MAXLINE = 100;
	char line[MAXLINE];

	while ( in.getline(line, MAXLINE) )
	{
		std::stringstream ss;
		ss.str(line);

		double t, corr;
		ss >> t;
		ss >> corr;

		time.push_back( t );
		correlation.push_back( corr );
	}

	in.close();
}

int main()
{
	std::cout << "-------------------------------------------------" << std::endl;
	std::cout << "Program for spline approximation of correlation function." << std::endl;
	std::cout << "Input file: [by default: input.txt]" << std::endl;

	std::string inputFile;

	getline(std::cin, inputFile);
	if ( inputFile == "" )
		inputFile = "input.txt";

	std::cout << "Output file: [by default: output.txt]" << std::endl;

	std::string outputFile;

	getline(std::cin, outputFile);
	if ( outputFile == "" )
		outputFile = "output.txt";

	std::cout << "Approximation step size: [by default: 5 A.T.U.]" << std::endl;

	double step;
	std::string step_;
	getline(std::cin, step_);
	if ( step_ == "" )
		step = 5.0 * ATU;
	else
	{
		std::stringstream ss;
		ss.str(step_);
		ss >> step;	
		step *= ATU;
	}

	std::cout << "Multply time axis by A.T.U.: (y/n) [by default: y]" << std::endl;

	std::string mult;
	getline(std::cin, mult);
	if ( mult == "" )
		mult = "y";

	std::vector<double> time, correlation;
	readFile( inputFile, time, correlation );
	
	if ( mult == "y" )
	{
		for ( size_t i = 0; i < time.size(); ++i )
			time[i] *= ATU;	
	}
		
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc( gsl_interp_cspline, time.size() );
	gsl_spline_init( spline, &time[0], &correlation[0], time.size() );

	std::ofstream out(outputFile);
	for ( double t = time[0]; t <= time.back(); t += step )
		out << t << " " << gsl_spline_eval( spline, t, acc) << std::endl; 

	out.close();

	return 0;
}
