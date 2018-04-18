#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <sstream>

namespace constants
{
    // vacuum permittivity
    const double EPSILON0 = 8.854187817 * 1E-12;

    // atomic unit of dipole moment == e * a0
    const double ADIPMOMU = 8.478353552 * 1E-30;
}

void reader( std::ifstream & inFile, std::vector<std::pair<double, double>> & contents )
{
    if ( !inFile )
        std::invalid_argument("Can't open the file");

    const int MAXLINE = 100;
    char buffer[MAXLINE];
    std::stringstream ss;
    double time, corr;

    while( inFile.getline(buffer, MAXLINE ))
    {
        ss << buffer;
        ss >> time >> corr;
        contents.emplace_back(time, corr);
        ss.clear();
        ss.str("");
    }
}

// симметризация корреляционной функции и частичное преобразование размерности
void symmetrize( std::vector<std::pair<double, double>> const & input, std::vector<std::pair<double, double>> & res )
{
    res.resize( 2 * input.size() - 1 );
    std::cout << "input.size(): " << input.size() << std::endl;
    std::cout << "res.size(): " << res.size() << std::endl;

    int j = 0;
    for ( int i = (int) input.size() - 1; i >= 0; --i, ++j )
        res[j] = std::make_pair( -input[i].first, constants::ADIPMOMU * constants::ADIPMOMU / \
                                 (4.0 * M_PI * constants::EPSILON0) * input[i].second );

    for ( int i = 1; i < (int) input.size(); ++i, ++j )
        res[j] = std::make_pair( input[i].first, constants::ADIPMOMU * constants::ADIPMOMU / \
                                 (4.0 * M_PI * constants::EPSILON0) * input[i].second );
}

int main()
{
	std::cout << "---------------------------------------------" << std::endl;
	std::cout << "Symmetrizer started." << std::endl;
	std::cout << "Excepted input: one-sided non-normalized correlation function multiplied by Volume" << std::endl;
	std::cout << "Opening input file: input.txt" << std::endl; 

	std::vector<std::pair<double, double>> input;
	std::ifstream inFile("./input.txt");
	reader( inFile, input );

	std::vector<std::pair<double, double>> symm_input;
	symmetrize( input, symm_input );

	std::cout << "One-sided correlation function symmetrized!" << std::endl;
	std::cout << "Saving result to: symm.txt" << std::endl;
	
	std::ofstream outFile("./symm.txt");
	for ( size_t k = 0; k < symm_input.size(); ++k )
		outFile << symm_input[k].first << " " << symm_input[k].second << std::endl;
	outFile.close();
	
	std::cout << "Proceed with DTFT procedure" << std::endl;
	std::cout << "---------------------------------------------" << std::endl;

	return 0;
}
