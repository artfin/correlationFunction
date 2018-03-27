#pragma once

#include "parameters.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

// throwing exception
#include <stdexcept>

class FileReader
{
private:
    Parameters* parameters;
    const int MAXLINE = 100;

public:
    FileReader( std::string filename, Parameters* parameters );
    ~FileReader();

   void parse_file( std::ifstream& infile );

   double string_to_double( std::string const & value, int & line );
   int string_to_int( std::string const & value, int & line );
   bool string_to_bool( std::string const & value, int & line );
   std::vector<double> string_to_vector( std::string const & value, int & line );

   void parse_string( std::string curr_str,
                      std::string& keyword,
                      std::string& variable,
                      std::string& value,
                      bool& is_doolar,
                      bool& is_empty,
                      bool& is_assignment,
                      int& line
                    );

   void analyse_grid_group_line( std::string const & variable,
                                 std::string const & value,
                                 int& line
                                );

   void analyse_mcparameters_group_line( std::string const & variable,
                                         std::string const & value,
                                         int& line
                                        );

   void analyse_files_group_line( std::string const & variable,
                                  std::string const & value,
                                  int& line
                                );

   void analyse_trajectory_group_line( std::string const & variable,
                                       std::string const & value,
                                       int& line
                                      );

   void analyse_conditions_group_line( std::string const & variable,
                                       std::string const & value,
                                       int& line
                                     );
};

