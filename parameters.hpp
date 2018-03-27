#pragma once
#include <vector>
#include <string>

class Parameters
{
public:
    Parameters () { }
    ~Parameters() { }

    double B_MIN;
    double B_MAX;
    int B_PARTS;

    double V0_MIN;
    double V0_MAX;
    int V0_PARTS;

    int DIM;
    int NPOINTS;
    double alpha;
    int subchain_length;
    std::vector<double> initial_point;
    double gunsight_upper_bound;

    double RDIST;
    double sampling_time;
    int MaxTrajectoryLength;
    double FREQ_MAX;

    double Temperature;

    std::string output_directory = "";
    std::string specfunc_filename = "";
    std::string spectrum_filename = "";
    std::string m2_filename = "";
};
