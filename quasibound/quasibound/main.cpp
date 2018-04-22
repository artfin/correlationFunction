#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>

#include "ar_he_pot.hpp"
#include "ar_he_pot_derivative.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_spline.h>

#define SHOW_ITERATIONS 0

const double MU = 6632.09;

// positive effective potential
double effective_potential_pos( const gsl_vector *v, void * params )
{
    double x = gsl_vector_get(v, 0);
    double J = *(double*) params;
    return ar_he_pot(x) + std::pow(J, 2.0) / (2.0 * MU * std::pow(x, 2.0));
}

// negative effective potential
double effective_potential_neg( const gsl_vector *v, void * params )
{
    double x = gsl_vector_get(v, 0);
    double J = *(double*) params;
    return - ar_he_pot(x) - std::pow(J, 2.0) / (2.0 * MU * std::pow(x, 2.0));
}

// derivative of negative effective potential
void effective_potential_neg_der( const gsl_vector *v, void *params, gsl_vector *df )
{
    double x = gsl_vector_get(v, 0);
    double J = *(double*) params;

    gsl_vector_set(df, 0, - ar_he_pot_derivative(x) + std::pow(J, 2.0) / (MU * std::pow(x, 3.0)));
}

// derivative of positive effective potential
void effective_potential_pos_der( const gsl_vector *v, void *params, gsl_vector *df )
{
    double x = gsl_vector_get(v, 0);
    double J = *(double*) params;

    gsl_vector_set(df, 0, ar_he_pot_derivative(x) - std::pow(J, 2.0) / (MU * std::pow(x, 3.0)));
}

// function to compute both f and df
void fdf_neg( const gsl_vector *x, void * params, double *f, gsl_vector *df )
{
    *f = effective_potential_neg(x, params);
    effective_potential_neg_der(x, params, df);
}

void fdf_pos( const gsl_vector *x, void * params, double *f, gsl_vector *df )
{
    *f = effective_potential_pos(x, params);
    effective_potential_pos_der(x, params, df);
}

void write_effective_potential( std::string const & filename, double J, std::vector<double> const & grid = std::vector<double>() )
{
    std::ofstream out(filename);

    if ( grid.size() == 0 )
    {
        for ( double xmin = 3.0; xmin < 25.0; xmin += 0.1 )
        {
            gsl_vector * v = gsl_vector_alloc(1);
            gsl_vector_set(v, 0, xmin);
            out << xmin << " " << effective_potential_pos(v, &J) << std::endl;
            gsl_vector_free(v);
        }
    }
    else
    {
        for ( size_t k = 0; k < grid.size(); ++k )
        {
            gsl_vector * v = gsl_vector_alloc(1);
            gsl_vector_set(v, 0, grid[k]);
            out << grid[k] << " " << effective_potential_pos(v, &J) << std::endl;
            gsl_vector_free(v);
        }
   }

   out.close();
}

double optimize( double J, double a, double b, std::string const& type )
{
    int status;
    size_t iter = 0;
    size_t maxiter = 100;

    double par[1] = { J }; // list of parameters for minimization problem

    gsl_multimin_function_fdf my_func;

    my_func.n = 1; // dimension

    if ( type == "negative" )
    {
        my_func.f = effective_potential_neg;
        my_func.df = effective_potential_neg_der;
        my_func.fdf = fdf_neg;
    }

    if ( type == "positive" )
    {
        my_func.f = effective_potential_pos;
        my_func.df = effective_potential_pos_der;
        my_func.fdf = fdf_pos;
    }

    my_func.params = par;

    gsl_vector *x = gsl_vector_alloc(1); // initial guess
    gsl_vector_set(x, 0, 0.5 * (a + b));

    const gsl_multimin_fdfminimizer_type *T;
    T = gsl_multimin_fdfminimizer_conjugate_fr;

    gsl_multimin_fdfminimizer * s;
    s = gsl_multimin_fdfminimizer_alloc(T, 1);

    if ( type == "negative" )
        gsl_multimin_fdfminimizer_set(s, &my_func, x, 0.01, 1e-6);
    else
        gsl_multimin_fdfminimizer_set(s, &my_func, x, 0.01, 1e-6);

    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);

        if ( status )
            break;

        status = gsl_multimin_test_gradient(s->gradient, 1e-9);

        //if ( status == GSL_SUCCESS )
        //    std::cout << "Minimum found at: " << std::endl;

        //std::cout << "iter: " << iter << " " << gsl_vector_get(s->x, 0) << " " << s->f << std::endl;
    } while ( status == GSL_CONTINUE && iter < maxiter );

    double res = gsl_vector_get(s->x, 0);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);

    return res;
}

int main()
{

    double jmin = 1.0;
    double jmax = 9.9;
    double jstep = 0.1;

    std::vector<double> j_values;
    std::vector<double> xmax_values;

    double maximum;

    for ( double j = jmin; j <= jmax; j += jstep )
    {
        j_values.push_back( j );

        double minimum = optimize(j, 6.0, 10.0, "positive");
        std::cout << "j: " << j << "; min: " << minimum;

        maximum = optimize(j, minimum, 30.0, "negative");
        std::cout << "; max; " << maximum << std::endl;

        xmax_values.push_back( maximum );

        write_effective_potential("../pot/effpot_j=" + std::to_string(j) + ".txt", j);
    }

    std::ofstream out("../result.txt");
    out << std::fixed << std::setprecision(8);
    for ( size_t k = 0; k < xmax_values.size(); ++k )
        out << j_values[k] << " " << xmax_values[k] << std::endl;
    out.close();

    // spline interpolation
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc( gsl_interp_cspline, xmax_values.size() );

    gsl_spline_init(spline, &j_values[0], &xmax_values[0], xmax_values.size() );

    std::ofstream spline_out("../spline.txt");
    spline_out << std::fixed << std::setprecision(8);

    double j_step = 1e-3;
    for ( double j = j_values[0]; j < j_values.back(); j += j_step )
        spline_out << j << " " << gsl_spline_eval(spline, j, acc) << std::endl;
    spline_out.close();

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return 0;
}
