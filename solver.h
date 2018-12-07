#include <stdio.h>
#include <gsl/gsl_linalg.h>

using namespace std;

void solve(double a_data[], double b_data[], double x_sol[], const int dim)
{
	//double a_data[] = { 0.18, 0.60, 0.57, 0.96,
	//					0.41, 0.24, 0.99, 0.58,
	//					0.14, 0.30, 0.97, 0.66,
	//					0.51, 0.13, 0.19, 0.85 };
	//
	//double b_data[] = { 1.0, 2.0, 3.0, 4.0 };

	gsl_matrix_view m 
		= gsl_matrix_view_array (a_data, dim, dim);

	gsl_vector_view b
		= gsl_vector_view_array (b_data, dim);

	gsl_vector *x = gsl_vector_alloc (dim);
  
	int s;

	gsl_permutation * p = gsl_permutation_alloc (dim);

	gsl_linalg_LU_decomp (&m.matrix, p, &s);

	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

	//printf ("x = \n");
	//gsl_vector_fprintf (stdout, x, "%g");
	
	for (int i = 0; i < dim; ++i)
		x_sol[i] = gsl_vector_get (x, i);

	gsl_permutation_free (p);
	gsl_vector_free (x);

	return;
}

