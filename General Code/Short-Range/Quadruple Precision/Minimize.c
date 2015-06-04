#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>


double ederivwrapper_(int *Omega, int *NumTerms, int *qmax, int *pmax, int *IsTriplet, int *Ordering, int *Method, int *EigenRoutine, double *Alpha, double *Beta, double *Gamma, int *Iter, double *Energy);


void testfortran(int *x)
{
	printf("%d\n", *x);
}


double test(const gsl_vector *v, void *params)
{
	double a1, b1, g1, Energy;
	int *p = (int*)params;
	int Omega = p[0], NumTerms = p[1], qmax = p[2], pmax = p[3], IsTriplet = p[4], Ordering = p[5], Method = p[6], Iter = p[7], EigenRoutine = p[8];

	a1 = gsl_vector_get(v, 0);
	b1 = gsl_vector_get(v, 1);
	g1 = gsl_vector_get(v, 2);

	//@TODO: REMOVE THIS!
	//printf("Parameters: %d %d %d %d %d %d\n", p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
	
	//@TODO: Make this into a subroutine?
	ederivwrapper_(&Omega, &NumTerms, &qmax, &pmax, &IsTriplet, &Ordering, &Method, &EigenRoutine, &a1, &b1, &g1, &Iter, &Energy);
	//printf("Test: %f %f %f\n", a1, b1, g1);
	printf("Energy: %f\n", Energy);

	return Energy;
}


int optimizewavefn(int MaxIter, int Omega, int NumTerms, int qmax, int pmax, int IsTriplet, int Ordering, int Method, int EigenRoutine, double alpha, double beta, double gamma)
{
	int par[9] = { Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, 0, EigenRoutine };
	//printf("Test: %f %f %f\n", alpha, beta, gamma);

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;  //gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	 
	size_t iter = 0;
	int status;
	double size;
	
	x = gsl_vector_alloc(3);
	gsl_vector_set(x, 0, alpha);
	gsl_vector_set(x, 1, beta);
	gsl_vector_set(x, 2, gamma);
	 
	/* Set initial step sizes to 1 */
	//ss = gsl_vector_alloc (2);
	ss = gsl_vector_alloc(3);
	gsl_vector_set_all(ss, 0.05);
	 
	/* Initialize method and iterate */
	//minex_func.n = 2;
	minex_func.n = 3;
	minex_func.f = (void*)test;
	minex_func.params = par;
	 
	//s = gsl_multimin_fminimizer_alloc (T, 2);
	s = gsl_multimin_fminimizer_alloc(T, 3);
	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
	 
	do
	{
		iter++;
		par[7] = iter;  // Lets us keep track internally of the iteration number
		status = gsl_multimin_fminimizer_iterate(s);
		   
		if (status) 
			break;
	 
		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-5);
	 
		if (status == GSL_SUCCESS) {
			printf ("converged to minimum at\n");
		}
	 
		printf ("%5d %10.5e f() = %12.7f size = %.5f\n", 
				iter,
				gsl_vector_get (s->x, 0), 
				gsl_vector_get (s->x, 1), 
				gsl_vector_get (s->x, 2), 
				s->fval, size);
	} while (status == GSL_CONTINUE && iter < MaxIter);
	   
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);
	 
	return 0;  //@TODO: How do I make this into a void subroutine?
}
