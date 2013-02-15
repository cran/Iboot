//Last modified 2013/02/13 by Nicola Lunardon.

#include<R.h>
#include<R_ext/Applic.h>

//structure to pass additional arguments to f_min_beta and der_f_min_beta
typedef struct args_func
{
	int n, *sample;
	double *gradient;
}args_func, *ArgsStruct;


//objective function to obtain the Lagrange multiplier (see, Owen 1990)
static double f_min_lambda(int p, double *x, void *ArgsIn)
{
	register int i, j;
	int n;
	ArgsStruct args = (ArgsStruct) ArgsIn;
	double z, logeps, res=0.0, dbln= (double) args->n, nsq;

	logeps = -log(dbln); 
	nsq = dbln*dbln;
	n = args->n;

	for(i=0; i<n; i++)
	{
		z = 1.0;

		for(j=0; j<p; j++) 	z += *(x+j)**(args->gradient + *(args->sample + i) + j*n);

		if(z*dbln<1.0)
			res += logeps - 1.5 + 2.0*dbln*z - 0.5*nsq*(z*z);
		else
			res += log(z);
	}
	return -res;
}	

//first order derivative of "f_min_lambda"
static void der_f_min_lambda(int p, double *x, double *gr, void *ArgsIn)
{
	register int i, j;
	int n;
	ArgsStruct args = (ArgsStruct) ArgsIn;
	double z, dbln=(double)args->n, nsq;

	nsq = dbln*dbln;
	n = args->n;

	for(j=0; j<p; j++) *(gr+j) = 0.0;

	for(i=0; i<n; i++)
	{
		z = 1.0;

		for(j=0; j<p; j++) 	z += *(x+j)**(args->gradient + *(args->sample + i) + j*n);

		if(z*dbln<1.0)
			for(j=0; j<p; j++) *(gr+j) += 2.0**(args->gradient + *(args->sample + i) + j*n)*dbln - z**(args->gradient + *(args->sample + i) + j*n)*nsq;
		else
			for(j=0; j<p; j++) *(gr+j) += *(args->gradient + *(args->sample + i) + j*n)/z;
	}
	for(j=0; j<p; j++) *(gr+j) *= -1.0;
}



//find lagrange multipliers
static void min_lambda(double *lambda, double *gradient, int *sample, int *n, int *p, int *fail, double *pgtol, double *factr, int *lmm, int *trace, int *maxit, int *nREPORT)
{
	register int i;
	int *nbd, fncount=0, grcount=0;
	double *l, *u, Fmin=0.0;
	char msg[60];
	
	args_func args;
	ArgsStruct args_ptr;

	args.n = *n;
	args.gradient = &gradient[0];
	args.sample = &sample[0];
	args_ptr = &args;

	l = (double *) R_alloc(*p, sizeof(double));	
	u = (double *) R_alloc(*p, sizeof(double));	
	nbd = (int *)  R_alloc(*p, sizeof(int));	

	for(i=0; i<*p; i++) *(nbd+i)=0;

	lbfgsb(*p, *lmm, lambda, l, u, nbd, &Fmin, f_min_lambda, der_f_min_lambda, fail, args_ptr, *factr, *pgtol, &fncount, &grcount, *maxit, msg, *trace, *nREPORT);

}
	
//find weights that satisfies the null hypothesis
static void bstr_w0(double *res, double *gradient, int *sample, int *n, int *p, int *fault, double *pgtol, double *sum_w0, double *factr, int *lmm, int *trace, int *maxit, int *nREPORT)
{
	register int i, j;
	double *lambda, tmp, tmp1=0.0, dbln=(double)*n;
	lambda = (double *) R_alloc(*p, sizeof(double));
	
	for(i=0; i<*p; i++) *(lambda+i)=0.0;

	//optimization step
	min_lambda(lambda, gradient, sample, n, p, fault, pgtol, factr, lmm, trace, maxit, nREPORT);

	//compute Owen's weights
	for(i=0; i<*n; i++)	
	{
		tmp = 1.0;
		for(j=0; j<*p; j++) tmp += *(lambda+j)**(gradient + *(sample+i) + j**n);
		*(res+i) = 1.0/( dbln*tmp );
		tmp1 += *(res+i); 
	}
	*sum_w0 = tmp1;
}

static void input_ProbSampleReplace(int *n, int *ids, double *probs)
{
	//in input ids: identities --> in output: will contain the permutation of the ids
	//in input probs: probabilities --> in output: cumulative probabilities

	register int i;

   /* sort the probabilities into descending order */
 	revsort(probs, ids, *n);
	/* compute cumulative probabilities */
	for(i = 1; i < *n; i++)	*(probs+i) += *(probs+i-1);

}

//function to sample with replacement (unequal probabilities). 
//Modified version of ProbSampleReplace that receives in input quantities 
//that need to be computed only once in resamplings (by using "input_ProbSampleReplace)
static void bstr_ProbSampleReplace(int *n, double *cum_probs, int *ids, int *nans, int *ans)
{
	register int i, j;
	double rU;

	/* compute the sample */
	for(i = 0; i < *nans; i++)
	{
		rU = unif_rand();
		for(j = 0; j < *n-1; j++)
		{
			if(rU <= *(cum_probs+j))
			break;
		}
		*(ans+i) = *(ids+j);
	}
}

//computes quadratic form
static double bstr_quad_form(double *gradient, int *sample, int *n, int *p)
{
	register int k, j;
	double tmp=0.0, tmp1 = 0.0;

	for(k=0; k<*p; k++)
	{
		tmp = 0.0;
			for(j=0; j<*n; j++) tmp += *(gradient + *(sample+j) + k**n);
		tmp1 += tmp*tmp;
	}
	return tmp1;
}


//Prepivoting with deterministic stopping rule
void bstr_iboot_stop(double *pwus_oss, double *res_boot, double *res_prep, double *gradient, int *n, int *p, int *B, int *M, int *R, int *kB, double *fails_outer, int *count_good, int *overall_failure, double *pgtol, double *factr, int *lmm, int *trace, int *maxit, int *nREPORT)
{
	register int i, j;
	int *sample, *index, ite2=0, ind=0, *id, *ind_perm, *ind_permR, *ids, fail=1;
	double *w0, *w0_boot, *W0, *res_prep_in, res_prep_tmp=0.0, dbl_M=(double)*M, boot_tmp=0.0, minPrep=0.0, sum_w0=0.0;
	
	ind_perm = ( int* ) R_alloc(*B, sizeof(int));
	for(i=0; i<*B; i++) *(ind_perm+i) = i;
	ind_permR = ( int* ) R_alloc(*R, sizeof(int));
	for(i=0; i<*R; i++) *(ind_perm+i) = i;

	id = (int*) R_alloc(*n, sizeof(int));
	ids = (int*) R_alloc(*n, sizeof(int));

		for(i=0; i<*n; i++) *(ids+i) = i;

	index = ( int* ) R_alloc(*B**n, sizeof(int));
	sample = (int *) R_alloc(*n, sizeof(int));
	W0 = ( double* ) R_alloc(*B**n, sizeof(double));
	w0 = ( double* ) R_alloc(*n, sizeof(double));
	w0_boot = (double *) R_alloc(*n, sizeof(double));

	res_prep_in = (double *) R_alloc(*R, sizeof(double));
	
	//observed value of the statistic
	*pwus_oss = bstr_quad_form(gradient, ids, n, p);

	//compute weights under theta
	bstr_w0(w0, gradient, ids, n, p, &fail, pgtol, &sum_w0, factr, lmm, trace, maxit, nREPORT);

	//compute quantities needed for resampling in the outer level
	input_ProbSampleReplace(n, ids, w0);

	//check for convergence. If not return immediately
	if( fail==0 && sum_w0>1.0-1e-4 && sum_w0<1.0+1e-4 )
	{	
		GetRNGstate();

		//outer level of B bootstrap replications
		while(ind<*B && ite2<*kB+*B)
		{
			R_CheckUserInterrupt();

			ite2 += 1;
			bstr_ProbSampleReplace(n, w0, ids, n, sample);

			//compute the weights
			bstr_w0(w0_boot, gradient, sample, n, p, &fail, pgtol, &sum_w0, factr, lmm, trace, maxit, nREPORT);

			//check convergence
			if( fail==0 && sum_w0>1.0-1e-4 && sum_w0<1.0+1e-4 )
			{
				//store the bootstrap value
				*(res_boot+ind) = bstr_quad_form(gradient, sample, n, p);

				//store the corresponding weights and indices (prepared for further steps of resampling by "input_ProbSampleReplace")
				input_ProbSampleReplace(n, sample, w0_boot);
				for(j=0; j<*n; j++)
				{ 
					*(index + ind**n + j) = *(sample+j);
					*(W0 + ind**n + j) = *(w0_boot+j);
				}
				ind += 1;
				*count_good = ind;
			}
		}
		*fails_outer = (double) (ite2-*B)/ *B;
		
		//check wheter maximum number of iterations have been exceeded. If so, exit
		if(ite2==*kB+*B)
		{
			*overall_failure = 2;
			*R=*B=0;
		}		

		//order bootstrapped values in descending order and store the indices of permutations in ind_perm
		revsort(res_boot, ind_perm, *B);

		//inner loop of bootstrap for the largest min(alpha)*B values (full blown)
		for(i=0; i<*R; i++)
		{
			R_CheckUserInterrupt();		

			*(res_prep_in + i) = 0.0;
			ind = 0;

			w0 = W0 + *( ind_perm + i )**n;
			id = index + *( ind_perm + i )**n;

			while(ind<*M)
			{
				R_CheckUserInterrupt();

				bstr_ProbSampleReplace(n, w0, id, n, sample);

				//bootstrap inner
				boot_tmp = bstr_quad_form(gradient, sample, n, p);

					ind += 1;		
					if( boot_tmp <= *(res_boot + i) ) *(res_prep_in + i) += 1.0;
			}
		}

		//sort the map in descending order 
		revsort(res_prep_in, ind_permR, *R);

		minPrep = *( res_prep_in + *R-1 );

		//start the last loop (stopping rule)
		for(i=*R; i<*B; i++)
		{
			R_CheckUserInterrupt();

			w0 = W0 + *( ind_perm + i )**n;
			id = index + *( ind_perm + i )**n;

			res_prep_tmp = 0.0;
			ind = 0;

			while(ind<*M)
			{
				R_CheckUserInterrupt();

				ind += 1;
				bstr_ProbSampleReplace(n, w0, id, n, sample);

				//bootstrap inner
				boot_tmp = bstr_quad_form(gradient, sample, n, p);

					if( boot_tmp <= *(res_boot + i) ) res_prep_tmp += 1.0;

					//check whether it is possible to stop
					if( res_prep_tmp + dbl_M - (double)(ind+1) <= minPrep ) break;
			}

			if(ind==*M)
			{
				if( res_prep_tmp > *( res_prep_in + *R - 1 ) )
				{
					*( res_prep_in + *R - 1 ) = res_prep_tmp;
					revsort(res_prep_in, ind_permR, *R);
					minPrep = *( res_prep_in + *R-1 );
				}
			}
		}
		PutRNGstate();

		for(i=0; i<*R; i++) *( res_prep + i ) = *( res_prep_in + i );
	}
	else
	{
		*overall_failure = 1;
	}
}

