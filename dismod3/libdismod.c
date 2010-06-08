#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_randist.h>

// dS/da = -(i + m) * S + r * C
// dC/da = i * S - (r + m + f) * C

int func (double t, const double y[], double fun[], void *params)
{
    double i = ((double *)params)[0];
    double m = ((double *)params)[1];
    double r = ((double *)params)[2];
    double f = ((double *)params)[3];
    fun[0] = -(i + m) * y[0] + r * y[1];
    fun[1] = i * y[0] - (r + m + f) * y[1];
    return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    double i = ((double *)params)[0];
    double m = ((double *)params)[1];
    double r = ((double *)params)[2];
    double f = ((double *)params)[3];
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * mat = &dfdy_mat.matrix; 
    gsl_matrix_set (mat, 0, 0, -(i + m));
    gsl_matrix_set (mat, 0, 1, r);
    gsl_matrix_set (mat, 1, 0, i);
    gsl_matrix_set (mat, 1, 1, -(i + m + f));
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}

double trim(double x, double a, double b)
{
    if(x < a)
        x = a;
    if(x > b)
        x = b;
    return x;
}

/**
Name              Type         Size
----------------------------------------------------------------
SC_0              double[]     2
i                 double[]     Length of estimate age mesh
r                 double[]     Length of estimate age mesh
f                 double[]     Length of estimate age mesh
m_all_cause       double[]     Length of estimate age mesh
age_mash          int[]        Length of parameter age mesh
len               int          NA
step              double       NA
NEARLY_ZERO       double       NA
double            SCpm[]       Length of parameter age mesh * 4 
*/
double* scpm(double SC_0[], double i[], double r[], double f[], double m_all_cause[],
             int age_mesh[], int len, double step, double NEARLY_ZERO, double SCpm[])
{
    int ii, len2 = len + len, len3 = len2 + len;
    double params[4];
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk4;
     
    double t, t1;
    double y[2] = {SC_0[0], SC_0[1]}, y_err[2];
    double dydt_in[2], dydt_out[2];

    SCpm[0]    = SC_0[0];
    SCpm[len]  = SC_0[1];
    SCpm[len2] = SC_0[1] / (SC_0[0] + SC_0[1]);
    SCpm[len3] = trim(m_all_cause[age_mesh[0]] - f[age_mesh[0]] * SCpm[len2],
                      .1 * m_all_cause[age_mesh[0]], 1 - NEARLY_ZERO);

    for(ii = 0; ii < len - 1; ii++)
    {
        gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 2);
        params[0] = i[age_mesh[ii]];
        params[1] = SCpm[len3 + ii];

        params[2] = r[age_mesh[ii]];
        params[3] = f[age_mesh[ii]];
        gsl_odeiv_system sys = {func, jac, 2, &params};
        GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);

        t = (double)age_mesh[ii];
        t1 = (double)age_mesh[ii + 1];
        //printf ("time        S            C\n");
        while (t < t1)
        {
            int status = gsl_odeiv_step_apply (s, t, step, y, y_err, 
                                               dydt_in, 
                                               dydt_out, 
                                               &sys);
            if (status != GSL_SUCCESS)
                break;
     
            dydt_in[0] = dydt_out[0];
            dydt_in[1] = dydt_out[1];
     
            t += step;
            //printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
        }
     
        SCpm[ii + 1]        = y[0];
        SCpm[len + ii + 1]  = y[1];
        SCpm[len2 + ii + 1] = trim(SCpm[len + ii + 1] / (SCpm[ii + 1] + SCpm[len + ii + 1]),
                                   NEARLY_ZERO, 1 - NEARLY_ZERO);
        SCpm[len3 + ii + 1] = trim(m_all_cause[age_mesh[ii + 1]] - f[age_mesh[ii + 1]] * SCpm[len2 + ii + 1],
                                   .1 * m_all_cause[age_mesh[ii + 1]], 1 - NEARLY_ZERO);

        gsl_odeiv_step_free (s);
    }

    return SCpm;
}

/**
Name              Type         Size
----------------------------------------------------------------
value             double[]     No. of data points
N                 double[]     No. of data points
Xa                double[]     No. of data points x No. of alpha
Xb                double[]     No. of data points x No. of Beta
alpha             double[]     No. of alpha
beta              double[]     No. of beta
gamma             double[]     Length of estimate age mesh
delta             double       NA
age_indices       int[]        Sum of No. of ages in all studies 
age_weights       double[]     Sum of No. of ages in all studies
ndata             int          NA
nalpha            int          NA
nbeta             int          NA
ngamma            int          NA
ages              int[]        No. of data points
lb                int          NA
logp              double       NA
*/
double obs(double value[], double N[], double Xa[], double Xb[], double alpha[],
           double beta[], double gamma[], double delta, int age_indices[],
           double age_weights[], int ndata, int nalpha, int nbeta, int ngamma,
           int ages[], int lb)
{
    int i = 0, j = 0, k = 0, l;
    double logp = 0;
    double* shifts = (double*) malloc(ndata * sizeof(double));
    double* mu_i = (double*) malloc(ndata * sizeof(double));
    double* exp_gamma = (double*) malloc(ngamma * sizeof(double));

    // exp_gamma = exp(gamma)
    for(i = 0; i < ngamma; i++)
        exp_gamma[i] = exp(gamma[i]);

    // shifts = exp((dot(Xa * alpha) + dot(Xb * beta))
    // mu_i = [dot(w, s * exp_gamma[a]) for s, a, w in zip(shifts, age_indices, age_weights)] * N
    if(lb == 0)
        // logp = mc.negative_binomial_like(value, mu_i*N, delta)
        for(i = 0; i < ndata; i++)
        {
            shifts[i] = 0;
            l = i * nalpha;
            for(j = 0; j < nalpha; j++)
                shifts[i] += Xa[l + j] * alpha[j];
            l = i * nbeta;
            for(j = 0; j < nbeta; j++)
                shifts[i] += Xb[l + j] * beta[j];
            mu_i[i] = 0;
            for(j = 0; j < ages[i]; j++)
                mu_i[i] += age_weights[k + j] * exp_gamma[age_indices[k + j]];
            k += ages[i];
            mu_i[i] *= exp(shifts[i]) * N[i];

            // double gsl_ran_negative_binomial_pdf(unsigned int k, double p, double n)
            // p(k) = {\Gamma(n + k) \over \Gamma(k + 1) \Gamma(n)} p^n (1 - p)^k
            // k = value[i], n = delta, pdelta / (mu_i[i] + delta) = delta / (mu_i[i] * N[i] + delta)
            logp += log(gsl_ran_negative_binomial_pdf(value[i], delta / (mu_i[i] + delta), delta));
        }
    else
        for(i = 0; i < ndata; i++)
        {
            shifts[i] = 0;
            l = i * nalpha;
            for(j = 0; j < nalpha; j++)
                shifts[i] += Xa[l + j] * alpha[j];
            l = i * nbeta;
            for(j = 0; j < nbeta; j++)
                shifts[i] += Xb[l + j] * beta[j];
            mu_i[i] = 0;
            for(j = 0; j < ages[i]; j++)
                mu_i[i] += age_weights[k + j] * exp_gamma[age_indices[k + j]];
            k += ages[i];
            mu_i[i] *= exp(shifts[i]) * N[i];
            if(mu_i[i] < value[i])
                logp += log(gsl_ran_negative_binomial_pdf(value[i], delta / (mu_i[i] + delta), delta));
        }

    free(shifts);
    free(mu_i);
    free(exp_gamma);

    return logp;
}

int main (void)
{
    double SC_0[2] = { 100, 0 };
    double i[11] = { .01, .01, .01, .01, .01, .01, .01, .01, .01, .01, .01 } ;
    double r[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    double f[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    double m_all_cause[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int age_mesh[3] = { 0, 5, 10 };
    int len = 3;
    double step = 1, NEARLY_ZERO = 1e-7;
    double buffer[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    double* SCpm = scpm(SC_0, i, r, f, m_all_cause, age_mesh, len, step, NEARLY_ZERO, buffer);
    printf ("   age = 0     age = 5     age = 10\n");
    printf ("S  %.5e %.5e %.5e\n", SCpm[0],       SCpm[1],           SCpm[2]);
    printf ("C  %.5e %.5e %.5e\n", SCpm[len],     SCpm[len + 1],     SCpm[len + 2]);
    printf ("p  %.5e %.5e %.5e\n", SCpm[len * 2], SCpm[len * 2 + 1], SCpm[len * 2 + 2]);
    printf ("m  %.5e %.5e %.5e\n", SCpm[len * 3], SCpm[len * 3 + 1], SCpm[len * 3 + 2]);

    double value[3] = { 100, 100, 100 };
    double N[3] = { 1000, 1000, 1000 };
    double Xa[6] = { 0, 0, 0, 0, 0, 0 };
    double Xb[6] = { 0, 0, 0, 0, 0, 0 };
    double alpha[2] = { 0, 0 };
    double beta[2] = { 0, 0 };
    double gamma[11] = { -10., -10., -10., -10., -10., -10., -10., -10., -10., -10., -10. };
    double delta = 2.0;
    int age_indices[7] = { 1, 2, 3, 4, 5, 6, 7 };
    double age_weights[7] = { .5, .5, .33, .33, .33, .5, .5 };
    int ndata = 3;
    int nalpha = 2;
    int nbeta = 2;
    int ngamma = 11;
    int ages[3] = { 2, 3, 2 };
    int lb = 0;

    double logp = obs(value, N, Xa, Xb, alpha, beta, gamma, delta, age_indices, age_weights,
                      ndata, nalpha, nbeta, ngamma, ages, lb);
    printf ("logp = %.5e", logp);
    
    return 0;
}
// gcc libdismod.c -o libdismod -I/usr/local/include/gsl/ -lgsl -lgslcblas
// ./libdismod
// gcc -shared -fPIC -DNDEBUG -O2 libdismod.c -o libdismod.so -I/usr/local/include/gsl/ -lgsl -lgslcblas
