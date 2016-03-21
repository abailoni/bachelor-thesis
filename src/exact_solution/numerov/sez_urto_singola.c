#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

#include <gsl/gsl_sf.h>
#include "../../../new_libraries/numerov.c"

#define alpha 0.035171
#define h_planck 6.626e-34
#define h_planck_tagliato 1.054571e-34
#define k_B 1.38e-23
#define uma_kg 1.6605e-27

// This is V*(r*)=V(r) divided for epsilon:
double reduced_pot(double r, void *input)
{
	return(4*(pow(r,-12)-pow(r,-6)));
}

// This is the effective potential (reduced) and divided for epsilon:
double eff_pot(double r, void *input)
{
	double null[1], *params = (double*)input;
	int l=(int) params[2];
	return(reduced_pot(r,null)+alpha*l*(l+1)/pow(r,2));
}

double F_schr(double r, void *input)
{
	double *params = (double*)input;
	double E_rid=params[1];
	return(1/alpha*(eff_pot(r,params)-E_rid));
}

int check_fun(struct numerov_check_struct *ck_val)
{
	double *check_params=ck_val->ck_par;
	double *output_vars=ck_val->output;
	double x = ck_val->x, h = *ck_val->h;
	if ( (x+h)<check_params[0] )
		return(0);
	else {
		output_vars[0]=x,output_vars[1]=x+h,output_vars[2]=ck_val->psi[0],output_vars[3]=ck_val->psi[1];
		return(1);
	}
}

double psi_start(double r)
{
	return(exp(-sqrt(4.0/(25.0*alpha))*pow(r,-5)));
}

int main(int argc, char const *argv[])
{
	//Dati vari:
	// double m_H=1.00794*uma_kg, m_Kr=83.798*uma_kg;
	// double m=m_Kr*m_H/(m_Kr+m_H);
	double sigma=3.18e-10, eps_rid=68.505;
	// double alpha2=pow(h_planck_tagliato,2)/(2*m*eps_rid*k_B*pow(sigma,2));
	
	//Variabili:
	double eVolt=1.602176e-19;
	double E_rid, k_rid;
	// printf("%g mev\n", 1*eps_rid*k_B/eVolt*1e3);
	// printf("%g rid\n", 5.90053*1e-3*eVolt/(eps_rid*k_B));

	//Calcolo sezione d'urto totale in funzione di E:
	double r_max=20,r_start=0.6, h=0.01;
	double start_val_fun[2]={psi_start(r_start),psi_start(r_start+h)};
	double F_params[3],check_params[]={r_max},schr[4];
	double r1,r2,psi1,psi2,K;
	double sigma_tot,phase_shift;

	int l = 4;

	F_params[0]=0;
	for (double E=0.001;E<=5;E+=0.01) 
	{
		E_rid=E*1e-3*eVolt;
		E_rid=E_rid/(eps_rid*k_B);
		// printf("#Ridotto: %g,  %g meV\n", E_rid,E);
		// return(0);
		// E_rid=0.1;
		F_params[1]=E_rid;
		k_rid=sqrt(E_rid/alpha);
		F_params[2]=l;
		numerov(r_start,h,start_val_fun,schr,F_params,F_schr,check_params,check_fun);
		r1=schr[0],r2=schr[1],psi1=schr[2],psi2=schr[3];
		//Calcolo phase-shift:
		K=(psi2/r2)*(r1/psi1);
		phase_shift=atan( (gsl_sf_bessel_jl(l,k_rid*r2)-K*gsl_sf_bessel_jl(l,k_rid*r1)) / (gsl_sf_bessel_yl(l,k_rid*r2)-K*gsl_sf_bessel_yl(l,k_rid*r1)) );
		sigma_tot=4*M_PI/pow(k_rid/sigma,2)*(2.0*l+1)*pow(sin(phase_shift),2);
		// printf("#E_rid=%g l=%d phase_shift=%g; delta=%g\n",E_rid,l,phase_shift, sigma_tot*1e20 );

		// printf("sigma_tot=%g; sigma_tot_old=%g\n", sigma_tot,sigma_tot_old);
		printf("%g\t%g\n",E_rid,sigma_tot*1e20);
	}


	
	return 0;
}
