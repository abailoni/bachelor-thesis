#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

#include <gsl/gsl_sf.h>
#include "../../functions/numerov.c"

#define h_planck 6.626e-34
#define h_planck_tagliato 1.054571e-34
#define k_B 1.38e-23
#define uma_kg 1.6605e-27

double F_schr(double r, double *params)
{
	double alpha=params[0],E_rid=params[1],l=params[2];
	return(1/alpha*(4*(pow(r,-12)-pow(r,-6))-E_rid)+l*(l+1)*pow(r,-2));
}

int check_fun(double x,double h, double *psi, double *output_vars, double *check_params)
{
	if ( (x+h)<check_params[0] )
		return(0);
	else {
		output_vars[0]=x,output_vars[1]=x+h,output_vars[2]=psi[0],output_vars[3]=psi[1];
		return(1);
	}
}

double psi_start(double r, double alpha)
{
	return(exp(-sqrt(4.0/(25.0*alpha))*pow(r,-5)));
}

int main(int argc, char const *argv[])
{
	//Dati vari:
	double m_H=1.00794*uma_kg, m_Kr=83.798*uma_kg;
	double m=m_Kr*m_H/(m_Kr+m_H);
	double sigma=3.18e-10, eps_rid=68.505;
	double alpha=pow(h_planck_tagliato,2)/(2*m*eps_rid*k_B*pow(sigma,2));
	//Variabili:
	double eVolt=1.602176e-19;
	double E_rid, k_rid;
	// printf("%g mev\n", 1*eps_rid*k_B/eVolt*1e3);
	// printf("%g rid\n", 5.90053*1e-3*eVolt/(eps_rid*k_B));

	//Calcolo sezione d'urto totale in funzione di E:
	double r_max=20,r_start=0.6, h=0.01,precisione =1e-8;
	double start_val_fun[2]={psi_start(r_start,alpha),psi_start(r_start+h,alpha)};
	double F_params[3],check_params[]={r_max},schr[4];
	double r1,r2,psi1,psi2,K;
	double sigma_tot,sigma_tot_old,phase_shift;

	F_params[0]=alpha;
	for (double E=0.001;E<=5;E+=0.01) 
	{
		E_rid=E*1e-3*eVolt;
		E_rid=E_rid/(eps_rid*k_B);
		printf("#Ridotto: %g,  %g meV\n", E_rid,E);
		// return(0);
		// E_rid=0.1;
		F_params[1]=E_rid;
		k_rid=sqrt(E_rid/alpha);
		sigma_tot=0;
		for (int l=0;1;l++){
			F_params[2]=l;
			numerov(r_start,h,start_val_fun,schr,F_params,F_schr,check_params,check_fun);
			r1=schr[0],r2=schr[1],psi1=schr[2],psi2=schr[3];
			//Calcolo phase-shift:
			K=(psi2/r2)*(r1/psi1);
			phase_shift=atan( (gsl_sf_bessel_jl(l,k_rid*r2)-K*gsl_sf_bessel_jl(l,k_rid*r1)) / (gsl_sf_bessel_yl(l,k_rid*r2)-K*gsl_sf_bessel_yl(l,k_rid*r1)) );
			sigma_tot_old=sigma_tot;
			sigma_tot+=4*M_PI/pow(k_rid/sigma,2)*(2.0*l+1)*pow(sin(phase_shift),2);
			printf("#E_rid=%g l=%d phase_shift=%g; delta=%g\n",E_rid,l,phase_shift, sigma_tot*1e20 );

			// printf("sigma_tot=%g; sigma_tot_old=%g\n", sigma_tot,sigma_tot_old);
			if ( fabs((sigma_tot-sigma_tot_old)/sigma_tot_old)<precisione )
				break;
		}
		printf("%g\t%g\n",E_rid,sigma_tot*1e20);
	}


	
	return 0;
}
