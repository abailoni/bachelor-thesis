#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

#include <gsl/gsl_sf.h>
#include "../../functions/numerov.c"
#include "../../new_libraries/fun_max.c"
#include "../../new_libraries/fun_roots.c"


#define h_planck 6.626e-34
#define h_planck_tagliato 1.054571e-34
#define k_B 1.38e-23
#define uma_kg 1.6605e-27

#define alpha 0.035171
#define sigma 3.18
#define eps_rid 5.9033

double F_schr(double r, double *params)
{
	double E_rid=params[1],l=params[2];
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

double psi_start(double r)
{
	return(exp(-sqrt(4.0/(25.0*alpha))*pow(r,-5)));
}

double sez_urto(double E_rid, void *input)
{
	//Variabili:
	// double eVolt=1.602176e-19;
	double k_rid;
	double *input_par=(double*)input;

	//Calcolo sezione d'urto totale in funzione di E:
	double r_max=20,r_start=0.6, h=0.01,precisione =1e-8;
	double start_val_fun[2]={psi_start(r_start),psi_start(r_start+h)};
	double F_params[3],check_params[]={r_max},schr[4];
	double r1,r2,psi1,psi2,K;
	double sigma_tot,phase_shift;

	// E_rid=E*1e-3*eVolt;
	// E_rid=E_rid/(eps_rid*k_B);
	// printf("#Ridotto: %g,  %g meV\n", E_rid,E);
	// return(0);
	// E_rid=0.1;
	F_params[1]=E_rid;
	k_rid=sqrt(E_rid/alpha);
	
	int l=(int) input_par[0];

	F_params[2]=l;
	numerov(r_start,h,start_val_fun,schr,F_params,F_schr,check_params,check_fun);
	r1=schr[0],r2=schr[1],psi1=schr[2],psi2=schr[3];
	//Calcolo phase-shift:
	K=(psi2/r2)*(r1/psi1);
	phase_shift=atan( (gsl_sf_bessel_jl(l,k_rid*r2)-K*gsl_sf_bessel_jl(l,k_rid*r1)) / (gsl_sf_bessel_yl(l,k_rid*r2)-K*gsl_sf_bessel_yl(l,k_rid*r1)) );
	sigma_tot=4*M_PI/pow(k_rid/sigma,2)*(2.0*l+1)*pow(sin(phase_shift),2);
	return(sigma_tot);
}

double diff(double E_rid, void *input)
{
	double *params=(void*) input;
	double meta_altezza=params[1];
	return(sez_urto(E_rid,input)-meta_altezza);
}

int main(int argc, char const *argv[])
{	
	int l=4;

	double m_H=1.00794*uma_kg, m_Kr=83.798*uma_kg;
	double m=m_Kr*m_H/(m_Kr+m_H);
	fprintf(stderr, "Massa: %g\n", m);

	// // Stampo sezione d'urto:
	// double sez_params[2]={l,0};
	// for(double E=0.001; E<1; E+=0.001)
	// 	printf("%g\t%g\n", E,sez_urto(E,sez_params)/pow(sigma,2));
	// return(0);

	// Calcolo valore risonanza: (trovo massimo)
	double params[2]={l,0};
	double estremi_max[2]={0.075,0.1}, precisione_max=1e-8;
	fun_max(estremi_max,precisione_max,params,sez_urto);
	double E_max = estremi_max[0], sez_max = sez_urto(estremi_max[0],params);
	printf("#---> Massimo ad energia: %g || Sezione d'urto: %g <---#\n", E_max, sez_max);
	
	// Calcolo larghezza risonanza:
	double estremi_width[2]={0.06,0.11}, precisione_width=1e-8, zeros[4];
	params[1]=sez_max/2.;
	int out=fun_roots(estremi_width,2,precisione_width,zeros,1,params,diff);
	if (out<=0)
	{
		fprintf(stderr, "Error fun_roots!!\n");
		return(0);
	}
	double larghezza=fabs(zeros[0]-zeros[2]);
	printf("## Larghezza: %g ##\n",larghezza);
	return 0;
	
}
