#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <complex.h>
#include <gsl/gsl_sf.h>
#include "../../new_libraries/zeri_newton_complex.c"
#include "../../new_libraries/numerov.c"
#include "../../new_libraries/int_trapezi.c"
#include "../../new_libraries/int_trapezi_complex.c"
#include "../../new_libraries/fun_roots.c"
#include "../../new_libraries/varie.c"


#define h_planck 6.626e-34
#define h_planck_tagliato 1.054571e-34
#define k_B 1.38e-23
#define uma_kg 1.6605e-27

#define alpha 0.035171
#define sigma 3.18
#define eps 5.9033

// #define r1 1
// #define r2 2
// #define r3 3


// Potenziale efficace e Co.
double pot_eff(double r, void *input)
{
	double *params=(double*)(input);
	double l=params[1];
	int ck2 = (int)(params[0]); //E' zero nel caso in cui non abbia potenziale 
	return(4*ck2*(1/pow(r,12)-1/pow(r,6)) + alpha*pow(l+1/2.,2)/pow(r,2));
}
complex double diff_E_pot(double r, void *input)
{
	// check = +1 ---> E - potenziale
	// check = -1 ---> potenziale - E
	complex double *params=(complex double*)(input);
	int check = (int) params[3];
	complex double E = params[2];
	if (check==1)
		return(E-pot_eff(r,input));
	else
		return(pot_eff(r,input)-E);
}
double diff_E_pot_real(double r, void *input)
{
	// check = +1 ---> E - potenziale
	// check = -1 ---> potenziale - E
	double *params=(double*)(input);
	int check = (int) params[3];
	double E = params[2];
	if (check==1)
		return(E-pot_eff(r,input));
	else
		return(pot_eff(r,input)-E);
}
complex double rad_diff(double r, void *input)
{
	// check = +1 ---> k_eff (ridotto)
	// check = -1 ---> beta_eff (ridotto)
	return(csqrt(diff_E_pot(r,input)));
}


complex double energie_meta(complex double E, void *input)
{
	// Variabili:
	double *params=(double*)(input);
	int n = (int)(params[0]), l=(int)(params[1]);

	// fprintf(stderr, "%g+i%g\n", creal(E), cimag(E));
	// getchar();

	// ------------------------------------------------
	// Calcolo punti di raccordo:
	// ------------------------------------------------
	// Elementi array pot_params:
	// - int che è 1 se ho V(r) e 0 se considero il caso libero
	// - momento l
	// - energia E (complessa)
	// - int che è +1 o -1 a seconda di cosa voglia calcolare (vedere funzione diff_E_pot)
	double pot_params_real[4]={1,l,creal(E),-1};
	double r_start=0.001, r_stop=20;
	// Trovo zeri funzione diff_E_pot:
	int count=0;
	double zeri[6];
	double last = diff_E_pot_real(r_start,pot_params_real), new, dr_start=1e-2;
	for (double r=r_start+dr_start; r<=r_stop; r+=dr_start) {
		new = diff_E_pot_real(r,pot_params_real);
		if (new*last<0) 
		{
			double estremi_r[2]={r-dr_start,r}, roots_array[2];
			int out=fun_roots(estremi_r,1,1e-8,roots_array,1,pot_params_real,diff_E_pot_real);
			if (out==0){
				fprintf(stderr, "## ERRORE: fun_roots non trova lo zero ##\n");
				return(0);
			}
		else if (roots_array[1]<0){ //ho trovato esattamente lo zero 
				roots_array[1]=roots_array[0];
				// fprintf(stderr, "Attenzione1...\n");
			}
			zeri[2*count]=roots_array[0], zeri[2*count+1]=roots_array[1];
			count++;
		}	
		last=new;
	}
	if (count==0) {
		fprintf(stderr, "## ERRORE: nessuno zero trovato (E=%g) ##\n",creal(E));
		return(0);
	}
	if (count==1)
	{
		fprintf(stderr, "Uscito dal range delle energie compatibili: E=%g\n",creal(E));
		getchar();
		return(0);
	}

	// Scambio ordine nel caso con la buca:
	// (nell'array in prima posizione ho lo zero più grande=più esterno)
	if (count==3) {
		double temp1=zeri[0], temp2=zeri[1];
		zeri[0]=zeri[4], zeri[1]=zeri[5];
		zeri[4]=temp1, zeri[5]=temp2;
	}

	// Calcolo J:
	double estremi[2]={zeri[5],zeri[3]}, dr=1e-3;
	complex double pot_params[]={1,l,E,+1};
	complex double J = int_trapezi_complex(estremi,dr,pot_params,rad_diff);
	J/=sqrt(alpha);
	// Calcolo K:
	estremi[0]=zeri[2], estremi[1]=zeri[0];
	pot_params[3]=-1;
	complex double K = int_trapezi_complex(estremi,dr,pot_params,rad_diff);
	K/=sqrt(alpha);
	
	// Calcolo equazione: 
	return( M_PI*(n+0.5) + I/2. * clog( (1-1/4.*cexp(-2*K))/(1+1/4.*cexp(-2*K)) ) - J );
}

// PROGRAMMA:
// int n (numero stato metastabile e parte da zero)
// int l

int main(int argc, char const *argv[])
{
	int n = atoi(argv[1]);
	int l = atoi(argv[2]);
	double precisione_energia=1e-3, params[]={n,l};
	complex double energia_partenza = +0.1 + I*0.01;
	complex double energia_metastabile = newton_complex(energia_partenza,precisione_energia,params,energie_meta);
	printf("Energia trovata: E0=%g; parte complessa=%g\n",creal(energia_metastabile),cimag(energia_metastabile) );
	return 0;
}


















