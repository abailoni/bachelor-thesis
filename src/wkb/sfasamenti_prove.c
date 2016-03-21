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


// Potenziale efficace e Co.
double pot_eff(double r, void *input)
{
	double *params=(double*)(input);
	double l=params[1];
	int ck2 = (int)(params[0]); //E' zero nel caso in cui non abbia potenziale 
	return(4*ck2*(1/pow(r,12)-1/pow(r,6)) + alpha*pow(l+1/2.,2)/pow(r,2));
}
double diff_E_pot(double r, void *input)
{
	// check = +1 ---> E - potenziale
	// check = -1 ---> potenziale - E
	double *params=(double*)(input);
	int check = (int) params[3];
	if (check==1)
		return(params[2]-pot_eff(r,input));
	else
		return(pot_eff(r,input)-params[2]);
}
double rad_diff(double r, void *input)
{
	// check = +1 ---> k_eff (ridotto)
	// check = -1 ---> beta_eff (ridotto)
	return(sqrt(diff_E_pot(r,input)));
}


complex double funzione_onda3(double r, double r_0, double J, double K, double *pot_params)
{
	// Primo integrale:
	double estremi[2]={r_0,r}, precisione_int=1e-3;
	pot_params[3]= +1;
	double integr_1=int_trapezi(estremi,precisione_int,pot_params,rad_diff,1); 
	integr_1/=sqrt(alpha);

	complex double term1 = cexp(-I*M_PI/4.) * (cos(J)*exp(K)+I/4.*sin(J)*exp(-K)) * cexp(I*integr_1); 
	complex double term2 = cexp(I*M_PI/4.) * (cos(J)*exp(K)-I/4.*sin(J)*exp(-K)) * cexp(-I*integr_1);
	return( (term1+term2)/sqrt(rad_diff(r,pot_params)) ); 

}

double funzione_onda1(double r, double r_0, double *pot_params)
{
	// Primo integrale:
	double estremi[2]={r_0,r}, precisione_int=1e-4;
	pot_params[3]= +1;
	double integr_1=int_trapezi(estremi,precisione_int,pot_params,rad_diff,1); 
	integr_1/=sqrt(alpha);

	return(sin(integr_1+M_PI/4)/sqrt(rad_diff(r,pot_params)));
}

double phase_shift(double E, int l)
{
	// VARIABILI:
	double r_start=1e-5, r_stop=100;

	// ------------------------------------------------
	// Trovo se potenziale efficace ha una buca o meno:
	// ------------------------------------------------
	// Elementi array pot_params:
	// - int che è 1 se ho V(r) e 0 se considero il caso libero
	// - momento l
	// - energia E
	// - int che è +1 o -1 a seconda di cosa voglia calcolare (vedere funzione diff_E_pot)
	double pot_params[4]={1,l,E,-1};
	// Trovo zeri funzione diff_E_pot:
	int count=0;
	double zeri[6];
	double last = diff_E_pot(r_start,pot_params), new, dr_start=1e-2;
	for (double r=r_start+dr_start; r<=r_stop; r+=dr_start) {
		new = diff_E_pot(r,pot_params);
		if (new*last<0) 
		{
			double estremi_r[2]={r-dr_start-0.001,r+0.001}, roots_array[2];
			int out=fun_roots(estremi_r,1,1e-8,roots_array,5,pot_params,diff_E_pot);
			if (out==0){
				fprintf(stderr, "## ERRORE: fun_roots non trova lo zero ##\n");
				fprintf(stderr, "## Dati errore: E=%g l=%d\n", E,l);
				fprintf(stderr, "## Estremi r: %g e %g\n",estremi_r[0],estremi_r[1] );
				fprintf(stderr, "Potenziale eff: %g\n", pot_eff(estremi_r[0],pot_params));
				fprintf(stderr, "## Premere per continuare.\n");
				getchar();
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
		fprintf(stderr, "## ERRORE: nessuno zero trovato ##\n");
		fprintf(stderr, "## Dati di errore: E=%g l=%d\n", E,l);
		fprintf(stderr, "## l troppo alto e r_max non sufficiente!\n");
		fprintf(stderr, "## Premere per continuare.\n");
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

	////////////////////////////////////
	// CASO 1: una sola intersezione //
	////////////////////////////////////
	r_stop=zeri[1]+100;
	double r1=r_stop-1, r2=r1+1e-2;
	complex double psi1=0, psi2=0;
	if (count==1)
	{
		printf("#(zero=%g)#",zeri[1]);
		psi1=funzione_onda1(r1,zeri[1],pot_params), psi2=funzione_onda1(r2,zeri[1],pot_params);
		// for (double r = zeri[1]; r < r_stop; r+=0.001)
		// 	printf("%g\t%g\n",r,funzione_onda1(r,zeri[1],pot_params));
		// return(0);
	}
	///////////////////////////////
	// CASO 2: tre intersezioni //
	///////////////////////////////
	else if (count==3)
	{
		// Trovo J:
		double estremi[2], precisione_int=1e-3;
		estremi[0]=zeri[5], estremi[1]=zeri[3], pot_params[3]=+1;
		double J = int_trapezi(estremi,precisione_int,pot_params,rad_diff,1);
		J /= sqrt(alpha);
		// Trovo K:
		estremi[0]=zeri[2], estremi[1]=zeri[0], pot_params[3]=-1;
		double K = int_trapezi(estremi,precisione_int,pot_params,rad_diff,1);
		K /= sqrt(alpha);
		pot_params[3]=+1;
		psi1=funzione_onda3(r1,zeri[1],J,K,pot_params), psi2=funzione_onda3(r2,zeri[1],J,K,pot_params);
		// printf("l=%d, J=%g K=%g ps1=%g\n",l,J,K,creal(psi1));
		// getchar();
		// for (double r = zeri[1]; r < r_stop; r+=0.001)
		// 	printf("%g\t%g\n",r,creal(funzione_onda3(r,zeri[1],J,K,pot_params)));
		// return(0);
	}
	///////////////////////////
	// TROVO SFASAMENTO:     //
	///////////////////////////
	complex double KONST=(psi2/r2)*(r1/psi1); //COSTANTE COMPLESSA?!?!?
	// printf("real=%g imag=%g\n",creal(KONST),cimag(KONST) );
	// return(0);
	double k_rid=sqrt(E/alpha);
	KONST=creal(KONST); //la parte immaginaria è nulla...
	return(atan( (gsl_sf_bessel_jl(l,k_rid*r2)-KONST*gsl_sf_bessel_jl(l,k_rid*r1)) / (gsl_sf_bessel_yl(l,k_rid*r2)-KONST*gsl_sf_bessel_yl(l,k_rid*r1)) ));
}


// PROGRAMMA:

int main(int argc, char const *argv[])
{
	// Calcolo phase-shift:
	// int l=atoi(argv[1]);
	double E = atof(argv[1]);

	// Stampo:
	for (int l = 0; l < 75; ++l)
	{
		double shift =phase_shift(E, l);
		printf("#l=%d, Shift=%g\n",l, shift);
	}
	return 0;
}










