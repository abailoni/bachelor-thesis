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
#include "../../new_libraries/zero_bisezione.c"
#include "../../new_libraries/varie.c"


#define alpha 0.1


// #define r1 1
// #define r2 2
// #define r3 3


// double analitic_potential(double r, void *input)
// {
// 	if (r>=1 && r<=2)
// 		return(1 + alpha*pow(l+1/2.,2)/pow(r,2));
// 	else
// 		return(0 + alpha*pow(l+1/2.,2)/pow(r,2));
// }

// Potenziale efficace e Co.
double pot_eff(double r, void *input)
{
	double *params=(double*)(input);
	double l=params[1];

	if (r>=1 && r<=2)
		return(1 + alpha*pow(l+1/2.,2)/pow(r,2));
	else
		return(0 + alpha*pow(l+1/2.,2)/pow(r,2));

	// return(analitic_potential(r,input));

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
	double zeri[6];

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
	double r_start=1e-2, r_stop=3;
	// Trovo zeri funzione diff_E_pot:
	int count=0;
	double last = diff_E_pot_real(r_start,pot_params_real), new, dr_start=1e-4;
	// fprintf(stderr, "###last=%g e %g\n",last,zeri[3]);
	// getchar();
	for (double r=r_start+dr_start; r<=r_stop; r+=dr_start) {
		new = diff_E_pot_real(r,pot_params_real);
		if (new*last<0) 
		{
			// double estremi_r[2]={r-dr_start,r}, roots_array[2];
			// int out=fun_roots(estremi_r,1,1e-8,roots_array,1,pot_params_real,diff_E_pot_real);
			// if (out==0){
			// 	fprintf(stderr, "## ERRORE: fun_roots non trova lo zero tra %g e %g ##\n", r-dr_start, r);
			// 	return(0);
			// }
			// else if (roots_array[1]<0){ //ho trovato esattamente lo zero 
			// 	roots_array[1]=roots_array[0];
			// }
			// zeri[2*count]=roots_array[0], zeri[2*count+1]=roots_array[1];
			zeri[2*count]=r-dr_start, zeri[2*count+1]=r;
			count++;
			if (count==3)
				break;
		}	
		last=new;
	}

	if (count==0) {
		fprintf(stderr, "## ERRORE: nessuno zero trovato (E=%g) ##\n",creal(E));
		getchar();
		return(0);
	}
	if (count==1)
	{
		fprintf(stderr, "Uscito dal range delle energie compatibili: E=%g\n",creal(E));
		getchar();
		return(0);
	} 
	else if (count==2) {
		fprintf(stderr, "Non trovo terzo zero perché troppo lontano...\n");
		getchar();
	}

	// Scambio ordine nel caso con la buca:
	// (nell'array in prima posizione ho lo zero più grande=più esterno)
	if (count==3) {
		double temp1=zeri[0], temp2=zeri[1];
		zeri[0]=zeri[4], zeri[1]=zeri[5];
		zeri[4]=temp1, zeri[5]=temp2;
	}

	double dr=1e-3;
	// zeri[0]=2-dr, zeri[1]=2+dr, zeri[2]=1+dr, zeri[3]=1-dr, zeri[4]=0, zeri[5]=0;

	// Calcolo J:
	double estremi[2]={zeri[5],zeri[3]};
	complex double pot_params[]={1,l,E,+1};
	complex double J = int_trapezi_complex(estremi,dr,pot_params,rad_diff);
	J/=sqrt(alpha);


	// Calcolo K:
	estremi[0]=zeri[2], estremi[1]=zeri[0];
	pot_params[3]=-1;
	complex double K = int_trapezi_complex(estremi,dr,pot_params,rad_diff);
	K/=sqrt(alpha);
	
	// Calcolo equazione: 
	// return( M_PI*(n+0.5) + I/2. * clog( (1-1/4.*cexp(-2*K))/(1+1/4.*cexp(-2*K)) ) - J );
	// fprintf(stderr, "-");
	return( M_PI*(n+0.5) - I/4. *cexp(-2*K) - J );

}

// PROGRAMMA:
// int n (numero stato metastabile e parte da zero)
// int l

int main(int argc, char const *argv[])
{
	// int n = atoi(argv[1]);
	// int l = atoi(argv[2]);
	// double precisione_energia=1e-7, params[]={n,l};
	// complex double energia_partenza = +0.3 + I*0.01;
	// complex double energia_metastabile = newton_complex(energia_partenza,precisione_energia,params,energie_meta);
	// printf("Energia trovata: E0=%g; parte complessa=%g\n",creal(energia_metastabile),cimag(energia_metastabile) );

	int n = atoi(argv[1]);
	int l = atoi(argv[2]);
	double params[]={n,l};
	// complex double energia_partenza = +0.1 + I*0.01;

	double estremi_real[]={0.1,0.99}, estremi_imag[]={-0.05,0.05};
	double d_real=1e-2, d_imag=1e-2;

	printf("#z_Real;z_Imag;modulo;real;imag\n");
	for (double real = estremi_real[0]; real < estremi_real[1]; real+=d_real)
		for (double imag = estremi_imag[0]; imag < estremi_imag[1]; imag+=d_imag)
		{
			complex double E = real + I*imag;
			complex double equation = energie_meta(E, params);
			printf("%g\t%g\t%g\t%g\t%g\n",real,imag,cabs(equation),creal(equation),cimag(equation));
		}




	// complex double energia_metastabile = newton_complex(energia_partenza,precisione_energia,params,energie_meta);
	// printf("Energia trovata: E0=%g; parte complessa=%g\n",creal(energia_metastabile),cimag(energia_metastabile) );
	return 0;
}


















