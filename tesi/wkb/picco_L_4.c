#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <complex.h>
#include <gsl/gsl_sf.h>
#include "../../new_libraries/numerov.c"
#include "../../new_libraries/int_trapezi.c"
#include "../../new_libraries/fun_roots.c"
#include "../../new_libraries/varie.c"
#include "../../new_libraries/fun_max.c"


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

void phase_shift(int l, double E, double *shift)
{
	// ------------------------------------------------
	// Trovo se potenziale efficace ha una buca o meno:
	// ------------------------------------------------
	// Elementi array pot_params:
	// - int che è 1 se ho V(r) e 0 se considero il caso libero
	// - momento l
	// - energia E
	// - int che è +1 o -1 a seconda di cosa voglia calcolare (vedere funzione diff_E_pot)
	double pot_params[4]={1,l,E,-1};
	double r_start=0.001, r_stop=100;
	// Trovo zeri funzione diff_E_pot:
	int count=0;
	double zeri[6];
	double last = diff_E_pot(r_start,pot_params), new, dr_start=1e-2;
	for (double r=r_start+dr_start; r<=r_stop; r+=dr_start) {
		new = diff_E_pot(r,pot_params);
		if (new*last<0) 
		{
			double estremi_r[2]={r-dr_start,r}, roots_array[2];
			int out=fun_roots(estremi_r,1,1e-8,roots_array,1,pot_params,diff_E_pot);
			if (out==0){
				fprintf(stderr, "## ERRORE: fun_roots non trova lo zero ##\n");
				fprintf(stderr, "## Dati errore: E=%g l=%d\n", E,l);
				fprintf(stderr, "## Estremi r: %g e %g\n",estremi_r[0],estremi_r[1] );
				fprintf(stderr, "Potenziale eff: %g\n", pot_eff(estremi_r[0],pot_params));
				fprintf(stderr, "## Premere per continuare.\n");
				getchar();
				return;
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
		return;
	}
	// Scambio ordine nel caso con la buca:
	// (nell'array in prima posizione ho lo zero più grande=più esterno)
	if (count==3) {
		double temp1=zeri[0], temp2=zeri[1];
		zeri[0]=zeri[4], zeri[1]=zeri[5];
		zeri[4]=temp1, zeri[5]=temp2;
	}

	// ------------------------------------------------
	// Calcolo parte regolare del phase-shift: (in entrambi i casi)
	// ------------------------------------------------
	r_stop=zeri[1]+30;

	double prec_integrali=1e-5;
	double estremi_int[2]={zeri[1],r_stop};
	pot_params[3]=+1; //mi interessa k_eff
	// printf("r1=%g r2=%g r3=%g // l=%d; E=%g\n",zeri[1],zeri[2],zeri[4],l,E);
	double val_int = int_trapezi(estremi_int,prec_integrali,pot_params,rad_diff,1);
	// Sfasamento senza potenziale:
	pot_params[0]=0; //spengo potenziale
	pot_params[3]=-1;
	double estremi_r[2]={r_start,r_stop}, roots_array[2];
	if (diff_E_pot(r_start,pot_params)*diff_E_pot(r_stop,pot_params)<0)
		fun_roots(estremi_r,1,1e-8,roots_array,0.1,pot_params,diff_E_pot);
	else
		roots_array[1]=r_start;
	if(roots_array[1]<0){
		// fprintf(stderr, "Attenzione2...\n");
		roots_array[1]=roots_array[0];
	}
	pot_params[3]=+1;
	estremi_int[0]=roots_array[1], estremi_int[1]=r_stop;
	double val_int2 = int_trapezi(estremi_int,prec_integrali,pot_params,rad_diff,1);
	shift[0]=(val_int-val_int2)/sqrt(alpha);

	// ------------------------------------------------
	// Calcolo parte risonante del phase-shift: (solo con buca)
	// ------------------------------------------------
	if (count==1){
		shift[1]=0;
		return;
	}
	pot_params[0]=1; //riaccendo potenziale
	// Trovo J:
	estremi_int[0]=zeri[5], estremi_int[1]=zeri[3], pot_params[3]=+1;
	val_int = int_trapezi(estremi_int,prec_integrali,pot_params,rad_diff,1);
	double J = val_int/sqrt(alpha);
	// Trovo K:
	estremi_int[0]=zeri[2], estremi_int[1]=zeri[0], pot_params[3]=-1;
	val_int = int_trapezi(estremi_int,prec_integrali,pot_params,rad_diff,1);
	double K = val_int/sqrt(alpha);

	complex double shift_complesso;
	shift_complesso = 1/tan(J) + I * exp(-2*K)/4;
	shift[1]=carg(shift_complesso);

	return;
}

double sezione_urto(double E, void *input)
{
	double *params=(void*) input;
	int l=(int) params[0];
	double k_2=E/alpha;
	double shift[2];
	phase_shift(l,E,shift);
	double shift_tot=+shift[0]+shift[1];
	double sez_tot = (2*l+1.)*pow(sin(shift_tot),2);
	return(sez_tot*4*M_PI/k_2*sigma*sigma);
}

double diff(double E, void *input)
{
	double *params=(void*) input;
	double meta_altezza=params[1];
	return(sezione_urto(E,input)-meta_altezza);
}


double buca(double x, void *input)
{
	return(exp(-pow(x-1.8,2))+exp(-pow(x+1.8,2)));
}

// PROGRAMMA:

int main(int argc, char const *argv[])
{
	int l= 4;
	
	// // Stampa potenziale:
	// double pot_params[]={1,l};
	// for(double r=0.01; r<10; r+=0.001)
	// 	printf("%g\t%g\n", r,pot_eff(r,pot_params));

	// // Buca:
	// double buca_params[1];
	// for(double r=-7; r<7; r+=0.001)
	// 	printf("%g\t%g\n", r,buca(r,buca_params));

	// // Stampo sezione d'urto:
	// double sez_params[2]={l,0};
	// for(double E=0.02; E<1; E+=0.001)
	// 	printf("%g\t%g\n", E,sezione_urto(E,sez_params)/pow(sigma,2));
	// return(0);

	// Calcolo valore risonanza: (trovo massimo)
	double params[2]={l,0};
	double estremi_max[2]={0.075,0.1}, precisione_max=1e-8;
	fun_max(estremi_max,precisione_max,params,sezione_urto);
	double E_max = estremi_max[0], sez_max = sezione_urto(estremi_max[0],params);
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
