#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <complex.h>
#include <gsl/gsl_sf.h>
#include "../../new_libraries/numerov.c"
#include "../../new_libraries/int_trapezi.c"
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

// PROGRAMMA:
// calcola la sezione d'urto totale al variare di E

int main(int argc, char const *argv[])
{
	// Calcolo sezione d'urto:
	double precisione_sez=1e-8, shift[2];
	for(double E=0.01; E<1; E+=0.001)
	{
		double k_2=E/alpha;
		double sez_tot=0;
		for (int l = 0; 1; ++l)
		{
			double sez_tot_old=sez_tot;
			phase_shift(l,E,shift);
			double shift_tot=+shift[0]+shift[1];
			sez_tot+=(2*l+1.)*pow(sin(shift_tot),2);
			if (fabs((sez_tot-sez_tot_old)/sez_tot_old)<precisione_sez || l>40){
				fprintf(stderr, ".%d",l);
				break;
			}
		}
		printf("%g\t%g\n",E,sez_tot*4*M_PI/k_2*sigma*sigma);
	}

	fprintf(stderr, "\n");
	return 0;
}
