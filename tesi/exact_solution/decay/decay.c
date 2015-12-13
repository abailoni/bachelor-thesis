#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_sf.h>
#include "../../../new_libraries/varie.c"

#define alpha 0.035171

double 	h 			= 1e-3,
		tempo_max 	= 10; //tempo massimo di attesa per la ricerca degli zeri

// This is V*(r*)=V(r) divided for epsilon:
double reduced_lenjon_pot(double r, void *input)
{
	return(4*(pow(r,-12)-pow(r,-6)));
}

double potenziale(double r, void *input)
{
	// return(0);
	
	if (r!=0)
	{
		double null[1];
		int l = 4; //choose l=4 for resonances
		double potential = reduced_lenjon_pot(r,null)+alpha*l*(l+1)/pow(r,2);
		
		if (potential>1000)
			return(1000);
		else
			return(potential);
	} else {
		return(0);
	}
}

double normalizzazione(Fun_1D *stato_legato)
{
	int len=stato_legato->len;
	double integrale = 1/2.*(stato_legato->y[0]*stato_legato->y[0]+stato_legato->y[len-1]*stato_legato->y[len-1]); 
	for (int i = 1; i < len-1; ++i)
		integrale += stato_legato->y[i]*stato_legato->y[i];
	integrale *= h;
	for (int i = 0; i < len; ++i)
		stato_legato->y[i]/=sqrt(integrale);
	return(integrale);
}

//////////////////////////////////////
// FUNZIONI NECESSARIE A numerov(): //
//////////////////////////////////////
#include "../../../new_libraries/numerov.c"

struct variables
{
	double A;
	double E_rid;
	int check;
	Fun_1D *stato_legato;
} vars;

// funzione F(x) dell'equazione di schrodinger:
double F_schr(double x, void *input)
{
	double null[1];
	return(1/alpha*(potenziale(x,null) - vars.E_rid));
}

// check_fun() per numerov in avanti:
int check_fun_avanti(struct numerov_check_struct *ck_val)
{
	int num_passi=ck_val->contatore;
	int num_max_passi=(int) (ck_val->ck_par[0]);
	Fun_1D *stato_legato=(Fun_1D*) ck_val->output;

	if ( num_passi>num_max_passi )
		return(1);
	else {
		stato_legato->x[num_passi]=ck_val->x, stato_legato->y[num_passi]=ck_val->psi[0];
		return(0);
	}
}

// check_fun() per numerov indietro:
int check_fun_indietro(struct numerov_check_struct *ck_val)
{
	int num_passi=ck_val->contatore, num_passi_tot= ck_val->ck_par[0]+ck_val->ck_par[1];
	int num_passi_ind=ck_val->ck_par[1];
	Fun_1D *stato_legato=(Fun_1D*) ck_val->output;

	if ( num_passi>num_passi_ind)
		return(1);
	else {
		stato_legato->x[num_passi_tot-num_passi]=ck_val->x, stato_legato->y[num_passi_tot-num_passi]=ck_val->psi[0];
		return(0);
	}
}

//////////////////////////////////////////
// FUNZIONI NECESSARIE A fun_roots():   //
//////////////////////////////////////////

#include "../../../new_libraries/fun_roots_MOD.c" //per poter stampare le funzioni dei "doppietti" è stata modificata la libreria

// Tale funzione data un energia E_rid calcola il "gap" in x_0:
double gap(double E_rid, void *input)
{
	//Variabili per numerov: 
	double x_start_av=0.7, x_start_ind=4.1;
	double x_0=1.3;
	
	// printf("%g-",E_rid );

	// Determino inizio e fine precisi (e numero di passi)
	int num_passi_av= (int)((x_0-x_start_av)/h), num_passi_ind;
	x_0=x_start_av+num_passi_av*h;
	for (int i = 1; 1; ++i)
		if ( x_0+h*i >= x_start_ind ) 
		{
			num_passi_ind=i;
			x_start_ind=x_0+h*i;
			break;
		}

	// Importo variabili:
	Fun_1D *stato_legato = vars.stato_legato;

	// Numerov in avanti:
	double start_val_fun[2]={1e-10,2e-10}, null[1];
	double check_params[]={num_passi_av, num_passi_ind};
	vars.E_rid=E_rid;
	numerov(x_start_av,h,start_val_fun,stato_legato,null,F_schr,check_params,check_fun_avanti);
	double u_avanti_r0=stato_legato->y[num_passi_av];


	// Numerov indietro:
	numerov_inverse(x_start_ind,h,start_val_fun,stato_legato,null,F_schr,check_params,check_fun_indietro);
	double u_indietro_r0=stato_legato->y[num_passi_av];

	
	// Raccordo le due funzioni:
	int passi_tot=num_passi_av+num_passi_ind;
	stato_legato->len=passi_tot+1;
	for (int j=num_passi_av; j<=passi_tot; j++)
		stato_legato->y[j]*=u_avanti_r0/u_indietro_r0;
	
	// Calcolo gap derivate prime in x_0:
	double F_r0=F_schr(x_0,null);
	double gap= (stato_legato->y[num_passi_av-1]+stato_legato->y[num_passi_av+1]-(2+pow(h,2)*F_r0)*u_avanti_r0) /h;
	return( gap );
}

////////////////
// PROGRAMMA: //
////////////////

// Input:
//  - numero di stati da trovare
//  - con secondo input evito la stampa delle funzioni d'onda e mi limito alla stampa delle energie

int main(int argc, char const *argv[])
{
	// Inizializzo variabili:
	int STATI=atoi(argv[1]);
	Fun_1D stati_legati[STATI+1];
	for (int i = 0; i < STATI+1; ++i)
		init_fun_2(1e8,stati_legati+i);

	// // Input:
	// vars.A=atof(argv[1]);
	vars.A =1;

	// Calcolo stati:
	double prec_NR=1e-4, estremi_E[2], zeros[(STATI+1)*2], null[1];
	vars.stato_legato=stati_legati;
	double dE = 1e-3, E_start = -0.4;
	estremi_E[0]=E_start;
	int n=0;
	for (int i = 1; 1; ++i)
	{
		estremi_E[0]= E_start+dE*(i-1), estremi_E[1]= E_start+dE*i;
		double first=gap(estremi_E[0],null), second=gap(estremi_E[1],null);
		// printf("#Energia: %g --> %g e %g\n",estremi_E[0],first,second);

		if (first*second<0) 
		{		
			printf("#Energy = %g\n", (estremi_E[0]+estremi_E[1])/2 );
			normalizzazione(vars.stato_legato);
			vars.stato_legato++;
			n++;
			if (n>STATI-1)
				break;
		}
		// if (first*second<0) 
		// {
		// 	int out=fun_roots(estremi_E,1,prec_NR,zeros,tempo_max,null,gap);
		// 	if (out!=0) //ho trovato qualche stato 
		// 	{
		// 		printf("#Energia=%g (n=%d)\n",zeros[0],n);
		// 		n++;
		// 		// Controllo se ne ho trovati 2 (ho doppietto)
		// 		if (out>0 && n<=STATI-1) {
		// 			printf("#Energia=%g (n=%d)\n",zeros[2],n);
		// 			n++;
		// 		}

		// 		// Controllo se ho finito:
		// 		if (n>STATI-1)
		// 			break;
		// 	}
		// }


	}

	// Controllo se stampare anche funz. d'onda:
	if(argc>2)
		return(0);

	// Salvo potenziale per output:
	double pot[stati_legati->len];
	for (int i = 0; i < stati_legati->len; ++i)
		pot[i]=potenziale(stati_legati->x[i],null);
	

	// Stampo funzioni d'onda:
	for (int i = 0; i < stati_legati->len; ++i) 
	{
		printf("%g\t", stati_legati[0].x[i]);
		for (int j = 0; j < STATI; ++j)
			printf("%g\t", stati_legati[j].y[i]);
		printf("%g\n", pot[i]);
	}

	return 0;
}

