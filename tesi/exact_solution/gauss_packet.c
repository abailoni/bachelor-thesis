#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_fft_complex.h>

#define zeta 1

#include "../../new_libraries/varie.c"
#include "../../new_libraries/FFT_library.c"
#include "../../new_libraries/split_operator.c"
#include "../../new_libraries/int_trapezi.c"

#define h_tagl 1
#define m 1
#define a 16

#define N 2048

#define dt_anim 25
#define dt 0.5
#define T_MAX 25000

double dx=1;

double potenziale(double x, void *input)
{
	double pot_start=N*dx/2, *params=(double*)input;
	if (x>pot_start && x<pot_start+a)
		return(1/(10.*potenza(params[0],2)));
	else
		return(0);
}

int check_split_op(check_struct_split_op *ck_val)
{
	FFT_vect *pacchetto=ck_val->funz_onda;
	int count = (int) ck_val->check_params[0];
	// Check if print:
	if (  ck_val->t>=count*dt_anim)
	{
		printf("%g",MODULO(pacchetto->main,0));
		for(int n=1;n<N;n++)
			printf(";%g",MODULO(pacchetto->main,n));
		printf("\n");
		ck_val->check_params[0]+=1;
	}

	// Check when stop:
	return( ck_val->t<T_MAX ? 0 : 1 );
}

void pacc_gauss(FFT_var *vars)
{
	double correzione=N*dx/2;
	double sigma=vars->params[0], k=vars->params[1];
	complex double pack_val;
	// pack_val = 20*cexp(I*k*vars->x)+cexp(-I*k*vars->x);
	pack_val = exp(- potenza((vars->x - correzione+2*sigma)/sigma, 2)) * cexp(I*k*vars->x);
	vars->f_r=creal(pack_val), vars->f_i=cimag(pack_val);	
}

int main(int argc, char const *argv[])
{
	// Variabili:
	double sigma=N*dx/16, k=M_PI/(dx*16);
	k *=3;

	// Salvo potenziale:
	double pot_params[]={dx};
	Fun_1D potential = create_fun(0,dx*(N-1),dx,pot_params,potenziale);

	// Stampo potenziale?
	if (argc>1)
	{
		print_fun(potential);
		return(0);
	}

	// Creo pacchetto:
	FFT_vect pacchetto, pacchetto_finale;
	double params[]={sigma,k};
	create_FFT(N,dx,&pacchetto,params,pacc_gauss);
	copy_FFT(pacchetto, &pacchetto_finale);


	// Evolvo nel tempo:
	double check_params[1]={0}, null[1];
	split_operator(N,dt,dx*N,&pacchetto_finale,null,&potential,return_val_ordered,check_params,check_split_op);

	// // Ricopio tutto in vettori Fun1D:
	// Fun_1D pacc_start=init_fun(1e3), pacc_stop=init_fun(1e3);
	// for (int i = 0; i < N; ++i)
	// {
	// 	pacc_start.x[i]=pacchetto.x[i], pacc_start.y[i]=potenza(MODULO(pacchetto.main,i),2), pacc_start.len=N;
	// 	pacc_stop.x[i]=pacchetto.x[i], pacc_stop.y[i]=potenza(MODULO(pacchetto_finale.main,i),2), pacc_stop.len=N;
	// }

	// // Calcolo coefficienti trasmissione e riflessione:
	// double estremi[]={0,N*dx/2.+a/2.};
	// double rifless=int_trapezi(estremi,dx,&pacc_stop,return_val_ordered,2);
	// estremi[0]=N*dx/2.+a/2.+dx, estremi[1]=(N-1)*dx;
	// double trasm=int_trapezi(estremi,dx,&pacc_stop,return_val_ordered,2);
	// estremi[0]=0, estremi[1]=(N-1)*dx;
	// double totale=int_trapezi(estremi,dx,&pacc_start,return_val_ordered,2);

	// // Stampo pacchetto:
	// printf("#Trasmesso: %g \tRiflesso: %g \n",trasm/totale,rifless/totale);
	// for(int n=0;n<N;n++)
 //  		printf("%g\t%g\t%g\t%g\n",pacchetto.x[n],MODULO(pacchetto.main,n),MODULO(pacchetto_finale.main,n),potential.y[n]);
	return 0;
}
