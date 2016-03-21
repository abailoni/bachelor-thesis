#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_fft_complex.h>

#define zeta 0.035171
#define alpha 0.035171

#include "../../new_libraries/varie.c"
#include "../../new_libraries/FFT_library.c"
#include "../../new_libraries/split_operator_MOD_without_0.5.c"
#include "../../new_libraries/int_trapezi.c"
#include "../../new_libraries/derivata_FFT.c"

#define N 2048

#define dt_anim 0.02
#define dt 0.01
#define T_MAX 10

double dr = 12./N;

// This is V*(r*)=V(r) divided for epsilon:
double reduced_lenjon_pot(double r, void *input)
{
	return(4*(pow(r,-12)-pow(r,-6)));
}

double potential_fun(double r, void *input)
{
	// return(0);
	
	if (r!=0)
	{
		double null[1];
		int l = 4; //choose l=4 for resonances
		return(reduced_lenjon_pot(r,null)+alpha*l*(l+1)/pow(r,2));
	} else {
		return(0);
	}

}

//--------------------------------------------
// Routines necessary to mean_energy:

void coniugato_wave_fun(FFT_var *vars)
{
	vars->f_r=vars->old_r, vars->f_i=-vars->old_i;
}

double normalizzazione_FFT_vect(FFT_vect *funzione)
{
	// int len=stato_legato->len;
	double integrale = 1/2. * ( potenza(MODULO(funzione->main,0), 2) + 
		potenza(MODULO(funzione->main,N-1), 2)); 
	for (int i = 1; i < N-1; ++i)
		integrale += potenza(MODULO(funzione->main,i), 2);
	integrale *= dr;
	for (int i = 0; i < N; ++i)
		REAL(funzione->main,i)/=sqrt(integrale), IMAG(funzione->main,i)/=sqrt(integrale);
	return(integrale);
}

// The packet must be normalized!
double mean_energy_fun(FFT_vect *wave_function, Fun_1D potential)
{
	// Calculate the derivata:
	double L = dr*N;
	FFT_vect derivata1, derivata2;
	derivata_FFT(N,L,*wave_function,&derivata1);
	derivata_FFT(N,L,derivata1,&derivata2);


	// Coniugato:
	FFT_vect coniug_wave_function;
	copy_FFT(*wave_function,&coniug_wave_function);
	double null[1];
	update_FFT(&coniug_wave_function,null,coniugato_wave_fun);

	// print_FFT(3, derivata1);
	// return(0);


	// Prodotto primo integrale: (modifico derivata2)
	product(&derivata2,coniug_wave_function);



	// Prepare integrands:
	Fun_1D integrand_1 = init_fun(1e6), integrand_2 = init_fun(1e6);
	integrand_1.len=N, integrand_2.len=N;
	for (int i=0; i < N; i++)
	{
		double r = i*dr;
		integrand_1.x[i]=r, integrand_2.x[i]=r;
		integrand_1.y[i]=-alpha*REAL(derivata2.main,i);
		integrand_2.y[i]=pow(MODULO(wave_function->main,i),2)*return_val_ordered(r,&potential);
	}

	

	// Do integrals:
	double estremi[2]={0,L};
	double integral1 = int_trapezi(estremi,dr,&integrand_1,return_val_ordered,2);
	double integral2 = int_trapezi(estremi,dr,&integrand_2,return_val_ordered,2);

	// // Print integrals:
	// fprintf(stderr, "Integrali: %g e %g\n", integral1, integral2);

	return(integral1+integral2);
}
//--------------------------------------------

int check_split_op(check_struct_split_op *ck_val)
{
	FFT_vect *pacchetto=ck_val->funz_onda;
	int count = (int) ck_val->check_params[0];
	Fun_1D *potential = (Fun_1D*) ck_val->output;

	//Calculate mean energy:
	double mean_energy = mean_energy_fun(pacchetto,*potential);

	// Check if print:
	if (  ck_val->t>=count*dt_anim)
	{
		// Al primo posto salvo l'energia:
		printf("%g",mean_energy);
		for(int n=0;n<N;n++)
			printf(";%g",MODULO(pacchetto->main,n));
		printf("\n");
		ck_val->check_params[0]+=1;
	}

	// Check when stop:
	return( ck_val->t<T_MAX ? 0 : 1 );
}

void pacc_gauss(FFT_var *vars)
{
	double correzione=N*dr/2;
	double sigma=vars->params[0], k=vars->params[1];
	complex double pack_val;
	// pack_val = 20*cexp(I*k*vars->x)+cexp(-I*k*vars->x);
	pack_val = exp(- potenza((vars->x - correzione+2*sigma)/sigma, 2))* cexp(-I*k*vars->x);
	vars->f_r=creal(pack_val), vars->f_i=cimag(pack_val);	
}

int main(int argc, char const *argv[])
{
	// Variabili:
	double sigma=N*dr/16, k=M_PI/(dr*16);
	k /=3;

	// Salvo potenziale:
	double pot_params[]={dr};
	Fun_1D potential = create_fun(0,dr*(N-1),dr,pot_params,potential_fun);

	// Stampo potenziale?
	if (argc>1)
	{
		print_fun(potential);
		return(0);
	}

	// Creo pacchetto:
	FFT_vect pacchetto, pacchetto_finale;
	double params[]={sigma,k};
	create_FFT(N,dr,&pacchetto,params,pacc_gauss);
	normalizzazione_FFT_vect(&pacchetto);
	copy_FFT(pacchetto, &pacchetto_finale);

	// // Print packet:
	// print_FFT(1,pacchetto);
	// return(0);

	// // Calculate mean energy:
	// double average_energy=mean_energy_fun(&pacchetto,potential);
	// printf("#Energia: %g\n",average_energy);
	// printf("k_iniziale = %g; k_energia_cinetica=%g\n", k, sqrt(average_energy/alpha));
	// return(0);

	// Evolvo nel tempo:
	double check_params[1]={0};
	//Al posto di output inserisco il potenziale che mi serve per la media:
	split_operator(N,dt,dr*N,&pacchetto_finale,&potential,&potential,return_val_ordered,check_params,check_split_op);

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
