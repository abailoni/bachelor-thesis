// Try to derivate the complex function in two different parts: before the real part and only then the imaginary part.

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
#include "../../new_libraries/split_operator_complex_potential.c"
#include "../../new_libraries/int_trapezi.c"
#include "../../new_libraries/derivata_FFT.c"

#define N 8192
// #define N 4096

#define delta_E_res 0.0054827 // gamma/2

double E_kin = 0.081446;
// double E_kin = 0.8;
#define dt_anim 2.2
#define dt 0.05 //Di più di 0.05 scazza...
#define T_MAX 1500

double dr = 160./N;
double r_start = 20;
double width = 3.2;
double k;

FILE *fp;
char *nome_file = "matplot/dati/well_probability8.csv";



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
		double potential = reduced_lenjon_pot(r,null)+alpha*l*(l+1)/pow(r,2);
		
		if (potential>100)
			return(100);
		else
			return(potential);
	} else {
		return(0);
	}
}

double complex_potential_fun(double r, void *input)
{
	// Parabola:
	double a = 0.0005, r_start = 55;
	if (r>r_start)
	{
		double b = -2*a*r_start;
		double c = b*b/(4*a);
		return(a*r*r + b*r + c);
	} else
	{
		return(0);
	}

	// // Retta:
	// double a = 0.01, b = 55.;
	// if (r>b)
	// 	return(a*(r-b));
	// else
	// 	return(0.);
}

//--------------------------------------------
// Routines necessary to mean_energy:

void real_part_update(FFT_var *vars)
{
	vars->f_r=vars->old_r, vars->f_i=0;
}
void imaginary_part_update(FFT_var *vars)
{
	vars->f_r=vars->old_i, vars->f_i=0;
}

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
	// DERIVATIVE:
	double L = dr*N;
	FFT_vect derivata1_real, derivata1_imag, derivata2_real, derivata2_imag;
	FFT_vect imag_wave_function, real_wave_function;
	copy_FFT(*wave_function,&imag_wave_function);
	copy_FFT(*wave_function,&real_wave_function);
	double null[1];
	update_FFT(&imag_wave_function,null,imaginary_part_update);
	update_FFT(&real_wave_function,null,real_part_update);

	// Real part:
	derivata_FFT(N,L,real_wave_function,&derivata1_real);
	derivata_FFT(N,L,derivata1_real,&derivata2_real);

	// Imag part:
	derivata_FFT(N,L,imag_wave_function,&derivata1_imag);
	derivata_FFT(N,L,derivata1_imag,&derivata2_imag);

	// Merge:
	FFT_vect final_derivative2;
	copy_FFT(*wave_function,&final_derivative2);
	for (int i = 0; i < N; ++i)
	{
		REAL(final_derivative2.main,i)=REAL(derivata2_real.main,i);
		IMAG(final_derivative2.main,i)=REAL(derivata2_imag.main,i);
		// IMAG(final_derivative2.main,i)=0;
	}

	// Coniugato:
	FFT_vect coniug_wave_function;
	copy_FFT(*wave_function,&coniug_wave_function);
	update_FFT(&coniug_wave_function,null,coniugato_wave_fun);


	// Prodotto primo integrale: (modifico derivata2)
	product(&final_derivative2,coniug_wave_function);

	// Prepare integrands:
	Fun_1D integrand_1 = init_fun(1e6), integrand_2 = init_fun(1e6), integrand_3 = init_fun(1e6);
	integrand_1.len=integrand_2.len=integrand_3.len = N;
	for (int i=0; i < N; i++)
	{
		double r = i*dr;
		integrand_1.x[i]=integrand_2.x[i]=integrand_3.x[i]=r;
		integrand_1.y[i]=-alpha*REAL(final_derivative2.main,i);
		integrand_3.y[i]=-alpha*IMAG(final_derivative2.main,i);
		integrand_2.y[i]=pow(MODULO(wave_function->main,i),2)*return_val_ordered(r,&potential);
	}

	// Do integrals:
	double estremi[2]={0,L};
	double integral1 = int_trapezi(estremi,dr,&integrand_1,return_val_ordered,2);
	double integral3 = int_trapezi(estremi,dr,&integrand_3,return_val_ordered,2);
	if (integral3>0.01) {
		fprintf(stderr,"REALE = %g;  COMPLESSO = %g\n",integral1,integral3);
		getchar();
	}
	double integral2 = int_trapezi(estremi,dr,&integrand_2,return_val_ordered,2);

	// Free space:
	FFT_vect * vect_FFT[8] = {&coniug_wave_function,&derivata1_real,&derivata1_imag,&derivata2_real,&derivata2_imag,&imag_wave_function,&real_wave_function,&final_derivative2};
	free_FFTs(vect_FFT,8);
	Fun_1D * vect_Fun[2] = {&integrand_1,&integrand_2};
	free_funs(vect_Fun,2);

	// fprintf(stderr, "%g e %g\n", integral1, integral2);
	return(integral1+integral2);
}
//--------------------------------------------

double total_probability_FFT_vect(FFT_vect *funzione)
{
	// int len=stato_legato->len;
	double integrale = 1/2. * ( potenza(MODULO(funzione->main,0), 2) + 
		potenza(MODULO(funzione->main,N-1), 2)); 
	for (int i = 1; i < N-1; ++i)
		integrale += potenza(MODULO(funzione->main,i), 2);
	integrale *= dr;
	return(integrale);
}

double well_probability(FFT_vect *funzione)
{
	// int len=stato_legato->len;
	double integrale = 1/2. * ( potenza(MODULO(funzione->main,0), 2) + 
		potenza(MODULO(funzione->main,N-1), 2)); 
	for (int i = 1; i < N-1; ++i) {
		if (i*dr<=2)
			integrale += potenza(MODULO(funzione->main,i), 2);
		else
			break;
	}
	integrale *= dr;
	return(integrale);
}

void print_reciproc_packet(FFT_vect *packet, double anim_time)
{
	FFT_vect trasf;
	forward_FFT(*packet,&trasf);
	normalizzazione_FFT_vect(&trasf);

	FILE *f_trasf;
	f_trasf=fopen("matplot/dati/trasm_packet_res.csv", "w");	

	fprintf(f_trasf, "#WIDTH: %g and %g\n", width, 1/width);
	fprintf(f_trasf, "#k=%g\n", k);

	fprintf(f_trasf, "#Time: %g\n", anim_time);
	for (int i = 0; i < N; ++i)
	{
		fprintf(f_trasf,"%g\t%g\n",trasf.x_sym[i],MODULO(trasf.main,i));
	}
	fclose(f_trasf);

}

int check_split_op(check_struct_split_op *ck_val)
{
	FFT_vect *pacchetto=ck_val->funz_onda;
	int count = (int) ck_val->check_params[0];
	Fun_1D *potential = (Fun_1D*) ck_val->output;



	// Check if print:
	if (  ck_val->t>=count*dt_anim)
	{
		//Calculate mean energy:
		double mean_energy = mean_energy_fun(pacchetto,*potential);

		// Calculate normalization and prob.:
		double norma = total_probability_FFT_vect(pacchetto);
		double prob = well_probability(pacchetto);

		// Log:
		if (count==0) {
			fprintf(stderr, "\nStart average energy: %g (center: %g)",mean_energy, E_kin);
			fprintf(fp, "#Start average energy: %g (kinetic: %g)\n",mean_energy, E_kin);
			fprintf(fp, "#Width: %g; r_start: %g;\n",width,r_start);
			fprintf(fp, "#N=%d; dr=%g; R_MAX=%g\n", N,dr,dr*N);
		}
		if (count%30==0)
			fprintf(stderr, "\nTime analysed = %gs\t",ck_val->t);
		fprintf(stderr, ".");

		
		// Print reciproc:
		if (count==70) 
			print_reciproc_packet(pacchetto, count*dt_anim);


		//In the first places save the energy and normalization:
		printf("%g;%g;%g",mean_energy, norma, ck_val->t);
		for(int n=0;n<N;n++)
			printf(";%g",MODULO(pacchetto->main,n));
		printf("\n");

		// Save dates for probabality:
		fprintf(fp, "%g\t%g\n", ck_val->t, prob);

		ck_val->check_params[0]+=1;
	}

	// Check when stop:
	return( ck_val->t<T_MAX ? 0 : 1 );
}

void pacc_gauss(FFT_var *vars)
{
	// double correzione=N*dr/2;
	double k=vars->params[0];
	complex double pack_val;
	// pack_val = 20*cexp(I*k*vars->x)+cexp(-I*k*vars->x);
	pack_val = exp(-0.5* potenza((vars->x - r_start)/width, 2))* cexp(-I*k*vars->x);
	vars->f_r=creal(pack_val), vars->f_i=cimag(pack_val);	
}

int main(int argc, char const *argv[])
{
	// Variabili:
	k = sqrt(E_kin/alpha - 1/(width*width));
	double k_plane_ris_wave = sqrt(E_kin/alpha);
	// double width_2 = sqrt( 1/(E_kin/alpha-pow(k_plane_ris_wave,2)) );
	double min_width = sqrt(alpha/delta_E_res);
	// fprintf(stderr, "k=%g; \nk onda piana risonante=%g\nwidth2=%g\n",k, k_plane_ris_wave,width_2 );
	// return(0);
	
	k = k_plane_ris_wave;
	fprintf(stderr, "WIDTH: %g and %g\n", width, 1/width);
	fprintf(stderr, "k=%g\n", k);

	// Salvo potenziale:
	double pot_params[]={dr};
	Fun_1D potential = create_fun(0,dr*(N-1),dr,pot_params,potential_fun);

	// Stampo potenziale reale?
	if (argc==2)
	{
		print_fun(potential);
		return(0);
	}

	// Save imaginary potential:
	double null[1];
	Fun_1D complex_pot = create_fun(0,dr*(N-1),dr,null,complex_potential_fun);

	// Stampo potenziale immaginario?
	if (argc==3)
	{
		print_fun(complex_pot);
		return(0);
	}	

	fp=fopen(nome_file, "w");

	// Creo pacchetto:
	FFT_vect pacchetto, pacchetto_finale;
	double params[]={k};
	create_FFT(N,dr,&pacchetto,params,pacc_gauss);
	normalizzazione_FFT_vect(&pacchetto);
	copy_FFT(pacchetto, &pacchetto_finale);

	// // Print reciproc:
	print_reciproc_packet(&pacchetto,0);
	return(0);

	// Log:
	fprintf(stderr, "\n############### START ##############");

	// Evolvo nel tempo:
	double check_params[1]={0};
	//Al posto di output inserisco il potenziale che mi serve per la media:
	split_operator_complex_potential(N,dt,dr*N,&pacchetto_finale,&potential,&potential,return_val_ordered,&complex_pot,return_val_ordered,check_params,check_split_op);

	// Log:
	fprintf(stderr, "\n");
	fprintf(stderr, "################ END ###############\n\n");
	fclose(fp);
	return 0;
}
