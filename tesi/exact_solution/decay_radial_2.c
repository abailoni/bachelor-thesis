// Try to derivate the complex function in two different parts: before the real part and only then the imaginary part.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_fft_complex.h>
#include <time.h>

#define zeta 0.035171
#define alpha 0.035171

#include "../../new_libraries/varie.c"
#include "../../new_libraries/FFT_library.c"
#include "../../new_libraries/split_operator_complex_potential.c"
#include "../../new_libraries/int_trapezi.c"
#include "../../new_libraries/derivata_FFT.c"
#include "../../new_libraries/zero_bisezione.c"



////////////////////////////
/// VARIOUS COSTANTS:
////////////////////////////

#define N 8192
// #define N 4096
// #define E_resonance 0.081446
#define E_resonance 0.09145
#define delta_E_res 0.00560675
// #define E_resonance 0.07


double E_rid = 0.;
#define dt_anim 0.6
#define dt 0.05 //Di più di 0.05 scazza...
#define T_MAX 250

double wanted_energy = E_resonance - 0.003;

double x_start_av=0.7, x_start_ind=1.4; //Per stato legato
double dr = 80./N;
double h = 80./N;
double r_start = 1.8;
double width = 0.66;
int packet_power = 6;
double k;
FILE *F_PROB, *F_ANIM, *F_STATE;
Fun_1D *well_prob;
char *nome_file = "matplot/dati/decay_time.csv";


////////////////////////////
/// ROUTINES:
////////////////////////////

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
		
		if (potential>1000)
			return(1000);
		else
			return(potential);
	} else {
		return(0);
	}
}

double complex_potential_fun(double r, void *input)
{
	// Parabola:
	double a = 0.005, r_start = 55;
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


////////////////////////////////////////////////
// CREAZIONE PACCHETTO MEDIANTE STATO LEGATO: //
////////////////////////////////////////////////
#include "decay/decay_bounded_state_lib.c"


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
	}

	// Coniugato:
	FFT_vect coniug_wave_function;
	copy_FFT(*wave_function,&coniug_wave_function);
	update_FFT(&coniug_wave_function,null,coniugato_wave_fun);


	// Prodotto primo integrale: (modifico derivata2)
	product(&final_derivative2,coniug_wave_function);

	// Prepare integrands:
	Fun_1D integrand_1 = init_fun(1e6), integrand_2 = init_fun(1e6);
	integrand_1.len=N, integrand_2.len=N;
	for (int i=0; i < N; i++)
	{
		double r = i*dr;
		integrand_1.x[i]=r, integrand_2.x[i]=r;
		integrand_1.y[i]=-alpha*REAL(final_derivative2.main,i);
		integrand_2.y[i]=pow(MODULO(wave_function->main,i),2)*return_val_ordered(r,&potential);
	}

	// Do integrals:
	double estremi[2]={0,L};
	double integral1 = int_trapezi(estremi,dr,&integrand_1,return_val_ordered,2);
	double integral2 = int_trapezi(estremi,dr,&integrand_2,return_val_ordered,2);

	// Free space:
	FFT_vect * vect_FFT[8] = {&coniug_wave_function,&derivata1_real,&derivata1_imag,&derivata2_real,&derivata2_imag,&imag_wave_function,&real_wave_function,&final_derivative2};
	free_FFTs(vect_FFT,8);
	Fun_1D * vect_Fun[2] = {&integrand_1,&integrand_2};
	free_funs(vect_Fun,2);

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

double decay_probability_FFT_vect(FFT_vect *funzione)
{
	// int len=stato_legato->len;
	double integrale = 1/2. * ( potenza(MODULO(funzione->main,0), 2) + 
		potenza(MODULO(funzione->main,N-1), 2)); 
	for (int i = 1; i < N-1; ++i) {
		if (i*dr<=3.7)
			integrale += potenza(MODULO(funzione->main,i), 2);
		else
			break;
	}
	integrale *= dr;
	return(integrale);
}

//--------------------------------------------
// WELL PROBABILITY AND INTEGRALS:
//--------------------------------------------

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

//--------------------------------------------
// CHECK SPLIT OPERATOR:
//--------------------------------------------

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

		// Calculate normalization and decay:
		double norma = total_probability_FFT_vect(pacchetto);

		// Log:
		if (count==0) {
			fprintf(stderr, "\nStart average energy: %g (res: %g)",mean_energy, E_resonance);
			// Prob. file: 
			fprintf(F_PROB, "#Start average energy: %g (kinetic: %g)\n",mean_energy, E_resonance);
			fprintf(F_PROB, "#Width: %g; k_start: %g; r_start: %g;\n",width,k,r_start);
			fprintf(F_PROB, "#N=%d; dr=%g; R_MAX=%g\n", N,dr,dr*N);
			fprintf(F_PROB, "%g\t%g\n",width,k);
			
			// Anim. file: 
			fprintf(F_ANIM, "#Start average energy: %g (kinetic: %g)\n",mean_energy, E_resonance);
			fprintf(F_ANIM, "#Width: %g; k_start: %g; r_start: %g;\n",width,k,r_start);
			fprintf(F_ANIM, "#N=%d; dr=%g; R_MAX=%g\n", N,dr,dr*N);
			// fprintf(F_ANIM, "%g;%g;%g\n",width,k);

		}
		if (count%30==0)
			fprintf(stderr, "\nTime analysed = %gs\t",ck_val->t);
		fprintf(stderr, ".");


		//In the first places save the energy and normalization:
		fprintf(F_ANIM,"%g;%g;%g",mean_energy, norma, ck_val->t);
		for(int n=0;n<N;n++)
			fprintf(F_ANIM,";%g",MODULO(pacchetto->main,n));
		fprintf(F_ANIM,"\n");

		// // Save dates for decay:
		// double decay = decay_probability_FFT_vect(pacchetto);
		// fprintf(F_PROB, "%g\t%g\n", ck_val->t, decay);
		
		// Save probabality:
		double prob = well_probability(pacchetto);
		fprintf(F_PROB, "%g\t%g\n", ck_val->t, prob);
		well_prob->x[count]=ck_val->t, well_prob->y[count]=prob;
		well_prob->len=count+1;

		// Check if the well is empty:
		if (ck_val->t>500 && prob<1e-4)
			return(1);
		
		ck_val->check_params[0]+=1;
	}

	// Check when stop:
	return( ck_val->t<T_MAX ? 0 : 1 );
}


void create_packet(FFT_var *vars)
{
	// double correzione=N*dr/2;
	// double k=vars->params[0];
	complex double pack_val;
	// pack_val = 20*cexp(I*k*vars->x)+cexp(-I*k*vars->x);
	// pack_val = exp(- 0.5*potenza((vars->x - r_start)/width, packet_power))* cexp(-I*k*vars->x);
	if (vars->x > x_start_av+0.001 && vars->x < x_start_ind-0.001) {
		// printf("Ciao! x=%g\n",vars->x);
		pack_val = return_val_ordered(vars->x,global_vars.stato_legato);
	}
	else
		pack_val = 0;
	vars->f_r=creal(pack_val), vars->f_i=cimag(pack_val);	
}

double energia_fondamentale(double x_max, void *input)
{
	x_start_ind=x_max;

	// Calcolo stati:
	double estremi_E[2];
	double dE = 1e-4, E_start = 0.02;
	estremi_E[0]=E_start;
	double null[1];
	for (int i = 1; 1; ++i)
	{
		estremi_E[0]= E_start+dE*(i-1), estremi_E[1]= E_start+dE*i;
		double first=gap(estremi_E[0],null), second=gap(estremi_E[1],null);
		// printf("#Energia: %g --> %g e %g\n",estremi_E[0],first,second);

		if (first*second<0) 
		{		
			// fprintf(stderr,"#Energy = %g\n", (estremi_E[0]+estremi_E[1])/2 );
			normalizzazione(global_vars.stato_legato);
			break;
		}
	}
	return((estremi_E[0]+estremi_E[1])/2 - wanted_energy);
}

double fit_function(double x, void *input)
{
	double *params = (void*) input;
	double width=params[0], x_0=params[1], norm=params[2];
	return(norm*exp(- pow((x-x_0),2)/(2*width*width) -1/pow(x,7)));
}

void fit_packet(FFT_var *vars)
{
	double fit_params[3]={1.59476924886, 0.358531619177, 1.2400583615};
	// double fit_params[3]={1.1897596719, 1.05080845552, 1.34927301079};
	
	vars->f_r=fit_function(vars->x,fit_params), vars->f_i=0;	
}

//--------------------------------------------
//--------------------------------------------
// MAIN PROGRAM:
//--------------------------------------------
//--------------------------------------------

int main(int argc, char const *argv[])
{
	////////////////////////////////
	// POTENTIAL INITIALIZATION: //
	////////////////////////////////

	// Save real potential:
	double pot_params[]={dr};
	Fun_1D potential = create_fun(0,dr*(N-1),dr,pot_params,potential_fun);
	
	// Save imaginary potential:
	double null[1];
	Fun_1D complex_pot = create_fun(0,dr*(N-1),dr,null,complex_potential_fun);

	// Print real potential?
	if (argc==2) {
		print_fun(potential);
		return(0);
	}
	// Print imaginary potential?
	if (argc==3) {
		print_fun(complex_pot);
		return(0);
	}	

	/////////////////////////
	// INITIAL VARIABLES: //
	/////////////////////////

	double min_width = sqrt(alpha/delta_E_res);

	// Random number:
	srand(time(NULL));
	double random_num = ((double) rand() / (RAND_MAX));
	random_num*=1e4;
	// random_num=9999;
	fprintf(stderr,"## NUMBER: %4.4g ##\n",random_num );

	/////////////////////////////////////
	// COMPUTE TIME-LIFE decay-state: //
	/////////////////////////////////////
	for (int i = 0; i<1; ++i)
	{
		// Variables:
		Fun_1D local_probabilty=init_fun(1e5);
		well_prob = &local_probabilty;

		// Open files:
		char name_file_prob[80], name_file_anim[80], name_file_decay[80];
		sprintf(name_file_prob,"matplot/dati/prob_decay/prob_%4.4g_decay.csv",random_num); 
		sprintf(name_file_anim,"matplot/dati/anim_decay/anim_%4.4g_decay.csv",random_num); 
		F_PROB=fopen(name_file_prob, "w");
		F_ANIM=fopen(name_file_anim, "w");
		sprintf(name_file_decay,"matplot/dati/decay_states/decay_%4.4g.csv",random_num); 
		F_STATE=fopen(name_file_decay, "w");
		fprintf(F_PROB, "#NUMBER: %4.4g\n",random_num);
		fprintf(F_ANIM, "#NUMBER: %4.4g\n",random_num);

		////////////////////////////
		// Create gauss packet: //
		////////////////////////////
		
		fprintf(stderr, "Computing minimum energy...\n");
		// Inizializzo variabili:
		Fun_1D stati_legati = init_fun(1e6);
		global_vars.stato_legato=&stati_legati;


		// Cerco energia minima:
		double estremi_min[]={1.5,7}, prec_min=1e-5;
		double x_min = bisezione(estremi_min,prec_min,null,energia_fondamentale);

		// Calcolo stato iniziale:
		double E_min = energia_fondamentale(x_min,null) + wanted_energy;
		fprintf(stderr, "#Minimum energy: %g; start: %g\n",E_min, x_start_ind);



		FFT_vect pacchetto, pacchetto_finale;
		double params[]={k};
		create_FFT(N,dr,&pacchetto,params,create_packet);
		normalizzazione_FFT_vect(&pacchetto);
		copy_FFT(pacchetto, &pacchetto_finale);

		for (int i = 0; i < N; ++i)
		{
			fprintf(F_STATE, "%g\t%g\n",pacchetto.x[i],MODULO(pacchetto.main,i));
		}


		// Evolution in time:
		fprintf(stderr, "\n########### Decay state: x_min=%g ###########",x_min);
		double check_params[1]={0};
		split_operator_complex_potential(N,dt,dr*N,&pacchetto_finale,&potential,&potential,return_val_ordered,&complex_pot,return_val_ordered,check_params,check_split_op);
		// (al posto di output inserisco il potenziale che mi serve per la media)
		fprintf(stderr, "\n");
		
		// Compute life-time:
		int len = well_prob->len;
		double life_time = 1/2. * (well_prob->y[0] + well_prob->y[len-1]); 
		for (int j = 1; j < len; ++j)
			life_time += well_prob->y[j];
		life_time *= fabs(well_prob->x[2]-well_prob->x[1]);
		printf("%g;%g;%g;%g\n",width,life_time,k,E_resonance);
		fprintf(stderr, "Time of life: %g\n",life_time);
		fprintf(stderr, "####### TIME OF LIFE: %g #########\n\n",life_time);

		free_FFT(&pacchetto), free_FFT(&pacchetto_finale);
		free_fun(&local_probabilty);
		fclose(F_ANIM), fclose(F_PROB), fclose(F_STATE);
	}



	///////////////////////////////////
	// COMPUTE TIME-LIFE fit-state: //
	///////////////////////////////////

	for (int i = 0; i<1; ++i)
	{
		// Variables:
		Fun_1D local_probabilty=init_fun(1e5);
		well_prob = &local_probabilty;

		// Open files:
		char name_file_prob[80], name_file_anim[80];
		sprintf(name_file_prob,"matplot/dati/prob_decay/prob_fit_%4.4g_fit.csv",random_num); 
		sprintf(name_file_anim,"matplot/dati/anim_decay/anim_%4.4g_fit.csv",random_num); 
		F_PROB=fopen(name_file_prob, "w");
		F_ANIM=fopen(name_file_anim, "w");
		// F_STATE=fopen("matplot/dati/DECAY_STATE.csv", "w");
		fprintf(F_PROB, "#NUMBER: %4.4g\n",random_num);
		fprintf(F_ANIM, "#NUMBER: %4.4g\n",random_num);

		////////////////////////////
		// Create fit state:      //
		////////////////////////////
		
		FFT_vect pacchetto, pacchetto_finale;
		double params[]={k};
		create_FFT(N,dr,&pacchetto,params,fit_packet);
		normalizzazione_FFT_vect(&pacchetto);
		copy_FFT(pacchetto, &pacchetto_finale);


		// Evolution in time:
		fprintf(stderr, "\n########### Fit state ###########");
		double check_params[1]={0};
		split_operator_complex_potential(N,dt,dr*N,&pacchetto_finale,&potential,&potential,return_val_ordered,&complex_pot,return_val_ordered,check_params,check_split_op);
		// (al posto di output inserisco il potenziale che mi serve per la media)
		fprintf(stderr, "\n");
		
		// Compute life-time:
		int len = well_prob->len;
		double life_time = 1/2. * (well_prob->y[0] + well_prob->y[len-1]); 
		for (int j = 1; j < len; ++j)
			life_time += well_prob->y[j];
		life_time *= fabs(well_prob->x[2]-well_prob->x[1]);
		printf("%g;%g;%g;%g\n",width,life_time,k,E_resonance);
		fprintf(stderr, "Time of life: %g\n",life_time);
		fprintf(stderr, "####### TIME OF LIFE: %g #########\n\n",life_time);

		free_FFT(&pacchetto), free_FFT(&pacchetto_finale);
		free_fun(&local_probabilty);
		fclose(F_ANIM), fclose(F_PROB);
	}

	// Log:
	fprintf(stderr, "\n");
	fprintf(stderr, "################ END ###############\n\n");

	return 0;
}
