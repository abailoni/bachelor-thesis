
//+++++++++++++++++++++++++++++++++++++++++++
// DEFINIZIONI UTILI PER TRASFORMATE: 
//+++++++++++++++++++++++++++++++++++++++++++

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define MODULO(z,i) (sqrt( potenza((z)[2*(i)],2) + potenza((z)[2*(i)+1],2) ))

//------------------------------------------------
// Struttura FFT_vect:
// contiene tutti i dati per le trasformate.
//------------------------------------------------
typedef struct {
	double *main; //Vettore lungo 2N da dare in pasto a FFT
	double *x; //vettore [N] con x oppure k
	double *x_sym; //vettore [N] con x anche negative centrate in zero
	int N;
	double step;
} FFT_vect;

//------------------------------------------------
// Struttura FFT_var:
// quella che solitamente hanno a disposizione le funzioni update.
//------------------------------------------------
typedef struct {
	int n;
	int N_tot;
	double x;
	double old_r;
	double old_i;
	double f_r;
	double f_i;
	double *params;
	void *input; //possibili funzioni o altro
} FFT_var;




//+++++++++++++++++++++++++++++++++++++++++++
// FUNZIONI UTILI: 
//+++++++++++++++++++++++++++++++++++++++++++

//------------------------------------------------
// Libera vettori allocati con malloc:
//------------------------------------------------
void free_FFT(FFT_vect *fft)
{
	free(fft->main),free(fft->x),free(fft->x_sym);
}

void free_FFTs(FFT_vect * vectorOfFFT[], int lenght)
{
	for (int i = 0; i < lenght; ++i)
	{
		FFT_vect *fft = vectorOfFFT[i];
		free(fft->main),free(fft->x),free(fft->x_sym);
	}
}


//------------------------------------------------
// Crea FFT_vect data una certa funzione:
//------------------------------------------------
void create_FFT(int N, double dx, FFT_vect *vector, double *params, void (*fun)(FFT_var*))
{
	// N=numero di elementi (non indice massimo)
	vector->main=(double*) malloc(2*N*sizeof(double));
	vector->x=(double*) malloc(N*sizeof(double));
	vector->x_sym=(double*) malloc(N*sizeof(double));
	FFT_var vars;
	vars.params=params, vector->step=dx, vector->N=N, vars.N_tot=N;
	for(vars.n=0;vars.n<N;vars.n++)
	{
		vars.x=vars.n*dx;
		fun(&vars);
		vector->x[vars.n]=vars.n*dx, vector->main[2*vars.n]=vars.f_r, vector->main[2*vars.n+1]=vars.f_i;
		vector->x_sym[vars.n]= dx * (vars.n <= N/2 ? vars.n : -(N-vars.n));
		// printf("%d %g\n",vars.n,main[2*vars.n]);
	}
}

//------------------------------------------------
// Modifica un FFT_vect applicandoci una funzione update che ha a disposizione la struttura FFT_var:
//------------------------------------------------
void update_FFT(FFT_vect *vector, double *params, void (*fun)(FFT_var*))
{
	FFT_var vars;
	vars.params=params;
	int N=vector->N;
	vars.N_tot=N;
	for(int k=0;k<N;k++)
	{
		vars.x=vector->x[k], vars.n=k, vars.old_r=vector->main[2*k], vars.old_i=vector->main[2*k+1];
		fun(&vars);
		vector->main[2*k]=vars.f_r, vector->main[2*k+1]=vars.f_i;
	}
}

//------------------------------------------------
// Modifica vector1 moltpilicandolo per vector2 (con prodotto complesso):
//------------------------------------------------
void product(FFT_vect *vector1, FFT_vect vector2)
{
	for(int i=0; i<vector1->N; i++)
	{
		double real=REAL(vector1->main,i), imag=IMAG(vector1->main,i);
		REAL(vector1->main,i)=real*REAL(vector2.main,i) - imag*IMAG(vector2.main,i);
		IMAG(vector1->main,i)=real*IMAG(vector2.main,i) + imag*REAL(vector2.main,i);
	}
}

//------------------------------------------------
// Modifica vector1 moltpilicandolo per un vettore (lungo 2*N):
//------------------------------------------------
void product_array(FFT_vect *vector1, double *array)
{
	for(int i=0; i<vector1->N; i++)
	{
		double real=REAL(vector1->main,i), imag=IMAG(vector1->main,i);
		REAL(vector1->main,i)=real*REAL(array,i) - imag*IMAG(array,i);
		IMAG(vector1->main,i)=real*IMAG(array,i) + imag*REAL(array,i);
	}
}

//------------------------------------------------
// Copia semplicemente il primo FFT_vect nel secondo FFT_vect: 
//------------------------------------------------
void copy_FFT(FFT_vect old, FFT_vect *new)
{
	int N_tot = old.N;
	new->main=(double*) malloc(2*N_tot*sizeof(double));
	new->x=(double*) malloc(N_tot*sizeof(double));
	new->x_sym=(double*) malloc(N_tot*sizeof(double));
	new->N=old.N, new->step=old.step;
	for (int i = 0; i < old.N; ++i){
		IMAG(new->main,i)=IMAG(old.main,i), REAL(new->main,i)=REAL(old.main,i), new->x[i]=old.x[i], new->x_sym[i]=old.x_sym[i];	
	}
}

//------------------------------------------------
// Stampa modulo (mod=1), parte reale (mod=2) o immaginaria (mod=3): 
//------------------------------------------------
void print_FFT(int mod, FFT_vect funzione)
{
	for (int i = 0; i < funzione.N; ++i)
	{
		switch(mod) {
		case 1: 
			printf("%g\t%g\n",funzione.x[i],MODULO(funzione.main,i));
			break;
		case 2:
			printf("%g\t%g\n",funzione.x[i],REAL(funzione.main,i));
			break;
		case 3:
			printf("%g\t%g\n",funzione.x[i],IMAG(funzione.main,i));
			break;
		}

	}
}
//------------------------------------------------
// Stampa modulo (mod=1), parte reale (mod=2) o immaginaria (mod=3): 
//------------------------------------------------
void print_FFT_sym(int mod, FFT_vect funzione)
{
	for (int i = 0; i < funzione.N; ++i)
	{
		switch(mod) {
		case 1: 
			printf("%g\t%g\n",funzione.x_sym[i],MODULO(funzione.main,i));
			break;
		case 2:
			printf("%g\t%g\n",funzione.x_sym[i],REAL(funzione.main,i));
			break;
		case 3:
			printf("%g\t%g\n",funzione.x_sym[i],IMAG(funzione.main,i));
			break;
		}

	}
}
//------------------------------------------------
// Stampa modulo su una riga (in genere per animazioni): 
//------------------------------------------------
void print_FFT_anim(FFT_vect funzione)
{
	printf("%g",MODULO(funzione.main,0) );
	for (int i = 1; i < funzione.N; ++i)
		printf(";%g", MODULO(funzione.main,i));
	printf("\n");
}


//------------------------------------------------
//  ### TRASFORMATE ###
//------------------------------------------------
// per problemi di allocamento bisogna dichiarare il secondo FFT_vect in cui verrà memorizzata la trasformata ma non è necessario inserire nessun valore, verrà riempito da se.


// ---> DIRETTA <---
void forward_FFT(FFT_vect vector_x, FFT_vect *vector_k)
{
	int N=vector_x.N;
	vector_k->main=(double*) malloc(2*N*sizeof(double));
	vector_k->x=(double*) malloc(N*sizeof(double));
	vector_k->x_sym=(double*) malloc(N*sizeof(double));

	double dK=2*M_PI/(N*vector_x.step);
	copy_array(vector_x.main,vector_k->main,2*N);
	vector_k->N=N;
	vector_k->step=dK;

	// Calcolo le trasformate:
	gsl_fft_complex_radix2_forward(vector_k->main, 1, N);

	// Calcolo le posizioni k:
	for (int i = 0; i < N; ++i)
		vector_k->x[i]=dK*i, vector_k->x_sym[i]=dK * (i <= N/2 ? i : -(N-i));
}

// ---> INVERSA <---
void inverse_FFT(FFT_vect vector_k, FFT_vect *vector_x)
{
	int N=vector_k.N;
	vector_x->main=(double*) malloc(2*N*sizeof(double));
	vector_x->x=(double*) malloc(N*sizeof(double));
	vector_x->x_sym=(double*) malloc(N*sizeof(double));

	double dx=2*M_PI/(N*vector_k.step);
	copy_array(vector_k.main,vector_x->main,2*N);
	vector_x->N=N;
	vector_x->step=dx;

	// Calcolo le trasformate:
	gsl_fft_complex_radix2_inverse(vector_x->main, 1, N);

	// Calcolo le posizioni k:
	for (int i = 0; i < N; ++i)
		vector_x->x[i]=dx*i, vector_x->x_sym[i]=dx * (i <= N/2 ? i : -(N-i));
}

// ---> INVERSA <---
// Nel caso in cui il secondo vettore non sia da inizializzare:
void inverse_FFT_2(FFT_vect vector_k, FFT_vect *vector_x)
{
	int N=vector_k.N;

	copy_array(vector_k.main,vector_x->main,2*N);
	vector_x->N=N;

	// Calcolo le trasformate:
	gsl_fft_complex_radix2_inverse(vector_x->main, 1, N);

}

