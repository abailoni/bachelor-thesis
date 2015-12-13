// INPUT:
// - N
// - dt
// - L
// - funzione d'onda iniziale (verrà poi aggiornata mano mano ad ogni passo)
// - void *output
// - parametri potenziale (void)
// - potenziale
// - parametri check (double)
// - funzione check
// 
// ATTENZIONE: PER IL MOMENTO E' STATO REALIZZATO SUPPONENDO H_TAGL=1 E M=1

typedef struct {
	int N_tot;
	double dt; 
	double L;
	double t;
	FFT_vect *funz_onda;
	void *output; //Vettore dove inserire cose da esportare (double, funzione, FFT_vect, ecc...)
	double *check_params;
} check_struct_split_op;

void split_operator(int N, double dt, double L, FFT_vect *funz_onda, void *output, void *pot_params, double (*potenziale)(double,void*), double *check_params, int (*check_fun)(check_struct_split_op*) )
{
	double dx=L/N;
	check_struct_split_op ck_val={N,dt,L,0,funz_onda,output,check_params};

	while (check_fun(&ck_val)==0) {
		// Moltiplico per ro:
		double ro[2*N];

		for (int i = 0; i < N; ++i)
		{
			REAL(ro,i)= cos( dt*potenziale(funz_onda->x[i],pot_params) );
			IMAG(ro,i)= -sin( dt*potenziale(funz_onda->x[i],pot_params) );
		}
		product_array(funz_onda,ro);

		// Trasformo:
		FFT_vect trasformata;
		forward_FFT(*(funz_onda),&trasformata);

		// Moltiplico antitrasformata per xi:
		double xi[2*N];
		for (int i = 0; i < N; ++i)
		{
			double delta_K = 2*M_PI/(N*dx) * (i <= N/2 ? i : N - i);
			REAL(xi,i) =   cos(1/2. * dt * delta_K*delta_K);
			IMAG(xi,i) = - sin(1/2. * dt * delta_K*delta_K);
		}
		product_array(&trasformata,xi);

		// Antitrasformo:
		inverse_FFT_2(trasformata,funz_onda);

		// Avanzo dati:
		ck_val.t+=dt;
	}
}
