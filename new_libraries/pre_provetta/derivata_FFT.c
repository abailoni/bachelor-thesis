// INPUT:
// - N
// - L (periodicità)
// - funzione FFT_vect da derivare
// - puntatore a FFT_vect dove salvare la derivata

void update_trasf(FFT_var *vars)
{
	int k=vars->n, N_tot=vars->N_tot;
	double K_imag=2*M_PI/N_tot*k + (k <= N_tot/2 ? 0 : -2*M_PI);
	vars->f_r = -vars->old_i * K_imag;
	vars->f_i = vars->old_r * K_imag ;
}

void derivata_FFT(int N, double L, FFT_vect input, FFT_vect *output)
{
	// Trasformo:
	FFT_vect trasformata;
	forward_FFT(input,&trasformata);

	// Moltiplico:
	double null[1];
	update_FFT(&trasformata,null,update_trasf);

	// Antitrasformo:
	inverse_FFT(trasformata,output);
	for (int i = 0; i < N; ++i)
		REAL(output->main,i)*=N/L;
}
