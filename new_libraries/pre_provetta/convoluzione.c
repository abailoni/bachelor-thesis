// conv basta dichiararla (non è necessario allocare...)
void convoluzione(FFT_vect funzione1, FFT_vect funzione2, FFT_vect *conv)
{
	// Variabili:
	int N_tot= funzione2.N;
	double L  = funzione1.x[N_tot-1];

	// Trasformo:
	FFT_vect trasformata1, trasformata2;
	forward_FFT(funzione1,&trasformata1), forward_FFT(funzione2,&trasformata2);

	// Eseguo prodotto:
	product(&trasformata1,trasformata2);

	// Ritrasformo il prodotto:
	inverse_FFT(trasformata1,conv);
	for (int i = 0; i < N_tot; ++i)
		REAL(conv->main,i)*=L/N_tot;
}
