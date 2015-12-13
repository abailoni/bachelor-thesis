/////////////////////
// Varie generali: //
/////////////////////


double potenza(double x, int N)
{
	double pot=x;
	for (int i = 2; i <= N; ++i)
		pot*=x;
	return(pot);
}

void copy_array(double *array, double *new_array, int N)
{
	for (int i = 0; i < N; ++i)
		new_array[i]=array[i];
}


////////////////////////////////
// FUNZIONI DI UNA VARIABILE: //
////////////////////////////////

typedef struct {
	double *x;
	double *y;
	int len; //è il numero di elementi nel vettore e non il massimo indice!!
	double *params;
} Fun_1D;

void init_fun_2(double malloc_sz, Fun_1D *funzione)
{	
	funzione->x=(double*) malloc(malloc_sz*sizeof(double)), funzione->y=(double*) malloc(malloc_sz*sizeof(double)); 
}

Fun_1D init_fun(double malloc_sz)
{	
	Fun_1D funzione;
	funzione.x=(double*) malloc(malloc_sz*sizeof(double)), funzione.y=(double*) malloc(malloc_sz*sizeof(double)); 
	return(funzione);
}

void free_fun(Fun_1D *funzione)
{	free(funzione->x), free(funzione->y); }

void free_funs(Fun_1D * vectOfFuns[], int lenght)
{	
	for (int i = 0; i < lenght; ++i)
	{
		Fun_1D *funzione = vectOfFuns[i];
		free(funzione->x), free(funzione->y); 
	}
}

// Crea una struttura funzione con i valori di una funzione passata:
Fun_1D create_fun(double x_start, double x_end, double dx, void *params, double (*fun)(double,void*) )
{
	int len = (int) ( (x_end-x_start)/dx +1 );
	Fun_1D funzione;
	funzione.x=(double*) malloc(len*sizeof(double)), funzione.y=(double*) malloc(len*sizeof(double)), funzione.len=len;
	for (int i = 0; i < len; ++i)
		funzione.x[i]=x_start+i*dx, funzione.y[i]=fun(x_start+i*dx,params);
	return(funzione);
}

// La funzione restituisce la y corrispondente al valore x passato in input. Essa funziona solo se le x sono ordinate nel vettore funz.x[]! Inoltre è ottimizzata per il caso in cui le x siano equidistanziate (gran parte dei casi) e in tal caso trova subito il valore senza scorrerli tutti. 
double return_val_ordered(double x, void *input)
{
	Fun_1D *funzione=(Fun_1D*) input;
	int len = funzione->len;

	// Controllo prima se incremento è costante e in caso trovo subito il valore:
	
	double dx = funzione->x[1]-funzione->x[0];
	int stimato= (int) ( fabs((x-funzione->x[0])/dx) );
	if (stimato < len)
	{
		if ( fabs(x-funzione->x[stimato])<dx/2. )
			return( funzione->y[stimato] );
	}

	// Altrimenti divido l'intervallo a due a due:
	int start=0, end=len-1, new, diff;
	do
	{
		new=(end+start)/2;
		if (funzione->x[new]>x) 
			end=new;
		else
			start=new;
		diff=end-start;
	} while (diff>1);

	// Scelgo tra gli ultimi due valori trovati:
	if (  fabs(funzione->x[start]-x) < fabs(funzione->x[end]-x) )
		return(funzione->y[start]);
	else 
		return(funzione->y[end]);
}

double return_val_ordered_2(double x, void *input)
{
	Fun_1D *funzione=(Fun_1D*) input;
	// int len = funzione->len;

	for (int i = 0; 1; ++i)
	{
		if (funzione->x[i]>x)
			return(funzione->y[i]);
	}

	// // Altrimenti divido l'intervallo a due a due:
	// int start=0, end=len-1, new, diff;
	// do
	// {
	// 	new=(end+start)/2;
	// 	if (funzione->x[new]>x) 
	// 		end=new;
	// 	else
	// 		start=new;
	// 	diff=end-start;
	// } while (diff>1);

	// // Scelgo tra gli ultimi due valori trovati:
	// if (  fabs(funzione->x[start]-x) < fabs(funzione->x[end]-x) )
	// 	return(funzione->y[start]);
	// else 
	// 	return(funzione->y[end]);
}


void print_fun(Fun_1D funzione)
{
	if (funzione.len==0)
		fprintf(stderr, "ATTENZIONE: non è presente la lunghezza della funzione.\n");
	for(int i=0; i<funzione.len; i++)
		printf("%g\t%g\n",funzione.x[i], funzione.y[i]);
}

void print_fun_2(Fun_1D funzione)
{
	int i=0;
	while (funzione.x[i]==funzione.x[i] && funzione.y[i]==funzione.y[i]) {
		printf("%g\t%g\n",funzione.x[i], funzione.y[i]);
		i++;
	}
}

void copy_fun(Fun_1D funz1, Fun_1D *funz2)
{ // Attenzione: la nuova funzione deve essere già stata inizializzata!
	funz2->len=funz1.len;
	for (int i = 0; i < funz1.len; ++i)
		funz2->x[i]=funz1.x[i], funz2->y[i]=funz1.y[i];
}
