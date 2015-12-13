/*

//+++++++++++++++++++
// RK4()
//+++++++++++++++++++
Essa prende in input i seguenti argomenti:

	  - numero N di variabili
	  - passo dt
	  - t_start
	  - condizioni iniziali (array double)
	  - output (puntatore void)
	  - parametri F (array double)
	  - funzione F()
	  - parametri funzione check
	  - funzione check_fun()



//+++++++++++++++++++++
// check_fun():
//+++++++++++++++++++++
La funzione check_fun() ha il ruolo di controllare quando fermare l'algoritmo RK4 ed eventualemente fare dei printf o salvare i dati nella variabile output (puntatore a qualsiasi cosa).

Argomenti che si hanno a disposizione:
	- tutte le variabili presenti nella struttura check_struct_RK4 (visibile sotto)

Output da far restituire:
	- INT che se è zero dice a RK4 di continuare, altrimenti lo ferma.



//+++++++++++++++++++++
// F():
//+++++++++++++++++++++

Input della funzione:

		- vettore di output; [double]
		- t (tempo a cui calcolare);
		- vettore y (contiene i valori delle n funzioni); [double]
		- vettore parametri; [double]


 */ 

// NOTA:
// volendo penso si possa togliere come input anche y_start tanto le condizioni iniziali possono essere passate mediante l'uso della funzione F controllando quando si ha t=0...

struct check_struct_RK4
{
	int contatore;
	int N;
	double *dt;	
	double t;
	double *y;
	double *y_start;
	double t_start;
	double *params;
	void *output;
} ck_val;

void RK4(int N, double dt, double t_start, double *y_start, void *output, double *F_params, void (*F)(double*,double,double*,double*), double *check_params, int (*check_fun) (struct check_struct_RK4*))
{
	//Variabili:
	double K1[N], K2[N], K3[N], K4[N];
	double t=t_start, y[N], y_temp[N];
	for (int i = 0; i < N; ++i)
		y[i]=y_start[i];
	ck_val.dt=&dt, ck_val.y_start=y_start, ck_val.N=N, ck_val.t_start=t_start, ck_val.params=check_params, ck_val.contatore=0, ck_val.t=t_start, ck_val.y=y_start, ck_val.output=output;
	check_fun(&ck_val);

	//Inizio ad avanzare:
	do 
	{
		F(K1,t,y,F_params); //Trovo K1
		for (int i = 0; i < N; ++i)
			y_temp[i]=y[i]+dt/2.*K1[i];
		F(K2,t+dt/2.,y_temp,F_params); //K2
		for (int i = 0; i < N; ++i)
			y_temp[i]=y[i]+dt/2.*K2[i];
		F(K3,t+dt/2.,y_temp,F_params); //K3
		for (int i = 0; i < N; ++i)
			y_temp[i]=y[i]+dt*K3[i];
		F(K4,t+dt,y_temp,F_params); //K4
		// Trovo il nuovo dato:
		for (int i = 0; i < N; ++i)
			y[i]=y[i]+dt*(1./6*K1[i]+1./3*K2[i]+1./3*K3[i]+1./6*K4[i]);
		t+=dt;
		// Preparo la struttura per check_fun:
		ck_val.t=t, ck_val.y=y, ck_val.contatore++;
	} while ( check_fun(&ck_val)==0 );
	return;
}
