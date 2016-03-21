/*

//+++++++++++++++++++
// numerov()
//+++++++++++++++++++
Essa prende in input i seguenti argomenti:

	- x da cui partire; [double]
	- passo h; [double]
	- array di dimensione 2 con primi due valori della funzione; [double]
	- array di double dove restituire eventuali dati di output;
	- parametri F; [puntatore void]
	- funzione F;
	- parametri check_fun(); [array double]
	- funzione check_fun();

//+++++++++++++++++++++
// check_fun():
//+++++++++++++++++++++
La funzione check_fun() ha il ruolo di controllare quando fermare l'algoritmo di numerov ed eventualemente fare dei printf o salvare i dati nell'array output_vars.

Argomenti che si hanno a disposizione:
	- penultimo x in cui è stata calcolata psi;
	- passo h; [puntatore in modo che possa essere modificato in corsa]
	- array con ultimi due valori di psi calcolati;
	- array output_vars;
	- array di eventuali altri parametri esterni;

Output da far restituire:
	- INT che se è zero dice a verlet di continuare, altrimenti lo ferma.

 */ 

struct numerov_check_struct
{
	int contatore; //mi dice che al momento x=x_start+contatore*h (a che passo sono)
	double x;
	double *h; //passo che può essere modificato
	double *ck_par; //parametri passati (che possono essere modificati)
	double *psi; //valori in x e x+h
	void *output; //può puntare a qualsiasi cosa (anche una struttura Fun1D dove salvare i dati)
} ck_val;

void numerov(double x_start, double h, double start_val_fun[2], void *output, void *F_params, double (*F)(double,void*), double *check_params, int (*check_fun)(struct numerov_check_struct*) )
{
	double psi[2]={start_val_fun[0],start_val_fun[1]};
	double F_value[2]={F(x_start,F_params),F(x_start+h,F_params)};
	double x=x_start, new_F_value, new_psi;
	// Preparo dati e check iniziale:
	ck_val.contatore=0, ck_val.x=x, ck_val.h=&h, ck_val.ck_par=check_params, ck_val.psi=psi, ck_val.output=output;
	check_fun(&ck_val);
	do
	{
		//calcolo nuovo punto:
		new_F_value=F(x+2*h,F_params);
		new_psi=( (2+5.0/6.0*pow(h,2)*F_value[1])*psi[1]-(1-pow(h,2)/12.0*F_value[0])*psi[0] ) / (1-pow(h,2)/12.0*new_F_value);
		//aggiorno valori:
		psi[0]=psi[1], psi[1]=new_psi;
		F_value[0]=F_value[1], F_value[1]=new_F_value;
		x+=h;
		ck_val.x=x, ck_val.psi=psi, (ck_val.contatore)++;
	} while( check_fun(&ck_val)==0 );
	return;
}


//+++++++++++++++++++
// numerov_inverse()
//+++++++++++++++++++
// Fa le stesse cose di numerov solo che procede facendo dei passi all'indietro
void numerov_inverse(double x_start, double h, double start_val_fun[2], void *output, void *F_params, double (*F)(double,void*), double *check_params, int (*check_fun)(struct numerov_check_struct*) )
{
	double psi[2]={start_val_fun[0],start_val_fun[1]};
	double F_value[2]={F(x_start,F_params),F(x_start+h,F_params)};
	double x=x_start, new_F_value, new_psi;
	// Preparo dati e check iniziale:
	ck_val.contatore=0, ck_val.x=x, ck_val.h=&h, ck_val.ck_par=check_params, ck_val.psi=psi, ck_val.output=output;
	check_fun(&ck_val);
	do
	{
		//calcolo nuovo punto:
		new_F_value=F(x-2*h,F_params);
		new_psi=( (2+5.0/6.0*pow(h,2)*F_value[1])*psi[1]-(1-pow(h,2)/12.0*F_value[0])*psi[0] ) / (1-pow(h,2)/12.0*new_F_value);
		//aggiorno valori:
		psi[0]=psi[1], psi[1]=new_psi;
		F_value[0]=F_value[1], F_value[1]=new_F_value;
		x-=h;
		ck_val.x=x, ck_val.psi=psi, (ck_val.contatore)++;
	} while( check_fun(&ck_val)==0 );
	return;
}
