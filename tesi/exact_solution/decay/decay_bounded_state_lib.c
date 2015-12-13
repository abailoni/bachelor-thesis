///
///
////////////////////////////////////////////////
// CREAZIONE PACCHETTO MEDIANTE STATO LEGATO: //
////////////////////////////////////////////////
///
///


double normalizzazione(Fun_1D *stato_legato)
{
	int len=stato_legato->len;
	double integrale = 1/2.*(stato_legato->y[0]*stato_legato->y[0]+stato_legato->y[len-1]*stato_legato->y[len-1]); 
	for (int i = 1; i < len-1; ++i)
		integrale += stato_legato->y[i]*stato_legato->y[i];
	integrale *= h;
	for (int i = 0; i < len; ++i)
		stato_legato->y[i]/=sqrt(integrale);
	return(integrale);
}

//////////////////////////////////////
// FUNZIONI NECESSARIE A numerov(): //
//////////////////////////////////////
#include "../../../new_libraries/numerov.c"

struct variables
{
	double E_rid;
	int check;
	Fun_1D *stato_legato;
} global_vars;

// funzione F(x) dell'equazione di schrodinger:
double F_schr(double x, void *input)
{
	double null[1];
	return(1/alpha*(potential_fun(x,null) - global_vars.E_rid));
}

// check_fun() per numerov in avanti:
int check_fun_avanti(struct numerov_check_struct *ck_val)
{
	int num_passi=ck_val->contatore;
	int num_max_passi=(int) (ck_val->ck_par[0]);
	Fun_1D *stato_legato=(Fun_1D*) ck_val->output;

	if ( num_passi>num_max_passi )
		return(1);
	else {
		stato_legato->x[num_passi]=ck_val->x, stato_legato->y[num_passi]=ck_val->psi[0];
		return(0);
	}
}

// check_fun() per numerov indietro:
int check_fun_indietro(struct numerov_check_struct *ck_val)
{
	int num_passi=ck_val->contatore, num_passi_tot= ck_val->ck_par[0]+ck_val->ck_par[1];
	int num_passi_ind=ck_val->ck_par[1];
	Fun_1D *stato_legato=(Fun_1D*) ck_val->output;

	if ( num_passi>num_passi_ind)
		return(1);
	else {
		stato_legato->x[num_passi_tot-num_passi]=ck_val->x, stato_legato->y[num_passi_tot-num_passi]=ck_val->psi[0];
		return(0);
	}
}

//////////////////////////////////////////
// FUNZIONI NECESSARIE A fun_roots():   //
//////////////////////////////////////////

// #include "../../new_libraries/fun_roots_MOD.c" //per poter stampare le funzioni dei "doppietti" è stata modificata la libreria

// Tale funzione data un energia E_rid calcola il "gap" in x_0:
double gap(double E_rid, void *input)
{
	//Variabili per numerov: 
	
	double x_0=1.3;
	
	// printf("%g-",E_rid );

	// Determino inizio e fine precisi (e numero di passi)
	int num_passi_av= (int)((x_0-x_start_av)/h), num_passi_ind;
	x_0=x_start_av+num_passi_av*h;
	for (int i = 1; 1; ++i)
		if ( x_0+h*i >= x_start_ind ) 
		{
			num_passi_ind=i;
			x_start_ind=x_0+h*i;
			break;
		}

	// Importo variabili:
	Fun_1D *stato_legato = global_vars.stato_legato;

	// Numerov in avanti:
	double start_val_fun[2]={1e-10,2e-10}, null[1];
	double check_params[]={num_passi_av, num_passi_ind};
	global_vars.E_rid=E_rid;
	numerov(x_start_av,h,start_val_fun,stato_legato,null,F_schr,check_params,check_fun_avanti);
	double u_avanti_r0=stato_legato->y[num_passi_av];


	// Numerov indietro:
	numerov_inverse(x_start_ind,h,start_val_fun,stato_legato,null,F_schr,check_params,check_fun_indietro);
	double u_indietro_r0=stato_legato->y[num_passi_av];

	
	// Raccordo le due funzioni:
	int passi_tot=num_passi_av+num_passi_ind;
	stato_legato->len=passi_tot+1;
	for (int j=num_passi_av; j<=passi_tot; j++)
		stato_legato->y[j]*=u_avanti_r0/u_indietro_r0;
	
	// Calcolo gap derivate prime in x_0:
	double F_r0=F_schr(x_0,null);
	double gap= (stato_legato->y[num_passi_av-1]+stato_legato->y[num_passi_av+1]-(2+pow(h,2)*F_r0)*u_avanti_r0) /h;
	return( gap );
}

///
///
//////////
// FINE //
//////////
///
///
