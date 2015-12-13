
/////////////////////////////////////////////////////////////////////
// TROVA N RADICI DI UNA FUNZIONE USANDO SOLO IL METODO DI NEWTON:
/////////////////////////////////////////////////////////////////////

/*
// *****************************
// INPUT fun_roots(): 
// *****************************
 - estremi in cui cercare gli zeri (double[2])
 - numero di zeri da cercare (int)
 - precisione con cui calcolare gli zeri (double)
 - array double di 2N elementi in cui verranno inseriti i dati relativi agli N zeri trovati
 - tempo massimo di ricerca per zero (double)
 - array double con eventuali parametri da passare alla funzione
 - puntatore alla funzione di cui calcolare le radici


// ******************************
// OUTPUT fun_roots(): 
// ******************************
 
 - int che in caso di successo restituisce il numero di tentativi in cui gli zeri sono stati trovati; in caso di insuccesso (fuori tempo massimo) restituisce un numero negativo che simboleggia l'opposto del numero di zeri trovati (sarà minore di N in valore assoluto)

I dati relativi agli zeri trovati vengono memorizzati nell'array di double passato in input. Per ogni zero sono salvati due valori che dipendono a seconda del caso:

 	- nel caso in cui venga trovato esattamente lo zero, allora come primo valore viene salvato x_0 e come secondo -989898.0;

 	- nel caso non si riesca a trovare esattamente lo zero ma solo un intervallo che lo comprende, come primo dato viene salvato l'estremo in cui f(x)>0 e come secondo l'altro estremo dove invece f(x)<0;

 	- nei precedenti casi si hanno sempre zeri in cui la funzione cambia segno; se invece se ne ha uno dove il segno della funzione rimane costante prima e dopo lo zero, tipo f(x)=x^2, allora il programma restituisce solo un valore approssimato dello zero e poi come secondo dato restituisce -999999.0

(I dati relativi al primo zero sono salvati nelle prime due posizioni dell'array, quelli del secondo nella terza e quarta posizione, ecc..)


// ********************
// FUNZIONAMENTO:	  
// ********************

Gli zeri vengono trovati mediante il metodo di Newton che si è scoperto essere molto più veloce ed efficace di quello della bisezione. Per non parlare del caso con più zeri dove quello della bisezione sarebbe di difficile implementazione e molto lento. 
Inoltre con tale metodo è possibile anche trovare zeri in cui la funzione non cambia segno.

Nell'ordine vengono richiamate le funzioni:

 	- fun_roots() seleziona gli x_0 (compresi tra gli estremi) da cui far partire l'algoritmo
 	
 	- newton() ricerca zeri mediante il METODO DI NEWTON partendo dagli x_0 che gli passa fun_roots()
 	
 	- infine la funzione zero_append() ogni volta che viene trovato uno zero accettabile controlla che non sia già stato trovato in precedenza e in caso lo memorizza.

*/

// *************
// CODICE:	  
// *************

#include <time.h>

#define N_PRECISIONI 20 //Numero di precisioni da compiere nel caso di zeri "complicati" tipo x^2 o x^4 che toccano soltanto l'asse ma non fanno cambiare il segno a f(x). Più tale numero è alto più aumenta il tempo richiesto all'algoritmo.


//++++++++++++++++++++++++++++
// newton():
//++++++++++++++++++++++++++++
// Se tutto va bene restituisce l'x in cui vi è lo zero; in caso di singolarità restituisce -999998; in caso di altri errori vari (punti usciti dall'intervallo, funzioni con asintoti o con derivate circa nulle, ecc..) restituisce -999999.
double newton(double x_0, double precisione, double* estremi, int* contatore, double TEMPO_MAX, void *params, double (*fun)(double,void*), clock_t *start)
{
	//Controllo di non essere uscito da estremi:
	if ((x_0<estremi[0]) || (x_0>estremi[1]))
		return(-999999);
	(*contatore)++;

	//--------------------------------------------
	// Trovo il primo punto "manualmente" con derivata approssimata:
	//--------------------------------------------
	double h=fabs(1e-4), fun_value[2];
	if (x_0==0)
		h=1e-5;
	fun_value[0]=fun(x_0,params);
	if (fabs(fun_value[0])==INFINITY) //controllo che non abbia preso una singolarità
		return(-999998);
	double derivata=(fun(x_0+h,params)-fun(x_0-h,params))/h;
	//Controllo su derivata necessario solo in casi estremi:
	if ( (derivata==0) || (fabs(derivata)<1e-3) )
		return(-999999);
	double x_1=x_0-fun_value[0]*(2/derivata);
	double x_values[]={x_0,x_1};
	//Controllo di non essere uscito da estremi:
	if ((x_1<estremi[0]) || (x_1>estremi[1]))
		return(-999999);
	fun_value[1]=fun(x_1,params);

	//--------------------------------------------
	// VADO AVANTI FINCHE' NON SUPERO PRECISIONE:
	//--------------------------------------------
	double new_x, new_fun_value;
	int contatore_ric=0;
	while (contatore_ric++<1e3)//Evito troppe ricorsioni (in casi complicati)
	{
		//Controllo tempo:
		if ( ((float)(clock()-(*start))/CLOCKS_PER_SEC)>TEMPO_MAX )
			return(-999999);
		// Controllo di non aver trovato esattamente lo zero:
		if (x_values[0]==x_values[1])
			return(x_values[0]);
		// Calcolo il nuovo x:
		new_x=x_values[1]-fun_value[1]*(x_values[1]-x_values[0])/(fun_value[1]-fun_value[0]);
		//Controllo di non essere uscito da estremi:
		if ((new_x<estremi[0]) || (new_x>estremi[1]))
			return(-999999);
		new_fun_value=fun(new_x,params);
		if (fabs(new_x-x_values[1])<precisione)
			break;
		//Non avendo raggiunto la precisione dimezzo l'intervallo e continuo:
		x_values[0]=x_values[1];
		x_values[1]=new_x;
		fun_value[0]=fun_value[1];
		fun_value[1]=new_fun_value;
	}

	//-------------------
	//ESPORTO IL VALORE:
	//-------------------
	if (fabs(new_fun_value)>0.1) 
		return(-999999); //evito possibili errori...
	//Controllo se ho trovato lo zero esatto:
	if (new_fun_value==0) 
		return(new_x);
	//Controllo se entro una precisione trovo lo zero:
	if ( fabs(new_fun_value+fun_value[1])!=(fabs(new_fun_value)+fabs(fun_value[1])))
		return(new_x);
	// in caso contrario mi avvicino ulteriormente allo zero:
	int sgn=1;
	if ((new_x-x_values[1])<0)
		sgn=-1;
	//Oltre le N_PRECISIONI concludo di avere uno zero che tocca l'asse tipo f(x)=x^2:
	for(int i=1;fabs(i)<=N_PRECISIONI;i+=sgn) 
	{
		double fun_value_better_x=fun(new_x+i*precisione,params);
		if (fabs(new_fun_value+fun_value_better_x)!=(fabs(new_fun_value)+fabs(fun_value_better_x)))
			return(new_x+i*precisione);
	}
	return(new_x); //sono nel caso tipo x^2
}


//++++++++++++++++++++++++++++
// zero_append():
//++++++++++++++++++++++++++++
// Controlla se zero1 è un nuovo zero e in caso lo inserisce in zeros[].
// In tal caso salva i due dati del nuovo zero.
void zero_append(double zero1, double *zeros, int *num_zeri, double *estremi, double precisione,void *params, double (*fun)(double,void*), clock_t *start)
{ 
	//--------------------------------------------
	// Controllo che lo zero non sia già stato trovato:
	//--------------------------------------------
	for(int i=0;i<2*(*num_zeri);i+=2) {
		int sgn=1;
		if ((zero1-zeros[i])<0)
			sgn=-1;
		if (  (sgn*(zeros[i]+sgn*precisione))>=(sgn*(zero1-sgn*precisione))  )
			return;
		if ((zeros[i+1]==-999999) &&  ((zeros[i]+N_PRECISIONI*precisione)>zero1) && ((zeros[i]-N_PRECISIONI*precisione)<zero1) )
			return;
	}
	
	//------------------------------
	// Caso in cui lo zero è nuovo:
	//------------------------------
	*start = clock(); //azzero il tempo
	double zero2=zero1+precisione;
	double val_fun[2]={fun(zero1,params),fun(zero2,params)};
	if (val_fun[0]==0){ //ho trovato esattamente lo zero
		zeros[(*num_zeri)*2]=zero1;
		zeros[(*num_zeri)*2+1]=-989898;
		*num_zeri+=1;
		return;
	}
	if (fabs(val_fun[0]+val_fun[1])==(fabs(val_fun[0])+fabs(val_fun[1]))) zero2=zero1-precisione;
	val_fun[1]=fun(zero2,params);
	if (fabs(val_fun[0]+val_fun[1])==(fabs(val_fun[0])+fabs(val_fun[1]))) 
	{ //sono nel caso tipo f(x)=x^2
		zeros[(*num_zeri)*2]=zero1;
		zeros[(*num_zeri)*2+1]=-999999;
		*num_zeri+=1;
		return;
	}

	//Ordino gli zeri nel modo giusto:
	if (val_fun[0]>0) {
		zeros[(*num_zeri)*2]=zero1;
		zeros[(*num_zeri)*2+1]=zero2;
	} else {
		zeros[(*num_zeri)*2+1]=zero1;
		zeros[(*num_zeri)*2]=zero2;
	}
	*num_zeri+=1;
	return;
}


//++++++++++++++++++++++++++++
// fun_roots():
//++++++++++++++++++++++++++++
int fun_roots(double *estremi, int N, double precisione, double *zeros, double TEMPO_MAX,  void *params, double (*fun)(double,void*))
{
	int num_zeri=0;
	int contatore=0;
	clock_t start = clock();
	// Uso gli estremi come punti di partenza:
	double newton_temp=newton(estremi[0],precisione,estremi,&contatore,TEMPO_MAX,params,fun,&start);
	if ((newton_temp!=-999999) && (newton_temp!=-999998))
		zero_append(newton_temp,zeros,&num_zeri,estremi,precisione,params,fun,&start);
	if (num_zeri==N) //Se N=1 controllo se ho già concluso:
		return(contatore);
	newton_temp=newton(estremi[1],precisione,estremi,&contatore,TEMPO_MAX,params,fun,&start);
	if ((newton_temp!=-999999) && (newton_temp!=-999998))
		zero_append(newton_temp,zeros,&num_zeri,estremi,precisione,params,fun,&start);
	// Se N=2 controllo se ho già concluso usando gli estremi:
	if (num_zeri==N)
		return(contatore);

	//----------------------------------
	// Continuo finché non ne trovo N:
	//----------------------------------
	double dx=estremi[1]-estremi[0];
	while (num_zeri<N)
	{
		dx=dx/2;
		for(double x=estremi[0]+dx;x<estremi[1];x+=(2*dx))
		{
			newton_temp=newton(x,precisione,estremi,&contatore,TEMPO_MAX,params,fun,&start);
			if ((newton_temp!=-999999) && (newton_temp!=-999998)) 
			{
					zero_append(newton_temp,zeros,&num_zeri,estremi,precisione,params,fun,&start);
				if (num_zeri==N)
					break;
			} 
			//Controllo tempo:
			if ( ((float)(clock()-start)/CLOCKS_PER_SEC)>TEMPO_MAX )
				return(-num_zeri);
		}
	}
	return(contatore);
}
