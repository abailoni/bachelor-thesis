
///////////////////////////////////////////////////
// TROVA ZERO FUNZIONE CON METODO DELLA BISEZIONE:
///////////////////////////////////////////////////
//
// INPUT bisezione():
// - due estremi da cui partire
// - intervallo minimo di precisione che si vuole raggiungere
// - eventuali parametri
// - funzione

double bisezione_ric(double estremi[2], double valori_fun[2], double precisione, void *params, double (*fun)(double,void*))
{
	double dx=estremi[1]-estremi[0];
	double new_value=fun(estremi[0]+dx/2,params);
	
	///// SELEZIONO IL NUOVO INTERVALLO:
	// prima maggiore di zero:
	int k=0; 
	int z=1;
	if (valori_fun[0]<0) { //dopo maggiore di zero
		k=1;
		z=0;
	}
	if (new_value>0) { //new_value maggiore di zero
		estremi[k]=estremi[0]+dx/2;
		valori_fun[k]=new_value;
	} else { //new_value minore di zero
		estremi[z]=estremi[0]+dx/2;
		valori_fun[z]=new_value;
	}
	
	///// CONTROLLO PRECISIONE RAGGIUNTA:
	if ((dx/2)<precisione)
		return(estremi[0]+dx/4);
	else
		return(bisezione_ric(estremi,valori_fun,precisione,params,fun));
}

double bisezione(double estremi[2], double precisione, void *params, double (*fun)(double,void*))
{
	double valori_fun[]={fun(estremi[0],params),fun(estremi[1],params)};
	return(bisezione_ric(estremi,valori_fun,precisione,params,fun));
}

///////////////////////////////////////////////////////////////////
