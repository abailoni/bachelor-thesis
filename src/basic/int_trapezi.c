
/////////////////////////////////////////
// INTEGRALE FINITO CON METODO TRAPEZI //
/////////////////////////////////////////

/*
INPUT int_trapezi():
 - array con estremi
 - precisione relativa per modalità 1 // passo dx per modalità 2
 - eventuali parametri funzione [puntatore void]
 - funzione 
 - modalità: 1 o 2 [int]
*/

double int_trapezi(double *estremi, double precisione, void *params, double (*fun) (double,void*), int mod) 
{ 
	double x_min=estremi[0], x_max=estremi[1];
	
	//CASO IN CUI SI VUOLE UNA CERTA PRECISIONE:
	if (mod==1) { 
		// Calcolo f(x_max) e f(x_min) in modo da dare la prima stima dell'integrale semplicemente come area di un trapezio di altezza x_max-x_min:
		double sum=(x_max-x_min)*(fun(x_max,params)+fun(x_min,params))/2;

		//Continuo finché non raggiungo la precisione:
		double new_sum, total_sum, old_sum=sum, h=(x_max-x_min)/2;
		for(int contatore=0;1;contatore++) {
			if (contatore>1e3)
				return(-999999);
			// Calcolo i nuovi termini della sommatoria:
			new_sum=0;
			for(double x=estremi[0]+h;x<estremi[1];x+=(2*h))
				new_sum+=fun(x,params);
			new_sum=new_sum*h;
			// Sommo i nuovi contributi a quelli vecchi (dato che h è dimezzato quelli vecchi dovranno essere divisi per 2):
			total_sum=new_sum+old_sum/2;
			h=h/2;
			if ((fabs(total_sum-old_sum)/total_sum)<=precisione)
				return(total_sum);
			old_sum=total_sum;
		}
	//CASO IN CUI h E' BEN DEFINITO:
	} else { 
		double sum=(fun(x_max,params)+fun(x_min,params))/2;
		double h=precisione;
		for(double x=estremi[0]+h;x<estremi[1];x+=(h)) {
			sum+=fun(x,params);
		}
		return(sum*=h);
	}

}

