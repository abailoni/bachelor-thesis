/////////////////////////////////////////////
// Richiamare la funzione fun_max() con INPUT:
//   * estremi in cui cercare massimo
//   * precisione con cui trovarlo
//   * eventuali parametri da passare a funzione
//   * funzione di cui calcolare il massimo
/////////////////////////////////////////////

#define sez_aurea 0.38197

void fun_max(double *estremi, double precisione, void *params, double (*fun)(double,void*))
{
	double dx=estremi[1]-estremi[0];
	int memory;
	double new_val_fun[]={fun(estremi[0]+dx*sez_aurea,params),fun(estremi[1]-dx*sez_aurea,params)};
	if (new_val_fun[0]<new_val_fun[1]) { //tengo l'intervallo a dx
		memory=1;
		estremi[0]=estremi[0]+dx*sez_aurea;
	} else { //tengo l'intervallo a sx
		memory=0;
		estremi[1]=estremi[1]-dx*sez_aurea;
	}
	//VADO AVANTI FINCHE' NON SUPERO PRECISIONE:
	dx=estremi[1]-estremi[0];
	double old_val_fun[2],new_points[2];
	while (fabs(dx)>=precisione)
	{
		old_val_fun[0]=new_val_fun[0],old_val_fun[1]=new_val_fun[1];
		new_points[0]=estremi[0]+dx*sez_aurea,new_points[1]=estremi[1]-dx*sez_aurea;
		//Controllo cosa ho fatto il passo precedente:
		if (memory==0){
			new_val_fun[0]=fun(new_points[0], params);
			new_val_fun[1]=old_val_fun[0];
		} else {
			new_val_fun[0]=old_val_fun[1];
			new_val_fun[1]=fun(new_points[1], params);
		}
		//Controllo quale intervallo tenere:
		if (new_val_fun[0]<new_val_fun[1]) {
			memory=1;
			estremi[0]=new_points[0];
		} else {
			memory=0;
			estremi[1]=new_points[1];
		}
		dx=estremi[1]-estremi[0];
	}
}

