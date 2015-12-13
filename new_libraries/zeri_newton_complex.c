
//////////////////////////////////////////////
// TROVA ZERO FUNZIONE CON METODO DI NEWTON:
//////////////////////////////////////////////
// Richiede libreria complex.h


// INPUT newton():
// - x_0 di partenza (complesso)
// - precisione
// - eventuali parametri [void]
// - funzione (deve restituire e richiedere un complesso)

// double newton_ric_complex(complex double x_values[2], complex double fun_value[2], double precisione, void *params, complex double (*fun)(complex double,void*))
// {
// 	// Calcolo il nuovo x:
// 	complex double new_x=x_values[1]-fun_value[1]*(x_values[1]-x_values[0])/(fun_value[1]-fun_value[0]);
// 	// Controllo la precisione:
// 	if (fabs(new_x-x_values[1])<precisione)
// 		return(new_x);
// 	else {
// 		x_values[0]=x_values[1];
// 		x_values[1]=new_x;
// 		fun_value[0]=fun_value[1];
// 		fun_value[1]=fun(new_x,params);
// 		return(newton_ric_complex(x_values,fun_value,precisione,params,fun));
// 	}
// }

complex double newton_complex(complex double x_0, double precisione, void *params, complex double (*fun)(complex double,void*))
{
	// Trovo il primo punto "manualmente":
	double h=1e-10;
	complex double fun_value[2];
	fun_value[0]=fun(x_0,params);
	// printf("%g\n", creal(fun_value[0]));
	// getchar();
	complex double x_1=x_0-fun_value[0]*(2*h/(fun(x_0+h,params)-fun(x_0-h,params)));
	complex double x_values[]={x_0,x_1};
	fun_value[1]=fun(x_1,params);
	// Richiamo la funzione ricorsiva:
	while(1)
	{
		complex double new_x=x_values[1]-fun_value[1]*(x_values[1]-x_values[0])/(fun_value[1]-fun_value[0]);
		if (cabs((new_x-x_values[1])/new_x)<precisione)
			return(new_x);
		x_values[0]=x_values[1];
		x_values[1]=new_x;
		fun_value[0]=fun_value[1];
		fun_value[1]=fun(new_x,params);
	}


	// return(newton_ric_complex(x_values,fun_value,precisione,params,fun));
}
