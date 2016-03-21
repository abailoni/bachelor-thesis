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
