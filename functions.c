#include "functions.h"

double exp_c(double x){

	if (-x>675)
		return 0;
	else if (-x<CN)	
		return 1.0 + x;

        return exp(x); 	

}

void init_parameters(param_t &p, size_t n){

  p.alfa   = (double *)malloc(n * sizeof(double));
  p.gama   = (double *)malloc(n * sizeof(double));
  p.Af     = (double *)malloc(n * sizeof(double));
  p.Aa     = (double *)malloc(n * sizeof(double));
  p.Ac     = (double *)malloc(n * sizeof(double));
  p.lambda = (double *)malloc(n * sizeof(double));

  p.beta   = (double *)malloc(n * sizeof(double));
  p.betaa  = (double *)malloc(n * sizeof(double));

  p.zaid   = (unsigned *)malloc(n * sizeof(unsigned));

  p.sum_flag   = (bool *)malloc(n * sizeof(bool));

}

void clean_parameters(param_t &p){

  free(p.alfa);     p.alfa     = NULL;
  free(p.gama);     p.gama     = NULL;
  free(p.Af);       p.Af       = NULL;
  free(p.Aa);       p.Aa       = NULL;
  free(p.Ac);       p.Ac       = NULL;
  free(p.lambda);   p.lambda   = NULL;
  free(p.zaid);     p.zaid     = NULL;
  free(p.sum_flag); p.sum_flag = NULL;
  free(p.beta);     p.beta     = NULL;
  free(p.betaa);    p.betaa    = NULL;

}

void calc_param(param_t &pr, double fluxo[], std::vector<chain_t>::iterator it, size_t n, double dt){

   for (int m=0; m != n ; ++m) {
  
      pr.Af[m]=pr.Ac[m]=pr.Aa[m]=pr.alfa[m]=pr.gama[m]=pr.lambda[m]=0.0; // inicia vetores para cadeia
  
      for(int g=0; g!=4; ++g){
         pr.Af[m] += (*it)[m]->sig_f[3-g] * fluxo[g] * 1.0E-24; // fissão
         pr.Ac[m] += (*it)[m]->sig_g[3-g] * fluxo[g] * 1.0E-24; // captura 
         pr.Aa[m] += (*it)[m]->sig_a[3-g] * fluxo[g] * 1.0E-24; // absorção 
      }   
  
      if ((*it)[m]->reac_id == 1){      //no_coupling
         pr.gama[m] = 0;
      }
      else if ((*it)[m]->reac_id == 2){ //decay
         pr.gama[m] = (*it)[m]->lamb;
      }
      else if ((*it)[m]->reac_id == 3){ //absorption
         pr.gama[m] = pr.Aa[m];
      }
      else if ((*it)[m]->reac_id == 4){ //capture
         pr.gama[m] = pr.Aa[m]-pr.Af[m];
      }
  
      pr.alfa[m]   = (*it)[m]->alpha;
      pr.lambda[m] = (*it)[m]->lamb; 
  
      pr.zaid[m]     = (*it)[m]->zaid;
      pr.sum_flag[m] = (*it)[m]->sum_flag;
   }  

//--   for
	for (int i=0; i!=n; ++i) {
		double cmi = 2.0 * (i+1) * CM;
        	pr.beta[i] = pr.lambda[i]+pr.Aa[i];
        
        	double trr  = pr.beta[i];
        	double trrt = trr;

		if (trr*dt <= cmi)
			trr=cmi/dt;
 
        	bool flag_i = false;

		if(i!=0){
			for(int j=0; j!=(i);++j){
				double cmm = cmi;
				if(cmm<cmi/dt)
					cmm=cmi/dt;
				if(abs(trr-pr.beta[j])<cmm){
//				}
//				else{
					trr=trr+trr*cmm;
					flag_i = true;
				}
			}	

		}

        	pr.beta[i]=trr;

  		if((trr-trrt)*dt>=CM)
			flag_i = true;

		if (flag_i == true){
		//	std::cerr<<"aviso :  !!!";
		}
	}
//--
//--   for
	for (int i=n-1; i!=-1; --i) {
		double cmi = 2.0 * (i+1) * CM;
        	pr.betaa[i] = pr.lambda[i]+pr.Aa[i];
        
        	double trr  = pr.betaa[i];
        	double trrt = trr;

		if (trr*dt <= cmi)
			trr=cmi/dt;
 
        	bool flag_i = false;

		if(i!=n-1){
			for(int j=i-1; j!=-1;--j){
				double cmm = cmi;
				if(cmm<cmi/dt)
					cmm=cmi/dt;
				if(abs(trr-pr.betaa[j])>=cmm){
				}
				else{
					trr=trr+trr*cmm;
					flag_i = true;
				}
			}	

		}

        	pr.betaa[i]=trr;

  		if((trr -trrt)*dt>=CM)
			flag_i = true;

		if (flag_i == true){
		//	std::cerr<<"aviso :  !!!";
		}
	}
//--

}

void calc_chains_cp(storage &nuclids, storage &nuclids_t, param_t &pr, size_t n, double dt){

	for (int ni=0; (ni!=n) ; ++ni) {
          
		double N_ind = 0.0;
	 
		for (int m=0; (m<=ni) ; ++m) {
	
			double prod_1 = 1;
	
			// primeiro produtório alfas e gamas
			if (ni!=0)
			for (int k=m; k<=(ni-1);++k){ // (n-1) na soma e o vetor vai de 0 a (n-1)
				prod_1 = prod_1*pr.gama[k]*pr.alfa[k];
			}
	
			double sum_2 = 0;
	
			// soma sem fissão
			double yield = nuclids.GetYield( pr.zaid[m] ); 

			bool vlt = false;
			double te = 0.0;
			
			if(yield > CUTOFF){
	    
				for(int j = m; j<=(ni); ++j){
	                        
					double prod_int = 1;
	                      		if (ni!=0)
					for(int i = m; i<=(ni); ++i){
						if (i!=j){
							double bibj = pr.beta[i]-pr.beta[j];
							prod_int = prod_int*(bibj);

                                                        if ( abs(pr.beta[i]*dt) < 1.0E-6 || abs(bibj*dt) < 1.0E-6)
                                                                vlt = true;
						}
					}
					if ( vlt != true)
						sum_2 += exp( - pr.beta[j] * dt )/prod_int;
				}
			}
	    
			N_ind += yield*sum_2*prod_1; // adiciona no somatório principal
			
		}      
	
		if ((pr.sum_flag[ni] == 1)&&(N_ind>CUTOFF)&&(N_ind>0)){ // se esta na soma e o montante é numericamente expressivo add ao yield
			if (ni!=0)nuclids_t.SumYield(pr.zaid[ni],N_ind, abs(pr.gama[ni-1]*pr.alfa[ni-1]*nuclids.GetYield( pr.zaid[ni-1] )-pr.beta[ni]*nuclids.GetYield( pr.zaid[ni] )));
			else nuclids_t.SumYield(pr.zaid[ni],N_ind,abs(-pr.beta[ni]*nuclids.GetYield( pr.zaid[ni] )));
		}
		else if ((pr.sum_flag[ni] == 1)&&(N_ind>CUTOFF)&&(N_ind<0)){ // se a concentração é negativa emite um aviso 
			nuclids_t.SumYield(pr.zaid[ni],N_ind);
			std::cerr<<"aviso : negative concentration calculated !!!";
		}
	} 
}

void calc_chains_cp_ad(storage &adjunts, storage &adjunts_t, param_t &pr, size_t n, double dt, double t_f){

	for (int ni=0; (ni!=n) ; ++ni) {
   
		double N_ind = 0.0;
 
		int i_n = (n-1) - ni;

		for (int k=(n-1); (k>=i_n) ; --k) {


			double prod_1 = 1;

			// primeiro produtório alfas e gamas

			if (i_n!=n-1)
			for (int m=k-1; m>=i_n;--m){ // (n-1) na soma e o vetor vai de 0 a (n-1)
				prod_1 = prod_1*pr.gama[m]*pr.alfa[m];
			}

			double sum_2 = 0;

			// soma sem fissão
			double yield = adjunts.GetYield( pr.zaid[k] ); 
    
			if(yield > CUTOFF){
    
				for(int l = k; l>=(i_n); --l){
                        
					double prod_int = 1;
                       
					if (i_n!=n-1)
					for(int j = k; j>=(i_n); --j){
						if (j!=l){
                       
							double bjbl = pr.betaa[j]-pr.betaa[l];
							prod_int = prod_int*(bjbl);

                                                        if ( pr.beta[j]*(t_f-dt) <  1.0E-9 && (t_f!=dt) )
                                                                prod_int=0.0;
                           			}
                       			}
					sum_2 += exp( -pr.betaa[l] * (t_f-dt) )/prod_int; 
                    		}
                	}

   
                	N_ind += yield*sum_2*prod_1; // adiciona no somatório principal
            	}
            
            	if ((pr.sum_flag[i_n] == 1)&&(N_ind>CUTOFF)&&(N_ind>0)){ // se esta na soma e o montante é numericamente expressivo add ao yield
                	adjunts_t.SumYield(pr.zaid[i_n],N_ind);
		}
            	else if ((pr.sum_flag[i_n] == 1)&&(N_ind>CUTOFF)&&(N_ind<0)){ // se a concentração é negativa emite um aviso 
                	adjunts_t.SumYield(pr.zaid[i_n],N_ind);
                	std::cerr<<"aviso : negative concentration calculated !!!";
            	}

       	}

}


double int_chains_cp_ad(storage &nuclids, param_t &pr, int N, int p, double t_f){

	double prod_1 = 1;

	for(int n = N-1; n>=(p); --n){
		if (p!=N-1)
			prod_1 = prod_1 * pr.gama[n] * pr.alfa[n]; 
	} 

	double soma_3 = 0;
	for(int m = 0; m <= (p); ++m){

		double prod_2 = 1;
		for(int k = m; k<=p-1; ++k){
			if (m!=p)
				prod_2 = prod_2 * pr.gama[k] * pr.alfa[k];
		}

		double yield = nuclids.GetYield( pr.zaid[m] ); 

		double soma_2=0;
		for(int h = m; h <= (p); ++h){

			double prod_3 = 1;
			for(int i = m; i<=p ; ++i){
				if (i!=h)
					prod_3 = prod_3 * (pr.beta[i] - pr.beta[h]); 
			}

                        double soma_1=0;
			for(int l = N; l >= (p); --l){

				double prod_4 = 1;
				for(int j = N; j>=p ; --j){
					if (j!=l)
						prod_4 = prod_4 * (pr.betaa[j] - pr.betaa[l]); 
				}
				if (h!=l)
 					soma_1 += ( exp(-pr.betaa[l]*t_f) - exp(-pr.beta[h]*t_f) ) / ( pr.beta[h] - pr.betaa[l] )/prod_4;
				else
 					soma_1 += ( exp(-pr.betaa[l]*t_f) * t_f ) / prod_4;

			}

			soma_2 = soma_2 + soma_1/prod_3; 	
		}

		soma_3 = soma_3 + yield * prod_2 * soma_2; 

	}

	return soma_3*prod_1; 
}


double int_chains_cp_ad2(storage &nuclids, param_t &pr, int N, int p, double t_f){

	double prod_1 = 1;

	for(int n = N-1; n>=(p+1); --n){
		if (p!=N-1)
			prod_1 = prod_1 * pr.gama[n] * pr.alfa[n]; 
	} 

	double soma_3 = 0;
	for(int m = 0; m <= (p); ++m){

		double prod_2 = 1;
		for(int k = m; k<=p-1; ++k){
			if (m!=p)
				prod_2 = prod_2 * pr.gama[k] * pr.alfa[k];
		}

		double yield = nuclids.GetYield( pr.zaid[m] ); 

		double soma_2=0;
		for(int h = m; h <= (p); ++h){

			double prod_3 = 1;
			for(int i = m; i<=p ; ++i){
				if (i!=h)
					prod_3 = prod_3 * (pr.beta[i] - pr.beta[h]); 
			}

                        double soma_1=0;
			for(int l = N; l >= (p+1); --l){

				double prod_4 = 1;
				for(int j = N; j>=(p+1) ; --j){
					if (j!=l)
						prod_4 = prod_4 * (pr.betaa[j] - pr.betaa[l]); 
				}
				if (h!=l)
 				soma_1 += ( exp(-pr.betaa[l]*t_f) - exp(-pr.beta[h]*t_f) ) / ( pr.beta[h] - pr.betaa[l] )/prod_4;

			}

			soma_2 = soma_2 + soma_1/prod_3; 	
		}

		soma_3 = soma_3 + yield * prod_2 * soma_2; 

	}

	return soma_3*prod_1; 
}


int  find_el(int zaid, param_t &pr, size_t n ){

	int pos = -1;

	for(int i=0; i!=n; ++i){

		if (zaid == pr.zaid[i])
			pos = i;
	}

	return pos;

}

//  manipulação de string

int get_substr_i(std::string str_){
	return std::stoi(str_.substr(0,6));
}


















