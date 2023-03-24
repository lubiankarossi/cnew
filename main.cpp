
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>

#include <map>
#include <vector>
#include <stdlib.h>

#include "structures.h"
#include "storage.h"
#include "functions.h"


int main() {

	struct geral_data data_c;

	std::string file_data ="in.dat";
	read_data(file_data, data_c);

//	print_data(data_c);

	std::string file_data_r ="in_rep.dat";
	replace_data(file_data_r, data_c);

//	print_data(data_c);


	double dt_  = 1.0 * (60*60*24);

// ====================================================================
// fluxos e fatores de auto-blindagem

	double g_1 = 0.44414E+00,
	       g_2 = 1.0,
	       g_3 = 1.0,
	       g_4 = 1.0;

	double p_14 = 0.23640E01, //termico(1)rapido(4)
	       p_24 = 0.32283E01, //epitermico(2)rapido(4)
	       p_34 = 0.24089E01; //epitemico(3)rapido(4)

	double p4 =   0.5E14;  //rapido

// ====================================================================
// passos de tempo

	std::vector<double> ts; // passos de tempo

// for (int k = 1; k!=121; ++k)
//         ts.push_back(1);

        ts.push_back(1);
        ts.push_back(1);
        ts.push_back(3);
        ts.push_back(5);
       ts.push_back(10);
       ts.push_back(10);
       ts.push_back(10);
       ts.push_back(10);
       ts.push_back(10);
       ts.push_back(30);
       ts.push_back(30);


// ====================================================================
// tabela de nuclídeos no combustível
	storage nuclids;

	nuclids.AddYield(922340, 5.84302E-06);
	nuclids.AddYield(922350, 6.89307E-04);
	nuclids.AddYield(922380, 2.22817E-02);
	nuclids.AddYield( 80160, 4.59538E-02);
      
// ====================================================================
// parametros de sensibilidade


	int id_ref = 942390; // alvo

	std::vector<int>    zl_s; // relação dos nuclídeos

// dados nucleares

	zl_s.push_back(922380);
	zl_s.push_back(922350);
	zl_s.push_back(922360);
	zl_s.push_back(922370);
	zl_s.push_back(922390);
	zl_s.push_back(932370);
	zl_s.push_back(932380);
	zl_s.push_back(932390);
	zl_s.push_back(942380);
	zl_s.push_back(942390);
	zl_s.push_back(942400);
	zl_s.push_back(942410);

	storage adjunts;
	adjunts.AddYield(id_ref, 1.0);

// ====================================================================
// lista para gerção das tabelas
//Calcula as densidades

	map<int, std::vector<double> > result;
	map<int, std::vector<double> > result_ad;

	std::vector<int>    zl; // lista de Zaids a serem guardados


	// add na lista

	zl.push_back(922350);
	zl.push_back(922360);
	zl.push_back(922370);
	zl.push_back(922380);
	zl.push_back(922390);
	zl.push_back(932370);
	zl.push_back(932380);
	zl.push_back(932390);
	zl.push_back(942380);
	zl.push_back(942390);
	zl.push_back(942400);
	zl.push_back(942410);

	for (std::vector<int>::iterator it = (zl).begin() ; it != (zl).end(); ++it){

		std::vector<double> vec;
		vec.push_back( nuclids.GetYield(*it) );		

		std::vector<double> vec_ad;
		vec_ad.push_back( adjunts.GetYield(*it) );		

		result.insert(map<int,std::vector<double> >::value_type( (*it) , vec));
		result_ad.insert(map<int,std::vector<double> >::value_type( (*it) , vec_ad));
	}

// ====================================================================

	double fluxo[4];

	fluxo[0] = 1.96868E+14;           //p_14*p4*g_1  termico 
	fluxo[1] = 3.82630E+14;           //p_24*p4*g_2
	fluxo[2] = 5.25830E+14;           //p_34*p4*g_3
	fluxo[3] = 3.62514E+14;           //p4*g_4      rapido

// ====================================================================

	storage nuclids_c;
	nuclids_c = nuclids; // cópia dos iniciais

	storage nuclids_t;

	storage adjunts_t;

// ====================================================================

// cálculo de t_f

	double  t_f = 0 ;
	for (std::vector<double>::iterator it = (ts).begin() ; it != (ts).end(); ++it){
		t_f += (*it);
	}
	t_f *= dt_;

// ====================================================================
// cálculo das cadeis 
// ====================================================================

	double time_ = 0; // acc de tempo

	for(int ps = 0 ; ps!=(ts).size() ; ++ps){

   		time_    += ts[ps]*dt_;
                double dt = ts[ps]*dt_;

		// calculo das cadeias principais
		for (std::vector<chain_t>::iterator it = (data_c.cp).begin() ; it != (data_c.cp).end(); ++it){

			int n = (*it).size(); // numero de elementos na cadeia

			param_t pr;
			init_parameters(pr,n);
			calc_param(pr, fluxo, it, n, dt);

			calc_chains_cp(nuclids, nuclids_t, pr, n, dt);
			calc_chains_cp_ad(adjunts, adjunts_t, pr, n, time_, t_f);

			clean_parameters(pr);
		}


		//std::cout << "\n\n\t\t Instante[seg] : " << time_ <<"\n" ;
		//std::cout <<     "\t\t ----------------------- "  <<"\n\n" ;

		//std::cout << nuclids_t ;
		//std::cout << "\n" ;
		//std::cout << adjunts_t ;
		//std::cout << "\n" ;
   
		nuclids=nuclids_t;

		for (std::vector<int>::iterator it = (zl).begin() ; it != (zl).end(); ++it){
			result[(*it)].push_back( nuclids_t.GetYield(*it) ); 
			result_ad[(*it)].push_back( adjunts_t.GetYield(*it) ); 
		}

		nuclids_t.reset();
		adjunts_t.reset();
	}

// ====================================================================

	std::cout << "####################################### \n";
	std::cout << "#       Concentrações Diretas         # \n";
	std::cout << "####################################### \n";

	for (std::vector<int>::iterator it = (zl).begin() ; it != (zl).end(); ++it){

		std::cout << "\n\t\t"   << "Zaid "     << "\t" << (*it)                 << "\n";
		std::cout <<   "\t\t"   << "----------------"  << "\n\n";
		std::cout << "\t\t"   << "tempo[dias]" << "\t" << "c[at./barns/cm]"  << "\n\n";
		std::cout << "\t\t\t" << double(0.0)   << "\t" << (result[(*it)])[0] << "\n"  ;

		double tem = 0;		
		for(int ps = 0 ; ps!=(ts).size() ; ++ps){
			tem += ts[ps];
			std::cout << "\t\t\t" << tem << "\t" << (result[(*it)])[ps+1] << "\n";
		}
		std::cout << "\n"; 
	}

	std::cout << "####################################### \n";
	std::cout << "#       Concentrações Adjunta         # \n";
	std::cout << "####################################### \n";


	for (std::vector<int>::iterator it = (zl).begin() ; it != (zl).end(); ++it){

		std::cout << "\n\t\t"   << "Zaid "     << "\t" << (*it)                 << "\n";
		std::cout <<   "\t\t"   << "----------------"  << "\n\n";
		std::cout << "\t\t"   << "tempo[dias]" << "\t" << "c[adjunctors]"       << "\n\n";

		double tem = 0;		
		for(int ps = 0 ; ps!=(ts).size() ; ++ps){
			tem += ts[ps];
			std::cout << "\t\t\t" << tem << "\t" << (result_ad[(*it)])[ps+1] << "\n";
		}
		std::cout << "\n"; 
	}

// ====================================================================
// cálculo dos coeficientes  
// ====================================================================
 
	std::cout << "####################################### \n";
	std::cout << "#    Coeficientes de Sensibilidade    # \n";
	std::cout << "####################################### \n";

	std::cout << "\n\t\t" << "Referência : " <<  id_ref << "\n";
	std::cout <<   "\t\t" << "--------------------" << "\n\n";


	for (std::vector<int>::iterator itr = (zl_s).begin() ; itr != (zl_s).end(); ++itr){

		std::cout << "\n\t\t" << "nuclideo : " <<  *itr << "\n";
		std::cout <<   "\t\t" << "-------------------------" << "\n\n";

		double s_abs_m[4],s_cap_m[4],s_fis_m[4], count;
		double s_abs_max[4],s_cap_max[4],s_fis_max[4];
		count = 0;

		for(int g=0; g!=4; ++g){
			s_abs_m[g]=0;
			s_cap_m[g]=0;
			s_fis_m[g]=0;
			s_abs_max[g]=0;
			s_cap_max[g]=0;
			s_fis_max[g]=0;
		}   


		int chain_number = 0;

		for (std::vector<chain_t>::iterator it = (data_c.cp).begin() ; it != (data_c.cp).end(); ++it){
		
			++chain_number;

			int n = (*it).size(); // numero de elementos na cadeia

			param_t pr;
			init_parameters(pr,n);
			calc_param(pr, fluxo, it, n, t_f);


			int id_se = *itr;
	
			int i_par  = find_el( id_se, pr, n );
			int i_prod = find_el( id_ref, pr, n );
	
			double int_1 = 0,
			       int_2 = 0;

			if ((i_par!=-1)&&(i_prod!=-1)&&(pr.sum_flag[i_prod]!=0)){
				int_1 =  int_chains_cp_ad( nuclids, pr, i_prod, i_par, t_f);
				int_2 = int_chains_cp_ad2( nuclids, pr, i_prod, i_par, t_f);

				count = count + 1.0;

				// abs

				double s_abs[4],s_cap[4],s_fis[4];

				double sa = int_1;
				double sc = int_1;
				double sf = int_1;

				if ( ((*it)[i_par]->reac_id == 3)||((*it)[i_par]->reac_id == 4) ) {
					sa += int_2 * (*it)[i_par]->alpha;
					sc += int_2 * (*it)[i_par]->alpha;
				}

				if ( ((*it)[i_par]->reac_id == 4) ) {
					sf -= int_2 * (*it)[i_par]->alpha;
				}
				else if ( ((*it)[i_par]->reac_id == 3) ) {
					sf += int_2 * (*it)[i_par]->alpha;
				}

				for(int g=0; g!=4; ++g){

					s_abs[g]=1.E-24*((*it)[i_par]->sig_a[3-g])*sa*fluxo[g]/nuclids.GetYield(id_se);
					s_cap[g]=1.E-24*((*it)[i_par]->sig_g[3-g])*sc*fluxo[g]/nuclids.GetYield(id_se);
					s_fis[g]=1.E-24*((*it)[i_par]->sig_f[3-g])*sf*fluxo[g]/nuclids.GetYield(id_se);				

					s_abs_m[g]+=s_abs[g];
					s_cap_m[g]+=s_cap[g];
					s_fis_m[g]+=s_fis[g];

					s_abs_max[g]=(abs(s_abs_max[g])>abs(s_abs[g]))?s_abs_max[g]:s_abs[g];
					s_cap_max[g]=(abs(s_cap_max[g])>abs(s_cap[g]))?s_cap_max[g]:s_cap[g];
					s_fis_max[g]=(abs(s_fis_max[g])>abs(s_fis[g]))?s_fis_max[g]:s_fis[g];

				}   

				std::cout << "\n\t\t" << "Número da cadeia : " <<  chain_number << "\n\n";
				std::cout << "\t\t" << "group" << "\t" << "S[abs]" <<"\t" << "S[cap]" <<"\t" << "S[fis]" << "\n";
				for(int g=0; g!=4; ++g){
					std::cout << "\t\t" << g+1 
					          << "\t" << s_abs[g]
					          << "\t" << s_cap[g]
					          << "\t" << s_fis[g] 
					          << std::scientific <<"\n";
				}   
				std::cout << "\n";
  
			}

			clean_parameters(pr);

		}
		std::cout << "\n\t\t" << "Média nas cadeias " <<  "\n\n";
		std::cout << "\t\t" << "group" << "\t" << "S[abs]" <<"\t" << "S[cap]" <<"\t" << "S[fis]" << "\n";
		for(int g=0; g!=4; ++g){
			std::cout << "\t\t" << g+1 
			          << "\t" << s_abs_m[g]/count
			          << "\t" << s_cap_m[g]/count
			          << "\t" << s_fis_m[g]/count
			          << std::scientific <<"\n";
		}   
		std::cout << "\n";
		std::cout << "\n\t\t" << "Maximo nas cadeias (conservativo)" <<  "\n\n";
		std::cout << "\t\t" << "group" << "\t" << "S[abs]" <<"\t" << "S[cap]" <<"\t" << "S[fis]" << "\n";
		for(int g=0; g!=4; ++g){
			std::cout << "\t\t" << g+1 
			          << "\t" << s_abs_max[g]
			          << "\t" << s_cap_max[g]
			          << "\t" << s_fis_max[g]
			          << std::scientific <<"\n";
		}   
		std::cout << "\n";

  

	}

// ====================================================================
// ====================================================================
// cálculo das derivativas  
// ====================================================================
 
	std::cout << "####################################### \n";
	std::cout << "#            Derivativas              # \n";
	std::cout << "####################################### \n";

	std::cout << "\n\t\t" << "Referência : " <<  id_ref << "\n";
	std::cout <<   "\t\t" << "--------------------" << "\n\n";


	for (std::vector<int>::iterator itr = (zl_s).begin() ; itr != (zl_s).end(); ++itr){

		std::cout << "\n\t\t" << "nuclideo : " <<  *itr << "\n";
		std::cout <<   "\t\t" << "-------------------------" << "\n\n";

		double d_abs_m[4],d_cap_m[4],d_fis_m[4], count;
		double d_abs_max[4],d_cap_max[4],d_fis_max[4];
		count = 0;

		for(int g=0; g!=4; ++g){
			d_abs_m[g]=0;
			d_cap_m[g]=0;
			d_fis_m[g]=0;
			d_abs_max[g]=0;
			d_cap_max[g]=0;
			d_fis_max[g]=0;
		}   


		int chain_number = 0;

		for (std::vector<chain_t>::iterator it = (data_c.cp).begin() ; it != (data_c.cp).end(); ++it){
		
			++chain_number;

			int n = (*it).size(); // numero de elementos na cadeia

			param_t pr;
			init_parameters(pr,n);
			calc_param(pr, fluxo, it, n, t_f);


			int id_se = *itr;
	
			int i_par  = find_el( id_se, pr, n );
			int i_prod = find_el( id_ref, pr, n );
	
			double int_1 = 0,
			       int_2 = 0;

			if ((i_par!=-1)&&(i_prod!=-1)&&(pr.sum_flag[i_prod]!=0)){
				int_1 =  int_chains_cp_ad( nuclids, pr, i_prod, i_par, t_f);
				int_2 = int_chains_cp_ad2( nuclids, pr, i_prod, i_par, t_f);

				count = count + 1.0;

				// abs

				double d_abs[4],d_cap[4],d_fis[4];

				double sa = int_1;
				double sc = int_1;
				double sf = int_1;

				if ( ((*it)[i_par]->reac_id == 3)||((*it)[i_par]->reac_id == 4) ) {
					sa += int_2 * (*it)[i_par]->alpha;
					sc += int_2 * (*it)[i_par]->alpha;
				}

				if ( ((*it)[i_par]->reac_id == 4) ) {
					sf -= int_2 * (*it)[i_par]->alpha;
				}
				else if ( ((*it)[i_par]->reac_id == 3) ) {
					sf += int_2 * (*it)[i_par]->alpha;
				}

				for(int g=0; g!=4; ++g){

					d_abs[g]=1.E-24*sa*fluxo[g];
					d_cap[g]=1.E-24*sc*fluxo[g];
					d_fis[g]=1.E-24*sf*fluxo[g];				

					d_abs_m[g]+=d_abs[g];
					d_cap_m[g]+=d_cap[g];
					d_fis_m[g]+=d_fis[g];

					d_abs_max[g]=(abs(d_abs_max[g])>abs(d_abs[g]))?d_abs_max[g]:d_abs[g];
					d_cap_max[g]=(abs(d_cap_max[g])>abs(d_cap[g]))?d_cap_max[g]:d_cap[g];
					d_fis_max[g]=(abs(d_fis_max[g])>abs(d_fis[g]))?d_fis_max[g]:d_fis[g];

				}   

				std::cout << "\n\t\t" << "Número da cadeia : " <<  chain_number << "\n\n";
				std::cout << "\t\t" << "group" << "\t" << "dN/ds[abs]" <<"\t" << "dN/ds[cap]" <<"\t" << "dN/ds[fis]" << "\n";
				for(int g=0; g!=4; ++g){
					std::cout << "\t\t" << g+1 
					          << "\t" << d_abs[g]
					          << "\t" << d_cap[g]
					          << "\t" << d_fis[g] 
					          << std::scientific <<"\n";
				}   
				std::cout << "\n";
  
			}

			clean_parameters(pr);

		}
		std::cout << "\n\t\t" << "Média nas cadeias " <<  "\n\n";
		std::cout << "\t\t" << "group" << "\t" << "dN/ds[abs]" <<"\t" << "dN/ds[cap]" <<"\t" << "dN/ds[fis]" << "\n";
		for(int g=0; g!=4; ++g){
			std::cout << "\t\t" << g+1 
			          << "\t" << d_abs_m[g]/count
			          << "\t" << d_cap_m[g]/count
			          << "\t" << d_fis_m[g]/count
			          << std::scientific <<"\n";
		}   
		std::cout << "\n";
		std::cout << "\n\t\t" << "Maximo nas cadeias (conservativo)" <<  "\n\n";
		std::cout << "\t\t" << "group" << "\t" << "dN/ds[abs]" <<"\t" << "dN/ds[cap]" <<"\t" << "dN/ds[fis]" << "\n";
		for(int g=0; g!=4; ++g){
			std::cout << "\t\t" << g+1 
			          << "\t" << d_abs_max[g]
			          << "\t" << d_cap_max[g]
			          << "\t" << d_fis_max[g]
			          << std::scientific <<"\n";
		}   
		std::cout << "\n";

  

	}

// ====================================================================


	clear_data(data_c); // limpa os dados 

	return 0;
}

