
#include "structures.h"

#include <stdbool.h>
#include <stdio.h>

static inline bool isspace(char c){
	return (c == ' ' || c == '\t' || c == '\n' || c == '\12');
}

double get_double(unsigned pos, unsigned length, std::string line){

   std::string temp;
   double result=0;

   for(unsigned i=pos; i!=pos+length; ++i)
      if (!isspace(line[i]))
         temp += line[i];

   sscanf(temp.c_str(), "%lf", &result);

   return result;
}

int get_int(unsigned pos, unsigned length, std::string line){

   std::string temp;
   int result=0;

   for(unsigned i=pos; i!=pos+length; ++i)
      if (!isspace(line[i]))
         temp += line[i];

   sscanf(temp.c_str(), "%d", &result);

   return result;
}

std::string get_string(unsigned pos, unsigned length, std::string line){

   std::string temp;
   int result=0;

   for(unsigned i=pos; i!=pos+length; ++i)
      if (!isspace(line[i]))
         temp += line[i];
   return temp;
}

template<typename T>
struct deleter {
    void operator()(T* ptr) {
        if(ptr!=NULL)
           free(ptr);
        ptr=NULL;
    }
};

template<typename T>
struct clear_d {
    void operator()(T ptr) {
        ptr.clear();
    }
};


void clear_data_vector(chain_t &data){

   std::for_each(data.begin(), data.end(), deleter<base_data_t>());
   data.clear();
   return;
}

void clear_data(struct geral_data &data){

   for (std::vector<chain_t>::iterator it = (data.cp).begin() ; it != (data.cp).end(); ++it)
       clear_data_vector(*it);

   data.cp.clear();

   for (std::vector<chain_t>::iterator it = data.ca.begin() ; it != data.ca.end(); ++it)
       clear_data_vector(*it);

   data.ca.clear();
   data.fp.clear();

   return;
}

void read_data(std::string file, geral_data_t &data_c){

    std::ifstream fortranfile;
    fortranfile.open(file.c_str());

    if (!fortranfile.is_open()) {
        exit(1);
    }

    std::string line;

    getline(fortranfile, line); // first line is comment

    unsigned       cp    = 0, 
                   ca    = 0, 
                   fp    = 0, 
                   dummy = 0;

    std::vector<int> nel_cp;
    std::vector<int> nel_ca;

    fortranfile >> cp;
    fortranfile >> dummy;
    fortranfile >> ca;
    fortranfile >> dummy;
    fortranfile >> fp;

    for(int i=0;i!=cp;++i){
       fortranfile >> dummy;
       nel_cp.push_back(dummy);
    }
       for(int i=0;i!=ca;++i){
         fortranfile >> dummy;
         nel_ca.push_back(dummy);
    }
    getline(fortranfile, line); // line break 

    for(unsigned m=0; m!=cp; ++m){

       chain_t local_data; 

       for(int n=0; n!=nel_cp[m]; ++n){

          getline(fortranfile, line);

          struct base_data_t* el = (struct base_data_t*) malloc(sizeof(struct base_data_t));

          if(NULL == el){
             std::cerr << "\n Node creation failed \n";
             exit(1);
          }
        
          bool flag_alpha = false;

          el->alpha = 1.0;

          if (get_int(67, 1, line) == 1)
              flag_alpha = true;

          if (get_int(68, 1, line) == 1)
              el->sum_flag = true;
          else
              el->sum_flag = false;

          for(int i=0; i!=NGPR; ++i)
             el->sig_a[i] = get_double((18 + i*12), 12, line);
          el->zaid = get_int(72, 6, line);

          getline(fortranfile, line);

          el->node_number  = get_int(3, 2, line);
          el->reac_id      = get_int(5, 1, line);    
 
          el->lamb         = get_double(6, 12, line);
          for(int i=0; i!=NGPR; ++i){
             el->sig_f[i] = get_double((18 + i*12), 12, line);
             el->sig_g[i] = el->sig_a[i] - el->sig_f[i];
          }

          if (flag_alpha){
             getline(fortranfile, line);
             el->alpha = get_double(6, 12, line);
          }
          local_data.push_back(el);
       }
     
       data_c.cp.push_back(local_data);
       local_data.clear();

    }

//-----------------------------------------------------

    for(unsigned m=0; m!=ca; ++m){

       chain_t local_data; 

       for(int n=0; n!=nel_ca[m]; ++n){

          getline(fortranfile, line);

          struct base_data_t* el = (struct base_data_t*) malloc(sizeof(struct base_data_t));

          if(NULL == el){
             std::cerr << "\n Node creation failed \n";
             exit(1);
          }
        
          bool flag_alpha = false;

          el->alpha = 1.0;

          if (get_int(67, 1, line) == 1)
              flag_alpha = true;

          for(int i=0; i!=NGPR; ++i)
             el->sig_a[i] = get_double((18 + i*12), 12, line);
          el->zaid = get_int(72, 6, line);

          getline(fortranfile, line);

          el->node_number  = get_int(3, 2, line);
          el->reac_id      = get_int(5, 1, line);    
 
          el->lamb         = get_double(6, 12, line);
          for(int i=0; i!=NGPR; ++i){
             el->sig_f[i] = get_double((18 + i*12), 12, line);
             el->sig_g[i] = el->sig_a[i] - el->sig_f[i];
          }

          if (get_int(68, 1, line) == 1)
              el->fis_flag = true;
          else
              el->fis_flag = false;

          if (flag_alpha){
             getline(fortranfile, line);
             el->alpha = get_double(6, 12, line);
          }
          local_data.push_back(el);
       }
       data_c.ca.push_back(local_data);
       local_data.clear();
    }

    struct base_yield_t byl;

    for (int i=0; i!=(fp/5); ++i){

       getline(fortranfile, line); 
       for (int j=0; j!=5; ++j){
          byl.id = get_string((6+12*j), 12, line);
          data_c.fp.push_back(byl);
       }
    }

    if ((fp%5)!=0){
       getline(fortranfile, line); 
       for (int j=0; j!=(fp%5); ++j){
          byl.id = get_string((6+12*j), 12, line);
          data_c.fp.push_back(byl);
       }
    }

//

    for (unsigned i=0; i!=(fp/5); ++i){

       getline(fortranfile, line); 
       for (int j=0; j!=5; ++j){
          data_c.fp[(5*i+j)].Q = get_double((6+12*j), 12, line);
       }
    }

    if ((fp%5)!=0){
       getline(fortranfile, line); 
       for (int j=0; j!=(fp%5); ++j){
          data_c.fp[(fp+j)].Q = get_double((6+12*j), 12, line);
       }
    }

//

    std::vector< std::vector<double> > temp(fp);

    for (int k=0; k!=ca; ++k){

       for (int n=0; n!=nel_ca[k]; ++n){

          for (unsigned i=0; i!=(fp/5); ++i){
             getline(fortranfile, line); 
             for (int j=0; j!=5; ++j)
                temp[(5*i+j)].push_back(get_double((6+12*j), 12, line));
          }

          if ((fp%5)!=0){
             getline(fortranfile, line); 
             for (int j=0; j!=(fp%5); ++j){
                temp[(fp+j)].push_back(get_double((6+12*j), 12, line));
             }
          }
       }
    }

    int count = 0;
    for (int k=0; k!=fp; ++k){

       for (int j=0; j!=ca; ++j){

          std::vector<double> temp_2; 
          for (int i=0; i!=nel_ca[j]; ++i){
 
             temp_2.push_back(temp[k][count]);
             count++;
          }
          data_c.fp[k].yield.push_back(temp_2);
          temp_2.clear();     
       }
       count =0;
    }
    temp.clear();

    fortranfile.close();

    return;
}

void replace_data(std::string file, geral_data_t &data_c ){

   	double   sig_a[NGPR], /* abs. cross section     */
	         sig_f[NGPR], /* fis. cross section     */
	         sig_g[NGPR]; /* cap. cross section     */

	unsigned zaid;        /* nuclide id             */

	std::ifstream fortranfile;
	fortranfile.open(file.c_str());

	if (!fortranfile.is_open()) {
		exit(1);
	}

	while (fortranfile>>zaid){

		for(int i=0;i!=NGPR;++i)
			fortranfile>>sig_a[i];
		for(int i=0;i!=NGPR;++i)
			fortranfile>>sig_f[i];
		for(int i=0;i!=NGPR;++i)
			fortranfile>>sig_g[i];

		for (std::vector<chain_t>::iterator it = (data_c.cp).begin() ; it != (data_c.cp).end(); ++it){

			int n = (*it).size(); // numero de elementos na cadeia

   			for (int m=0; m != n ; ++m) {
				if ( (*it)[m]->zaid ==  zaid){
					for(int i=0;i!=NGPR;++i){
						(*it)[m]->sig_a[i] = sig_a[i]; 
						(*it)[m]->sig_f[i] = sig_f[i]; 
						(*it)[m]->sig_g[i] = sig_g[i]; 
					}
				}
			}
		}



	}

	fortranfile.close();
	return;  
}


void print_data( geral_data_t &data_c ){

	std::cout <<"\n PRINT BEG *****************************************  \n\n"; 


	for (std::vector<chain_t>::iterator it = (data_c.cp).begin() ; it != (data_c.cp).end(); ++it){

		int n = (*it).size(); // numero de elementos na cadeia

		for (int m=0; m != n ; ++m) {
			std::cout <<"\n\n\n\n" << (*it)[m]->zaid << "  \n\n"; 

			for(int i=0;i!=NGPR;++i){
				std::cout << (*it)[m]->sig_f[i] << "  "; 
				std::cout << (*it)[m]->sig_a[i] << "  "; 
				std::cout << (*it)[m]->sig_g[i] << "\n"; 
			}
		}
	}

	std::cout <<"\n PRINT END *****************************************  \n\n"; 

	return;  

}




