#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#define CUTOFF  1.0E-18
#define LIM_DIV 1.0E-18

#define CN 1.0E-6
#define CM 1.0E-6


#include <stdlib.h>
#include <stdbool.h>

#include <vector>
#include <string>

#include "storage.h"
#include "structures.h"

struct parameters {

   double  *alfa ,
           *gama ,
             *Af ,
             *Aa ,
             *Ac ,
         *lambda ,
           *beta ,
          *betaa ;

   unsigned *zaid;

   bool *sum_flag;
};

typedef parameters param_t;

double exp_c(double x);

void init_parameters(param_t &p, size_t n);

void clean_parameters(param_t &p);

void calc_param(param_t &pr, double fluxo[], std::vector<chain_t>::iterator it, size_t n, double dt);

void calc_chains_cp(storage &nuclids, storage &nuclids_t, param_t &pr, size_t n, double dt);

void calc_chains_cp_ad(storage &adjunts, storage &adjunts_t, param_t &pr, size_t n, double dt, double t_f);

// calc int cof

double int_chains_cp_ad(storage &nuclids, param_t &pr, int N, int p, double t_f);

double int_chains_cp_ad2(storage &nuclids, param_t &pr, int N, int p, double t_f);

int  find_el(int zaid, param_t &pr, size_t n );

//  manipulação de string

int get_substr_i(std::string str_);

#endif /* end guards FUNCTIONS_H */

















