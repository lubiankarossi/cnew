#ifndef STRUCTURES_H
#define STRUCTURES_H
#define NGPR 4

#include <stdbool.h>

#include <string>
#include <algorithm>
#include <vector>


#include <cstdio>
#include <iostream>
#include <fstream>

#include <stdlib.h>


enum {
        no_coupling = 1,
        decay       = 2,
        absorption  = 3,
        capture     = 4
};

struct base_data_t {

   double   sig_a[NGPR], /* abs. cross section     */
            sig_f[NGPR], /* fis. cross section     */
            sig_g[NGPR]; /* cap. cross section     */

   double   lamb,        /* decay constant         */
            alpha;       /* branching factor       */ 

   unsigned reac_id,     /* reaction id (in chain) */
            zaid,        /* nuclide id             */
            node_number;

   bool sum_flag,
        fis_flag;
};

struct base_yield_t {

   std::vector< std::vector<double> > yield;
   std::string id;
   double   Q;       
};

typedef std::vector< base_data_t*  > chain_t;
typedef std::vector< base_yield_t  > yield_t;

struct geral_data {
   std::vector< chain_t > cp; /* vector with all principal chains */
   std::vector< chain_t > ca; /* vector with all actinides chains */
   yield_t                fp; /* vector with all fission products */
};

typedef geral_data geral_data_t;

void clear_data_vector(chain_t &data);
void clear_data(geral_data_t &data);

void read_data(std::string file, geral_data_t &data_c);

/* replace mtd */

void replace_data(std::string file, geral_data_t &data_c );


void print_data( geral_data_t &data_c );


extern "C"{

double      get_double(unsigned pos, unsigned length, std::string line);
int         get_int   (unsigned pos, unsigned length, std::string line);
std::string get_string(unsigned pos, unsigned length, std::string line);

};

#endif /* end guards STRUCTURES_H */


