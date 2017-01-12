#include <vector>
#include <stdio.h>      /* printf */
#include <math.h>       /* log */

class Data {
  std::vector <double> swab_data_v; // nrow=total number of vaccinated, ncol=swabs
  std::vector <double> swab_data_nv; // nrow=total number of vaccinated, ncol=swabs
  std::vector <double> swab_times; // timing of swabs
  std::vector <double> ab_data; // nrow=total number of vaccinated, ncol=serotypes in a vaccine
  int n_vacc;
  int n_nvacc;
  int n_swabs;
  int n_vtypes;
  int n_nvtypes;
  int n_tot;
public:
  Data (std::vector<double>, std::vector<double>, std::vector<double>,
        std::vector <double>, int, int);
};
