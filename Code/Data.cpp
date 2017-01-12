#include "Data.h"

Data::Data (std::vector<double> vdata, std::vector<double>nvdata, std::vector<double>ab_data,
      std::vector <double>timings, int vtypes, int nvtypes) {
  swab_times = timings;
  n_vtypes = vtypes;
  n_nvtypes = nvtypes;
  n_tot = n_vtypes + n_nvtypes;
  swab_data_v = vdata;
  swab_data_nv = nvdata;
  ab_data = ab_data;
}