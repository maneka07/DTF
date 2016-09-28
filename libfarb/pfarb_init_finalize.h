#ifndef PFARB_INIT_FINALIZE_H_INCLUDED
#define PFARB_INIT_FINALIZE_H_INCLUDED

int load_config(const char *ini_name, const char *service_name);
void clean_config();
int init_comp_comm();
void finalize_comp_comm();
int init_data_distr();

#endif // PFARB_INIT_FINALIZE_H_INCLUDED
