#ifndef DTF_INIT_FINALIZE_H_INCLUDED
#define DTF_INIT_FINALIZE_H_INCLUDED

int load_config(const char *ini_name, const char *service_name);
void clean_config();
int init_comp_comm();
void finalize_comp_comm();
void finalize_files();
#endif // dtf_INIT_FINALIZE_H_INCLUDED
