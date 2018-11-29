/* This function was automatically generated by CasADi */
#ifdef __cplusplus
extern "C" {
#endif

#ifndef real_t
#define real_t double
#endif /* real_t */

int model_mpc_casadi(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
void model_mpc_casadi_incref(void);
void model_mpc_casadi_decref(void);
int model_mpc_casadi_n_in(void);
int model_mpc_casadi_n_out(void);
const char* model_mpc_casadi_name_in(int i);
const char* model_mpc_casadi_name_out(int i);
const int* model_mpc_casadi_sparsity_in(int i);
const int* model_mpc_casadi_sparsity_out(int i);
int model_mpc_casadi_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
#ifdef __cplusplus
} /* extern "C" */
#endif