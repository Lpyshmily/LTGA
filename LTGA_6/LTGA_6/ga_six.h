#ifndef GA_SIX_H_
#define GA_SIX_H_

double ga_six_obj_nlopt(unsigned n, const double* x, double* grad, void* para);
void ga_six_obj_nlopt_test();

double ga_six_obj_PSO(const double* x, const double* para);
void ga_six_PSO();

void ga_six_nlopt();

#endif