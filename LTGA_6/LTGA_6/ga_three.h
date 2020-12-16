#ifndef GA_THREE_H_
#define GA_THREE_H_

double ga_three_obj_nlopt(unsigned n, const double* x, double* grad, void* para);
void ga_three_obj_nlopt_test();
void ga_three_obj_nlopt_list();
void ga_three_obj_nlopt_denselist();

double ga_three_obj_PSO(const double* x, const double* para);
void ga_three_PSO();

void ga_three_nlopt();

#endif