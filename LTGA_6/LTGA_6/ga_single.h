#ifndef GA_SINGLE_H_
#define GA_SINGLE_H_

double ga_single1_obj(double x);
void ga_single1_obj_test();
void ga_single1_list();

double ga_single2_obj(double x);
void ga_single2_list();

double ga_single1_obj_nlopt(unsigned n, const double* x, double* grad, void* para);
void ga_single1_nlopt();

#endif