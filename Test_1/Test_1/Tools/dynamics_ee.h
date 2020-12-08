#ifndef DYNAMICS_EE_H_
#define DYNAMICS_EE_H_

int dynamics_ee_top(double t, const double* x, double* dx, const double* dfpara);
int dynamics_ee_fop(double t, const double* x, double* dx, const double* dfpara);

#endif