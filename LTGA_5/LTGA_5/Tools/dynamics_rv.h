#ifndef DYNAMICS_RV_H_
#define DYNAMICS_RV_H_

int dynamics_rv_top(double t, const double* x, double* dx, const double* dfpara);
int dynamics_rv_fop(double t, const double* x, double* dx, const double* dfpara);

#endif