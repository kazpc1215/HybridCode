#include "hybrid.h"
#include "func.h"


#if FRAGMENTATION
void MassFlux(int i, CONST struct orbital_elements *ele_p, struct fragmentation *frag_p, CONST struct parameter *para_p){
  double F;
  double a_inv = 1.0/((ele_p+i)->axis);
  double alpha = (para_p->alpha);
  double v = ((frag_p+i)->v_ave);
  double sigma = ((frag_p+i)->sigma);
  double s_1 = (para_p->s_1);
  double s_2 = (para_p->s_2);
  double s_3 = (para_p->s_3);
  double Q_D = (para_p->Q_D);
  double h_0 = (para_p->h_0);

#if !defined(G) && !defined(M_0)
  F = - (2.0 - alpha) * (2.0 - alpha) / cbrt(M_MAX) * sigma * sigma * sqrt(a_inv * a_inv * a_inv) * h_0;
#else
  F = - (2.0 - alpha) * (2.0 - alpha) / cbrt(M_MAX) * sigma * sigma * sqrt(G * M_0 * a_inv * a_inv * a_inv) * h_0;
#endif

  F *= pow(v * v * 0.5 / Q_D, alpha - 1.0);
  F *= ((- log(EPSILON_FRAG) + 1.0 / (2.0-B_FRAG)) * s_1 + s_2 + s_3);

  ((frag_p+i)->flux) = F;

  return;
}
#endif


#if FRAGMENTATION
double Depletion_Time(int i, CONST struct fragmentation *frag_p){
  return - XI * ((frag_p+i)->sigma) / ((frag_p+i)->flux);
}
#endif


#if FRAGMENTATION
double MassDepletion(int i, double mass, double t_dyn, CONST struct fragmentation *frag_p){
  double t_frag = ((frag_p+i)->t_frag);
  double tau_dep = - ((frag_p+i)->sigma) / ((frag_p+i)->flux);
  if(i > global_n_p && t_dyn > 1.0E-10){  //惑星でない かつ 初期でない.
    return mass / (1.0 + (t_dyn - t_frag) / tau_dep);
  }else
    return mass;  //惑星の場合は変化させない.
}
#endif


#if FRAGMENTATION
double s_1_FRAG_trapezoid(int n, double dx, double ini, CONST struct parameter *para_p){
  return 0.5 * dx * (s_1_FRAG_integrand(ini + n * dx, para_p) + s_1_FRAG_integrand(ini + (n + 1) * dx, para_p));
}
#endif


#if FRAGMENTATION
double s_2_FRAG_trapezoid(int n, double dx, double ini, CONST struct parameter *para_p){
  return 0.5 * dx * (s_2_FRAG_integrand(ini + n * dx, para_p) + s_2_FRAG_integrand(ini + (n + 1) * dx, para_p));
}
#endif


#if FRAGMENTATION
double s_3_FRAG_trapezoid(int n, double dx, double ini, CONST struct parameter *para_p){
  return 0.5 * dx * (s_3_FRAG_integrand(ini + n * dx, para_p) + s_3_FRAG_integrand(ini + (n + 1) * dx, para_p));
}
#endif


#if FRAGMENTATION
double s_1_FRAG_integrand(double x, CONST struct parameter *para_p){
  return exp((2.0 - (para_p->alpha)) * x) / (1.0 + exp(x));
}
#endif


#if FRAGMENTATION
double s_2_FRAG_integrand(double x, CONST struct parameter *para_p){
  return - exp((2.0 - (para_p->alpha)) * x) / (1.0 + exp(x)) * (x - 2.0 * log(1 + exp(x)));
}
#endif


#if FRAGMENTATION
double s_3_FRAG_integrand(double x, CONST struct parameter *para_p){
  return exp((1.0 - (para_p->alpha)) * x) / (1.0 + exp(x)) * log(1.0 + exp(x));
}
#endif


#if FRAGMENTATION
double s_1_FRAG(struct parameter *para_p){
  int n, n_max;
  double dx, sum=0.0, sum_pre=0.0;
  double ini=-36.0, fin=36.0;
  double eps=1.0E-7;

  n_max = 1;
  do{
    dx = (fin - ini) / (double)n_max;
    sum_pre = sum;
    sum=0;
    for(n=0;n<n_max;n++){
      sum += s_1_FRAG_trapezoid(n, dx, ini, para_p);
    }
    //fprintf(fplog,"n_max=%d\n",n_max);
    n_max *= 2;
  }while(fabs(sum_pre-sum)>eps);

  return sum;
}
#endif


#if FRAGMENTATION
double s_2_FRAG(struct parameter *para_p){
  int n, n_max;
  double dx, sum=0.0, sum_pre=0.0;
  double ini=-36.0, fin=36.0;
  double eps=1.0E-7;

  n_max = 1;
  do{
    dx = (fin - ini) / (double)n_max;
    sum_pre = sum;
    sum=0;
    for(n=0;n<n_max;n++){
      sum += s_2_FRAG_trapezoid(n, dx, ini, para_p);
    }
    //fprintf(fplog,"n_max=%d\n",n_max);
    n_max *= 2;
  }while(fabs(sum_pre-sum)>eps);

  return sum;
}
#endif


#if FRAGMENTATION
double s_3_FRAG(struct parameter *para_p){
  int n, n_max;
  double dx,sum=0.0, sum_pre=0.0;
  double ini=-36.0, fin=36.0;
  double eps=1.0E-7;

  n_max = 1;
  do{
    dx = (fin - ini) / (double)n_max;
    sum_pre = sum;
    sum=0;
    for(n=0;n<n_max;n++){
      sum += s_3_FRAG_trapezoid(n, dx, ini, para_p);
    }
    //fprintf(fplog,"n_max=%d\n",n_max);
    n_max *= 2;
  }while(fabs(sum_pre-sum)>eps);

  return sum;
}
#endif
