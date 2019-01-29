#include "hybrid.h"
#include "func.h"


/*初期 タイムステップ計算*/
double Timestep_i_0(int i, CONST double a_0[][4], CONST double adot_0[][4]){
  int k;
  double abs_a = 0.0;
  double abs_adot = 0.0;

  for(k=1;k<=3;++k){
    abs_a += a_0[i][k] * a_0[i][k];
    abs_adot += adot_0[i][k] * adot_0[i][k];
  }  //k loop

  abs_a = sqrt(abs_a);
  abs_adot = sqrt(abs_adot);

  //fprintf(fplog,"abs_a[%d]=%f\tabs_adot[%d]=%f\n",i,abs_a[i],i,abs_adot[i]);
  return ETA * abs_a / abs_adot;
}


/*i_sys のみのタイムステップ計算*/
double Timestep_i_sys(int i_sys, CONST double a[][4], CONST double adot[][4], CONST double adot2_dt2[][4], CONST double adot3_dt3[][4], CONST double dt_[]){

  int k;
  double dt_inv = 1.0 / dt_[i_sys];

  double abs_a = 0.0;
  double abs_adot = 0.0;
  double abs_adot2 = 0.0;
  double abs_adot3 = 0.0;
  for(k=1;k<=3;++k){
    abs_a += a[i_sys][k] * a[i_sys][k];
    abs_adot += adot[i_sys][k] * adot[i_sys][k];
    abs_adot2 += (adot2_dt2[i_sys][k] + adot3_dt3[i_sys][k])*(adot2_dt2[i_sys][k] + adot3_dt3[i_sys][k]) * dt_inv * dt_inv * dt_inv * dt_inv;
    abs_adot3 += adot3_dt3[i_sys][k] * adot3_dt3[i_sys][k] * dt_inv * dt_inv * dt_inv * dt_inv * dt_inv * dt_inv;
  }  //k loop
  abs_a = sqrt(abs_a);
  abs_adot = sqrt(abs_adot);
  abs_adot2 = sqrt(abs_adot2);
  abs_adot3 = sqrt(abs_adot3);

  return ETA * sqrt((abs_a * abs_adot2 + abs_adot * abs_adot) / (abs_adot * abs_adot3 + abs_adot2 * abs_adot2));
}
