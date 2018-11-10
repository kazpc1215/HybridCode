#include "hybrid.h"
#include "func.h"

#if N_tr !=0
#if FRAGMENTATION
void NeighborSearch(int i, double t_dyn, CONST struct orbital_elements *ele_p, struct fragmentation *frag_p, CONST double x_0[][4]){

  int j,l,m;
  double radius[global_n+1];
  double theta[global_n+1];
  double S;
  double M;
  double v2;
  double delta_theta=DELTA_THETA;


  radius[i] = sqrt(x_0[i][1]*x_0[i][1] + x_0[i][2]*x_0[i][2]);
  theta[i] = atan2(x_0[i][2],x_0[i][1]);  //[-pi:pi]

  if(i<=global_n_p){
    fprintf(fplog,"i = %d, global_n_p = %d\n",i,global_n_p);
  }

  m = 1;
  do{

    for(l=1;l<=NEIGHBOR_MAX;l++){
      ((frag_p+i)->neighborlist[l]) = 0;
    }
    ((frag_p+i)->neighbornumber) = 0;

    ((frag_p+i)->delta_r_out) = (double)m * DELTA_R;  //外側.
    ((frag_p+i)->delta_r_in) = (double)m * DELTA_R;  //内側.

    S = 2.0 * (((frag_p+i)->delta_r_out) + ((frag_p+i)->delta_r_in)) * radius[i] * delta_theta;

    //fprintf(fplog,"i=%d\tS=%e\n",i,S);
    //fprintf(fplog,"delta_r[%d]=%f\n",i,((frag_p+i)->delta_r_in));

    l = 1;
    for(j=global_n_p+1;j<=global_n;j++){  //惑星抜き.
      if(j!=i){
	radius[j] = sqrt(x_0[j][1]*x_0[j][1] + x_0[j][2]*x_0[j][2]);
	theta[j] = atan2(x_0[j][2],x_0[j][1]);  //[-pi:pi]
	if((radius[j]-radius[i]>=0.0 && radius[j]-radius[i]<=((frag_p+i)->delta_r_out)) || (radius[i]-radius[j]>=0.0 && radius[i]-radius[j]<=((frag_p+i)->delta_r_in))){  //動径方向.
	  if(fabs(theta[j] - theta[i])<=delta_theta || 2.0*M_PI - fabs(theta[j] - theta[i])<=delta_theta){  //角度方向.
	    ((frag_p+i)->neighborlist[l]) = j;
	    l++;
	  }
	}
      }
    }

    if(NEIGHBOR_MAX < l-1){
      fprintf(fplog,"neighbornumber (%d) > NEIGHBOR_MAX (%d). Ending this program.\n",l-1,NEIGHBOR_MAX);
      exit(-1);
    }

    ((frag_p+i)->neighbornumber) = l-1;

    m++;

    //}while(((frag_p+i)->neighbornumber)<10);  //近傍粒子が10個未満なら、10個以上になるまでdelta_rをm倍に広げる.
  }while(((frag_p+i)->neighbornumber)<10 && m<10);  //近傍粒子が10個未満 かつ delta_r < 1.0 なら、1個以上になるまでdelta_rをm倍に広げる.

  v2 = 0.0;

  M = MassDepletion(i,((ele_p+i)->mass),t_dyn,frag_p);  //ターゲットiの質量も含める. 予測もする.

  if(((frag_p+i)->neighbornumber)!=0){
    for(j=1;j<=((frag_p+i)->neighbornumber);j++){
      v2 += SquareRandomVelocity(i,((frag_p+i)->neighborlist[j]),ele_p);

      M += MassDepletion(j,((ele_p+((frag_p+i)->neighborlist[j]))->mass),t_dyn,frag_p);  //領域iの総質量. 周りのトレーサーjの質量を予測してから足す.

      if(isnan(SquareRandomVelocity(i,((frag_p+i)->neighborlist[j]),ele_p))){
	fprintf(fplog,"i=%d,j=%d\tvij is nan.\n",i,((frag_p+i)->neighborlist[j]));
	fprintf(fplog,"i=%d\taxis=%e\tecc=%e\tinc=%e\tu=%e\tOmega=%e\tomega=%e\n",i,((ele_p+i)->axis),((ele_p+i)->ecc),((ele_p+i)->inc),((ele_p+i)->u),((ele_p+i)->Omega),((ele_p+i)->omega));
	fprintf(fplog,"j=%d\taxis=%e\tecc=%e\tinc=%e\tu=%e\tOmega=%e\tomega=%e\n",((frag_p+i)->neighborlist[j]),((ele_p+((frag_p+i)->neighborlist[j]))->axis),((ele_p+((frag_p+i)->neighborlist[j]))->ecc),((ele_p+((frag_p+i)->neighborlist[j]))->inc),((ele_p+((frag_p+i)->neighborlist[j]))->u),((ele_p+((frag_p+i)->neighborlist[j]))->Omega),((ele_p+((frag_p+i)->neighborlist[j]))->omega));
      }

    }
    if(isnan(v2)){
      fprintf(fplog,"i=%d\tv_tot_2 is nan.\n",i);
    }
    ((frag_p+i)->v_ave) = sqrt(v2/(double)((frag_p+i)->neighbornumber));  //領域iの平均速度(RMS).

    //fprintf(fplog,"i=%d\tmass=%e\n",i,ele[i].mass);
    //fprintf(fplog,"i=%d\tM=%e\n",i,M);

    ((frag_p+i)->sigma) = M / S;  //領域iの表面密度.
    ((frag_p+i)->n_s) = ((frag_p+i)->neighbornumber) / S;  //領域iの個数密度.
  }else{
    fprintf(fplog,"t_sys = %e [yr]\ti = %d\tneighbor number = %d\tv_ave, sigma, & n_s does NOT change.\n",t_dyn/2.0/M_PI,i,((frag_p+i)->neighbornumber));
    //((frag_p+i)->v_ave) = 0.0;
    //((frag_p+i)->sigma) = 0.0;
    //((frag_p+i)->n_s) = 0.0;
  }

  return;
}
#endif
#endif


#if N_tr != 0
#if FRAGMENTATION
double SquareRandomVelocity(int i, int j, CONST struct orbital_elements *ele_p){
  double eij2;
  double iij2;

  eij2 = ((ele_p+i)->ecc) * ((ele_p+i)->ecc) + ((ele_p+j)->ecc) * ((ele_p+j)->ecc) - 2.0 * ((ele_p+i)->ecc) * ((ele_p+j)->ecc) * cos(((ele_p+i)->omega) + ((ele_p+i)->Omega) - ((ele_p+j)->omega) - ((ele_p+j)->Omega));

  iij2 = ((ele_p+i)->inc) * ((ele_p+i)->inc) + ((ele_p+j)->inc) * ((ele_p+j)->inc) - 2.0 * ((ele_p+i)->inc) * ((ele_p+j)->inc) * cos(((ele_p+i)->Omega) - ((ele_p+j)->Omega));


  /*
  if(isnan((ele_p+i)->ecc)){
    fprintf(fplog,"i=%d\tecc is nan. (in RandomVelocity)\n",i);
  }
  if(isnan((ele_p+j)->ecc)){
    fprintf(fplog,"j=%d\tecc is nan. (in RandomVelocity)\n",j);
  }
  if(isnan((ele_p+i)->inc)){
    fprintf(fplog,"i=%d\tinc is nan. (in RandomVelocity).\n",i);
  }
  if(isnan((ele_p+j)->inc)){
    fprintf(fplog,"j=%d\tinc is nan. (in RandomVelocity)\n",j);
  }
  if(isnan((ele_p+i)->omega)){
    fprintf(fplog,"i=%d\tomega is nan. (in RandomVelocity)\n",i);
  }
  if(isnan((ele_p+j)->omega)){
    fprintf(fplog,"j=%d\tomega is nan. (in RandomVelocity)\n",j);
  }
  if(isnan((ele_p+i)->Omega)){
    fprintf(fplog,"i=%d\tOmega is nan. (in RandomVelocity)\n",i);
  }
  if(isnan((ele_p+j)->Omega)){
    fprintf(fplog,"j=%d\tOmega is nan. (in RandomVelocity)\n",j);
  }
  if(isnan(eij2)){
    fprintf(fplog,"i=%d,j=%d\teij2 is nan. (in RandomVelocity)\n",i,j);
    return -1;
  }
  if(isnan(iij2)){
    fprintf(fplog,"i=%d,j=%d\tiij2 is nan. (in RandomVelocity)\n",i,j);
    return -1;
  }
  if(isnan((ele_p+i)->axis)){
    fprintf(fplog,"i=%d\taxis is nan. (in RandomVelocity)\taxis=%f\n",i,((ele_p+i)->axis));
  }
  */

  if(
#if !defined(G) && !defined(M_0)
     isnan(sqrt((eij2 + iij2)/((ele_p+i)->axis)))
#else
     isnan(sqrt((eij2 + iij2)*G*M_0/((ele_p+i)->axis)))
#endif
     ){
    fprintf(fplog,"i=%d,j=%d\tvij is nan.(in RandomVelocity)\n",i,j);
    fprintf(fplog,"eij2=%e\tiij2=%e\taxis=%e\n",eij2,iij2,((ele_p+i)->axis));
  }


#if !defined(G) && !defined(M_0)
  return (eij2 + iij2) / ((ele_p+i)->axis);
#else
  return (eij2 + iij2) * G * M_0 / ((ele_p+i)->axis);
#endif
}
#endif
#endif
