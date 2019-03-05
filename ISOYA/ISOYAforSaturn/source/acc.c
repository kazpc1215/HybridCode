#include "hybrid.h"
#include "func.h"


//inline関数

/*相互重力加速度*/
static inline ALWAYS_INLINE double Acceleration_ij(int i, int j, int k, CONST double m_j, CONST double x_0[][4], CONST double abs_r[]){
  double rij3;

#ifndef EPSILON
  rij3 = abs_r[j]*abs_r[j]*abs_r[j];
#else
  rij3 = (abs_r[j]*abs_r[j] + EPSILON*EPSILON)*sqrt(abs_r[j]*abs_r[j] + EPSILON*EPSILON);
#endif
  rij3 = 1.0/rij3;

#ifndef G
  return m_j * (x_0[j][k] - x_0[i][k]) * rij3;
#else
  return G * m_j * (x_0[j][k] - x_0[i][k]) * rij3;
#endif
}


/*相互重力加加速度*/
static inline ALWAYS_INLINE double dAcceleration_ij(int i, int j, int k, CONST double m_j, CONST double x_0[][4], CONST double v_0[][4], CONST double r_dot_v_ij[], CONST double abs_r[]){
  double rij3;
  double rij5;

#ifndef EPSILON
  rij3 = abs_r[j]*abs_r[j]*abs_r[j];
  rij5 = rij3*abs_r[j]*abs_r[j];
#else
  rij3 = (abs_r[j]*abs_r[j] + EPSILON*EPSILON)*sqrt(abs_r[j]*abs_r[j] + EPSILON*EPSILON);
  rij5 = (abs_r[j]*abs_r[j] + EPSILON*EPSILON)*(abs_r[j]*abs_r[j] + EPSILON*EPSILON)*sqrt(abs_r[j]*abs_r[j] + EPSILON*EPSILON);
#endif
  rij3 = 1.0/rij3;
  rij5 = 1.0/rij5;

#ifndef G
  return m_j * ((v_0[j][k] - v_0[i][k]) * rij3 - 3.0 * r_dot_v_ij[j] * (x_0[j][k] - x_0[i][k]) * rij5);
#else
  return G * m_j * ((v_0[j][k] - v_0[i][k]) * rij3 - 3.0 * r_dot_v_ij[j] * (x_0[j][k] - x_0[i][k]) * rij5);
#endif
}


/*加速度indirect項*/
static inline ALWAYS_INLINE double Acceleration_indirect(int i, int k, CONST double m_i, CONST double x_0[][4], CONST double r_0[]){
  double r3;
  r3 = r_0[i]*r_0[i]*r_0[i];
  r3 = 1.0/r3;

#ifndef G
  return - m_i * x_0[i][k] * r3;
#else
  return - G * m_i * x_0[i][k] * r3;
#endif
}


/*加加速度indirect項*/
static inline ALWAYS_INLINE double dAcceleration_indirect(int i, int k, CONST double m_i, CONST double x_0[][4], CONST double v_0[][4], CONST double r_0[], CONST double r_dot_v[]){
  double r3;
  double r5;
  r3 = r_0[i]*r_0[i]*r_0[i];
  r5 = r3*r_0[i]*r_0[i];
  r3 = 1.0/r3;
  r5 = 1.0/r5;

#ifndef G
  return - m_i * (v_0[i][k] * r3 - 3.0 * r_dot_v[i] * x_0[i][k] * r5);
#else
  return - G * m_i * (v_0[i][k] * r3 - 3.0 * r_dot_v[i] * x_0[i][k] * r5);
#endif
}


/*外力加速度*/
static inline ALWAYS_INLINE double External_Acceleration(int i, int k, CONST double x_0[][4], CONST double r_0[]){
  double r3;
  r3 = r_0[i]*r_0[i]*r_0[i];
  r3 = 1.0/r3;

#if !defined(G) && !defined(M_0)
  return - x_0[i][k] * r3;
#else
  return - G * M_0 * x_0[i][k] * r3;
#endif
}


/*外力加加速度*/
static inline ALWAYS_INLINE double External_dAcceleration(int i, int k, CONST double x_0[][4], CONST double v_0[][4], CONST double r_0[], CONST double r_dot_v[]){
  double r3;
  double r5;
  r3 = r_0[i]*r_0[i]*r_0[i];
  r5 = r3*r_0[i]*r_0[i];
  r3 = 1.0/r3;
  r5 = 1.0/r5;

#if !defined(G) && !defined(M_0)
  return - (v_0[i][k] * r3 - 3.0 * r_dot_v[i] * x_0[i][k] * r5);
#else
  return - G M_0 * (v_0[i][k] * r3 - 3.0 * r_dot_v[i] * x_0[i][k] * r5);
#endif
}


/*i_sys のみの修正子*/
void Corrector_sys(int n_ite, int i_sys, CONST struct orbital_elements *ele_p, CONST double x_p[][4], CONST double v_p[][4], CONST double r_p[], double x_c[][4], double v_c[][4], double r_c[], double v2_c[], double r_dot_v[], CONST double a_0[][4], CONST double adot_0[][4], double a[][4], double adot[][4], double adot2_dt2[][4], double adot3_dt3[][4], CONST double dt_[]
#if FRAGMENTATION
		   , double t_dyn
		   , CONST struct fragmentation *frag_p
#endif
		   ){

  int j, k;
  double dt = dt_[i_sys];
  double abs_r[
#if INTERACTION_ALL
	       global_n
#elif INTERACTION_PLANET_TRACER
	       global_n
#elif INTERACTION_TEST_PARTICLE
	       global_n_p
#endif
	       +1];
  double r_dot_v_ij[
#if INTERACTION_ALL
		    global_n
#elif INTERACTION_PLANET_TRACER
		    global_n
#elif INTERACTION_TEST_PARTICLE
		    global_n_p
#endif
		    +1];
  double a_tmp[4]={}, adot_tmp[4]={}, adot2_dt2_tmp[4]={}, adot3_dt3_tmp[4]={};
  double x_c_tmp[4]={}, v_c_tmp[4]={};


  for(j=1;j<=
#if INTERACTION_ALL
	global_n
#elif INTERACTION_PLANET_TRACER
	global_n
#elif INTERACTION_TEST_PARTICLE
	global_n_p
#endif
	;++j){
    if(i_sys!=j){
      if(n_ite == 0){
	abs_r[j] = RelativeDistance(i_sys,j,x_p);  //絶対値.
	r_dot_v_ij[j] = RelativeInnerProduct(i_sys,j,x_p,v_p);  //r_ij,v_ijの内積.
      }else if(n_ite <= ITE_MAX){  //iterationではx_c,v_cを使う.
	abs_r[j] = RelativeDistance(i_sys,j,x_c);  //絶対値.
	r_dot_v_ij[j] = RelativeInnerProduct(i_sys,j,x_c,v_c);  //r_ij,v_ijの内積.
      }
    }else{
      abs_r[i_sys] = 0.0;
      r_dot_v_ij[i_sys] = 0.0;
    }
  }


  for(k=1;k<=3;++k){
    if(n_ite == 0){
      a_tmp[k] = All_Acceleration(i_sys,k,ele_p,x_p,r_p,abs_r
#if FRAGMENTATION
				  ,t_dyn,frag_p
#endif
				  );
      adot_tmp[k] = All_dAcceleration(i_sys,k,ele_p,x_p,v_p,r_dot_v,r_dot_v_ij,r_p,abs_r
#if FRAGMENTATION
				  ,t_dyn,frag_p
#endif
				      );
    }else if(n_ite <= ITE_MAX){  //iterationではx_c,v_cを使う.
      a_tmp[k] = All_Acceleration(i_sys,k,ele_p,x_c,r_c,abs_r
#if FRAGMENTATION
				  ,t_dyn,frag_p
#endif
				  );
      adot_tmp[k] = All_dAcceleration(i_sys,k,ele_p,x_c,v_c,r_dot_v,r_dot_v_ij,r_c,abs_r
#if FRAGMENTATION
				  ,t_dyn,frag_p
#endif
				      );
    }
  }


  for(k=1;k<=3;++k){  //修正子.

    adot2_dt2_tmp[k] = - 6.0 * (a_0[i_sys][k] - a_tmp[k]) - (4.0 * adot_0[i_sys][k] + 2.0 * adot_tmp[k]) * dt; //第2次導関数.
    adot3_dt3_tmp[k] = 12.0 * (a_0[i_sys][k] - a_tmp[k]) + 6.0 * (adot_0[i_sys][k] + adot_tmp[k]) * dt; //第3次導関数.

    //x_c_tmp[k] = x_p[i_sys][k] + adot2_dt2_tmp[k]*dt_[i_sys]*dt_[i_sys]/24.0 + adot3_dt3_tmp[k]*dt_[i_sys]*dt_[i_sys]/120.0;
    //v_c_tmp[k] = v_p[i_sys][k] + adot2_dt2_tmp[k]*dt_[i_sys]/6.0 +adot3_dt3_tmp[k]*dt_[i_sys]/24.0;

    x_c_tmp[k] = x_p[i_sys][k] + dt * dt * 0.125 * INV_3 * (adot2_dt2_tmp[k] + adot3_dt3_tmp[k] * 0.2);
    v_c_tmp[k] = v_p[i_sys][k] + dt * 0.5 * INV_3 * (adot2_dt2_tmp[k] + adot3_dt3_tmp[k] * 0.25);
  }  //k loop


  r_c[i_sys] = sqrt(x_c_tmp[1]*x_c_tmp[1] + x_c_tmp[2]*x_c_tmp[2] + x_c_tmp[3]*x_c_tmp[3]);
  v2_c[i_sys] = v_c_tmp[1]*v_c_tmp[1] + v_c_tmp[2]*v_c_tmp[2] + v_c_tmp[3]*v_c_tmp[3];
  r_dot_v[i_sys] = x_c_tmp[1]*v_c_tmp[1] + x_c_tmp[2]*v_c_tmp[2] + x_c_tmp[3]*v_c_tmp[3];


  for(k=1;k<=3;++k){
    a[i_sys][k] = a_tmp[k];
    adot[i_sys][k] = adot_tmp[k];
    adot2_dt2[i_sys][k] = adot2_dt2_tmp[k];
    adot3_dt3[i_sys][k] = adot3_dt3_tmp[k];
    x_c[i_sys][k] = x_c_tmp[k];
    v_c[i_sys][k] = v_c_tmp[k];
  }

  return;
}



/*i_sys のみのiteration*/
/*
void Iteration_sys(int i_sys,CONST struct orbital_elements *ele_p,CONST double x_p[][4],CONST double v_p[][4],double x_c[][4],double v_c[][4],double r_c[],double v2_c[],double r_dot_v[],CONST double a_0[][4],CONST double adot_0[][4],double a[][4],double adot[][4],double adot2_dt2[][4],double adot3_dt3[][4],CONST double dt_[]){

  int j,k;
  double dt = dt_[i_sys];
  double abs_r[
#if INTERACTION_ALL
	       global_n
#elif INTERACTION_PLANET_TRACER
	       global_n
#elif INTERACTION_TEST_PARTICLE
	       global_n_p
#endif
	       +1];
  double r_dot_v_ij[
#if INTERACTION_ALL
		    global_n
#elif INTERACTION_PLANET_TRACER
		    global_n
#elif INTERACTION_TEST_PARTICLE
		    global_n_p
#endif
		    +1];
  double a_tmp[4]={},adot_tmp[4]={},adot2_dt2_tmp[4]={},adot3_dt3_tmp[4]={};
  double x_c_tmp[4]={},v_c_tmp[4]={};


  for(j=1;j<=
#if INTERACTION_ALL
	global_n
#elif INTERACTION_PLANET_TRACER
	global_n
#elif INTERACTION_TEST_PARTICLE
	global_n_p
#endif
	;++j){
    if(i_sys!=j){
      abs_r[j] = RelativeDistance(i_sys,j,x_c);  //絶対値.
      r_dot_v_ij[j] = RelativeInnerProduct(i_sys,j,x_c,v_c);  //r_ij,v_ijの内積.
    }else{
      abs_r[i_sys] = 0.0;
      r_dot_v_ij[i_sys] = 0.0;
    }
  }  //j loop


  for(k=1;k<=3;++k){
    a_tmp[k] = All_Acceleration(i_sys,k,ele_p,x_c,r_c,abs_r);
    adot_tmp[k] = All_dAcceleration(i_sys,k,ele_p,x_c,v_c,r_dot_v,r_dot_v_ij,r_c,abs_r);
  }


  //修正子 (iteration).
  for(k=1;k<=3;++k){
    adot2_dt2_tmp[k] = -6.0*(a_0[i_sys][k] - a_tmp[k]) - (4.0*adot_0[i_sys][k] + 2*adot_tmp[k])*dt; //第2次導関数.
    adot3_dt3_tmp[k] = 12.0*(a_0[i_sys][k] - a_tmp[k]) + 6.0*(adot_0[i_sys][k] + adot_tmp[k])*dt; //第3次導関数.

    //x_c_tmp[k] = x_p[i_sys][k] + adot2_dt2_tmp[k]*dt_[i_sys]*dt_[i_sys]/24.0 + adot3_dt3_tmp[k]*dt_[i_sys]*dt_[i_sys]/120.0;
    //v_c_tmp[k] = v_p[i_sys][k] + adot2_dt2_tmp[k]*dt_[i_sys]/6.0 +adot3_dt3_tmp[k]*dt_[i_sys]/24.0;

    x_c_tmp[k] = x_p[i_sys][k] + dt*dt*0.125*INV_3*(adot2_dt2_tmp[k] + adot3_dt3_tmp[k]*0.2);
    v_c_tmp[k] = v_p[i_sys][k] + dt*0.5*INV_3*(adot2_dt2_tmp[k] +adot3_dt3_tmp[k]*0.25);
  }  //k loop


  r_c[i_sys] = sqrt(x_c_tmp[1]*x_c_tmp[1] + x_c_tmp[2]*x_c_tmp[2] + x_c_tmp[3]*x_c_tmp[3]);
  v2_c[i_sys] = v_c_tmp[1]*v_c_tmp[1] + v_c_tmp[2]*v_c_tmp[2] + v_c_tmp[3]*v_c_tmp[3];
  r_dot_v[i_sys] = x_c_tmp[1]*v_c_tmp[1] + x_c_tmp[2]*v_c_tmp[2] + x_c_tmp[3]*v_c_tmp[3];


  for(k=1;k<=3;++k){
    a[i_sys][k] = a_tmp[k];
    adot[i_sys][k] = adot_tmp[k];
    adot2_dt2[i_sys][k] = adot2_dt2_tmp[k];
    adot3_dt3[i_sys][k] = adot3_dt3_tmp[k];
    x_c[i_sys][k] = x_c_tmp[k];
    v_c[i_sys][k] = v_c_tmp[k];
  }

  return;
}
*/

/*全加速度*/
double All_Acceleration(int i, int k, CONST struct orbital_elements *ele_p, CONST double x_0[][4], CONST double r_0[], CONST double abs_r[]
#if FRAGMENTATION
			, double t_dyn
			, CONST struct fragmentation *frag_p
#endif
			){
  int j;
  double a_0;
  double m_j;

  a_0 = External_Acceleration(i,k,x_0,r_0);

#if INTERACTION_PLANET_TRACER
  if(i>global_n_p){  //i_sysがトレーサーの場合，惑星からの影響を計算.
#endif

#if INTERACTION_ALL
#pragma omp parallel for reduction(+:a_0)
#endif
    for(j=1;j<=
#if INTERACTION_ALL
	  global_n
#elif INTERACTION_PLANET_TRACER
	  global_n_p
#elif INTERACTION_TEST_PARTICLE
	  global_n_p
#endif
	  ;++j){

#if FRAGMENTATION
      m_j = MassDepletion(j,((ele_p+j)->mass),t_dyn,frag_p);
#else
      m_j = ((ele_p+j)->mass);
#endif

#if INDIRECT_TERM
      a_0 += Acceleration_indirect(j,k,m_j,x_0,r_0);
#endif

      if(i!=j){
	a_0 += Acceleration_ij(i,j,k,m_j,x_0,abs_r);
      }
    }



#if INTERACTION_PLANET_TRACER
  }else if(i<=global_n_p){  //i_sysが惑星の場合，自身以外のすべてからの影響を計算.

#pragma omp parallel for reduction(+:a_0)
    for(j=1;j<=global_n;++j){

#if FRAGMENTATION
      m_j = MassDepletion(j,((ele_p+j)->mass),t_dyn,frag_p);
#else
      m_j = ((ele_p+j)->mass);
#endif

#if INDIRECT_TERM
      a_0 += Acceleration_indirect(j,k,m_j,x_0,r_0);
#endif

      if(i!=j){
	a_0 += Acceleration_ij(i,j,k,m_j,x_0,abs_r);
      }
    }

  }
#endif

  return a_0;
}


/*全加加速度*/
double All_dAcceleration(int i, int k, CONST struct orbital_elements *ele_p, CONST double x_0[][4], CONST double v_0[][4], CONST double r_dot_v[], CONST double r_dot_v_ij[], CONST double r_0[], CONST double abs_r[]
#if FRAGMENTATION
			 , double t_dyn
			 , CONST struct fragmentation *frag_p
#endif
			 ){
  int j;
  double adot_0;
  double m_j;

  adot_0 = External_dAcceleration(i,k,x_0,v_0,r_0,r_dot_v);

#if INTERACTION_PLANET_TRACER
  if(i>global_n_p){  //i_sysがトレーサーの場合，惑星からの影響を計算.
#endif

#if INTERACTION_ALL
#pragma omp parallel for reduction(+:adot_0)
#endif
    for(j=1;j<=
#if INTERACTION_ALL
	  global_n
#elif INTERACTION_PLANET_TRACER
	  global_n_p
#elif INTERACTION_TEST_PARTICLE
	  global_n_p
#endif
	  ;++j){

#if FRAGMENTATION
      m_j = MassDepletion(j,((ele_p+j)->mass),t_dyn,frag_p);
#else
      m_j = ((ele_p+j)->mass);
#endif

#if INDIRECT_TERM
      adot_0 += dAcceleration_indirect(j,k,m_j,x_0,v_0,r_0,r_dot_v);
#endif

      if(i!=j){
	adot_0 += dAcceleration_ij(i,j,k,m_j,x_0,v_0,r_dot_v_ij,abs_r);
      }
    }

#if INTERACTION_PLANET_TRACER
  }else if(i<=global_n_p){  //i_sysが惑星の場合，自身以外のすべてからの影響を計算.

#pragma omp parallel for reduction(+:adot_0)
    for(j=1;j<=global_n;++j){

#if FRAGMENTATION
      m_j = MassDepletion(j,((ele_p+j)->mass),t_dyn,frag_p);
#else
      m_j = ((ele_p+j)->mass);
#endif

#if INDIRECT_TERM
      adot_0 += dAcceleration_indirect(j,k,m_j,x_0,v_0,r_0,r_dot_v);
#endif

      if(i!=j){
	adot_0 += dAcceleration_ij(i,j,k,m_j,x_0,v_0,r_dot_v_ij,abs_r);
      }
    }

  }
#endif

  return adot_0;
}
