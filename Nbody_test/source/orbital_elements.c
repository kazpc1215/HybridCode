#include "hybrid.h"
#include "func.h"
#include "SFMTdir/SFMT.h"


/*P計算*/
double Calculate_P(int i, int k, CONST struct orbital_elements *ele_p){
  if(k==1){
    return cos(((ele_p+i)->omega))*cos(((ele_p+i)->Omega)) - sin(((ele_p+i)->omega))*sin(((ele_p+i)->Omega))*cos(((ele_p+i)->inc));
  }else if(k==2){
    return cos(((ele_p+i)->omega))*sin(((ele_p+i)->Omega)) + sin(((ele_p+i)->omega))*cos(((ele_p+i)->Omega))*cos(((ele_p+i)->inc));
  }else{
    return sin(((ele_p+i)->omega))*sin(((ele_p+i)->inc));
  }
}


/*Q計算*/
double Calculate_Q(int i, int k, CONST struct orbital_elements *ele_p){
  if(k==1){
    return -sin(((ele_p+i)->omega))*cos(((ele_p+i)->Omega)) - cos(((ele_p+i)->omega))*sin(((ele_p+i)->Omega))*cos(((ele_p+i)->inc));
  }else if(k==2){
    return -sin(((ele_p+i)->omega))*sin(((ele_p+i)->Omega)) + cos(((ele_p+i)->omega))*cos(((ele_p+i)->Omega))*cos(((ele_p+i)->inc));
  }else{
    return cos(((ele_p+i)->omega))*sin(((ele_p+i)->inc));
  }
}


/*惑星の初期軌道要素*/
void InitialOrbitalElements_Planet(int i,struct orbital_elements *ele_p){

  sprintf((ele_p+i)->name,"Planet%02d",i);
  (ele_p+i)->mass = PLANET_MASS;
  //(ele_p+i)->axis = PLANET_AXIS;  //軌道長半径axisはすでに求めてある.
  (ele_p+i)->ecc = PLANET_ECC;
  (ele_p+i)->inc = PLANET_INC;
  (ele_p+i)->u = rand_func() * 2.0 * M_PI;
  (ele_p+i)->omega = rand_func() * 2.0 * M_PI;
  (ele_p+i)->Omega = rand_func() * 2.0 * M_PI;

#ifndef M_0
  (ele_p+i)->r_h = ((ele_p+i)->axis)*cbrt(((ele_p+i)->mass)/3.0);
#else
  (ele_p+i)->r_h = ((ele_p+i)->axis)*cbrt(((ele_p+i)->mass)/M_0/3.0);
#endif

  (ele_p+i)->radius = cbrt(3.0/4.0/M_PI*((ele_p+i)->mass)*1.989E33/PLANET_DENSITY)/1.496E13;
  (ele_p+i)->orinum = i;

  return;
}


#if ORBITING_SMALL_PARTICLE
/*トレーサーの初期軌道要素*/
void InitialOrbitalElements_Tracer(int i, double x_0[][4], struct orbital_elements *ele_p){

  //惑星の位置x_0[][4]はすでに求めてあることが前提.
#if N_p == 3
  double orbital_r_min = ((ele_p+1)->axis) / MutualHillRadius_to_SemimajorAxis(0.5*DELTA_HILL);
  double orbital_r_max = ((ele_p+3)->axis) * MutualHillRadius_to_SemimajorAxis(0.5*DELTA_HILL);
#elif N_p == 1
  double orbital_r_min = ((ele_p+1)->axis) / MutualHillRadius_to_SemimajorAxis(0.5*DELTA_HILL);
  double orbital_r_max = ((ele_p+1)->axis) * MutualHillRadius_to_SemimajorAxis(0.5*DELTA_HILL);
#elif N_p == 0 && N_SMALL == 1000 && N_LARGE == 0
  double orbital_r_min = 0.9293;  //10g/cm^2 @1AU.
  double orbital_r_max = 1.0707;
#elif N_p == 0 && N_SMALL == 800 && N_LARGE == 200
  double orbital_r_min = 0.8868;  //10g/cm^2 @1AU.
  double orbital_r_max = 1.1132;
#endif

  int j=0,k=0,flag=0;
  //double peri=0.0,apo=0.0;
  double orbital_r=0.0;

  //double index = 0.5;  //f(x) = C x^{-p}.

  sprintf((ele_p+i)->name,"tracer%06d",i-global_n_p);

  if(i<=N_SMALL){
    (ele_p+i)->mass = M_SMALL;  //質量.
  }else if(i<=N_tr){
    (ele_p+i)->mass = M_LARGE;  //質量.
  }

  do{
    flag = 0;

    (ele_p+i)->axis = rand_func() * (orbital_r_max - orbital_r_min) + orbital_r_min;  //面密度がa^{-1}に比例する分布. ==> 面数密度が一様.
    //(ele_p+i)->axis = pow((pow(orbital_r_max,1.0-index) - pow(orbital_r_min,1.0-index)) * rand_func() + pow(orbital_r_min,1.0-index), 1.0/(1.0-index));  //べき分布.

#if RAYLEIGH_DISTRIBUTION
    (ele_p+i)->ecc = sqrt(-log(rand_func())) * ECC_RMS;  //離心率.  //Rayleigh分布.
    (ele_p+i)->inc = sqrt(-log(rand_func())) * INC_RMS;  //軌道傾斜角.  //Rayleigh分布.
#else
    (ele_p+i)->ecc = ECC_RMS;  //離心率.
    (ele_p+i)->inc = INC_RMS;  //軌道傾斜角.
#endif
    (ele_p+i)->u = rand_func() * 2.0 * M_PI;
    (ele_p+i)->omega = rand_func() * 2.0 * M_PI;
    (ele_p+i)->Omega = rand_func() * 2.0 * M_PI;

    for(k=1;k<=3;++k){
      x_0[i][k] = ((ele_p+i)->axis) * Calculate_P(i,k,ele_p) * (cos(((ele_p+i)->u)) - ((ele_p+i)->ecc)) + ((ele_p+i)->axis) * sqrt(1.0 - ((ele_p+i)->ecc) * ((ele_p+i)->ecc)) * Calculate_Q(i,k,ele_p) * sin(((ele_p+i)->u));
    }

    //peri = ((ele_p+i)->axis)*(1.0 - (ele_p+i)->ecc);
    //apo = ((ele_p+i)->axis)*(1.0 + (ele_p+i)->ecc);
    orbital_r = RadiusFromCenter(i,x_0);

    if(orbital_r > orbital_r_min && orbital_r < orbital_r_max){  //orbital_r_minからorbital_r_maxの範囲にいる場合.

      for(j=1;j<=global_n_p;++j){
	if(i!=j){
	  if(RelativeDistance(i,j,x_0)>SEPARATE_HILL*((ele_p+j)->r_h)){  //それぞれの惑星からSEPARATE_HILLヒル以上離れている場合.
	    flag += 1;
	  }
	}
      }
    }
  }while(flag < global_n_p);  //orbital_r_minからorbital_r_maxの範囲にいる場合，かつglobal_n_p個の全ての惑星からSEPARATE_HILLヒル以上離れている場合のみ抜け出せるループ.


#ifndef M_0
  (ele_p+i)->r_h = ((ele_p+i)->axis) * cbrt(((ele_p+i)->mass) / 3.0);
#else
  (ele_p+i)->r_h = ((ele_p+i)->axis) * cbrt(((ele_p+i)->mass) / M_0 / 3.0);
#endif

  (ele_p+i)->radius = cbrt(3.0 / 4.0 / M_PI * ((ele_p+i)->mass) * 1.989E33 / PLANET_DENSITY) / 1.496E13;
  (ele_p+i)->orinum = i;

  return;
}
#endif  /*ORBITING_SMALL_PARTICLE*/


#if EJECTION
void EjectionOfTracerFromPlanet(double x_0[][4], double v_0[][4], double v2_0[], double r_dot_v[], double r_0[], CONST struct orbital_elements *ele_p){
  static double x_eject[N_p+N_tr+1][4]={};
  static double v_eject[N_p+N_tr+1][4]={};
  int i,k;
  double tmp_x=0.0, tmp_y=0.0, tmp_r=0.0, tmp_v=0.0, tmp_theta=0.0, tmp_rand=0.0;
#ifndef G
  double ejection_velocity = sqrt(2.0 * ((ele_p+PLANET_OF_EJECTION)->mass) / ((ele_p+PLANET_OF_EJECTION)->radius));
#else
  double ejection_velocity = sqrt(2.0 * G * ((ele_p+PLANET_OF_EJECTION)->mass) / ((ele_p+PLANET_OF_EJECTION)->radius));
#endif

  for(i=global_n_p+1;i<=global_n;++i){

    sprintf((ele_p+i)->name,"tracer%06d",i-global_n_p);
    (ele_p+i)->mass = M_TOT / (double)N_tr;
    (ele_p+i)->orinum = i;


    //coneを作る.
    //位置.
    tmp_r = ((ele_p+PLANET_OF_EJECTION)->radius) * (1.0 + 0.1 * (int)((double)(i - global_n_p - 1) / 100.0));
    tmp_rand = 2.0 * M_PI / 100.0 * (double)(i - global_n_p - 1);
    x_eject[i][1] = tmp_r * cos(EJECTION_CONE_ANGLE);  //破片のx座標.
    x_eject[i][2] = tmp_r * sin(EJECTION_CONE_ANGLE) * cos(tmp_rand);  //破片のy座標.
    x_eject[i][3] = tmp_r * sin(EJECTION_CONE_ANGLE) * sin(tmp_rand);  //破片のz座標.
    //速度.
    tmp_v = ejection_velocity * (0.9 + 0.1 * tmp_r / ((ele_p+PLANET_OF_EJECTION)->radius));
    tmp_theta = EJECTION_CONE_ANGLE * tmp_r / ((ele_p+PLANET_OF_EJECTION)->radius);
    v_eject[i][1] = tmp_v * cos(tmp_theta);  //破片の速度x成分.
    v_eject[i][2] = tmp_v * sin(tmp_theta) * cos(tmp_rand);  //破片の速度y成分.
    v_eject[i][3] = tmp_v * sin(tmp_theta) * sin(tmp_rand);  //破片の速度z成分.



    //z軸周りに回転.
    //位置. 回転させないとouter.
    //Rotation_3D_zaxis(i,x_eject,M_PI);  //180度回転. inner.
    //Rotation_3D_zaxis(i,x_eject,M_PI/2.0);  //90度回転. front.
    Rotation_3D_zaxis(i,x_eject,-M_PI/2.0);  //-90度回転. back.
    //速度. 回転させないとouter.
    //Rotation_3D_zaxis(i,v_eject,M_PI);  //180度回転. inner.
    //Rotation_3D_zaxis(i,v_eject,M_PI/2.0);  //90度回転. front.
    Rotation_3D_zaxis(i,v_eject,-M_PI/2.0);  //-90度回転. back.



    //fprintf(fplog,"%s\tx_eject[%d][1]=%f\tx_eject[%d][2]=%f\tx_eject[%d][3]=%f\n",ele[i].name,i,x_eject[i][1],i,x_eject[i][2],i,x_eject[i][3]);
    //fprintf(fplog,"%s\tv_eject[%d][1]=%f\tv_eject[%d][2]=%f\tv_eject[%d][3]=%f\n",ele[i].name,i,v_eject[i][1],i,v_eject[i][2],i,v_eject[i][3]);


    //////////////////ここまでで、惑星中心、xyは軌道面上、x軸は太陽から惑星の方向/////////////////////


    //x軸が惑星の近日点を向くようにz軸周りに回転.
    //位置.
    tmp_x = x_eject[i][1];
    tmp_y = x_eject[i][2];

    x_eject[i][1] = cos((((ele_p+PLANET_OF_EJECTION)->u) - ((ele_p+PLANET_OF_EJECTION)->ecc)) / (1.0 - ((ele_p+PLANET_OF_EJECTION)->ecc) * cos(((ele_p+PLANET_OF_EJECTION)->u)))) * tmp_x - sin((sqrt(1.0 - ((ele_p+PLANET_OF_EJECTION)->ecc) * ((ele_p+PLANET_OF_EJECTION)->ecc)) * sin(((ele_p+PLANET_OF_EJECTION)->u))) / (1.0 - ((ele_p+PLANET_OF_EJECTION)->ecc) * cos(((ele_p+PLANET_OF_EJECTION)->u)))) * tmp_y;

    x_eject[i][2] = sin((sqrt(1.0 - ((ele_p+PLANET_OF_EJECTION)->ecc) * ((ele_p+PLANET_OF_EJECTION)->ecc)) * sin(((ele_p+PLANET_OF_EJECTION)->u))) / (1.0 - ((ele_p+PLANET_OF_EJECTION)->ecc) * cos(((ele_p+PLANET_OF_EJECTION)->u)))) * tmp_x + cos((((ele_p+PLANET_OF_EJECTION)->u) - ((ele_p+PLANET_OF_EJECTION)->ecc)) / (1.0 - ((ele_p+PLANET_OF_EJECTION)->ecc) * cos(((ele_p+PLANET_OF_EJECTION)->u)))) * tmp_y;

    //速度.
    tmp_x = v_eject[i][1];
    tmp_y = v_eject[i][2];

    v_eject[i][1] = cos((((ele_p+PLANET_OF_EJECTION)->u) - ((ele_p+PLANET_OF_EJECTION)->ecc)) / (1.0 - ((ele_p+PLANET_OF_EJECTION)->ecc) * cos(((ele_p+PLANET_OF_EJECTION)->u)))) * tmp_x - sin((sqrt(1.0 - ((ele_p+PLANET_OF_EJECTION)->ecc) * ((ele_p+PLANET_OF_EJECTION)->ecc)) * sin(((ele_p+PLANET_OF_EJECTION)->u))) / (1.0 - ((ele_p+PLANET_OF_EJECTION)->ecc) * cos(((ele_p+PLANET_OF_EJECTION)->u)))) * tmp_y;

    v_eject[i][2] = sin((sqrt(1.0 - ((ele_p+PLANET_OF_EJECTION)->ecc) * ((ele_p+PLANET_OF_EJECTION)->ecc)) * sin(((ele_p+PLANET_OF_EJECTION)->u))) / (1.0 - ((ele_p+PLANET_OF_EJECTION)->ecc) * cos(((ele_p+PLANET_OF_EJECTION)->u)))) * tmp_x + cos((((ele_p+PLANET_OF_EJECTION)->u) - ((ele_p+PLANET_OF_EJECTION)->ecc)) / (1.0 - ((ele_p+PLANET_OF_EJECTION)->ecc) * cos(((ele_p+PLANET_OF_EJECTION)->u)))) * tmp_y;


    //////////////////ここまでで、惑星中心、xyは軌道面上、x軸は惑星の近日点の方向/////////////////////


    //基準系と平行になるようにオイラー回転.
    //位置.
    Rotation_3D_zaxis(i,x_eject,((ele_p+PLANET_OF_EJECTION)->Omega));
    Rotation_3D_xaxis(i,x_eject,((ele_p+PLANET_OF_EJECTION)->inc));
    Rotation_3D_zaxis(i,x_eject,((ele_p+PLANET_OF_EJECTION)->omega));

    //速度.
    Rotation_3D_zaxis(i,v_eject,((ele_p+PLANET_OF_EJECTION)->Omega));
    Rotation_3D_xaxis(i,v_eject,((ele_p+PLANET_OF_EJECTION)->inc));
    Rotation_3D_zaxis(i,v_eject,((ele_p+PLANET_OF_EJECTION)->omega));


    //////////////////ここまでで、惑星中心に平行移動した基準系/////////////////////


    //太陽中心の基準系に平行移動.
    for(k=1;k<=3;++k){
      x_0[i][k] = x_0[PLANET_OF_EJECTION][k] + x_eject[i][k];
      v_0[i][k] = v_0[PLANET_OF_EJECTION][k] + v_eject[i][k];
    }


    //////////////////ここまでで、基準系/////////////////////


    r_0[i] = RadiusFromCenter(i,x_0);  //中心星からの距離.
    v2_0[i] = SquareOfVelocity(i,v_0);  //速度の2乗.
    r_dot_v[i] = InnerProduct(i,x_0,v_0);  //r_i,v_iの内積.

  }

  return;
}
#endif  /*EJECTION*/


/*軌道要素計算*/
void Calculate_OrbitalElements(int i, CONST double x_c[][4], CONST double v_c[][4], struct orbital_elements *ele_p, CONST double r_c[], CONST double v2_c[], CONST double r_dot_v[]
#if FRAGMENTATION
			       , double t_dyn
			       , CONST struct fragmentation *frag_p
#endif
			       ){


  double m_i;

#if FRAGMENTATION
  m_i = MassDepletion(i,((ele_p+i)->mass),t_dyn,frag_p);
#else
  m_i = ((ele_p+i)->mass);
#endif


#if INDIRECT_TERM

#if !defined(G) && !defined(M_0)
  double mu = 1.0 + m_i;
#else
  double mu = G * (M_0 + m_i);
#endif

#else

#if !defined(G) && !defined(M_0)
  double mu = 1.0;
#else
  double mu = G * M_0;
#endif

#endif

  int k;
  double esin_u;
  double ecos_u;
  double sin_inc;
  double cos_inc;
  double sin_omega;
  double cos_omega;
  double sin_Omega;
  double cos_Omega;
  double radian;
  static double P[N_p+N_tr+1][4]={}, Q[N_p+N_tr+1][4]={};


  ((ele_p+i)->axis) = 1.0 / (2.0 / r_c[i] - v2_c[i] / mu);


  if(isnan((ele_p+i)->axis)){
    fprintf(fplog,"i=%d\taxis is nan.\n",i);
  }


  if(((ele_p+i)->axis)<0.0){
    fprintf(fplog,"i=%d\taxis=%e < 0 hyperbola orbit.\n",i,((ele_p+i)->axis));
    /*
    ((ele_p+i)->ecc) = NAN;
    ((ele_p+i)->u) = NAN;
    ((ele_p+i)->inc) = NAN;
    ((ele_p+i)->omega) = NAN;
    ((ele_p+i)->Omega) = NAN;
    ((ele_p+i)->r_h) = NAN;
    */
  }


  ((ele_p+i)->ecc) = sqrt((1.0 - r_c[i] / ((ele_p+i)->axis)) * (1.0 - r_c[i] / ((ele_p+i)->axis)) + r_dot_v[i] * r_dot_v[i] / mu / ((ele_p+i)->axis));


  if(isnan((ele_p+i)->ecc)){
    fprintf(fplog,"i=%d\tecc is nan.\n",i);
  }


  if(((ele_p+i)->ecc)==0.0){
    ((ele_p+i)->u) = 0.0;
  }else{
    esin_u = r_dot_v[i] / sqrt(mu * ((ele_p+i)->axis));
    ecos_u = 1.0 - r_c[i] / ((ele_p+i)->axis);
    radian = atan2(esin_u,ecos_u);
    if(radian<0.0){
      ((ele_p+i)->u) = radian + 2.0 * M_PI;
    }else{
      ((ele_p+i)->u) = radian;
    }
  }


  for(k=1;k<=3;++k){

    P[i][k] = x_c[i][k] * cos(((ele_p+i)->u)) / r_c[i] - sqrt(((ele_p+i)->axis) / mu) * v_c[i][k] * sin(((ele_p+i)->u));
    Q[i][k] = (x_c[i][k] * sin(((ele_p+i)->u)) / r_c[i] + sqrt(((ele_p+i)->axis) / mu) * v_c[i][k] * (cos(((ele_p+i)->u)) - ((ele_p+i)->ecc))) / sqrt(1.0 - ((ele_p+i)->ecc));
  }


  if(isnan(P[i][1])||isnan(P[i][2])||isnan(P[i][3])){
    fprintf(fplog,"i=%d\tP is nan.\t[1]=%f\t[2]=%f\t[3]=%f\n",i,P[i][1],P[i][2],P[i][3]);
  }
  if(isnan(Q[i][1])||isnan(Q[i][2])||isnan(Q[i][3])){
    fprintf(fplog,"i=%d\tQ is nan.\t[1]=%f\t[2]=%f\t[3]=%f\n",i,Q[i][1],Q[i][2],Q[i][3]);
  }


  sin_inc = sqrt(P[i][3] * P[i][3] + Q[i][3] * Q[i][3]);
  cos_inc = P[i][1] * Q[i][2] - P[i][2] * Q[i][1];
  radian = atan2(sin_inc,cos_inc);


  if(radian<0.0){
    ((ele_p+i)->inc) = radian + 2.0 * M_PI;
  }else{
    ((ele_p+i)->inc) = radian;
  }


  sin_omega = P[i][3] / sin_inc;
  cos_omega = Q[i][3] / sin_inc;
  radian = atan2(sin_omega,cos_omega);
  if(radian<0.0){
    ((ele_p+i)->omega) = radian + 2.0*M_PI;
  }else{
    ((ele_p+i)->omega) = radian;
  }


  sin_Omega = (P[i][2] * Q[i][3] - Q[i][2] * P[i][3]) / sin_inc;
  cos_Omega = (P[i][1] * Q[i][3] - Q[i][1] * P[i][3]) / sin_inc;
  radian = atan2(sin_Omega,cos_Omega);
  if(radian<0.0){
    ((ele_p+i)->Omega) = radian + 2.0 * M_PI;
  }else{
    ((ele_p+i)->Omega) = radian;
  }


  if(sin_inc==0.0){
    ((ele_p+i)->omega) = 0.0;
    ((ele_p+i)->Omega) = 0.0;
  }


#ifndef M_0
  ((ele_p+i)->r_h) = ((ele_p+i)->axis) * cbrt(m_i / 3.0);
#else
  ((ele_p+i)->r_h) = ((ele_p+i)->axis) * cbrt(m_i / M_0 / 3.0);
#endif


  if(isnan((ele_p+i)->inc)){
    fprintf(fplog,"i=%d\tinc is nan.\n",i);
  }
  if(isnan((ele_p+i)->omega)){
    fprintf(fplog,"i=%d\tomega is nan.\n",i);
  }
  if(isnan((ele_p+i)->Omega)){
    fprintf(fplog,"i=%d\tOmega is nan.\n",i);
  }


  return;
}


void Calculate_RMS(CONST struct orbital_elements *ele_p, double *ecc_p_rms, double *ecc_tr_rms, double *inc_p_rms, double *inc_tr_rms){

  int i;
  double ecc_2, ecc_2_mean, inc_2, inc_2_mean;

  ecc_2 = 0.0;
  inc_2 = 0.0;
  for(i=1;i<=global_n_p;++i){
    ecc_2 += ((ele_p+i)->ecc) * ((ele_p+i)->ecc);
    inc_2 += ((ele_p+i)->inc) * ((ele_p+i)->inc);
  }
  ecc_2_mean = ecc_2 / ((double)global_n_p);
  inc_2_mean = inc_2 / ((double)global_n_p);
  *ecc_p_rms = sqrt(ecc_2_mean);
  *inc_p_rms = sqrt(inc_2_mean);


  ecc_2 = 0.0;
  inc_2 = 0.0;
  for(i=global_n_p+1;i<=global_n;++i){
    ecc_2 += ((ele_p+i)->ecc) * ((ele_p+i)->ecc);
    inc_2 += ((ele_p+i)->inc) * ((ele_p+i)->inc);
  }
  ecc_2_mean = ecc_2 / ((double)(global_n-global_n_p));
  inc_2_mean = inc_2 / ((double)(global_n-global_n_p));
  *ecc_tr_rms = sqrt(ecc_2_mean);
  *inc_tr_rms = sqrt(inc_2_mean);

  return;
}


/*初期位置、速度計算*/
void InitialCondition(int i, double x_0[][4], double v_0[][4], double v2_0[], double r_dot_v[], double r_0[], CONST struct orbital_elements *ele_p){

#if INDIRECT_TERM

#if !defined(G) && !defined(M_0)
  double mu = 1.0 + ((ele_p+i)->mass);
#else
  double mu = G * (M_0 + ((ele_p+i)->mass));
#endif

#else

#if !defined(G) && !defined(M_0)
  double mu = 1.0;
#else
  double mu = G * M_0;
#endif

#endif

  int k;
  static double P[4]={},Q[4]={};


  for(k=1;k<=3;k++){
    P[k] = Calculate_P(i,k,ele_p);
    Q[k] = Calculate_Q(i,k,ele_p);

    x_0[i][k] = ((ele_p+i)->axis) * P[k] * (cos(((ele_p+i)->u)) - ((ele_p+i)->ecc)) + ((ele_p+i)->axis) * sqrt(1.0 - ((ele_p+i)->ecc) * ((ele_p+i)->ecc)) * Q[k] * sin(((ele_p+i)->u));
  }
  //fprintf(fplog,"x=%f\ty=%f\tz=%f\n",x_0[i][1],x_0[i][2],x_0[i][3]);

  r_0[i] = RadiusFromCenter(i,x_0);  //中心星からの距離.


  for(k=1;k<=3;++k){
    v_0[i][k] = sqrt(mu / ((ele_p+i)->axis)) / r_0[i] * (- ((ele_p+i)->axis) * P[k] * sin(((ele_p+i)->u)) + ((ele_p+i)->axis) * sqrt(1.0 - ((ele_p+i)->ecc) * ((ele_p+i)->ecc)) * Q[k] * cos(((ele_p+i)->u)));
  }

  r_dot_v[i] = InnerProduct(i,x_0,v_0);  //r_i,v_iの内積.
  v2_0[i] = SquareOfVelocity(i,v_0);  //速度の2乗.
  //fprintf(fplog,"vx=%f\tvy=%f\tvz=%f\n",v_0[i][1],v_0[i][2],v_0[i][3]);

  return;
}
