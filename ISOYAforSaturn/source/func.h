#ifndef INCLUDED_FUNC_H  //include-guard
#define INCLUDED_FUNC_H  //include-guard


#include "hybrid.h"

double AngularMomentum(CONST struct orbital_elements *ele_p, CONST double x_0[][4], CONST double v_0[][4]
#if FRAGMENTATION
		       , double t_dyn
		       , CONST struct fragmentation *frag_p
#endif
		       );

double Timestep_i_0(int i, CONST double a_0[][4], CONST double adot_0[][4]);

double Timestep_i_sys(int i_sys, CONST double a[][4], CONST double adot[][4], CONST double adot2_dt2[][4], CONST double adot3_dt3[][4], CONST double dt_[]);

double All_Acceleration(int i, int k, CONST struct orbital_elements *ele_p, CONST double x_0[][4], CONST double r_0[], CONST double abs_r[]
#if FRAGMENTATION
			, double t_dyn
			, CONST struct fragmentation *frag_p
#endif
			);

double All_dAcceleration(int i, int k, CONST struct orbital_elements *ele_p, CONST double x_0[][4], CONST double v_0[][4], CONST double r_dot_v[], CONST double r_dot_v_ij[], CONST double r_0[], CONST double abs_r[]
#if FRAGMENTATION
			 , double t_dyn
			 , CONST struct fragmentation *frag_p
#endif
			 );

void Calculate_OrbitalElements(int i, CONST double x_c[][4], CONST double v_c[][4], struct orbital_elements *ele_p, CONST double r_c[], CONST double v2_c[], CONST double r_dot_v[]
#if FRAGMENTATION
			       , double t_dyn
			       , CONST struct fragmentation *frag_p
#endif
			       );

void Calculate_RMS(CONST struct orbital_elements *ele_p, double *ecc_p_rms, double *ecc_tr_rms, double *inc_p_rms, double *inc_tr_rms);

double Calculate_Energy(CONST struct orbital_elements *ele_p, CONST double x_c[][4]
#if INDIRECT_TERM
			, CONST double v_c[][4]
			, CONST double v_G[]
#else
			, CONST double v2_c[]
#endif
			, CONST double r_c[]
#if FRAGMENTATION
			, double t_dyn
			, CONST struct fragmentation *frag_p
#endif
			);

double MutualHillRadius_to_SemimajorAxis(double ratio);

//double MutualHillRadius_to_SemimajorAxis(double m_1,double m_2,double ratio);

//double IsolationMass(double axis,double ratio,double sigma_0,double alpha);

//double Iteration_of_InitialAxis(double axis_1,double m_1,double ratio,double sigma_0, double alpha);

double Escape_Velocity(double mass_p, double r);


void InitialOrbitalElements_Planet(int i, struct orbital_elements *ele_p);

void InitialOrbitalElements_Tracer(int i, double x_0[][4], struct orbital_elements *ele_p);

//void EjectionOfTracerFromPlanet(double x_0[][4], double v_0[][4], double v2_0[], double r_dot_v[], double r_0[], CONST struct orbital_elements *ele_p);

void EjectionOfTracerFromPlanet(double x_0[][4], double v_0[][4], double v2_0[], double r_dot_v[], double r_0[], struct orbital_elements *ele_p, double x_col, double y_col, double z_col,double mass_p, double r_h_p);

void InitialCondition(int i, double x_0[][4], double v_0[][4], double v2_0[], double r_dot_v[], double r_0[], CONST struct orbital_elements *ele_p);

void Corrector_sys(int n_ite, int i_sys, CONST struct orbital_elements *ele_p, CONST double x_p[][4], CONST double v_p[][4], CONST double r_p[], double x_c[][4], double v_c[][4], double r_c[], double v2_c[], double r_dot_v[], CONST double a_0[][4], CONST double adot_0[][4], double a[][4], double adot[][4], double adot2_dt2[][4], double adot3_dt3[][4], CONST double dt_[]
#if FRAGMENTATION
		   , double t_dyn
		   , CONST struct fragmentation *frag_p
#endif
		   );

//void Iteration_sys(int i_sys, CONST struct orbital_elements *ele_p,CONST double x_p[][4],CONST double v_p[][4],double x_c[][4],double v_c[][4],double r_c[],double v2_c[],double r_dot_v[],CONST double a_0[][4],CONST double adot_0[][4],double a[][4],double adot[][4],double adot2_dt2[][4],double adot3_dt3[][4],CONST double dt_[]);

double rand_func();

void HeapSort(CONST int NUM_DATA, CONST double t[], CONST double dt[], CONST double t_tmp, int index[]);

void DownHeap(double a[], int index[], int leaf, int root);

#if COLLISION
bool Collision_Judgement(int i_sys, CONST struct orbital_elements *ele_p, CONST double x_p[][4], double abs_r[], int *i_col, int *j_col);

void Energy_Correction(int i_col, int j_col, CONST double x_0[][4], CONST double v_0[][4], CONST struct orbital_elements *ele_p, double *dE_heat, double *dE_grav, double *dE_c, double *v_imp
#if FRAGMENTATION
		       , double t_dyn
		       , CONST struct fragmentation *frag_p
#endif
		       );

#if COALESCENCE
void Coalescence(int i_col, int j_col, double x_0[][4], double v_0[][4], struct orbital_elements *ele_p
#if FRAGMENTATION
		 , double t_dyn
		 , struct fragmentation *frag_p
#endif
		 );
#endif  /*COALESCENCE*/
#endif  /*COLLISION*/

#if FRAGMENTATION
void NeighborSearch(int i, double t_dyn, CONST struct orbital_elements *ele_p, struct fragmentation *frag_p, CONST double x_0[][4]);

double SquareRandomVelocity(int i, int j, CONST struct orbital_elements *ele_p);

void MassFlux(int i, CONST struct orbital_elements *ele_p, struct fragmentation *frag_p, CONST struct parameter *para_p);

double Depletion_Time(int i, CONST struct fragmentation *frag_p);

double MassDepletion(int i, double mass, double t_dyn, CONST struct fragmentation *frag_p);

double s_1_FRAG_trapezoid(int n, double dx, double ini, CONST struct parameter *para_p);

double s_2_FRAG_trapezoid(int n, double dx, double ini, CONST struct parameter *para_p);

double s_3_FRAG_trapezoid(int n, double dx, double ini, CONST struct parameter *para_p);

double s_1_FRAG_integrand(double x, CONST struct parameter *para_p);

double s_2_FRAG_integrand(double x, CONST struct parameter *para_p);

double s_3_FRAG_integrand(double x, CONST struct parameter *para_p);

double s_1_FRAG(struct parameter *para_p);

double s_2_FRAG(struct parameter *para_p);

double s_3_FRAG(struct parameter *para_p);
#endif

#if EJECTION
void Rotation_3D_xaxis(int i, double x_eject[][4], double theta);

void Rotation_3D_yaxis(int i, double x_eject[][4], double theta);

void Rotation_3D_zaxis(int i, double x_eject[][4], double theta);
#endif

#if INDIRECT_TERM
void CenterOfGravity(CONST double x_0[][4], CONST double v_0[][4], double x_G[], double v_G[], CONST struct orbital_elements *ele_p
#if FRAGMENTATION
		     , double t_dyn
		     , CONST struct fragmentation *frag_p
#endif
		     );
#endif

#if EXECUTION_TIME
void Sort_Exetime(struct timeval realtime_start_main, struct timeval realtime_end_main);
#endif

#endif //include-guard


//わざとインクルードガードから外し、static付きで定義
//inline関数

/*x_i,v_iの内積*/
static inline double InnerProduct(int i, CONST double x[][4], CONST double v[][4]){
  return x[i][1]*v[i][1] + x[i][2]*v[i][2] + x[i][3]*v[i][3];
}

/*中心星からの距離*/
static inline double RadiusFromCenter(int i, CONST double x[][4]){
  return sqrt(x[i][1]*x[i][1] + x[i][2]*x[i][2] + x[i][3]*x[i][3]);
}

/*速度の2乗*/
static inline double SquareOfVelocity(int i, CONST double v[][4]){
  return v[i][1]*v[i][1] + v[i][2]*v[i][2] + v[i][3]*v[i][3];
}

/*相対距離*/
static inline double RelativeDistance(int i, int j, CONST double x[][4]){
  return sqrt((x[j][1] - x[i][1])*(x[j][1] - x[i][1]) + (x[j][2] - x[i][2])*(x[j][2] - x[i][2]) + (x[j][3] - x[i][3])*(x[j][3] - x[i][3]));
}

/*相対速度の2乗*/
static inline double SquareOfRelativeVelocity(int i, int j, CONST double v[][4]){
  return (v[j][1] - v[i][1])*(v[j][1] - v[i][1]) + (v[j][2] - v[i][2])*(v[j][2] - v[i][2]) + (v[j][3] - v[i][3])*(v[j][3] - v[i][3]);
}

/*x_ij, v_ijの内積*/
static inline double RelativeInnerProduct(int i, int j, CONST double x[][4], CONST double v[][4]){
  return (x[j][1] - x[i][1])*(v[j][1] - v[i][1]) + (x[j][2] - x[i][2])*(v[j][2] - v[i][2]) + (x[j][3] - x[i][3])*(v[j][3] - v[i][3]);
}

/*Swap*/
static inline void Swap_double(double *a, double *b){
  double tmp;
  tmp = (*a);
  (*a) = (*b);
  (*b) = tmp;
  return;
}


static inline void Swap_int(int *a, int *b){
  int tmp;
  tmp = (*a);
  (*a) = (*b);
  (*b) = tmp;
  return;
}


static inline int Min_int(int a, int b){
  if(a < b){
    return a;
  }else{
    return b;
  }
}


static inline int Max_int(int a, int b){
  if(a > b){
    return a;
  }else{
    return b;
  }
}


static inline double Cal_time(struct timeval time1, struct timeval time2){
  return (double)(time2.tv_sec - time1.tv_sec) + (double)(time2.tv_usec - time1.tv_usec)*1.0E-6;
}


static inline void Predictor(int i, CONST double x_0[][4], CONST double v_0[][4], CONST double a_0[][4], CONST double adot_0[][4], double x_p[][4], double v_p[][4], double r_p[], double v2_p[], double r_dot_v[], CONST double Dt[]){

  int k;
  double dt = Dt[i];

  for(k=1;k<=3;++k){
    //x_p[i][k] = x_0[i][k] + v_0[i][k]*Dt[i] + a_0[i][k]*Dt[i]*Dt[i]/2.0 + adot_0[i][k]*Dt[i]*Dt[i]*Dt[i]/6.0;
    //v_p[i][k] = v_0[i][k] + a_0[i][k]*Dt[i] + adot_0[i][k]*Dt[i]*Dt[i]/2.0;

    x_p[i][k] = x_0[i][k] + dt * (v_0[i][k] + dt * 0.5 * (a_0[i][k] + adot_0[i][k] * dt * INV_3));
    v_p[i][k] = v_0[i][k] + dt * (a_0[i][k] + adot_0[i][k] * dt * 0.5);
  }


  r_p[i] = RadiusFromCenter(i,x_p);  //中心星からの距離.
  v2_p[i] = SquareOfVelocity(i,v_p); //速度の2乗.
  r_dot_v[i] = InnerProduct(i,x_p,v_p);  //r_i,v_iの内積.

  return;
}
