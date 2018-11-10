#include "hybrid.h"
#include "func.h"


/*エネルギー計算*/
double Calculate_Energy(CONST struct orbital_elements *ele_p, CONST double x_c[][4]
#if INDIRECT_TERM
			, CONST double v_c[][4], CONST double v_G[]
#else
			, CONST double v2_c[]
#endif
			, CONST double r_c[]
#if FRAGMENTATION
			, double t_dyn
			, CONST struct fragmentation *frag_p
#endif
			){

  int i,j;
  double abs_rij;
  double E[global_n+1];
  double E_tot;
  double m_i,m_j;

#if INDIRECT_TERM
#ifndef M_0
  E_tot = 0.5 * (v_G[1]*v_G[1] + v_G[2]*v_G[2] + v_G[3]*v_G[3]);
#else
  E_tot = 0.5 * M_0 * (v_G[1]*v_G[1] + v_G[2]*v_G[2] + v_G[3]*v_G[3]);
#endif

#else
  E_tot = 0.0;
#endif

  for(i=1;i<=global_n;++i){

#if FRAGMENTATION
    m_i = MassDepletion(i,((ele_p+i)->mass),t_dyn,frag_p);
#else
    m_i = ((ele_p+i)->mass);
#endif

#if INDIRECT_TERM
    E[i] = 0.5 * m_i * ((v_c[i][1]-v_G[1])*(v_c[i][1]-v_G[1]) + (v_c[i][2]-v_G[2])*(v_c[i][2]-v_G[2]) + (v_c[i][3]-v_G[3])*(v_c[i][3]-v_G[3]));
#else
    E[i] = 0.5 * m_i * v2_c[i];
#endif

    for(j=i+1;j<=global_n;++j){

#if FRAGMENTATION
      m_j = MassDepletion(j,((ele_p+j)->mass),t_dyn,frag_p);
#else
      m_j = ((ele_p+j)->mass);
#endif

      abs_rij = RelativeDistance(i,j,x_c); //絶対値.

#ifndef G
      E[i] += - m_i * m_j / abs_rij;  //エネルギー計算.
#else
      E[i] += - G * m_i * m_j / abs_rij;  //エネルギー計算.
#endif

    }  //j loop

#if !defined(G) && !defined(M_0)
    E_tot += - m_i / r_c[i] + E[i];
#else
    E_tot += - G * M_0 * m_i / r_c[i] + E[i];
#endif

  }  //i loop
  return E_tot;
}


/*角運動量*/
double AngularMomentum(int i, CONST struct orbital_elements *ele_p, CONST double x_0[][4], CONST double v_0[][4]
#if FRAGMENTATION
		       , double t_dyn
		       , CONST struct fragmentation *frag_p
#endif
		       ){
  int k;
  double L[global_n+1][4];
  double L_tot_0[4];
  double abs_L_0;
  double m_i;


  for(i=1;i<=global_n_p;++i){

#if FRAGMENTATION
    m_i = MassDepletion(i,((ele_p+i)->mass),t_dyn,frag_p);
#else
    m_i = ((ele_p+i)->mass);
#endif

    L[i][1] = m_i * (x_0[i][2]*v_0[i][3] - x_0[i][3]*v_0[i][2]);
    L[i][2] = m_i * (x_0[i][3]*v_0[i][1] - x_0[i][1]*v_0[i][3]);
    L[i][3] = m_i * (x_0[i][1]*v_0[i][2] - x_0[i][2]*v_0[i][1]);
    //fprintf(fplog,"i=%d\t(ele_p+i)->mass=%e\n",i,(ele_p+i)->mass);
  }

  for(k=1;k<=3;++k){
    L_tot_0[k] = 0.0;
    for(i=1;i<=global_n;++i){
      L_tot_0[k] += L[i][k];
    }
  }

  abs_L_0 = 0.0;

  for(k=1;k<=3;++k){
    abs_L_0 += L_tot_0[k]*L_tot_0[k];
  }

  abs_L_0 = sqrt(abs_L_0);

  return abs_L_0;
}
