#include "hybrid.h"
#include "func.h"


double rand_func(){
  return sfmt_genrand_real2(&sfmt);  //generates a random number on [0,1)-real-interval.
}


void Rotation_3D_xaxis(int i, double x_eject[][4], double theta){
  double tmp_y = x_eject[i][2];
  double tmp_z = x_eject[i][3];
  x_eject[i][2] = cos(theta) * tmp_y - sin(theta) * tmp_z;
  x_eject[i][3] = sin(theta) * tmp_y + cos(theta) * tmp_z;
  return;
}


void Rotation_3D_yaxis(int i, double x_eject[][4], double theta){
  double tmp_x = x_eject[i][1];
  double tmp_z = x_eject[i][3];
  x_eject[i][1] = cos(theta) * tmp_x + sin(theta) * tmp_z;
  x_eject[i][3] = - sin(theta) * tmp_x + cos(theta) * tmp_z;
  return;
}


void Rotation_3D_zaxis(int i,double x_eject[][4],double theta){
  double tmp_x = x_eject[i][1];
  double tmp_y = x_eject[i][2];
  x_eject[i][1] = cos(theta) * tmp_x - sin(theta) * tmp_y;
  x_eject[i][2] = sin(theta) * tmp_x + cos(theta) * tmp_y;
  return;
}


void CenterOfGravity(CONST double x_0[][4], CONST double v_0[][4], double x_G[], double v_G[], CONST struct orbital_elements *ele_p
#if FRAGMENTATION
			 , double t_dyn
			 , CONST struct fragmentation *frag_p
#endif
		     ){
  int i, k;
  double M;

#ifndef M_0
  M = 1.0;
#else
  M = M_0;
#endif
  for(i=1;i<=
#if INTERACTION_ALL
	global_n
#else
	global_n_p
#endif
	;++i){

#if FRAGMENTATION
    M += MassDepletion(i,((ele_p+i)->mass),t_dyn,frag_p);
#else
    M += ((ele_p+i)->mass);
#endif
  }

  for(k=1;k<=3;++k){
    x_G[k] = 0.0;
    v_G[k] = 0.0;
    for(i=1;i<=
#if INTERACTION_ALL
	global_n
#else
	global_n_p
#endif
	  ;++i){

#if FRAGMENTATION
      x_G[k] += MassDepletion(i,((ele_p+i)->mass),t_dyn,frag_p) * x_0[i][k];
      v_G[k] += MassDepletion(i,((ele_p+i)->mass),t_dyn,frag_p) * v_0[i][k];
#else
      x_G[k] += ((ele_p+i)->mass) * x_0[i][k];
      v_G[k] += ((ele_p+i)->mass) * v_0[i][k];
#endif
    }
    x_G[k] = x_G[k] / M;
    v_G[k] = v_G[k] / M;
  }

  return;
}


double MutualHillRadius_to_SemimajorAxis(double ratio){
  return (1.0 / ratio + 0.5 * cbrt(2.0 * PLANET_MASS / 3.0)) / (1.0 / ratio - 0.5 * cbrt(2.0 * PLANET_MASS / 3.0));
}


#if EXECUTION_TIME
void Sort_Exetime(struct timeval realtime_start_main, struct timeval realtime_end_main){

  int i,j;
  double exetime_main = Cal_time(realtime_start_main,realtime_end_main);
  int exetime_num[7]={0,1,2,3,4,5,6};

  double exetime_array[7]={
    exetime.Energy[0],
    exetime.Orbital_Elements[0],
    exetime.Predictor[0],
    exetime.Corrector[0],
    exetime.Iteration[0],
    exetime.Collision_Judgement[0],
    exetime.Fragmentation[0]
  };


#if EXECUTION_TIME_FUNC
  double exetime_others = 0.0;

  char exetime_name[7][30]={
    "Energy\t\t\t",
    "Orbital_Elements\t",
    "Predictor\t\t",
    "Corrector\t\t",
    "Iteration\t\t",
    "Collision_Judgement\t",
    "Fragmentation\t\t"
  };
#endif

  for(i=0;i<7;++i){
    for(j=i+1;j<7;++j){
      if(exetime_array[i] < exetime_array[j]){
	Swap_int(&exetime_num[i],&exetime_num[j]);
	Swap_double(&exetime_array[i],&exetime_array[j]);
      }
    }
  }

  fprintf(fplog,"Execution Time\t(total\t= %e [s])\n",exetime_main);

#if EXECUTION_TIME_FUNC
  for(i=0;i<7;++i){
    fprintf(fplog,"%s= %e [s]\t%5.2f [%%]\n",exetime_name[exetime_num[i]],exetime_array[i],exetime_array[i]/exetime_main*100.0);
    exetime_others += exetime_array[i];
  }
  exetime_others = exetime_main - exetime_others;
  fprintf(fplog,"Others\t\t\t= %e [s]\t%5.2f [%%]\n",exetime_others,exetime_others/exetime_main*100.0);
#endif

  return;
}
#endif
