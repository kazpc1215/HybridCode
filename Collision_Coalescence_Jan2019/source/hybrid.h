#ifndef INCLUDED_hybrid_H  //include-guard
#define INCLUDED_hybrid_H  //include-guard


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <time.h>
#include <stdbool.h>
#include <mpi.h>
#include "SFMTdir/SFMT.h"


#ifdef _OPENMP
#include <omp.h>
#define OMP_SCHEDULE static
#endif


#if __GNUC__ == 7
#define CONST const
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define CONST
#define ALWAYS_INLINE
#endif


#define DIRECTORY ../data/N25_t1E8_dtlog_10RHM_1MMSN_Miso_ecc1E-2/  //ディレクトリ.
#define SUBDIRECTORY rand  //子ディレクトリ. rand%02d

#define STR_(str) #str
#define STR(str) STR_(str)

#define INV_3 0.33333333333333333333

#ifndef EXTERN
#define EXTERN extern
#endif


//////////////////////////////////////////////////
#define N_tr 0  //初期のトレーサーの数.
#define N_p 25  //初期の原始惑星の数.
#define ECC_RATIO 10.0  //ecc=0.01の何倍か. inc=ecc/2.
#define STEP_INTERVAL 5.0E6  //何ステップごとに標準出力するか.
//#define BREAK_TIME 100.0  //4h = 14400sec, 12h = 43200sec.
//#define BREAK_TIME 14100.0  //4h = 14400sec, 12h = 43200sec.
//#define BREAK_TIME 42900.0  //4h = 14400sec, 12h = 43200sec.
#define BREAK_TIME 864000.0  //10days = 864000sec.

#define RAYLEIGH_DISTRIBUTION true  //離心率や軌道傾斜角の分布 true : Rayleigh, false : v_relが軌道長半径によらず一定.

#define FRAGMENTATION false  //破壊 近傍粒子探索と質量フラックス計算.
#define COLLISION true  //衝突.
#if COLLISION
#define COALESCENCE true  //衝突後に合体.
#else
//#define SOFTENING true  //衝突しないようソフトニング.
#endif


#define PLANET_INNER_AXIS 0.4  //[AU].
#define DELTA_HILL 10.0  //相互ヒル半径の何倍か.
#define SIGMA_SOLID 7.1  //[g/cm^2]. //固体面密度（MMSN->7.1）. 1.5MMSN=10.65, 2MMSN=14.2, 2.5MMSN=17.75, 3MMSN=21.3.
#define ALPHA_SOLID 1.5  //面密度分布での軌道長半径のべき（MMSN->1.5）



EXTERN int global_n;  //グローバル変数.
EXTERN int global_n_p;
EXTERN int global_myid;
EXTERN sfmt_t sfmt;
EXTERN FILE *fplog;

//////////////////////////////////////////////////


//////////////////////////////////////////////////
#define ENERGY_FILE true  //エネルギー計算&ファイル作成.
#define ORBITALELEMENTS_FILE true  //軌道要素計算&ファイル作成.
#define POSI_VELO_FILE true  //位置速度ファイル作成.
#define COLLISION_FILE true  //衝突直前の位置速度ファイル作成.
#define TRACERLIST_FILE false  //トレーサーリストのファイル作成.
#define EXECUTION_TIME true  //mainの実行時間測定.
#define EXECUTION_TIME_FUNC false  //mainかつ関数ごとの実行時間測定.
#if EXECUTION_TIME
#include <sys/time.h>
#include <sys/resource.h>
#endif
//////////////////////////////////////////////////


//////////////////////////////////////////////////
#define INTERACTION_ALL false  //全粒子同士の重力相互作用.
#define INTERACTION_PLANET_TRACER true  //惑星とトレーサー間の相互作用.
#define INTERACTION_TEST_PARTICLE false  //トレーサーをテスト粒子として扱う.
#define INDIRECT_TERM true  //中心星が動く効果を補正.
#define EJECTION false  //初期に破片（トレーサー）を放出する.
#define ORBITING_SMALL_PARTICLE false  //初期に微惑星をケプラー運動させておく.
#define ELIMINATE_PARTICLE false  //太陽に飲みこまれるか系外へ出て行くかで粒子を消す.
//////////////////////////////////////////////////


//////////////////////////////////////////////////
//#define G 1.0  //重力定数.
//#define M_0 1.0  //主星の質量.
#if SOFTENING
#define EPSILON 5.21495378928615e-05   //ソフトニングパラメーター.
#endif
#define ETA 1.0E-2  //刻み幅調整.
#define ITE_MAX 2  //イテレーション回数（修正子計算の回数はITE_MAX+1）.
//////////////////////////////////////////////////


#if ELIMINATE_PARTICLE
//////////////////////////////////////////////////
#define SOLAR_RADIUS 0.00465040106951872  //[AU] 6.957E10/1.496E13.
#define SOLAR_SYSTEM_LIMIT 100.0  //[AU]
//////////////////////////////////////////////////
#endif


//////////////////////////////////////////////////
#define PLANET_MASS 3.0E-6  //地球質量M_E.
#define PLANET_AXIS 1.0  //[AU].
#define PLANET_ECC (0.01*ECC_RATIO)
#define PLANET_INC (PLANET_ECC*0.5)
#define PLANET_DENSITY 3.0  //[g/cc].
/*
Earth Mean Orbital Elements (J2000)
Semimajor axis (AU)                  1.00000011
Orbital eccentricity                 0.01671022
Orbital inclination (deg)            0.00005
Longitude of ascending node (deg)  -11.26064
Longitude of perihelion (deg)      102.94719
Mean Longitude (deg)               100.46435
*/
//////////////////////////////////////////////////



#if N_tr != 0
//////////////////////////////////////////////////
#define M_TOT (3.0E-7*N_p)  //0.1M_E * N_p  //トレーサーの総質量.

#if EJECTION
#define PLANET_OF_EJECTION 1
#define EJECTION_CONE_ANGLE M_PI/180.0*30.0  //30度.
#endif

#if ORBITING_SMALL_PARTICLE
#define ECC_RMS (0.01*ECC_RATIO)  //トレーサーの離心率の二乗平均平方根.  //Rayleigh分布.
#define INC_RMS (ECC_RMS*0.5)  //トレーサーの軌道傾斜角の二乗平均平方根.  //Rayleigh分布.
#define DELTA_HILL 10.0  //惑星を「相互」ヒル半径の何倍離すか（軌道長半径）.
#define SEPARATE_HILL 5.0  //初期に惑星とトレーサーをヒル半径の何倍以上離すか（相対距離）.
#endif

#if FRAGMENTATION
/* t_fragcheck : 初項 DT_FRAGCHECK，公比 GEOMETRIC_RATIO_FRAG の等比数列 */
#define DT_FRAGCHECK 2.0*M_PI*0.1  //0.1yr
#define GEOMETRIC_RATIO_FRAG pow(10.0,1.0/8.0) //10**(1/8)
#define DELTA_R 0.01  //Hill 近傍粒子探索用.
#define DELTA_THETA 0.125*M_PI  //近傍粒子探索用.
//#define DELTA_THETA 1.0*M_PI  //近傍粒子探索用.
#define NEIGHBOR_MAX 200  //近傍粒子リスト配列の最大値.
#define RHO 3.0  // [g/cc]  微惑星の物質密度.
#define EPSILON_FRAG 0.2
#define B_FRAG (5.0/3.0)
#define Q_0_FRAG 9.5E8 // [erg/g]  Q_D = Q_0*(rho/3[g/cc])^0.55*(m/10^21[g])^p
#define P_FRAG 0.453
#define XI 0.01 //面密度減少タイムスケールを自身のXI倍間隔で更新する.
#define M_MAX 5.00E-15  //最大微惑星質量. 1E19 g ~10kmサイズ.
//#define M_MAX 5.00E-18  //最大微惑星質量. 1E16 g ~1kmサイズ.
#endif  /*FRAGMENTATION*/
//////////////////////////////////////////////////
#endif  /*N_tr != 0*/


//////////////////////////////////////////////////
#define T_MAX (2.0*M_PI*1.0E8)  //10^8yr 全計算時間.
#define DT_LOG true  //true: t_eneをlogでとる. false: t_eneをlinearでとる.

/* linear では 初項 DT_ENE，公差 DT_ENE の等差数列 */
#define DT_ENE (2.0*M_PI*1.0E-1)  //dt_ene = 0.1yr

#if DT_LOG
/* log では 初項 DT_ENE，公比 GEOMETRIC_RATIO の等比数列 */
#define GEOMETRIC_RATIO pow(10.0,1.0/8.0) //10**(1/8)
#define GEOMETRIC_RATIO_LONGTERM pow(10.0,1.0/128.0) //10**(1/128) 10^6yrを超える時用.
#endif
//////////////////////////////////////////////////



struct orbital_elements{
  char name[30];  //名前（番号.
  double ecc;  //離心率.
  double axis;  //軌道長半径.
  double u;  //離心近点離角.
  double inc;  //軌道傾斜角.
  double Omega;  //昇交点経度.
  double omega;  //近点引数.
  double r_h;  //ヒル半径.
  double radius;  //物理半径.
  double mass;  //質量.
  int orinum;  //初期の番号.
};


#if FRAGMENTATION
struct fragmentation{
  int neighborlist[NEIGHBOR_MAX+1];  //近傍粒子リスト.
  int neighbornumber;  //近傍粒子の個数
  double delta_r_out;  //扇型外側.
  double delta_r_in;  //扇型内側.
  double sigma;  //表面密度.
  double n_s;  //表面数密度.
  double v_ave;  //領域内での平均速度.
  double flux;  //質量フラックス.
  double dt_frag;  //統計的計算の粒子ごとのタイムステップ.
  double t_frag;  //統計的計算の粒子ごとの時間.
  int fragtimes;  //何回統計的計算をしているか.
};


struct parameter{
  double s_1;
  double s_2;
  double s_3;
  double alpha;
  double h_0;
  double Q_D;
};
#endif


#if EXECUTION_TIME
struct execution_time{
  double Energy[3];  //[0]: real time, [1]: user(cpu) time, [2]: system time
  double Orbital_Elements[3];
  double Predictor[3];
  double Corrector[3];
  double Iteration[3];
  double Collision_Judgement[3];
  double Fragmentation[3];
};
EXTERN struct execution_time exetime;  //グローバル変数
#endif


#endif //include-guard
