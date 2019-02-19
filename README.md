# HybridCode

衝突・破壊をとりいれたN体計算

作成者：磯谷和秀

# 1. ディレクトリ構成

## source/
ソースコードとMakefile

## data/
データ保管

## image/
画像保管

## log/
ログ保管（最新版では「data/」中にlog.txtとして保存している）


# 2.「source/」内のファイル概要


## Makefile
コンパイル用のMakefile

コンパイルがうまくいくと、実行ファイルの名前の入力待ちになる

天文台のXC50を使う際には、ジョブを投げる用のシェルスクリプトを生成することもできる

## hybrid.h
ヘッダーファイル

パラメータは基本ここでいじる


## func.h
関数の宣言用ヘッダー


## hybrid_main.c
メインファイル


## acc.c
加速度や加加速度計算


## collision.c
衝突判定、衝突前後の操作


## energy.c
エネルギー計算


## heapsort.c
ヒープソート（階層化タイムステップを導入する際に必要）



## massflux.c
衝突・破壊の際の定常質量フラックス（Kobayashi & Tanaka, 2010）を計算



## neighbor.c
近傍トレーサーの探索、面密度と平均相対速度の計算


## orbital_elements.c
初期軌道要素の設定（or 初期位置の設定）、軌道要素計算


## SFMT.c
メルセンヌ・ツイスタ法による乱数生成


## SFMTdir/
メルセンヌ・ツイスタ法用のヘッダーが入ったディレクトリ

## sub.c
一行でかけるような細々した関数


## timestep.c
タイムステップ計算


## qsub_depend.sh
天文台のXC50にて、ジョブを投げるシェルスクリプト




# 3.ソースコードの詳細

## hybrid.h


```c:hybrid.h
#ifndef INCLUDED_hybrid_H  //include-guard
#define INCLUDED_hybrid_H  //include-guard
.
.
.
#endif //include-guard
```

ヘッダーを複数のファイルで読み込む場合、一度読んだものは再度読まないようにする。

```c:hybrid.h
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
```

stdlib.h, stdio.h, math.hの説明は省略。

sys/stat.h : 時間計測に使う。

time.h : 時間計測に使う。

stdbool.h : Boolean型の関数やtrue(==1), false(==0)を使う。

mpi.h : MPI並列に使う。Macには標準では入っていないので、Macportsとかで入れ、コンパイルオプションでライブラリなどを指定する。天文台のXC50ではコンパイルオプションはいらない。

SFMTdir/SFMT.h : メルセンヌ・ツイスタ法を使う。

omp.h : OpenMP並列に使う。gccではコンパイルオプションに-fopenmpをつける。このオプションがない場合、\_OPENMPは定義されないため、インクルードしないようにしている。


```c
#if __GNUC__ == 7
#define CONST const
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define CONST
#define ALWAYS_INLINE
#endif
```

関数の引数としてポインタを渡しつつもそれを変更しない場合（配列とか）、型の前にconstをつけることが推奨されているが、gccがver4より古いときにエラーになるので、CONSTというマクロを変わりに書いている。この例ではgccのver7の場合のみCONSTをconstとして定義している。

さらに、gccでは関数のインライン展開を強制的に行いたいときに、__attribute__((always_inline))を関数名の前につけるが、これも古いversionでは使えないため、ALWAYS_INLINEというマクロを変わりに書いている。

```c
#define DIRECTORY ../data/Ntr1E2_t1E8_dtlog_Mtot3E-5_Mmax5E-15_ecc1E-1_nofrag_acc/  //ディレクトリ.
#define SUBDIRECTORY rand  //子ディレクトリ. rand%02d
```

データを書き出すファイルを置くディレクトリを指定する。

```c
#define STR_(str) #str
#define STR(str) STR_(str)
```

マクロを文字列に変換するマクロ。

```c
#define INV_3 0.33333333333333333333
```

1/3を定義。

```c
#ifndef EXTERN
#define EXTERN extern
#endif
```

グローバル変数をこのヘッダーで宣言する際に、externを型の前につける。そして、グローバル変数を定義するファイル（hybrid_main.c）にてこのヘッダーを読み込む前にEXTERNマクロを定義するようにする。その他のファイルではEXTERNを定義しないようにしているため、定義が重複しない。

```c
//////////////////////////////////////////////////
#define N_tr 100  //初期のトレーサーの数.
#define N_p 1  //初期の原始惑星の数.
#define ECC_RATIO 10.0  //ecc=0.01の何倍か. inc=ecc/2.
#define STEP_INTERVAL 5.0E5  //何ステップごとに標準出力するか.
//#define BREAK_TIME 100.0  //4h = 14400sec, 12h = 43200sec.
#define BREAK_TIME 14100.0  //4h = 14400sec, 12h = 43200sec.
//#define BREAK_TIME 42900.0  //4h = 14400sec, 12h = 43200sec.

#define RAYLEIGH_DISTRIBUTION true  //離心率や軌道傾斜角の分布 true : Rayleigh, false : v_relが軌道長半径によらず一定.

#define FRAGMENTATION false  //破壊 近傍粒子探索と質量フラックス計算.
#define COLLISION true  //衝突.
#if COLLISION
#define COALESCENCE true  //衝突後に合体.
#else
//#define SOFTENING true  //衝突しないようソフトニング.
#endif
```

計算パラメータ。true, falseはそれぞれ1, 0を意味しており、各関数の中で\#ifを使って場合分けしている。

BREAK_TIMEは実行する実時間（秒）。XC50では一回の実行時間が制限されているため。

```c
EXTERN int global_n;  //グローバル変数.
EXTERN int global_n_p;
EXTERN int global_myid;
EXTERN sfmt_t sfmt;
EXTERN FILE *fplog;

//////////////////////////////////////////////////
```

グローバル変数の宣言。

```c
//////////////////////////////////////////////////
#define ENERGY_FILE true  //エネルギー計算&ファイル作成.
#define ORBITALELEMENTS_FILE true  //軌道要素計算&ファイル作成.
#define POSI_VELO_FILE true  //位置速度ファイル作成.
#define COLLISION_FILE true  //衝突直前の位置速度ファイル作成.
#define TRACERLIST_FILE true  //トレーサーリストのファイル作成.
#define EXECUTION_TIME true  //mainの実行時間測定.
#define EXECUTION_TIME_FUNC false  //mainかつ関数ごとの実行時間測定.
#if EXECUTION_TIME
#include <sys/time.h>
#include <sys/resource.h>
#endif
//////////////////////////////////////////////////
```

true, falseはそれぞれ1, 0を意味しており、各関数の中で\#ifを使って場合分けしている。

実行時間を図る場合、実時間についてはsystem/time.h、CPU時間とシステム時間についてはsys/resource.hが必要。

```c
//////////////////////////////////////////////////
#define INTERACTION_ALL false  //全粒子同士の重力相互作用.
#define INTERACTION_PLANET_TRACER true  //惑星とトレーサー間の相互作用.
#define INTERACTION_TEST_PARTICLE false  //トレーサーをテスト粒子として扱う.
#define INDIRECT_TERM true  //中心星が動く効果を補正.
#define EJECTION false  //初期に破片（トレーサー）を放出する.
#define ORBITING_SMALL_PARTICLE true  //初期に微惑星をケプラー運動させておく.
#define ELIMINATE_PARTICLE false  //太陽に飲みこまれるか系外へ出て行くかで粒子を消す.
//////////////////////////////////////////////////
```

true, falseはそれぞれ1, 0を意味しており、各関数の中で\#ifを使って場合分けしている。

```c
//////////////////////////////////////////////////
//#define G 1.0  //重力定数.
//#define M_0 1.0  //主星の質量.
#if SOFTENING
#define EPSILON 5.21495378928615e-05   //ソフトニングパラメーター.
#endif
#define ETA 1.0E-2  //刻み幅調整.
#define ITE_MAX 2  //イテレーション回数（修正子計算の回数はITE_MAX+1）.
//////////////////////////////////////////////////
```

エルミート法で必要なパラメータ。

```c
#if ELIMINATE_PARTICLE
//////////////////////////////////////////////////
#define SOLAR_RADIUS 0.00465040106951872  //[AU] 6.957E10/1.496E13.
#define SOLAR_SYSTEM_LIMIT 100.0  //[AU]
//////////////////////////////////////////////////
#endif
```

太陽に飲みこまれるか系外へ出て行くかで粒子を消す場合に必要なパラメータ。

```c
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
```

惑星の初期パラメータ。

```c
#if N_tr != 0
//////////////////////////////////////////////////
#define M_TOT (3.0E-7*N_p)  //0.1M_E * N_p  //トレーサーの総質量.
```

トレーサーの総質量。

```c
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
```

トレーサーの初期パラメータ。EJECTION（破片として放出）かORBITATING_SMALL_PARTICLE（軌道上を運動）を上で選ぶ。

```c
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
```

破壊計算のパラメータ。

```c
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
```

何年分計算するかを設定。重力定数を1とし、質量を太陽質量、距離を1AUで規格化すると、時間は1/2$\pi$年に規格化される。そのため、「年」を設定する場合には2$\pi$を掛ける。

```c
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
```

各構造体の宣言。
orbital_elementsは軌道要素の構造体（粒子毎にもつので配列）。
fragmentationは破壊計算に必要なデータをもつ構造体（トレーサー毎にもつので配列）。
parameterは破壊計算のパラメータの構造体。
execution_timeは実行開始からの時間をもつ構造体。



## func.h

インクルードガードあり。

関数の宣言が主だが、複数のファイルにある関数からよく呼ばれるような簡単な関数を、インライン関数として定義する。

```c:func.h
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
```

インライン関数。

内積や距離の計算や、int型とdouble型の変数の入れ替え、2つの数字の大小、実行時間計算


```c:func.h
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
```

予測子計算（インライン関数）。

1. i
i 粒子。
2. x_0[][4]
初期位置。
3. v_0[][4]
初期速度。
4. a_0[][4]
加速度。
5. adot_0[][4]
加加速度。
6. x_p[][4]
位置の予測子。
7. v_p[][4]
速度の予測子。
8. r_p[]
中心星からの距離の予測子。
9. v2_p[]
速度の予測子の大きさの2乗。
10. r_dot_v[]
位置ベクトルと速度ベクトルの内積。
11. Dt[]
タイムステップ。（システム時間）ー（粒子毎の時間）

## acc.c

```c:acc.c
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
```

i, j粒子間の加速度計算。
$G m_j \frac{(\vec{x}_j - \vec{x}_i)}{({r_{ij}^2 + \varepsilon^2 )}^{3/2}}$

1. i
i 粒子。
2. j
j 粒子。
3. k
ベクトル3成分。
4. m_j
j 粒子の質量。
5. x_0[][4]
粒子の位置x, y, zの配列。
6. abs_r[]
i, j 粒子間の距離。


```c:acc.c
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
```

i, j粒子間の加加速度計算。
$G m_j  \left[ \frac { (\vec{v}_j - \vec{v}_i) } { ({r_{ij}^2 + \varepsilon^2 )}^{3/2} } - 3 \frac { \left[ (\vec{v}_j - \vec{v}_i) \cdot (\vec{x}_j - \vec{x}_i) \right] (\vec{x}_j - \vec{x}_i) } { ({r_{ij}^2 + \varepsilon^2 )}^{5/2} } \right]$

1. i
i 粒子。
2. j
j 粒子。
3. k
ベクトル3成分。
4. m_j
j 粒子の質量。
5. x_0[][4]
粒子の位置$x, y, z$の配列。
6. v_0[][4]
粒子の速度$v_x, v_y, v_z$の配列。
7. r_dot_v_ij[]
i, j 粒子の相対位置ベクトルと相対速度ベクトルの内積の配列。
8. abs_r[]
i, j 粒子間の距離。

```c:acc.c
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
```

加速度のindirect項計算。
$- G m_i \frac{\vec{x}_i}{r_i}$

1. i
i 粒子。
2. k
ベクトル3成分。
3. m_i
i 粒子の質量。
4. x_0[][4]
粒子の位置$x, y, z$の配列。
5. r_0[]
i 粒子の原点からの距離（原点は中心星を含めた重心）。


```c:acc.c
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
```

加加速度のindirect項計算。
$- G m_i \left( \frac{\vec{v}_i}{r_i^3} - 3 \frac{(\vec{x}_i \cdot \vec{v}_i) \vec{x}_i}{r_i^5} \right)$

1. i
i 粒子。
2. k
ベクトル3成分。
3. m_i
i 粒子の質量。
4. x_0[][4]
粒子の位置$x, y, z$の配列。
5. v_0[][4]
粒子の速度$v_x, v_y, v_z$の配列。
6. r_0[]
i 粒子の原点からの距離（原点は中心星を含めた重心）。
7. r_dot_v_ij[]
i, j 粒子の相対位置ベクトルと相対速度ベクトルの内積の配列。



```c:acc.c
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
```

加速度の外力項。
$- G M_{\ast} \frac{\vec{x}_i}{r_i^3}$

1. i
i 粒子。
2. k
ベクトル3成分。
3. x_0[][4]
粒子の位置$x, y, z$の配列。
4. r_0[]
i 粒子の原点からの距離（原点は中心星を含めた重心）。

```c:acc.c
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
```

加加速度の外力項。
$- G M_{\ast} \left( \frac{\vec{v}_i}{r_i^3} - 3 \frac{(\vec{x}_i \cdot \vec{v}_i) \vec{x}_i}{r_i^5} \right)$

1. i
i 粒子。
2. k
ベクトル3成分。
3. x_0[][4]
粒子の位置$x, y, z$の配列。
4. v_0[][4]
粒子の速度$v_x, v_y, v_z$の配列。
5. r_0[]
i 粒子の原点からの距離（原点は中心星を含めた重心）。
6. r_dot_v[]
i 粒子の位置ベクトルと速度ベクトルの内積の配列。

```c:acc.c
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
```

修正子計算。

a_tmp[]などを新たにローカル関数として定義している理由は、加速度の関数All_Accelerationと加加速度の関数All_dAcceleration内で並列化をする際に、別の粒子を抜け駆けで更新しないようにするため。

1. n_ite
イテレーション回数。
2. i_sys
システム粒子。
3. *ele_p
軌道要素の構造体のポインタ。
4. x_p[][4]
位置の予測子。
5. v_p[][4]
速度の予測子。
6. r_p[]
中心星からの距離の予測子。
7. x_c[][4]
位置の修正子。
8. v_c[][4]
速度の修正子。
9. r_c[]
中心星からの距離の修正子。
10. v2_c[]
速度の修正子の大きさの2乗。
11. r_dot_v[]
位置ベクトルと速度ベクトルの内積。
12. a_0[][4]
現時点の加速度。
13. adot_0[][4]
現時点の加加速度。
14. a[][4]
次のステップの加速度。
15. adot[][4]
次のステップの加加速度。
16. adot2_dt2[][4]
次のステップの加速度の第2次導関数 × タイムステップの2乗。
17. adot3_dt3[][4]
次のステップの加速度の第3次導関数 × タイムステップの3乗。
18. dt_[]
粒子毎のタイムステップ。
19. t_dyn
システム時間。破壊計算で質量は減少するのでその見積もり用。
20. *frag_p
破壊計算に必要なデータをもつ構造体ポインタ。破壊計算で質量は減少するのでその見積もり用。

```c:acc.c
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
```

i 粒子の加速度。

1. i
i 粒子。
2. k
ベクトル3成分。
3. *ele_p
軌道要素の構造体のポインタ。
4. x_0[][4]
粒子の位置$x, y, z$の配列。
5. r_0[]
i 粒子の原点からの距離（原点は中心星を含めた重心）。
6. abs_r[]
i, j 粒子間の距離。
7. t_dyn
システム時間。破壊計算で質量は減少するのでその見積もり用。
8. *frag_p
破壊計算に必要なデータをもつ構造体ポインタ。破壊計算で質量は減少するのでその見積もり用。

```c:acc.c
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
```

i 粒子の加加速度。

1. i
i 粒子。
2. k
ベクトル3成分。
3. *ele_p
軌道要素の構造体のポインタ。
4. x_0[][4]
粒子の位置$x, y, z$の配列。
5. v_0[][4]
粒子の速度$v_x, v_y, v_z$の配列。
6. r_dot_v[]
i 粒子の位置ベクトルと速度ベクトルの内積の配列。
7. r_dot_v_ij[]
i, j 粒子の相対位置ベクトルと相対速度ベクトルの内積の配列。
8. r_0[]
i 粒子の原点からの距離（原点は中心星を含めた重心）。
9. abs_r[]
i, j 粒子間の距離。
10. t_dyn
システム時間。破壊計算で質量は減少するのでその見積もり用。
11. *frag_p
破壊計算に必要なデータをもつ構造体ポインタ。破壊計算で質量は減少するのでその見積もり用。


## collision.c

```c:collision.c
bool Collision_Judgement(int i_sys, CONST struct orbital_elements *ele_p, CONST double x_p[][4], double abs_r[], int *i_col, int *j_col){


#if EXECUTION_TIME && EXECUTION_TIME_FUNC
  struct timeval realtime_start, realtime_end;
  struct rusage usage_start, usage_end;
  gettimeofday(&realtime_start,NULL);
  getrusage(RUSAGE_SELF,&usage_start);
#endif


  int j;

  for(j=1;j<=global_n;++j){
    if(i_sys!=j){
      abs_r[j] = RelativeDistance(i_sys,j,x_p);  //絶対値.
      if(abs_r[j] < ((ele_p+i_sys)->radius) + ((ele_p+j)->radius)){
		(*i_col) = Min_int(i_sys,j);
		(*j_col) = Max_int(i_sys,j);

#if EXECUTION_TIME && EXECUTION_TIME_FUNC
		gettimeofday(&realtime_end,NULL);
		getrusage(RUSAGE_SELF,&usage_end);
		exetime.Collision_Judgement[0] += Cal_time(realtime_start,realtime_end);
		exetime.Collision_Judgement[1] += Cal_time(usage_start.ru_utime,usage_end.ru_utime);
		exetime.Collision_Judgement[2] += Cal_time(usage_start.ru_stime,usage_end.ru_stime);
#endif

		return (true);  //衝突した場合.
      }
    }
  }


#if EXECUTION_TIME && EXECUTION_TIME_FUNC
  gettimeofday(&realtime_end,NULL);
  getrusage(RUSAGE_SELF,&usage_end);
  exetime.Collision_Judgement[0] += Cal_time(realtime_start,realtime_end);
  exetime.Collision_Judgement[1] += Cal_time(usage_start.ru_utime,usage_end.ru_utime);
  exetime.Collision_Judgement[2] += Cal_time(usage_start.ru_stime,usage_end.ru_stime);
#endif

  return (false);  //衝突しない場合.
}
```

衝突判定。

1. i_sys
システム粒子。
2. *ele_p
軌道要素の構造体ポインタ。
3. x_p[][4]
位置の予測子。
4. abs_r[]
i, j 粒子間の距離。
5. *i_col
衝突した粒子番号ポインタ。若い方。
6. *j_col
衝突した粒子番号ポインタ。もう一つ。


```c:collision.c
void Energy_Correction(int i_col, int j_col, CONST double x_0[][4], CONST double v_0[][4], CONST struct orbital_elements *ele_p, double *dE_heat, double *dE_grav, double *dE_c, double *v_imp
#if FRAGMENTATION
		       , double t_dyn
		       , CONST struct fragmentation *frag_p
#endif
		       ){

  double m_1, m_2;
  double abs_v2 = SquareOfRelativeVelocity(i_col,j_col,v_0);  //相対速度2乗.
  double r_p12 = RelativeDistance(i_col,j_col,x_0);  //2粒子間の距離.
  double r_g12 = RadiusFromCenter(0,x_0);  //2粒子の重心と中心星との距離.
  double r_1 = RadiusFromCenter(i_col,x_0);
  double r_2 = RadiusFromCenter(j_col,x_0);

#if FRAGMENTATION
  m_1 = MassDepletion(i_col,((ele_p+i_col)->mass),t_dyn,frag_p);
  m_2 = MassDepletion(j_col,((ele_p+j_col)->mass),t_dyn,frag_p);
#else
  m_1 = ((ele_p+i_col)->mass);
  m_2 = ((ele_p+j_col)->mass);
#endif

  (*dE_heat) = - 0.5 * m_1 * m_2 / (m_1 + m_2) * abs_v2;  //完全合体することで、相対速度分の運動エネルギーが熱エネルギーとなって散逸する.

#ifndef G
  (*dE_grav) = m_1 * m_2 / r_p12;  //2粒子間の距離に対応する相互重力エネルギーがなくなっている.
#else
  (*dE_grav) = G * m_1 * m_2 / r_p12;  //2粒子間の距離に対応する相互重力エネルギーがなくなっている.
#endif

#if !defined(G) && !defined(M_0)
  (*dE_c) = - (m_1 + m_2) / r_g12 + m_1 / r_1 + m_2 / r_2;  //中心星ポテンシャルエネルギーが変わる.
#else
  (*dE_c) = G * M_0 * (- (m_1 + m_2) / r_g12 + m_1 / r_1 + m_2 / r_2);  //中心星ポテンシャルエネルギーが変わる.
#endif

  (*v_imp) = sqrt(abs_v2);

  fprintf(fplog,"dE_heat=%e\tdE_grav=%e\tdE_c=%e\tv_imp=%e\n",(*dE_heat),(*dE_grav),(*dE_c),(*v_imp));

  return;
}
```

衝突時のエネルギー補正。

1. i_col
衝突した粒子。若い方。
2. j_col
衝突した粒子。もう一つ。
3. x_0[][4]
位置。
4. v_0[][4]
速度。
5. *ele_p
軌道要素の構造体ポインタ。
6. *dE_heat
合体時に相対速度分の運動エネルギーが熱として散逸する量を表すポインタ。
7. *dE_grav
衝突した2粒子間の相互重力エネルギーが消える量を表すポインタ。
8. *dE_c
中心星ポテンシャルエネルギーが変わる量を表すポインタ。
9. *v_imp
衝突速度を表すポインタ。
10. t_dyn
システム時間。破壊計算で質量は減少するのでその見積もり用。
11. *frag_p
破壊計算に必要なデータをもつ構造体ポインタ。破壊計算で質量は減少するのでその見積もり用。

```c:collision.c
void Coalescence(int i_col, int j_col, double x_0[][4], double v_0[][4], struct orbital_elements *ele_p
#if FRAGMENTATION
		 , double t_dyn
		 , struct fragmentation *frag_p
#endif
		 ){

  int k;

  //i_colを新しい合体粒子の番号にする.
#if FRAGMENTATION
  ((ele_p+i_col)->mass) = MassDepletion(i_col,((ele_p+i_col)->mass),t_dyn,frag_p) + MassDepletion(j_col,((ele_p+j_col)->mass),t_dyn,frag_p);
#else
  ((ele_p+i_col)->mass) = ((ele_p+i_col)->mass) + ((ele_p+j_col)->mass);
#endif
  ((ele_p+i_col)->radius) = cbrt(3.0/4.0/M_PI*((ele_p+i_col)->mass)*1.989E33/PLANET_DENSITY)/1.496E13;
  for(k=1;k<=3;++k){
    x_0[i_col][k] = x_0[0][k];
    v_0[i_col][k] = v_0[0][k];
  }

  if(j_col<=global_n_p){

    //惑星同士の衝突合体
    //j_colとglobal_n_pを入れ替える.
    *(ele_p+0) = *(ele_p+j_col);  //構造体のためSwap関数は使えない. 0番目の要素はコピーに使うだけ.
    *(ele_p+j_col) = *(ele_p+global_n_p);
    *(ele_p+global_n_p) = *(ele_p+0);
#if FRAGMENTATION
    *(frag_p+0) = *(frag_p+j_col);  //構造体のためSwap関数は使えない. 0番目の要素はコピーに使うだけ.
    *(frag_p+j_col) = *(frag_p+global_n_p);
    *(frag_p+global_n_p) = *(frag_p+0);
#endif
    for(k=1;k<=3;++k){
      Swap_double(&x_0[j_col][k],&x_0[global_n_p][k]);
      Swap_double(&v_0[j_col][k],&v_0[global_n_p][k]);
    }

#if N_tr != 0
    //global_n_pとglobal_nを入れ替える.
    *(ele_p+0) = *(ele_p+global_n_p);  //構造体のためSWAP関数は使えない. 0番目の要素はコピーに使うだけ.
    *(ele_p+global_n_p) = *(ele_p+global_n);
    *(ele_p+global_n) = *(ele_p+0);
#if FRAGMENTATION
    *(frag_p+0) = *(frag_p+global_n_p);  //構造体のためSWAP関数は使えない. 0番目の要素はコピーに使うだけ.
    *(frag_p+global_n_p) = *(frag_p+global_n);
    *(frag_p+global_n) = *(frag_p+0);
#endif
    for(k=1;k<=3;++k){
      Swap_double(&x_0[global_n_p][k],&x_0[global_n][k]);
      Swap_double(&v_0[global_n_p][k],&v_0[global_n][k]);
    }
#endif

    //global_n_p,global_nを1つ減らす.
    global_n_p--;
    global_n--;

  }else if(j_col>global_n_p){

    //惑星-微惑星の衝突合体
    //j_colとglobal_nを入れ替える.
    *(ele_p+0) = *(ele_p+j_col);  //構造体のためSwap関数は使えない. 0番目の要素はコピーに使うだけ.
    *(ele_p+j_col) = *(ele_p+global_n);
    *(ele_p+global_n) = *(ele_p+0);
#if FRAGMENTATION
    *(frag_p+0) = *(frag_p+j_col);  //構造体のためSwap関数は使えない. 0番目の要素はコピーに使うだけ.
    *(frag_p+j_col) = *(frag_p+global_n);
    *(frag_p+global_n) = *(frag_p+0);
#endif
    for(k=1;k<=3;++k){
      Swap_double(&x_0[j_col][k],&x_0[global_n][k]);
      Swap_double(&v_0[j_col][k],&v_0[global_n][k]);
    }

    //global_nを1つ減らす.
    global_n--;
  }

  return;
}
```

衝突・合体計算。

1. i_col
衝突した粒子。若い方。
2. j_col
衝突した粒子。もう一つ。
3. x_0[][4]
位置。
4. v_0[][4]
速度。
5. *ele_p
軌道要素の構造体ポインタ。
6. t_dyn
システム時間。破壊計算で質量は減少するのでその見積もり用。
7. *frag_p
破壊計算に必要なデータをもつ構造体ポインタ。破壊計算で質量は減少するのでその見積もり用。



## energy.c
エネルギー計算

```c:energy.c
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

#if !INTERACTION_TEST_PARTICLE

#if INTERACTION_PLANET_TRACER
    if(i<=global_n_p){
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

#if INTERACTION_PLANET_TRACER
    }
#endif

#endif  /*!INTERACTION_TEST_PARTICLE*/

#if !defined(G) && !defined(M_0)
    E_tot += - m_i / r_c[i] + E[i];
#else
    E_tot += - G * M_0 * m_i / r_c[i] + E[i];
#endif

  }  //i loop
  return E_tot;
}
```

全エネルギー計算。

1. *ele_p
軌道要素の構造体ポインタ。
2. x_c[][4]
位置の修正子。
3. v_c[][4]
速度の修正子。
4. v_G[]
重心速度。
5. v2_c[]
速度の修正子の大きさの2乗。
6. r_c[]
中心星からの距離。
7. t_dyn
システム時間。破壊計算で質量は減少するのでその見積もり用。
8. *frag_p
破壊計算に必要なデータをもつ構造体ポインタ。破壊計算で質量は減少するのでその見積もり用。

```c:energy.c
double AngularMomentum(CONST struct orbital_elements *ele_p, CONST double x_0[][4], CONST double v_0[][4]
#if FRAGMENTATION
		       , double t_dyn
		       , CONST struct fragmentation *frag_p
#endif
		       ){
  int i, k;
  double L[global_n+1][4];
  double L_tot_0[4];
  double m_i;


  for(i=1;i<=global_n;++i){

#if FRAGMENTATION
    m_i = MassDepletion(i,((ele_p+i)->mass),t_dyn,frag_p);
#else
    m_i = ((ele_p+i)->mass);
#endif

    L[i][1] = m_i * (x_0[i][2]*v_0[i][3] - x_0[i][3]*v_0[i][2]);
    L[i][2] = m_i * (x_0[i][3]*v_0[i][1] - x_0[i][1]*v_0[i][3]);
    L[i][3] = m_i * (x_0[i][1]*v_0[i][2] - x_0[i][2]*v_0[i][1]);
    //fprintf(fplog,"i=%d\t(ele_p+i)->mass=%e\n",i,(ele_p+i)->mass);
    //fprintf(fplog,"L[%d][1] = %e\t[2] = %e\t[3] = %e\n",i,L[i][1],L[i][2],L[i][3]);
  }

  for(k=1;k<=3;++k){
    L_tot_0[k] = 0.0;
    for(i=1;i<=global_n;++i){
      L_tot_0[k] += L[i][k];
    }
  }

  return sqrt(L_tot_0[1]*L_tot_0[1] + L_tot_0[2]*L_tot_0[2] + L_tot_0[3]*L_tot_0[3]);
}
```

全角運動量計算。

1. *ele_p
軌道要素の構造体ポインタ。
2. x_0[][4]
位置。
3. v_0[][4]
速度。
4. t_dyn
システム時間。破壊計算で質量は減少するのでその見積もり用。
5. *frag_p
破壊計算に必要なデータをもつ構造体ポインタ。破壊計算で質量は減少するのでその見積もり用。

## heapsort.c
ヒープソート（階層化タイムステップを導入する際に必要）



## massflux.c
衝突・破壊の際の定常質量フラックス（Kobayashi & Tanaka, 2010）を計算



## neighbor.c
近傍トレーサーの探索、面密度と平均相対速度の計算


## orbital_elements.c
初期軌道要素の設定（or 初期位置の設定）、軌道要素計算


## SFMT.c
メルセンヌ・ツイスタ法による乱数生成


## SFMTdir/
メルセンヌ・ツイスタ法用のヘッダーが入ったディレクトリ

## sub.c
一行でかけるような細々した関数


## timestep.c
タイムステップ計算


## qsub_depend.sh
天文台のXC50にて、ジョブを投げるシェルスクリプト




# 5.文字の修飾(イタリック、太字)

* イタリック
「_」または「*」で文字をくくります。

>
\_イタリック\_
_イタリック_
>
\*イタリック\*
*イタリック*

* 太字
「__」または「**」で文字をくくります。

>
\_\_太字\_\_
__太字__
>
\*\*太字\*\*
**太字**

# 6.リスト
リストの上下に空白を入れないと正しく表示されないので注意。
また、記号と文の間に半角スペースを入れること。

* 順序なしリスト

文頭に「*」「+」「-」のいずれかを入れる。
>
\* 順序なしリスト
>
* 順序なしリスト


* 順序つきリスト

文頭に「数字.」を入れる。
見た目はほぼ変わりません。

>
1. リスト1
2. リスト2

# 7. 水平線

「*」か「-」を3つ以上一行に書く。
以下は全て水平線となる。
>
\*\*\*
\* \* \*
\-\-\-
\- \- \-
>
全部以下の水平線
***

# 8. テーブル

以下のようにテーブルを組みます。
基本は「|」でくくっていくことです。
2行目がポイントで、2行目のコロンの位置によってセル内の文字の配置が変わります。

>
\|左揃え|中央揃え|右揃え|
|:---|:---:|--:|
|align-left|align-center|align-right|
|セルの左揃えです|セルの中央揃えです|セルの右揃えです|
>
>
|左揃え|中央揃え|右揃え|
|:---|:---:|--:|
|align-left|align-center|align-right|
|セルの左揃えです|セルの中央揃えです|セルの右揃えです|

# 9. マークダウンのエスケープ
「\」をMarkdownの前につけることでMarkdownを無効化出来ます。
この記事ではこれを多用しました。
>
\\#見出しh1
とすると
\#見出しh1
となります。

# 10. 補足:ページをMarkdownで見る
Qiitaを見ていると「これはどんな記法で書いてあるんだろう」ときになることがあるかもしれません。
そんな時はMarkdown記法で見たいURLの最後に.mdをつければ見ることが出来ます。

##参考

[Markdown記法チートシート](http://qiita.com/Qiita/items/c686397e4a0f4f11683d)
<!--stackedit_data:
eyJoaXN0b3J5IjpbNDcyMDI1MDgsLTEyNjQ1OTM1MjMsMTE3MD
IyMzAwOCwtMTE2NjUyNDc1LDEzNDI3MzkwMzEsNTE5Mzg3MDAx
LC0xNTI5NjczNTYsMjEyMzk0MDQ4MywtMTU2Nzk3MDQzNSw5MT
k5NTYzNjUsMTYwOTcwOTA2MSwtMTQyMjQ1NTQ5OCw5NTE5NTMw
NjEsLTE4MzUxOTg5NTYsMTczODg1NzAxMiwtMTc1NTUzNjI5LC
03ODY3ODA1NTAsLTE5NDI0NzY5NywtMTM0MDc5ODE3NSwtNTE5
NjU1MTgyXX0=
-->