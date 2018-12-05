#include "hybrid.h"
#include "func.h"


/* ヒープソートを行う */
void HeapSort(CONST int NUM_DATA, CONST double t[], CONST double dt[], CONST double t_tmp, int index[]){

  int i, n;
  int leaf, root;
  double a[NUM_DATA+1];

  n = NUM_DATA;

  for(i=1;i<=n;i++){
    index[i] = i;
    a[i] = t[i] + dt[i] + t_tmp;
  }

  leaf = n;                           /* 初期値は末尾の要素 */
  root = n/2;                         /* 初期値はその親 */

  while(root >= 1){                   /* 半順序木を構成 */
    DownHeap(a,index,leaf,root);
    root--;
  }

  while(leaf >= 1){
    Swap_double(&a[1],&a[leaf]);      /* 半順序木の根と末尾の要素を交換 */
    Swap_int(&index[1],&index[leaf]);
    leaf--;                           /* 末尾の要素を半順序木から外す */
    DownHeap(a,index,leaf,1);         /* 半順序木を再構成する */
  }

  return;
}


/* 半順序木 ( ヒープ ) を構成する */
void DownHeap(double a[], int index[], int leaf, int root){

  int i;

  if(root==1){
    i = 2;
  }
  else{
    i = root * 2;
  }
  while(i <= leaf){
    if(i < leaf && a[i+1] > a[i]){    /* a[i] と a[i+1] の大きい方と */
      i++;
    }
    if(a[root] >= a[i]){              /* a[root] と比較 */
      break;                          /* 子の方が大きければ */
    }

    Swap_double(&a[root],&a[i]);      /* 交換 */
    Swap_int(&index[root],&index[i]);
    root = i;                         /* 更にその子についても調べる */
    i = root * 2;
  }

  return;
}
