/*********************************************************************************
!--- Description:
!--- Author : Qing Hou, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : qhou@scu.edu.cn
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!**********************************************************************************/

#ifndef __DIFFUSORLIST_H
#define __DIFFUSORLIST_H

#ifdef __cplusplus
struct CDiffusorDef{
#else
typedef struct{
#endif
  char symbol[20];

  //In free matrix
  int DiffusorValueType_Free;
  double DiffuseCoefficient_Free_Value;

  // If the DiffuseCoefficient type is by Arrhenius or BCluster(bigger cluster),use this
  double PreFactor_Free;
  double ActEnergy_Free;

  int ECRValueType_Free;
  double ECR_Free;

  //In GB
  int DiffusorValueType_InGB;
  double DiffuseCoefficient_InGB_Value;

  // If the DiffuseCoefficient type is by Arrhenius or BCluster(bigger cluster),use this
  double PreFactor_InGB;
  double ActEnergy_InGB;

  int ECRValueType_InGB;
  double ECR_InGB;

#ifdef __cplusplus
};
#else
}CDiffusorDef;
#endif

#ifdef __cplusplus
struct CDiffusorList{
#else
typedef struct CDiffusorList{
#endif

  CDiffusorDef data;

  #ifdef __cplusplus
  CDiffusorList* next;
  #else
  struct CDiffusorList* next;
  #endif

  int size;

#ifdef __cplusplus
};
#else
}CDiffusorList;
#endif

void InitDiffusorList(CDiffusorList **list);

void AppendDiffusor(CDiffusorList *list, CDiffusorDef *element);

#endif
