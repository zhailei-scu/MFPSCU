/*********************************************************************************
!--- Description:
!--- Author : Qing Hou, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : qhou@scu.edu.cn
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!**********************************************************************************/

#ifndef __REACTIONSLIST_H
#define __REACTIONSLIST_H

#ifdef __cplusplus
struct CReactionDef{
#else
typedef struct{
#endif
  char SubjectSymbol[20];
  char ObjectSymbol[20];

  int ReactionCoefficientType;
  double ReactionCoefficient_Value;

  // If the DiffuseCoefficient type is by Arrhenius or BCluster(bigger cluster),use this
  double PreFactor;
  double ActEnergy;

  int ECRValueType;
  double ECR;
#ifdef __cplusplus
};
#else
}CReactionDef;
#endif

#ifdef __cplusplus
struct CReactionsList{
#else
typedef struct CReactionsList{
#endif

  CReactionDef data;

  #ifdef __cplusplus
  CReactionsList* next;
  #else
  struct CReactionsList* next;
  #endif

  int size;

#ifdef __cplusplus
};
#else
}CReactionsList;
#endif

void InitReactionsList(CReactionsList **list);

void AppendReaction(CReactionsList *list, CReactionDef *element);

#endif
