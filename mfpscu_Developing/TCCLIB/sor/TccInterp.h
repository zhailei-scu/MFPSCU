/*********************************************************************************
!--- Description:
!--- Author : Qing Hou, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : qhou@scu.edu.cn
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!**********************************************************************************/
#ifndef _TccInterp_H
#define _TccInterp_H

#include "DiffusorList.h"
#include "ReactionsList.h"

#ifdef __cplusplus
extern "C" {
#endif
int InterpCScript_DiffusorsDef(char* scriptStr);
#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
extern "C"{
#endif
void GetInterpedDiffusorsArray(CDiffusorDef* theArray);
#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
extern "C" {
#endif
int InterpCScript_ReactionsDef(char* scriptStr);
#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
extern "C"{
#endif
void GetInterpedReactionsArray(CReactionDef* theArray);
#ifdef __cplusplus
}
#endif

#endif
