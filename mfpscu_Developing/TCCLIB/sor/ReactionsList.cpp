/*********************************************************************************
!--- Description:
!--- Author : Qing Hou, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : qhou@scu.edu.cn
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!**********************************************************************************/

#include "ReactionsList.h"
#ifdef __cplusplus
#include <iostream>
#include <cstring>
#else
#include <string.h>
#include <stdlib.h>
#endif

#ifdef __cplusplus
using namespace std;
#endif

void copyCReactionDef(CReactionDef *Dest,CReactionDef *Source){

  memset(Dest->SubjectSymbol,0,20);
  strcpy(Dest->SubjectSymbol,Source->SubjectSymbol);
  memset(Dest->ObjectSymbol,0,20);
  strcpy(Dest->ObjectSymbol,Source->ObjectSymbol);
  Dest->ReactionCoefficientType = Source->ReactionCoefficientType;
  Dest->ReactionCoefficient_Value = Source->ReactionCoefficient_Value;
  Dest->PreFactor = Source->PreFactor;
  Dest->ActEnergy = Source->ActEnergy;
  Dest->ECRValueType = Source->ECRValueType;
  Dest->ECR = Source->ECR;
}

#ifdef __cplusplus
void InitReactionsList(CReactionsList **list){
  if(NULL == *list){
    *list = new CReactionsList;
  }

  (*list)->size = 0;
  (*list)->next = NULL;
}
#else
void InitReactionsList(CReactionsList **list){
  if(NULL == *list){
    *list = malloc(sizeof(CReactionsList));
  }

  (*list)->size = 0;
  (*list)->next = NULL;
}
#endif

#ifdef __cplusplus
void AppendReaction(CReactionsList *list, CReactionDef *element){
  CReactionsList *cursor = NULL;
  if(0==list->size){
    copyCReactionDef(&(list->data),element);
    list->size++;
  }else{
    cursor = list;
    while(NULL != cursor->next){
      cursor = cursor->next;
    }
    cursor->next = new CReactionsList[1];
    copyCReactionDef(&(cursor->next->data),element);
    cursor->next->next = NULL;
    list->size++;
    cursor->next->size = list->size;
  }

}
#else
void AppendReaction(CReactionsList *list, CReactionDef *element){
  CReactionsList *cursor = NULL;
  if(0==list->size){
    copyCReactionDef(&(list->data),element);
    list->size++;
  }else{
    cursor = list;
    while(NULL != cursor->next){
      cursor = cursor->next;
    }
    cursor->next = malloc(sizeof(CReactionsList));
    copyCReactionDef(&(cursor->next->data),element);
    cursor->next->next = NULL;
    list->size++;
    cursor->next->size = list->size;
  }

}
#endif
