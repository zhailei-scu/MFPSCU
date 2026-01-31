/*********************************************************************************
!--- Description:
!--- Author : Qing Hou, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : qhou@scu.edu.cn
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!**********************************************************************************/

#include "DiffusorList.h"
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

void copyCDiffusorDef(CDiffusorDef *Dest,CDiffusorDef *Source){

  memset(Dest->symbol,0,20);
  strcpy(Dest->symbol,Source->symbol);

  Dest->DiffusorValueType_Free = Source->DiffusorValueType_Free;
  Dest->DiffuseCoefficient_Free_Value = Source->DiffuseCoefficient_Free_Value;
  Dest->PreFactor_Free = Source->PreFactor_Free;
  Dest->ActEnergy_Free = Source->ActEnergy_Free;
  Dest->ECRValueType_Free = Source->ECRValueType_Free;
  Dest->ECR_Free = Source->ECR_Free;

  Dest->DiffusorValueType_InGB = Source->DiffusorValueType_InGB;
  Dest->DiffuseCoefficient_InGB_Value = Source->DiffuseCoefficient_InGB_Value;
  Dest->PreFactor_InGB = Source->PreFactor_InGB;
  Dest->ActEnergy_InGB = Source->ActEnergy_InGB;
  Dest->ECRValueType_InGB = Source->ECRValueType_InGB;
  Dest->ECR_InGB = Source->ECR_InGB;
}

#ifdef __cplusplus
void InitDiffusorList(CDiffusorList **list){
  if(NULL == *list){
    *list = new CDiffusorList;
  }

  (*list)->size = 0;
  (*list)->next = NULL;
}
#else
void InitDiffusorList(CDiffusorList **list){
  if(NULL == *list){
    *list = malloc(sizeof(CDiffusorList));
  }

  (*list)->size = 0;
  (*list)->next = NULL;
}
#endif

#ifdef __cplusplus
void AppendDiffusor(CDiffusorList *list, CDiffusorDef *element){
  CDiffusorList *cursor = NULL;
  if(0==list->size){
    copyCDiffusorDef(&(list->data),element);
    list->size++;
  }else{
    cursor = list;
    while(NULL != cursor->next){
      cursor = cursor->next;
    }
    cursor->next = new CDiffusorList[1];
    copyCDiffusorDef(&(cursor->next->data),element);
    cursor->next->next = NULL;
    list->size++;
    cursor->next->size = list->size;
  }

}
#else
void AppendDiffusor(CDiffusorList *list, CDiffusorDef *element){
  CDiffusorList *cursor = NULL;
  if(0==list->size){
    copyCDiffusorDef(&(list->data),element);
    list->size++;
  }else{
    cursor = list;
    while(NULL != cursor->next){
      cursor = cursor->next;
    }
    cursor->next = malloc(sizeof(CDiffusorList));
    copyCDiffusorDef(&(cursor->next->data),element);
    cursor->next->next = NULL;
    list->size++;
    cursor->next->size = list->size;
  }

}
#endif
