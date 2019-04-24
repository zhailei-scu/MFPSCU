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
