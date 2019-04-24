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
