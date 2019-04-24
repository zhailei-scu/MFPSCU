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
