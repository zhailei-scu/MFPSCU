//As the code::block cannot identify the macro string, so in code::block IDE
//we fix the following macro. While we public the package, these should be remove because when we
//use the makefile instead of code::block IDE, we can define these macro while installing and building
#ifndef LIBTCCPATH
#define LIBTCCPATH "/home/zhailei/Development/MeanField/tcc-0.9.27/"
#endif

#ifndef LIBTCCHEADERF
#define LIBTCCHEADERF "/home/zhailei/Development/MeanField/tcc-0.9.27/libtcc.h"
#endif

#ifndef DIFFUSORDEFPATH
#define DIFFUSORDEFPATH "/home/zhailei/Development/MeanField/mcpscu_Developing_New/TCCLIB/sor/"
#endif

#include LIBTCCHEADERF
#ifdef __cplusplus
#include <string>
#include <iostream>
#include <sstream>
#else
#include <string.h>
#endif
#include "TccInterp.h"

#ifdef __cplusplus
using namespace std;
#endif


typedef void(*SetDiffusor)(CDiffusorList *list);
typedef void(*SetReaction)(CReactionsList *list);

extern void copyCDiffusorDef(CDiffusorDef *Dest,CDiffusorDef *Source);
extern void copyCReactionDef(CReactionDef *Dest,CReactionDef *Source);

CDiffusorList* DiffusorsList = NULL;
CReactionsList* ReactionsList = NULL;
/*
--------------------------------------------
*/

#ifdef __cplusplus
int InterpCScript_DiffusorsDef(char* scriptStr){

  TCCState *s = NULL;

  stringstream pathStream;

  pathStream.clear();
  pathStream.str("");

  // @todo (zhail#1#):
  //The current version of codeblock is bug at the string macro,it can not reconigze the string macro during making command,
  //so the TCC path is putted in program tempory, when public these project, we should move it to make file.
  //The way to move it to makefile please note the TestTCC demo
  //pathStream<<LIBTCCPATH;
  pathStream<<LIBTCCPATH;
  cout<< pathStream.str().c_str()<<endl;

  s = tcc_new();

  tcc_add_include_path(s, pathStream.str().c_str());

  tcc_set_lib_path(s, pathStream.str().c_str());


  pathStream.clear();
  pathStream.str("");

  // @todo (zhail#1#):
  //The current version of codeblock is bug at the string macro,it can not reconigze the string macro during making command,
  //so the TCC path is putted in program tempory, when public these project, we should move it to make file.
  //The way to move it to makefile please note the TestTCC demo
  //pathStream<<DIFFUSORDEFPATH;
  pathStream<<DIFFUSORDEFPATH;
  tcc_add_include_path(s,pathStream.str().c_str());

  if(!s){
    cout<<"MCPSCUERROR: The TCC interperation init failed!"<<endl;
    return -1;
  }

  tcc_set_output_type(s, TCC_OUTPUT_MEMORY);


  cout<<scriptStr<<endl;

  if(tcc_compile_string(s, scriptStr) == -1){
    cout<<scriptStr<<endl;
    cout<<"MCPSCUERROR: The above script compille failed !"<<endl;
    return -1;
  }

  tcc_add_symbol(s,"InitDiffusorList",(void*)InitDiffusorList);
  tcc_add_symbol(s,"AppendDiffusor",(void*)AppendDiffusor);

  if(tcc_relocate(s, TCC_RELOCATE_AUTO) < 0)
     return 1;

  SetDiffusor p_SetDiffusor = (SetDiffusor)tcc_get_symbol(s,"SetDiffusor");

  InitDiffusorList(&DiffusorsList);

  if(!p_SetDiffusor){
    cout<<scriptStr<<endl;
    cout<<"MCPSCUERROR: The above script interperate failed!"<<endl;
    return -1;
  }else{
    p_SetDiffusor(DiffusorsList);
  }

  return DiffusorsList->size;
}
#else
int InterpCScript_DiffusorsDef(char* scriptStr){

  TCCState *s = NULL;
  char path[100];

  // @todo (zhail#1#):
  //The current version of codeblock is bug at the string macro,it can not reconigze the string macro during making command,
  //so the TCC path is putted in program tempory, when public these project, we should move it to make file.
  //The way to move it to makefile please note the TestTCC demo
  //pathStream<<LIBTCCPATH;
  strcpy(path,LIBTCCPATH);

  s = tcc_new();

  tcc_add_include_path(s, path);
  tcc_set_lib_path(s, path);

  memset(path,0,sizeof(path));
  strcpy(path,DIFFUSORDEFPATH);
  tcc_add_include_path(s, path);

  if(!s){
    printf("%s\n","MCPSCUERROR: The TCC interperation init failed!");
    return -1;
  }

  tcc_set_output_type(s, TCC_OUTPUT_MEMORY);

  if(tcc_compile_string(s, scriptStr) == -1){
    printf("%s\n",scriptStr);
    printf("%s\n","MCPSCUERROR: The above script compille failed !");
    return -1;
  }

  tcc_add_symbol(s,"InitDiffusorList",(void*)InitDiffusorList);
  tcc_add_symbol(s,"AppendDiffusor",(void*)AppendDiffusor);

  if(tcc_relocate(s, TCC_RELOCATE_AUTO) < 0)
     return 1;

  SetDiffusor p_SetDiffusor = (SetDiffusor)tcc_get_symbol(s,"SetDiffusor");

  InitDiffusorList(&DiffusorsList);

  if(!p_SetDiffusor){
    printf("%s\n",scriptStr);
    printf("%s\n","MCPSCUERROR: The above script interperate failed!");
    return -1;
  }else{
    p_SetDiffusor(DiffusorsList);
  }

  return DiffusorsList->size;
}
#endif

/*
--------------------------------------------------------
*/
void GetInterpedDiffusorsArray(CDiffusorDef* theArray){

  int I;

  int theSize = DiffusorsList->size;

  CDiffusorList* cursor;

  cursor = DiffusorsList;

  for(I=0;I<theSize;I++){
    copyCDiffusorDef(&(theArray[I]),&(cursor->data));
    cursor = cursor->next;
  }

  return;
}


#ifdef __cplusplus
int InterpCScript_ReactionsDef(char* scriptStr){

  TCCState *s = NULL;

  stringstream pathStream;

  pathStream.clear();
  pathStream.str("");

  // @todo (zhail#1#):
  //The current version of codeblock is bug at the string macro,it can not reconigze the string macro during making command,
  //so the TCC path is putted in program tempory, when public these project, we should move it to make file.
  //The way to move it to makefile please note the TestTCC demo
  //pathStream<<LIBTCCPATH;
  pathStream<<LIBTCCPATH;
  cout<< pathStream.str().c_str()<<endl;

  s = tcc_new();

  tcc_add_include_path(s, pathStream.str().c_str());

  tcc_set_lib_path(s, pathStream.str().c_str());


  pathStream.clear();
  pathStream.str("");

  // @todo (zhail#1#):
  //The current version of codeblock is bug at the string macro,it can not reconigze the string macro during making command,
  //so the TCC path is putted in program tempory, when public these project, we should move it to make file.
  //The way to move it to makefile please note the TestTCC demo
  //pathStream<<DIFFUSORDEFPATH;
  pathStream<<DIFFUSORDEFPATH;
  tcc_add_include_path(s,pathStream.str().c_str());

  if(!s){
    cout<<"MCPSCUERROR: The TCC interperation init failed!"<<endl;
    return -1;
  }

  tcc_set_output_type(s, TCC_OUTPUT_MEMORY);


  cout<<scriptStr<<endl;

  if(tcc_compile_string(s, scriptStr) == -1){
    cout<<scriptStr<<endl;
    cout<<"MCPSCUERROR: The above script compille failed !"<<endl;
    return -1;
  }

  tcc_add_symbol(s,"InitReactionsList",(void*)InitReactionsList);
  tcc_add_symbol(s,"AppendReaction",(void*)AppendReaction);

  if(tcc_relocate(s, TCC_RELOCATE_AUTO) < 0)
     return 1;

  SetReaction p_SetReaction = (SetReaction)tcc_get_symbol(s,"SetReaction");

  InitReactionsList(&ReactionsList);

  if(!p_SetReaction){
    cout<<scriptStr<<endl;
    cout<<"MCPSCUERROR: The above script interperate failed!"<<endl;
    return -1;
  }else{
    p_SetReaction(ReactionsList);
  }

  return ReactionsList->size;
}
#else
int InterpCScript_ReactionsDef(char* scriptStr){

  TCCState *s = NULL;
  char path[100];

  // @todo (zhail#1#):
  //The current version of codeblock is bug at the string macro,it can not reconigze the string macro during making command,
  //so the TCC path is putted in program tempory, when public these project, we should move it to make file.
  //The way to move it to makefile please note the TestTCC demo
  //pathStream<<LIBTCCPATH;
  strcpy(path,LIBTCCPATH);

  s = tcc_new();

  tcc_add_include_path(s, path);
  tcc_set_lib_path(s, path);

  memset(path,0,sizeof(path));
  strcpy(path,DIFFUSORDEFPATH);
  tcc_add_include_path(s, path);

  if(!s){
    printf("%s\n","MCPSCUERROR: The TCC interperation init failed!");
    return -1;
  }

  tcc_set_output_type(s, TCC_OUTPUT_MEMORY);

  if(tcc_compile_string(s, scriptStr) == -1){
    printf("%s\n",scriptStr);
    printf("%s\n","MCPSCUERROR: The above script compille failed !");
    return -1;
  }

  tcc_add_symbol(s,"InitReactionsList",(void*)InitReactionsList);
  tcc_add_symbol(s,"AppendReaction",(void*)AppendReaction);

  if(tcc_relocate(s, TCC_RELOCATE_AUTO) < 0)
     return 1;

  SetReaction p_SetReaction = (SetReaction)tcc_get_symbol(s,"SetReaction");

  InitReactionsList(&ReactionsList);

  if(!p_SetReaction){
    printf("%s\n",scriptStr);
    printf("%s\n","MCPSCUERROR: The above script interperate failed!");
    return -1;
  }else{
    p_SetReaction(ReactionsList);
  }

  return ReactionsList->size;
}
#endif

/*
--------------------------------------------------------
*/
void GetInterpedReactionsArray(CReactionDef* theArray){

  int I;

  int theSize = ReactionsList->size;

  CReactionsList* cursor;

  cursor = ReactionsList;

  for(I=0;I<theSize;I++){
    copyCReactionDef(&(theArray[I]),&(cursor->data));
    cursor = cursor->next;
  }

  return;
}
