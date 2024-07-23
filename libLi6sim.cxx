// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME libLi6sim
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "src/rootoutput.h"
#include "src/correlations.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *PFragS_Dictionary();
   static void PFragS_TClassManip(TClass*);
   static void *new_PFragS(void *p = nullptr);
   static void *newArray_PFragS(Long_t size, void *p);
   static void delete_PFragS(void *p);
   static void deleteArray_PFragS(void *p);
   static void destruct_PFragS(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PFragS*)
   {
      ::PFragS *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PFragS));
      static ::ROOT::TGenericClassInfo 
         instance("PFragS", "src/rootoutput.h", 23,
                  typeid(::PFragS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PFragS_Dictionary, isa_proxy, 4,
                  sizeof(::PFragS) );
      instance.SetNew(&new_PFragS);
      instance.SetNewArray(&newArray_PFragS);
      instance.SetDelete(&delete_PFragS);
      instance.SetDeleteArray(&deleteArray_PFragS);
      instance.SetDestructor(&destruct_PFragS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PFragS*)
   {
      return GenerateInitInstanceLocal((::PFragS*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PFragS*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PFragS_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::PFragS*)nullptr)->GetClass();
      PFragS_TClassManip(theClass);
   return theClass;
   }

   static void PFragS_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CFragS_Dictionary();
   static void CFragS_TClassManip(TClass*);
   static void *new_CFragS(void *p = nullptr);
   static void *newArray_CFragS(Long_t size, void *p);
   static void delete_CFragS(void *p);
   static void deleteArray_CFragS(void *p);
   static void destruct_CFragS(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CFragS*)
   {
      ::CFragS *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CFragS));
      static ::ROOT::TGenericClassInfo 
         instance("CFragS", "src/rootoutput.h", 47,
                  typeid(::CFragS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CFragS_Dictionary, isa_proxy, 4,
                  sizeof(::CFragS) );
      instance.SetNew(&new_CFragS);
      instance.SetNewArray(&newArray_CFragS);
      instance.SetDelete(&delete_CFragS);
      instance.SetDeleteArray(&deleteArray_CFragS);
      instance.SetDestructor(&destruct_CFragS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CFragS*)
   {
      return GenerateInitInstanceLocal((::CFragS*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::CFragS*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CFragS_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::CFragS*)nullptr)->GetClass();
      CFragS_TClassManip(theClass);
   return theClass;
   }

   static void CFragS_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *SampledValues_Dictionary();
   static void SampledValues_TClassManip(TClass*);
   static void *new_SampledValues(void *p = nullptr);
   static void *newArray_SampledValues(Long_t size, void *p);
   static void delete_SampledValues(void *p);
   static void deleteArray_SampledValues(void *p);
   static void destruct_SampledValues(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SampledValues*)
   {
      ::SampledValues *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SampledValues));
      static ::ROOT::TGenericClassInfo 
         instance("SampledValues", "src/correlations.h", 14,
                  typeid(::SampledValues), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SampledValues_Dictionary, isa_proxy, 4,
                  sizeof(::SampledValues) );
      instance.SetNew(&new_SampledValues);
      instance.SetNewArray(&newArray_SampledValues);
      instance.SetDelete(&delete_SampledValues);
      instance.SetDeleteArray(&deleteArray_SampledValues);
      instance.SetDestructor(&destruct_SampledValues);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SampledValues*)
   {
      return GenerateInitInstanceLocal((::SampledValues*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SampledValues*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SampledValues_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SampledValues*)nullptr)->GetClass();
      SampledValues_TClassManip(theClass);
   return theClass;
   }

   static void SampledValues_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_PFragS(void *p) {
      return  p ? new(p) ::PFragS : new ::PFragS;
   }
   static void *newArray_PFragS(Long_t nElements, void *p) {
      return p ? new(p) ::PFragS[nElements] : new ::PFragS[nElements];
   }
   // Wrapper around operator delete
   static void delete_PFragS(void *p) {
      delete ((::PFragS*)p);
   }
   static void deleteArray_PFragS(void *p) {
      delete [] ((::PFragS*)p);
   }
   static void destruct_PFragS(void *p) {
      typedef ::PFragS current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PFragS

namespace ROOT {
   // Wrappers around operator new
   static void *new_CFragS(void *p) {
      return  p ? new(p) ::CFragS : new ::CFragS;
   }
   static void *newArray_CFragS(Long_t nElements, void *p) {
      return p ? new(p) ::CFragS[nElements] : new ::CFragS[nElements];
   }
   // Wrapper around operator delete
   static void delete_CFragS(void *p) {
      delete ((::CFragS*)p);
   }
   static void deleteArray_CFragS(void *p) {
      delete [] ((::CFragS*)p);
   }
   static void destruct_CFragS(void *p) {
      typedef ::CFragS current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::CFragS

namespace ROOT {
   // Wrappers around operator new
   static void *new_SampledValues(void *p) {
      return  p ? new(p) ::SampledValues : new ::SampledValues;
   }
   static void *newArray_SampledValues(Long_t nElements, void *p) {
      return p ? new(p) ::SampledValues[nElements] : new ::SampledValues[nElements];
   }
   // Wrapper around operator delete
   static void delete_SampledValues(void *p) {
      delete ((::SampledValues*)p);
   }
   static void deleteArray_SampledValues(void *p) {
      delete [] ((::SampledValues*)p);
   }
   static void destruct_SampledValues(void *p) {
      typedef ::SampledValues current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SampledValues

namespace {
  void TriggerDictionaryInitialization_libLi6sim_Impl() {
    static const char* headers[] = {
"src/rootoutput.h",
"src/correlations.h",
nullptr
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/home/Li6Webb/Desktop/Li6Plus2IAS/li6sim/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libLi6sim dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
struct __attribute__((annotate("$clingAutoload$src/rootoutput.h")))  PFragS;
struct __attribute__((annotate("$clingAutoload$src/rootoutput.h")))  CFragS;
struct __attribute__((annotate("$clingAutoload$src/correlations.h")))  SampledValues;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libLi6sim dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "src/rootoutput.h"
#include "src/correlations.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"CFragS", payloadCode, "@",
"PFragS", payloadCode, "@",
"SampledValues", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libLi6sim",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libLi6sim_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libLi6sim_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libLi6sim() {
  TriggerDictionaryInitialization_libLi6sim_Impl();
}
