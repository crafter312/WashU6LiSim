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
#include "ROOT/RConfig.hxx"
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
   static TClass *KinematicValues_Dictionary();
   static void KinematicValues_TClassManip(TClass*);
   static void *new_KinematicValues(void *p = nullptr);
   static void *newArray_KinematicValues(Long_t size, void *p);
   static void delete_KinematicValues(void *p);
   static void deleteArray_KinematicValues(void *p);
   static void destruct_KinematicValues(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::KinematicValues*)
   {
      ::KinematicValues *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::KinematicValues));
      static ::ROOT::TGenericClassInfo 
         instance("KinematicValues", "src/frame.h", 14,
                  typeid(::KinematicValues), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &KinematicValues_Dictionary, isa_proxy, 4,
                  sizeof(::KinematicValues) );
      instance.SetNew(&new_KinematicValues);
      instance.SetNewArray(&newArray_KinematicValues);
      instance.SetDelete(&delete_KinematicValues);
      instance.SetDeleteArray(&deleteArray_KinematicValues);
      instance.SetDestructor(&destruct_KinematicValues);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::KinematicValues*)
   {
      return GenerateInitInstanceLocal(static_cast<::KinematicValues*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::KinematicValues*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *KinematicValues_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::KinematicValues*>(nullptr))->GetClass();
      KinematicValues_TClassManip(theClass);
   return theClass;
   }

   static void KinematicValues_TClassManip(TClass* ){
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
         instance("SampledValues", "src/correlations.h", 18,
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
      return GenerateInitInstanceLocal(static_cast<::SampledValues*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::SampledValues*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SampledValues_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::SampledValues*>(nullptr))->GetClass();
      SampledValues_TClassManip(theClass);
   return theClass;
   }

   static void SampledValues_TClassManip(TClass* ){
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
         instance("CFragS", "src/rootoutput.h", 26,
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
      return GenerateInitInstanceLocal(static_cast<::CFragS*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::CFragS*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CFragS_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::CFragS*>(nullptr))->GetClass();
      CFragS_TClassManip(theClass);
   return theClass;
   }

   static void CFragS_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_KinematicValues(void *p) {
      return  p ? new(p) ::KinematicValues : new ::KinematicValues;
   }
   static void *newArray_KinematicValues(Long_t nElements, void *p) {
      return p ? new(p) ::KinematicValues[nElements] : new ::KinematicValues[nElements];
   }
   // Wrapper around operator delete
   static void delete_KinematicValues(void *p) {
      delete (static_cast<::KinematicValues*>(p));
   }
   static void deleteArray_KinematicValues(void *p) {
      delete [] (static_cast<::KinematicValues*>(p));
   }
   static void destruct_KinematicValues(void *p) {
      typedef ::KinematicValues current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::KinematicValues

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
      delete (static_cast<::SampledValues*>(p));
   }
   static void deleteArray_SampledValues(void *p) {
      delete [] (static_cast<::SampledValues*>(p));
   }
   static void destruct_SampledValues(void *p) {
      typedef ::SampledValues current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::SampledValues

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
      delete (static_cast<::CFragS*>(p));
   }
   static void deleteArray_CFragS(void *p) {
      delete [] (static_cast<::CFragS*>(p));
   }
   static void destruct_CFragS(void *p) {
      typedef ::CFragS current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::CFragS

namespace ROOT {
   static TClass *vectorlECFragSgR_Dictionary();
   static void vectorlECFragSgR_TClassManip(TClass*);
   static void *new_vectorlECFragSgR(void *p = nullptr);
   static void *newArray_vectorlECFragSgR(Long_t size, void *p);
   static void delete_vectorlECFragSgR(void *p);
   static void deleteArray_vectorlECFragSgR(void *p);
   static void destruct_vectorlECFragSgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<CFragS>*)
   {
      vector<CFragS> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<CFragS>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<CFragS>", -2, "vector", 339,
                  typeid(vector<CFragS>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlECFragSgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<CFragS>) );
      instance.SetNew(&new_vectorlECFragSgR);
      instance.SetNewArray(&newArray_vectorlECFragSgR);
      instance.SetDelete(&delete_vectorlECFragSgR);
      instance.SetDeleteArray(&deleteArray_vectorlECFragSgR);
      instance.SetDestructor(&destruct_vectorlECFragSgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<CFragS> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<CFragS>","std::vector<CFragS, std::allocator<CFragS> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<CFragS>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlECFragSgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<CFragS>*>(nullptr))->GetClass();
      vectorlECFragSgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlECFragSgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlECFragSgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<CFragS> : new vector<CFragS>;
   }
   static void *newArray_vectorlECFragSgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<CFragS>[nElements] : new vector<CFragS>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlECFragSgR(void *p) {
      delete (static_cast<vector<CFragS>*>(p));
   }
   static void deleteArray_vectorlECFragSgR(void *p) {
      delete [] (static_cast<vector<CFragS>*>(p));
   }
   static void destruct_vectorlECFragSgR(void *p) {
      typedef vector<CFragS> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<CFragS>

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
struct __attribute__((annotate("$clingAutoload$src/rootoutput.h")))  CFragS;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$src/rootoutput.h")))  KinematicValues;
class __attribute__((annotate("$clingAutoload$src/rootoutput.h")))  SampledValues;
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
"KinematicValues", payloadCode, "@",
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
