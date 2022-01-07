#ifndef MEMORY_CONTROL_HH
#define MEMORY_CONTROL_HH

#include "cudaUtils.hh"

#include "qs_assert.hh"

namespace MemoryControl
{
   enum AllocationPolicy {HOST_MEM, UVM_MEM, UNDEFINED_POLICY};

   template <typename T>
   T* allocate(const int size, const AllocationPolicy policy)
   {
      if (size == 0) { return NULL;}
      T* tmp = NULL;
      
      switch (policy)
      {
        case AllocationPolicy::HOST_MEM:
         tmp = new T [size];
         break;
        default:
         qs_assert(false);
         break;
      }
      return tmp;
   }

   template <typename T>
   void deallocate(T* data, const int size, const AllocationPolicy policy)
   {
      switch (policy)
      {
        case MemoryControl::AllocationPolicy::HOST_MEM:
         delete[] data; 
         break;
        default:
         qs_assert(false);
         break;
      }
   }
}


#endif
