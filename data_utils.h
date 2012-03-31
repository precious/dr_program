#ifndef DATA_UTILS_H
#define DATA_UTILS_H

// methods of this class act similar to the stl algorithms however don't use
// iterators and work only with containers that have operator[] (and with usual C arrays)

class DataUtils {
public:
   // C - container type (vector or pointer to array)
   // T - type of values in container
   template <typename C,typename T>
   static T reduce(function<T& (T&,T&)> f,C container,int size) {
       T result = container[0];
       for(int i = 1;i < size;i++)
           result = f(result,container[i]);
       return result;
   }

   template <typename C,typename T>
   static void map(function<void (T&)> f,C container,int size) {
       for(int i = 0;i < size;i++)
           f(container[i]);
   }
};

typedef DataUtils DU;

#endif // DATA_UTILS_H
