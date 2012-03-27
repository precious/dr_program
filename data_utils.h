#ifndef DATA_UTILS_H
#define DATA_UTILS_H


class DataUtils {
public:
   // C - container type (vector or pointer to array)
   // T - type of values in container
   template <typename C,typename T>
   static T reduce(function<T& (T,T)> f,C container,int size) {
       T result = container[0];
       for(int i = 1;i < size;i++)
           result = f(result,container[i]);
       return result;
   }
};

#endif // DATA_UTILS_H
