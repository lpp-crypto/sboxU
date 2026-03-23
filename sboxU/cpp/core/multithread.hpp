#ifndef _MULTITHREAD_
#define _MULTITHREAD_
#include<omp.h>

inline void set_max_threads(int n)
{
    omp_set_num_threads(n);
}

inline int get_max_threads()
{
    return  omp_get_max_threads();
}

inline int threads_from_size(BinWord size)
{
    return std::min(omp_get_max_threads(), std::max(1, ((int) size) >> 4));
}

#endif
