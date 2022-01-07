//Determine which atomics to use based on platform being compiled for
//
//If compiling with CUDA

#ifdef HAVE_OPENMP
    #define USE_OPENMP_ATOMICS
#elif HAVE_OPENMP_TARGET
    #define USE_OPENMP_ATOMICS
#endif


//If in a OpenMP section use the OpenMP atomics
#if defined (USE_OPENMP_ATOMICS)
    #define ATOMIC_WRITE( x, v ) \
        _Pragma("omp atomic write") \
        x = v;

    #define ATOMIC_ADD( x, v ) \
        _Pragma("omp atomic") \
        x += v;

    #define ATOMIC_UPDATE( x ) \
        _Pragma("omp atomic update") \
        x++;

        #define ATOMIC_CAPTURE( x, v, p ) \
            _Pragma("omp atomic capture") \
            {p = x; x = x + v;}

//If in a serial section, no need to use atomics
#else
    #define ATOMIC_WRITE( x, v ) \
        x = v;

    #define ATOMIC_UPDATE( x ) \
        x++;

    #define ATOMIC_ADD( x, v ) \
        x += v;

    #define ATOMIC_CAPTURE( x, v, p ) \
        {p = x; x = x + v;}
#endif
