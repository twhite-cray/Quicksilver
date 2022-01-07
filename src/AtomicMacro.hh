#define ATOMIC_WRITE( x, v ) \
  x = v;

#define ATOMIC_UPDATE( x ) \
  x++;

#define ATOMIC_ADD( x, v ) \
  x += v;

#define ATOMIC_CAPTURE( x, v, p ) \
{p = x; x = x + v;}
