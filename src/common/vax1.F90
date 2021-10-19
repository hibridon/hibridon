!comdeck vax1
#if defined(HIB_UNIX_SUN)
      integer shiftl,shiftr,popcnt,poppar,compl
      logical mnpar
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_VAX) || defined(HIB_UNIX_CONVEX) || defined(HIB_UNIX_NEC) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_IBM) || defined(HIB_UNIX_DEC)
      integer xor,or,and,shiftl,shiftr,popcnt,poppar,compl
      logical mnpar
#endif
#if defined(HIB_UNIX_IRIS)
      integer xor,or,and,shiftl,shiftr,popcnt,poppar,compl
      logical mnpar
#endif
