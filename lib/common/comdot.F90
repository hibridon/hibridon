!comdeck comdot
#if defined(HIB_FPS) || defined(HIB_VAX) || defined(HIB_UNIX) || defined(HIB_CRAY) || defined(HIB_MAC)
data dot /'.'/
#endif
#if defined(HIB_UNIVAC) || defined(HIB_CRAY_COS)
data dot /'$'/
#endif
