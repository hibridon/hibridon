!comdeck parlbuf
#if defined(HIB_VAX) || defined(HIB_MAC) || defined(HIB_FPS)
parameter (lbuf=512)
#endif
#if defined(HIB_UNIX_DEC)
parameter (lbuf=1024)
#endif
