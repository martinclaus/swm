/* Two point interpolation routine */
#if defined(__SX__)
  #define interp2(v, in, i, j) (in%weight_vec(1, i, j) * v(in%iind(1, i, j), in%jind(1, i, j)) \
                              + in%weight_vec(1, i, j) * v(in%iind(1, i, j), in%jind(1, i, j)))
#else
  #define interp2 interpolate
#endif

/* Four point interpolation routine */
#if defined(__SX__)
  #define interp4(v, in, i, j) (in%weight_vec(1, i, j) * v(in%iind(1, i, j), in%jind(1, i, j)) \
                              + in%weight_vec(2, i, j) * v(in%iind(2, i, j), in%jind(2, i, j)) \
                              + in%weight_vec(3, i, j) * v(in%iind(3, i, j), in%jind(3, i, j)) \
                              + in%weight_vec(4, i, j) * v(in%iind(4, i, j), in%jind(4, i, j)))
#else
  #define interp4 interpolate
#endif
