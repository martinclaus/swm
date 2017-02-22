/* Two point interpolation routine */
#define interp2(v, in, i, j) (in%weight_vec(1, i, j) * v(in%iind(1, i, j), in%jind(1, i, j)) \
                              + in%weight_vec(2, i, j) * v(in%iind(2, i, j), in%jind(2, i, j)))
/* Four point interpolation routine */
#define interp4(v, in, i, j) (in%weight_vec(1, i, j) * v(in%iind(1, i, j), in%jind(1, i, j)) \
                              + in%weight_vec(2, i, j) * v(in%iind(2, i, j), in%jind(2, i, j)) \
                              + in%weight_vec(3, i, j) * v(in%iind(3, i, j), in%jind(3, i, j)) \
                              + in%weight_vec(4, i, j) * v(in%iind(4, i, j), in%jind(4, i, j)))

