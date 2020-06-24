/*
 * p_vec_ten.cc
 *
 *  Created on: Jun 25, 2015
 *      Author: hassasp9
 */
/*
 * p_vec_ten.cc
 *
 *  Created on: Jan 21, 2014
 *      Author: hassasp9
 */
#include "p_vec_ten.hh"

ostream &operator<< (ostream &ostr, const p_ten2 &tensor) {
  // cout operator

    ostr
      << tensor.xx << " " << tensor.xy << " " << tensor.xz << " "
      << tensor.yx << " " << tensor.yy << " " << tensor.yz << " "
      << tensor.zx << " " << tensor.zy << " " << tensor.zz;

    return ostr;
  }








