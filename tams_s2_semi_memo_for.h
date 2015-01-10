  /**************************************************************************
  **************************************************************************

                           S2kit 1.0

          A lite version of Spherical Harmonic Transform Kit

   Peter Kostelec, Dan Rockmore
   {geelong,rockmore}@cs.dartmouth.edu

   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu

   Copyright 2004 Peter Kostelec, Dan Rockmore

   This file is part of S2kit.

   S2kit is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   S2kit is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with S2kit; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

   See the accompanying LICENSE file for details.

  ************************************************************************
  ************************************************************************/
/***********************************************************************
  * Author: Bo Sun
  * Afflication: TAMS, University of Hamburg
  * E-Mail: bosun@informatik.uni-hamburg.de
  *         user_mail@QQ.com
  * "tams_s2_semi_memo_for" is modified from "test_s2_semi_memo_for.c"
  * in package s2kit10.
  *
  * 1. Rather than read the real and imaginary part of the input function,
  *    the new function just accept the real part (imaginary part is 0).
  * 2. Rather than write the spherical harmonic coefficients in a file,
  *    the real and imaginary part of spherical harmonic coefficients are
  *    stored in two vectors.
  ***********************************************************************/

#ifndef _TAMS_S2_SEMI_MEMO_H
#define _TAMS_S2_SEMI_MEMO_H

extern void tams_s2_semi_memo_for(  const Eigen::VectorXf &TAMS_sei_real,
                                    int TAMS_bandwidth_,
                                    std::vector<double> &TAMS_sh_real,
                                    std::vector<double> &TAMS_sh_imag);

#endif /* _TAMS_S2_SEMI_MEMO_H */
