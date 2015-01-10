/***************************************************************************
  * Modified from "test_s2_semi_memo_for.c"
  * in package s2kit10
  **************************************************************************
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
  ************************************************************************
  *
  * Modified by Bo Sun
  * From TAMS, University of Hamburg
  * "tams_s2_semi_memo_for" is modified from "test_s2_semi_memo_for.c"
  * in package s2kit10.
  * 1. Rather than read the real and imaginary part of the input function,
  *    the new function just accept the real part (imaginary part is 0)
  * 2. Rather than store the spherical harmonic coefficients in a file,
  *    the real and imaginary part of spherical harmonic coefficients are
  *    stored in two vectors.
  * Author: Bo Sun
  * Afflication: TAMS, University of Hamburg
  * E-Mail: bosun@informatik.uni-hamburg.de
  *         user_mail@QQ.com
  */


#include <eigen3/Eigen/Dense>
#include <vector>

#include <iostream>

extern "C"{
#include "fftw3.h"
#include "makeweights.h"
#include "cospmls.h"
#include "FST_semi_memo.h"

void tams_s2_semi_memo_for(  const Eigen::VectorXf TAMS_sei_real,
                             int TAMS_bandwidth_,
                             std::vector<double> &TAMS_sh_real,
                             std::vector<double> &TAMS_sh_imag)
{
  int i, bw, size ;
  int cutoff;
  int rank, howmany_rank ;
  double *rdata, *idata ;
  double *rcoeffs, *icoeffs ;
  double *weights ;
  double *seminaive_naive_tablespace, *workspace;
  double **seminaive_naive_table ;
  fftw_plan dctPlan, fftPlan ;
  fftw_iodim dims[1], howmany_dims[1];


  bw = TAMS_bandwidth_;

  /*** ASSUMING WILL SEMINAIVE ALL ORDERS ***/
  cutoff = bw ;
  size = 2*bw;

  /* allocate memory */
  rdata = (double *) malloc(sizeof(double) * (size * size));
  idata = (double *) malloc(sizeof(double) * (size * size));
  rcoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  icoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  weights = (double *) malloc(sizeof(double) * 4 * bw);
  seminaive_naive_tablespace =
    (double *) malloc(sizeof(double) *
              (Reduced_Naive_TableSize(bw,cutoff) +
               Reduced_SpharmonicTableSize(bw,cutoff)));
  workspace = (double *) malloc(sizeof(double) *
                ((8 * (bw*bw)) +
                 (7 * bw)));


  /****
       At this point, check to see if all the memory has been
       allocated. If it has not, there's no point in going further.
  ****/

  if ( (rdata == NULL) || (idata == NULL) ||
       (rcoeffs == NULL) || (icoeffs == NULL) ||
       (seminaive_naive_tablespace == NULL) ||
       (workspace == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }

  /* now precompute the Legendres */
  seminaive_naive_table = SemiNaive_Naive_Pml_Table(bw, cutoff,
                            seminaive_naive_tablespace,
                            workspace);

  /* construct fftw plans */

  /* make DCT plan -> note that I will be using the GURU
     interface to execute these plans within the routines*/

  /* forward DCT */
  dctPlan = fftw_plan_r2r_1d( 2*bw, weights, rdata,
                  FFTW_REDFT10, FFTW_ESTIMATE ) ;

  /*
    fftw "preamble" ;
    note that this plan places the output in a transposed array
  */
  rank = 1 ;
  dims[0].n = 2*bw ;
  dims[0].is = 1 ;
  dims[0].os = 2*bw ;
  howmany_rank = 1 ;
  howmany_dims[0].n = 2*bw ;
  howmany_dims[0].is = 2*bw ;
  howmany_dims[0].os = 1 ;

  /* forward fft */
  fftPlan = fftw_plan_guru_split_dft( rank, dims,
                      howmany_rank, howmany_dims,
                      rdata, idata,
                      workspace, workspace+(4*bw*bw),
                      FFTW_ESTIMATE );


  /* now make the weights */
  makeweights( bw, weights );


  /* now read in samples */
  for(i = 0 ; i < size*size ; i++ )
    {
      /* first the real part of the sample */
      (*(rdata + i)) = TAMS_sei_real(i);
      /* now the imaginary part */
      (*(idata + i)) = 0;
    }
  /* For DEBUG
  for (i=0; i<size*size; i++)
  {
      std::cout << *(rdata+i) << std::endl;
  }
  for (i=0; i<size*size; i++)
  {
      std::cout << *(idata+i) << std::endl;
  }
   * For DEBUG
   */

  /* now do the forward spherical transform */
  FST_semi_memo(rdata, idata,
        rcoeffs, icoeffs,
        bw,
        seminaive_naive_table,
        workspace,
        0,      /* datformat = 0 -> samples are complex, =1 -> samples real*/
        cutoff,
        &dctPlan,
        &fftPlan,
        weights );


   for( i = 0 ; i < bw*bw ; i ++ )
   {
       TAMS_sh_real.push_back(rcoeffs[i]);
       TAMS_sh_imag.push_back(icoeffs[i]);
   }

   if (TAMS_sh_real.size()!=bw*bw |
           TAMS_sh_imag.size()!=bw*bw)
   {
       std::cout << "Some thing wrong in computing SH!\n" << std::endl;
   }

  /* clean up */
  fftw_destroy_plan( fftPlan );
  fftw_destroy_plan( dctPlan );

  free(workspace);
  free(seminaive_naive_table);
  free(seminaive_naive_tablespace);
  free(weights);
  free(icoeffs);
  free(rcoeffs);
  free(idata);
  free(rdata);

  return ;

} // end of function
} // end of extern "C"
