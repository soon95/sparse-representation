
TwIST 2.0

October, 1, 2007

---------------------------------------------------------------------
Copyright (2010): José M. Bioucas-Dias and Mário Figueiredo

TwIST is distributed under the terms of the GNU General Public
License 2.0.

Permission to use, copy, modify, and distribute this software for
any purpose without fee is hereby granted, provided that this entire
notice is included in all copies of any software which is or includes
a copy or modification of this software and in all copies of the
supporting documentation for such software.
This software is being provided "as is", without any express or
implied warranty.  In particular, the authors do not make any
representation or warranty of any kind concerning the merchantability
of this software or its fitness for any particular purpose."
---------------------------------------------------------------------

This set of MATLAB files contain an implementation of the algorithms
described in the papers

J. Bioucas-Dias, M. Figueiredo, "A new TwIST: two-step iterative
shrinkage/thresholding  algorithms for image restoration”, IEEE
Transactions on Image Processing, December 2007.
and

J. Bioucas-Dias and M. Figueiredo , "Two-step algorithms for
linear inverse problems with non-quadratic regularization",
IEEE International Conference on Image Processing -ICIP'2007,
San Antonio, TX, USA, September 2007.

This algorithm solves the otimization problem

min_x \| y - A x \|_2^2 + lambda Phi(x)

and has applications in compressed sensing, image restoration
and reconstruction, sparse regression, and several other problems.

Both papers as well as the latest version of the code are
available at http://www.lx.it.pt/~bioucas/TwIST/

---------------------------------------------------------------------

The main file implementing the algorithm is TwIST.m

For usage details, type  "help TwIST"
at the MATLAB prompt.

This package includes 8 demos.

The demo "demo_Piecewise_cubic_polynomial.m" uses the
SPARCO toolbox, available at
http://www.cs.ubc.ca/labs/scl/sparco/
and the Rice wavelet toolbox, available at
http://www-dsp.rice.edu/software/rwt.shtml


The demo "demo_l2_TV" uses the total variation
denoising script tvdenoise coded by Pascal Getreuer 
and modified by us. 

Contacts: bioucas@lx.it.pt
          mario.figueiredo@lx.it.pt

This code is in development stage; any comments or bug reports
are very welcome.
