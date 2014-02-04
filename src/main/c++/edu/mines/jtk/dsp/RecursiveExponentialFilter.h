
/*

C++ port from Java by Timothy Yue

The original documentation for RecursiveExponentialFilter.

/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************
package edu.mines.jtk.dsp;

import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Recursive symmetric exponential smoothing filter. Except perhaps near
 * the edges of input and output arrays, the impulse response of this 
 * two-sided filter is symmetric and decays exponentially from its peak 
 * value at zero lag. Specifically, the impulse response has the form
 * h[n] = a^abs(n)*(1-a)/(1+a), where a is a parameter in the range
 * [0:1) derived from a specified half-width sigma.
 * <p>
 * Like the Gaussian filter, the impulse response of the exponential
 * filter is nowhere zero. The half-width sigma for the exponential 
 * filter is here defined so that, for low frequencies, the frequency 
 * response of the exponential filter approximates that for a Gaussian 
 * filter with the same specified half-width sigma. Specifically, the
 * value, slope and curvature of the frequency responses will be the
 * same for exponential and Gaussian filters if the same half-widths
 * are specified.
 * <p>
 * This smoothing filter is faster than a recursive Gaussian filter. 
 * This filter also provides a variety of boundary conditions that
 * can be used to control the filtering of samples near the edges
 * of arrays. For most (but not all) of these boundary conditions, 
 * this filter is symmetric and positive-definite (SPD). This means, 
 * for example, that it can be used as a preconditioner in 
 * conjugate-gradient solutions of SPD systems of equations.
 * <p>
 * Multidimensional filters are applied as a cascade of one-dimensional
 * filters applied for each dimension of multidimensional arrays. In 
 * contrast to the Gaussian filter, this cascade for the exponential 
 * filter does not have isotropic impulse or frequency responses.
 * <p>
 * All smoothing can be performed in place, so that input and output 
 * arrays can be the same array.
 *
 * @author Dave Hale &amp; Simon Luo, Colorado School of Mines
 * @version 2011.10.01
 */

#include <vector>

using namespace std;

class RecursiveExponentialFilter 
{
public:
  enum Edges 
  {
    INPUT_ZERO_VALUE,
    INPUT_ZERO_SLOPE,
    OUTPUT_ZERO_VALUE,
    OUTPUT_ZERO_SLOPE
  };

  RecursiveExponentialFilter(double sigma);
 ~RecursiveExponentialFilter();

  void apply                    (const vector<float>& x, vector<float>& y);

  void apply1                   (const vector<float>& x, vector<float>& y);

private:
  static float aFromSigma       (double sigma);
  static void  smooth1          (bool ei,
                                 bool zs,
                                 float a,
                                 const vector<float>& x,
                                       vector<float>& y);

  static void smooth1Ei         (bool zs,
                                 float a,
                                 const vector<float>& x,
                                       vector<float>& y);

  static void smooth1Eo          (bool zs,
                                  float a,
                                  const vector<float>& x,
                                        vector<float>& y);

private:
  float _sigma1;
  float _a1;
  bool  _ei;
  bool  _zs;

};

#endif
