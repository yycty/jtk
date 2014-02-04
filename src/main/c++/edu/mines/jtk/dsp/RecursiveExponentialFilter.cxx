/*! \file RecursiveExponentialFilter.cxx */

/*
C++ port from Java by Timothy Yue

Please refer to the header file for complete documentation and licensing information.
*/


/******************************************************************************
 *  
 ******************************************************************************
 *
 * Member Functions:
 * RecursiveExponentialFilter::RecursiveExponentialFilter
 * RecursiveExponentialFilter::~RecursiveExponentialFilter
 * RecursiveExponentialFilter::apply
 * RecursiveExponentialFilter::apply1
 * RecursiveExponentialFilter::aFromSigma
 * RecursiveExponentialFilter::smooth1
 * RecursiveExponentialFilter::smooth1Ei
 * RecursiveExponentialFilter::smooth1Eo
 *
 ******************************************************************************
*/

#include "RecursiveExponentialFilter.h"

/**
   * Constructs a filter with specified half-width.
   * The same half-width is used when applying the filter for all 
   * dimensions of multidimensional arrays.
   * @param sigma filter half-width.
   */
RecursiveExponentialFilter::RecursiveExponentialFilter(double sigma)
{
  _sigma1 = (float) sigma;
  _a1     = aFromSigma(sigma);
}

RecursiveExponentialFilter::~RecursiveExponentialFilter()
{
}

 /**
   * Applies this filter.
   * @param x input array.
   * @param y output array.
   */
void RecursiveExponentialFilter::apply (const vector<float>& x, vector<float>& y)
{
  apply1(x, y);
}

/**
  * Applies this filter along the 1st (only) array dimension.
  * Input and output arrays can be the same array.
  * @param x input array.
  * @param y output array.
  */
void RecursiveExponentialFilter::apply1(const vector<float>& x, vector<float>& y) 
{
  smooth1(_ei,_zs,_a1,x,y);
}

float RecursiveExponentialFilter::aFromSigma(double sigma) 
{
  if (sigma<=0.0f)
    return 0.0f;
  double ss = sigma*sigma;
  return (float)((1.0+ss-sqrt(1.0+2.0*ss))/ss);
}

// Smooth a 1D array.
void RecursiveExponentialFilter::smooth1
  (bool ei, bool zs, float a, const vector<float>& x, vector<float>& y) 
{
  if (a==0.0f) {
    y = x;
  } else if (ei) {
    smooth1Ei(zs,a,x,y);
  } else {
    smooth1Eo(zs,a,x,y);
  }
}

// Smooth a 1D array for input boundary conditions.
void RecursiveExponentialFilter::smooth1Ei(bool zs, float a, const vector<float>& x, vector<float>& y)
{
  int n1 = x.size();
  float b = 1.0f-a;
  float sx = zs?1.0f:b;
  float sy = a;
  float yi = y[0] = sx*x[0];
  for (int i1=1; i1<n1-1; ++i1)
    y[i1] = yi = a*yi+b*x[i1];
  sx /= 1.0f+a;
  sy /= 1.0f+a;
  y[n1-1] = yi = sy*yi+sx*x[n1-1];
  for (int i1=n1-2; i1>=0; --i1)
    y[i1] = yi = a*yi+b*y[i1];
}

// Smooth a 1D array for output boundary conditions.
// Adapted from Algorithm 4.1 in Boisvert, R.F., Algorithms for
// special tridiagonal systems: SIAM J. Sci. Stat. Comput., v. 12,
// no. 2, pp. 423-442.
void RecursiveExponentialFilter::smooth1Eo(bool zs, float a, const vector<float>& x, vector<float>& y)
{
  int n1 = x.size();
  float aa = a*a;
  float ss = zs?1.0f-a:1.0f;
  float gg = zs?aa-a:aa;
  float c = (1.0f-aa-ss)/ss;
  float d = 1.0f/(1.0f-aa+gg*(1.0f+c*pow(aa,n1-1)));
  float e = (1.0f-a)*(1.0f-a)*FLT_EPSILON/4.0f;

  // copy scaled input to output
  //mul((1.0f-a)*(1.0f-a),x,y);
  y.resize(x.size());
  for (int i = 0; i < x.size(); i++)
  {
    y[i] = x[i] * (1.0f-a)*(1.0f-a);
  }

  // reversed triangular factorization
  int k1 = min((int)ceil(log(e)/log(a)),2*n1-2); // 1 <= k1 <= 2*n1-2
  float ynm1 = 0.0f;
  int m1 = k1-n1+1; // 2-n1 <= m1 <= n1-1
  for (int i1=m1; i1>0; --i1)
    ynm1 = a*ynm1+y[i1];
  ynm1 *= c;
  if (n1-k1<1)
    ynm1 = a*ynm1+(1.0f+c)*y[0];
  m1 = max(n1-k1,1); // 1 <= m1 <= n1-1
  for (int i1=m1; i1<n1; ++i1)
    ynm1 = a*ynm1+y[i1];
  ynm1 *= d;

  // reverse substitution
  y[n1-1] -= gg*ynm1;
  for (int i1=n1-2; i1>=0; --i1)
    y[i1] += a*y[i1+1];
  y[0] /= ss;

  // forward substitution
  for (int i1=1; i1<n1-1; ++i1)
    y[i1] += a*y[i1-1];
  y[n1-1] = ynm1;
}


