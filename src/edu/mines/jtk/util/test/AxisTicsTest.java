/****************************************************************************
Copyright (c) 2004, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package edu.mines.jtk.util.test;

import static java.lang.Math.*;
import junit.framework.*;
import edu.mines.jtk.util.*;

/**
 * Tests {@link edu.mines.jtk.util.AxisTics}.
 * @author Dave Hale, Colorado School of Mines 
 * @version 2004.12.14
 */
public class AxisTicsTest extends TestCase {
  public static void main(String[] args) {
    TestSuite suite = new TestSuite(AxisTicsTest.class);
    junit.textui.TestRunner.run(suite);
  }

  public void testAll() {

    check(0.0, 10.0, 1.0,  11, 1.0,   0.0,  101, 0.1,   0.0);
    check(0.0, 10.0,  11,  11, 1.0,   0.0,  101, 0.1,   0.0);

    check(0.0,-10.0, 1.0,  11, 1.0, -10.0,  101, 0.1, -10.0);
    check(0.0,-10.0,  11,  11, 1.0, -10.0,  101, 0.1, -10.0);

    check(1.0, 10.0, 1.0,  10, 1.0,   1.0,   91, 0.1,   1.0);
    check(1.0, 10.0,  10,  10, 1.0,   1.0,   91, 0.1,   1.0);

    check(1.0, 10.0, 2.0,   5, 2.0,   2.0,   10, 1.0,   1.0);
    check(1.0, 10.0,   9,   5, 2.0,   2.0,   10, 1.0,   1.0);

    check(0.9, 10.1, 2.0,   5, 2.0,   2.0,   10, 1.0,   1.0);
    check(0.9, 10.1,   9,   5, 2.0,   2.0,   10, 1.0,   1.0);
  }

  public void testCount() {
  }

  private void check(
    double xmin, double xmax, double dtic,
    int nmajor, double dmajor, double fmajor,
    int nminor, double dminor, double fminor)
  {
    check(new AxisTics(xmin,xmax,dtic),
      nmajor,dmajor,fmajor,nminor,dminor,fminor);
  }

  private void check(
    double xmin, double xmax, int ntic,
    int nmajor, double dmajor, double fmajor,
    int nminor, double dminor, double fminor)
  {
    check(new AxisTics(xmin,xmax,ntic),
      nmajor,dmajor,fmajor,nminor,dminor,fminor);
  }

  private void check(AxisTics at, 
    int nmajor, double dmajor, double fmajor,
    int nminor, double dminor, double fminor)
  {
    assertEquals(nmajor,at.getCountMajor());
    assertEquals(dmajor,at.getDeltaMajor());
    assertEquals(fmajor,at.getFirstMajor());
    assertEquals(nminor,at.getCountMinor());
    assertEquals(dminor,at.getDeltaMinor());
    assertEquals(fminor,at.getFirstMinor());
  }

  private void assertEquals(double e, double a) {
    double tiny = max(abs(e),abs(a))*100.0*M.DBL_EPSILON;
    assertEquals(e,a,tiny);
  }
}