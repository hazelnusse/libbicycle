#include <cmath>
#include "whipple.h"
#include "gtest/gtest.h"

TEST(WhippleTest, DefaultConstructor)
{
  Whipple a;
  EXPECT_EQ(a.w, 1.02);
  EXPECT_EQ(a.c, 0.08);
  EXPECT_EQ(a.lambda, M_PI/10.0);
  EXPECT_EQ(a.g, 9.81);
  EXPECT_EQ(a.rR, 0.3);
  EXPECT_EQ(a.mR, 2.0);
  EXPECT_EQ(a.IRxx, 0.0603);
  EXPECT_EQ(a.IRyy, 0.12);
  EXPECT_EQ(a.xB, 0.3);
  EXPECT_EQ(a.zB, -0.9);
  EXPECT_EQ(a.mB, 85.0);
  EXPECT_EQ(a.IBxx, 9.2);
  EXPECT_EQ(a.IByy, 11.0);
  EXPECT_EQ(a.IBzz, 2.8);
  EXPECT_EQ(a.IBxz, 2.4);
  EXPECT_EQ(a.xH, 0.9);
  EXPECT_EQ(a.zH, -0.7);
  EXPECT_EQ(a.mH, 4.0);
  EXPECT_EQ(a.IHxx, 0.05892);
  EXPECT_EQ(a.IHyy, 0.06);
  EXPECT_EQ(a.IHzz, 0.00708);
  EXPECT_EQ(a.IHxz, -0.00756);
  EXPECT_EQ(a.rF, 0.35);
  EXPECT_EQ(a.mF, 3.0);
  EXPECT_EQ(a.IFxx, 0.1405);
  EXPECT_EQ(a.IFyy, 0.28);
}
