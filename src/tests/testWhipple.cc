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
}
