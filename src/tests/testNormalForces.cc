#include <cmath>
#include "gtest/gtest.h"
#include "bicycle.h"
#include "wheelassemblygyrostat.h"
#include "whipple.h"

TEST(NormalForcesTest, BenchmarkReferenceConfiguration)
{
  using namespace bicycle;
  Bicycle b;      // Initial state is zero, all parameters zero, 9.81 gravity
  Whipple w;
  b.set_parameters(w);
  b.solve_configuration_constraint_and_set_state();
  Vector cf = b.steady_no_slip_constraint_forces();

  double rear_normal = -w.g*(w.w * w.mR + (w.w - w.xB) * w.mB + (w.w - w.xH) * w.mH)/w.w;
  double front_normal = -w.g*(w.xB * w.mB + w.xH * w.mH + w.w * w.mF)/w.w;

  EXPECT_DOUBLE_EQ(rear_normal + front_normal, -w.g*(w.mR + w.mB + w.mH + w.mF));
  EXPECT_DOUBLE_EQ(rear_normal + front_normal, cf[2] + cf[5]);
  EXPECT_DOUBLE_EQ(rear_normal, cf[2]);
  EXPECT_DOUBLE_EQ(front_normal, cf[5]);

  EXPECT_NEAR(0.0, cf[0], 1e-14);     // Longitudinal force
  EXPECT_NEAR(0.0, cf[1], 1e-14);     // Lateral force
  EXPECT_NEAR(0.0, cf[3], 1e-14);     // Longitudinal force
  EXPECT_NEAR(0.0, cf[4], 1e-14);     // Lateral force
}

TEST(NormalForcesTest, VerticalSteerAxis)
{
  using namespace bicycle;
  Bicycle b;      // Initial state is zero, all parameters zero, 9.81 gravity
  WheelAssemblyGyrostat r, f;
  r.m = f.m = 1.0; // mass
  r.R = f.R = 1.0; // wheel radius
  r.c = 0.5;       // distance to steer axis
  f.c = -0.5;      // distance to steer axis
  b.set_parameters(r, f, 0.0, 1.0);
  b.solve_configuration_constraint_and_set_state();
  Vector cf = b.steady_no_slip_constraint_forces();
  EXPECT_DOUBLE_EQ(cf[2], -1.0);
  EXPECT_DOUBLE_EQ(cf[5], -1.0);
  EXPECT_NEAR(0.0, cf[0], 1e-14);     // Longitudinal force
  EXPECT_NEAR(0.0, cf[1], 1e-14);     // Lateral force
  EXPECT_NEAR(0.0, cf[3], 1e-14);     // Longitudinal force
  EXPECT_NEAR(0.0, cf[4], 1e-14);     // Lateral force

  r.b = 0.25;
  f.b = -0.25;
  b.set_parameters(r, f, 0.0, 1.0);
  b.solve_configuration_constraint_and_set_state();
  cf = b.steady_no_slip_constraint_forces();
  EXPECT_DOUBLE_EQ(cf[2], -1.0);
  EXPECT_DOUBLE_EQ(cf[5], -1.0);
  EXPECT_NEAR(0.0, cf[0], 1e-14);     // Longitudinal force
  EXPECT_NEAR(0.0, cf[1], 1e-14);     // Lateral force
  EXPECT_NEAR(0.0, cf[3], 1e-14);     // Longitudinal force
  EXPECT_NEAR(0.0, cf[4], 1e-14);     // Lateral force

  // Move the mass centers in opposite directions along x-axis
  r.b = f.b = 0;
  r.a = .25;
  f.a = -.25;
  b.set_parameters(r, f, 0.0, 1.0);
  b.solve_configuration_constraint_and_set_state();
  cf = b.steady_no_slip_constraint_forces();
  EXPECT_DOUBLE_EQ(cf[2], -1.0);
  EXPECT_DOUBLE_EQ(cf[5], -1.0);
  EXPECT_NEAR(0.0, cf[0], 1e-14);     // Longitudinal force
  EXPECT_NEAR(0.0, cf[1], 1e-14);     // Lateral force
  EXPECT_NEAR(0.0, cf[3], 1e-14);     // Longitudinal force
  EXPECT_NEAR(0.0, cf[4], 1e-14);     // Lateral force

  // Vertical steer axis with smaller front wheel
  r.R = 1.0;
  f.R = 0.5;
  double ls = 0.5;
  r.a = r.b = f.a = f.b = 0;
  r.c = 0.5;
  f.c = -0.5;
  b.set_parameters(r, f, ls, 1.0);
  b.solve_configuration_constraint_and_set_state();
  cf = b.steady_no_slip_constraint_forces();
  EXPECT_DOUBLE_EQ(cf[2], -1.0);
  EXPECT_DOUBLE_EQ(cf[5], -1.0);
  EXPECT_NEAR(0.0, cf[0], 1e-14);     // Longitudinal force
  EXPECT_NEAR(0.0, cf[1], 1e-14);     // Lateral force
  EXPECT_NEAR(0.0, cf[3], 1e-14);     // Longitudinal force
  EXPECT_NEAR(0.0, cf[4], 1e-14);     // Lateral force

  // 45 degree head tube
  r.R = 1.0;
  f.R = 1.0;
  ls = 1.0;
  r.a = r.b = f.a = f.b = 0;
  r.c = 0.5;
  f.c = -0.5;
  b.set_parameters(r, f, ls, 1.0);
  b.solve_configuration_constraint_and_set_state();
  cf = b.steady_no_slip_constraint_forces();
  EXPECT_DOUBLE_EQ(cf[2], -1.0);
  EXPECT_DOUBLE_EQ(cf[5], -1.0);
  EXPECT_NEAR(0.0, cf[0], 1e-14);     // Longitudinal force
  EXPECT_NEAR(0.0, cf[1], 1e-14);     // Lateral force
  EXPECT_NEAR(0.0, cf[3], 1e-14);     // Longitudinal force
  EXPECT_NEAR(0.0, cf[4], 1e-14);     // Lateral force
  EXPECT_DOUBLE_EQ(M_PI/4.0, b.coordinate(2)); // Pitch
}

