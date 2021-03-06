#include <cmath>
#include <Eigen/Dense>
#include "bicycle.h"
#include "whipple.h"
#include "gtest/gtest.h"
#include "benchmark_steady_turns.h"

TEST(KinematicsTest, LambdaEqualsZero)
{
  using namespace bicycle;
  Bicycle b;      // Initial state is zero, all parameters zero, 9.81 gravity
  Whipple w;
  w.lambda = 0.0;
  w.w = 1.1;
  w.c = -0.1;
  w.mR = 0.0;
  w.mF = 0.0;
  w.rR = 0.1;
  w.xB = 0.2;
  w.zB = -0.3;
  w.rF = 0.2;
  w.xH = 1.3;
  w.zH = -0.4;
  b.set_parameters(w);
  b.solve_configuration_constraint_and_set_state();
  WheelAssemblyGyrostat rear, front;
  rear = b.rear_parameters();
  front = b.front_parameters();
  Matrix p = b.points_of_interest();

  EXPECT_DOUBLE_EQ(rear.a, w.xB);
  EXPECT_DOUBLE_EQ(rear.b, w.zB + w.rR);
  EXPECT_DOUBLE_EQ(rear.c, w.w + w.c);
  EXPECT_DOUBLE_EQ(front.a, w.xH - w.w);
  EXPECT_DOUBLE_EQ(front.b, w.zH + w.rF);

  EXPECT_DOUBLE_EQ(p(0, 0), 0.0);    // x of rear wheel center
  EXPECT_DOUBLE_EQ(p(0, 1), 0.0);    // y of rear wheel center
  EXPECT_DOUBLE_EQ(p(0, 2), -w.rR);  // z of rear wheel center
  EXPECT_DOUBLE_EQ(p(1, 0), w.mB*w.xB/(w.mB + w.mR));      // x of rear assembly mass center
  EXPECT_DOUBLE_EQ(p(1, 1), 0.0);    // y of rear assembly mass center
  EXPECT_DOUBLE_EQ(p(1, 2), w.mB*w.zB/(w.mB + w.mR));   // z of rear assembly mass center
  EXPECT_DOUBLE_EQ(p(2, 0), rear.c);      // x of steer axis point
  EXPECT_DOUBLE_EQ(p(2, 0), w.w + w.c);   // x of steer axis point
  EXPECT_DOUBLE_EQ(p(2, 1), 0.0);         // y of steer axis point
  EXPECT_DOUBLE_EQ(p(2, 2), -rear.R);     // z of steer axis point
  EXPECT_DOUBLE_EQ(p(3, 0), w.w);         // x of front wheel center
  EXPECT_DOUBLE_EQ(p(3, 1), 0.0);         // y of front wheel center
  EXPECT_DOUBLE_EQ(p(3, 2), -front.R);    // z of front wheel center
  EXPECT_DOUBLE_EQ(p(3, 2), -w.rF);       // z of front wheel center
  EXPECT_DOUBLE_EQ(p(4, 0), w.xH);        // x of front mass center 
  EXPECT_DOUBLE_EQ(p(4, 1), 0.0);         // y of front mass center
  EXPECT_DOUBLE_EQ(p(4, 2), w.zH);        // z of front mass center
  EXPECT_DOUBLE_EQ(p(5, 0), w.w + w.c);   // x of front steer axis point
  EXPECT_DOUBLE_EQ(p(5, 1), 0.0);         // y of front steer axis point
  EXPECT_DOUBLE_EQ(p(5, 2), -w.rF);       // z of front steer axis point
  EXPECT_DOUBLE_EQ(p(6, 0), w.w);         // x of front contact point
  EXPECT_DOUBLE_EQ(p(6, 1), 0.0);         // y of front contact point
  EXPECT_DOUBLE_EQ(p(6, 2), 0.0);         // z of front contact point
}


TEST(KinematicsTest, BenchmarkMassCenterLocations)
{
  using namespace bicycle;
  Bicycle b;
  Whipple w;
  b.set_parameters(w);
  b.solve_configuration_constraint_and_set_state();
  Matrix p = b.points_of_interest();
  
  // Rear wheel center
  double wc_r_x = 0.0,
         wc_r_z = -(w.rR + w.tR);
  EXPECT_DOUBLE_EQ(wc_r_x, p(0, 0));
  EXPECT_DOUBLE_EQ(0.0, p(0, 1));
  EXPECT_DOUBLE_EQ(wc_r_z, p(0, 2));
  // Rear assembly mass center location
  double mc_r_x = w.mB*w.xB/(w.mB + w.mR),
         mc_r_z = (w.mR*(-w.rR - w.tR) + w.mB*w.zB)/(w.mB + w.mR);
  EXPECT_DOUBLE_EQ(mc_r_x, p(1, 0));
  EXPECT_DOUBLE_EQ(0.0, p(1, 1));
  EXPECT_DOUBLE_EQ(mc_r_z, p(1, 2));
  // Front wheel center
  double wc_f_x = w.w,
         wc_f_z = -(w.rF + w.tF);
  EXPECT_DOUBLE_EQ(wc_f_x, p(3, 0));
  EXPECT_DOUBLE_EQ(0.0, p(3, 1));
  EXPECT_DOUBLE_EQ(wc_f_z, p(3, 2));
  // Front assembly mass center location
  double mc_f_x = (w.mF*w.w + w.mH*w.xH)/(w.mF + w.mH),
         mc_f_z = (w.mF*(-w.rF - w.tF) + w.mH*w.zH)/(w.mF + w.mH);
  EXPECT_DOUBLE_EQ(mc_f_x, p(4, 0));
  EXPECT_DOUBLE_EQ(0.0, p(4, 1));
  EXPECT_DOUBLE_EQ(mc_f_z, p(4, 2));
  // Front wheel contact point
  double gc_f_x = w.w,
         gc_f_z = 0.0;
  EXPECT_DOUBLE_EQ(gc_f_x, p(6, 0));
  EXPECT_DOUBLE_EQ(0.0, p(6, 1));
  EXPECT_NEAR(gc_f_z, p(6, 2), 1e-14);

  // Verify reference pitch angle equals lambda
  EXPECT_DOUBLE_EQ(w.lambda, b.coordinate(2));
  // Verify masses
  WheelAssemblyGyrostat r, f;
  r = b.rear_parameters(); f = b.front_parameters();
  EXPECT_DOUBLE_EQ(w.mR + w.mB, r.m);
  EXPECT_DOUBLE_EQ(w.mF + w.mH, f.m);
}

TEST(KinematicsTest, BasuMandalKinematicsTest)
{
  using ::bicycle::Whipple; using ::bicycle::Bicycle;
  using ::Eigen::Map; using ::Eigen::Matrix;
  Whipple w;
  Bicycle b;
  b.set_parameters(w);
  b.set_coordinates_basu_mandal(Map<Matrix<double, 9, 1>>(q));
  b.solve_configuration_constraint_and_set_state(1e-16);
  b.set_speeds_basu_mandal(Map<Matrix<double, 9, 1>>(q_dot));
  b.set_dependent_speeds({0, 2, 4});  // yaw, pitch, rear wheel rates dependent
  b.solve_velocity_constraints_and_set_state();

  EXPECT_DOUBLE_EQ(b.coordinate(0), -q[3]);
  EXPECT_DOUBLE_EQ(b.coordinate(1), M_PI/2.0 - q[4]);
  EXPECT_NEAR(b.coordinate(2), M_PI - q[5] + b.reference_pitch(), 1e-13);
  EXPECT_DOUBLE_EQ(b.coordinate(3), - q[6]);
  EXPECT_DOUBLE_EQ(b.coordinate(4), - q[7]);
  EXPECT_DOUBLE_EQ(b.coordinate(5), - q[8]);
  double z[3], xy[2];
  z[0] = cos(b.coordinate(1));
  z[1] = sqrt(pow(z[0], 2));
  z[2] = b.rear_parameters().R*z[0]*sin(b.coordinate(1));
  xy[0] = (q[0]*z[1] + z[2]*sin(b.coordinate(0)))/z[1];
  xy[1] = (q[1]*z[1] - z[2]*cos(b.coordinate(0)))/z[1];
  EXPECT_DOUBLE_EQ(b.coordinate(6), xy[0]);
  EXPECT_DOUBLE_EQ(b.coordinate(7), xy[1]);

}
