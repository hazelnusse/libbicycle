#include "bicycle.h"
#include "whipple.h"

namespace bicycle {

void Bicycle::set_parameters(const Whipple & w)
{
  double * z = new double[32];

  z[0] = sin(w.lambda);
  z[1] = cos(w.lambda);
  z[2] = w.mB + w.mR;
  z[3] = w.mF + w.mH;
  z[4] = w.rF + w.tF;
  z[5] = w.rR + w.tR;
  z[6] = w.zB + z[5];
  z[7] = w.w - w.xH;
  z[8] = w.IBxx*z[1] - w.IBxz*z[0];
  z[9] = w.IBxz*z[1] - w.IBzz*z[0];
  z[10] = w.IHxx*z[1] - w.IHxz*z[0];
  z[11] = w.IHxz*z[1] - w.IHzz*z[0];
  z[12] = -w.zH - z[4];
  z[13] = w.xB*z[0] + z[1]*z[6];
  z[14] = w.mB*w.xB/z[2] - w.xB;
  z[15] = w.xB*z[1] - z[0]*z[6];
  z[16] = -z[0]*z[12] + z[1]*z[7];
  z[17] = -w.mB*z[6]/z[2] + z[6];
  z[18] = z[0]*z[7] + z[1]*z[12];
  z[19] = w.mH*z[7]/z[3] - z[7];
  z[20] = pow(-z[0]*z[17] - z[1]*z[14], 2);
  z[21] = pow(-z[0]*z[14] + z[1]*z[17], 2);
  z[22] = z[0]*z[19] + z[1]*(w.mH*z[12]/z[3] - z[12]);
  z[23] = -z[0]*(w.mH*z[12]/z[3] - z[12]) + z[1]*z[19];
  z[24] = pow(z[22], 2);
  z[25] = pow(z[23], 2);
  z[26] = -z[0];
  z[27] = -z[1];
  z[28] = w.mB/z[2];
  z[29] = w.mH/z[3];
  z[30] = pow(w.mB, 2)*w.mR/pow(z[2], 2);
  z[31] = w.mF*pow(w.mH, 2)/pow(z[3], 2);

  rear_.Ixx = w.IRxx + w.mB*z[21] + pow(z[13], 2)*z[30] + z[26]*z[9] - z[27]*z[8];
  rear_.Iyy = w.IByy + w.IRyy + w.mB*(z[20] + z[21]) + w.mR*(pow(w.xB*z[0]*z[28] - z[27]*z[28]*z[6], 2) + pow(w.xB*z[1]*z[28] + z[26]*z[28]*z[6], 2));
  rear_.Izz = w.IRxx + w.mB*z[20] + pow(z[15], 2)*z[30] - z[26]*(-w.IBxx*z[26] + w.IBxz*z[1]) - z[27]*(w.IBxz*z[0] - w.IBzz*z[27]);
  rear_.Ixz = -w.mB*(z[14]*z[26] - z[17]*z[27])*(z[14]*z[27] + z[17]*z[26]) - z[13]*z[15]*z[30] - z[26]*z[8] - z[27]*z[9];
  rear_.J = w.IRyy;
  rear_.m = z[2];
  rear_.R = w.rR;
  rear_.r = w.tR;
  rear_.a = z[15]*z[28];
  rear_.b = z[13]*z[28];
  rear_.c = z[26]*z[5] - z[27]*(w.c + w.w);
  front_.Ixx = w.IFxx + w.mH*z[24] - z[10]*z[27] + z[11]*z[26] + pow(z[18], 2)*z[31];
  front_.Iyy = w.IFyy + w.IHyy + w.mF*(pow(z[12]*z[26]*z[29] - z[27]*z[29]*z[7], 2) + pow(-z[12]*z[27]*z[29] - z[26]*z[29]*z[7], 2)) + w.mH*(z[24] + z[25]);
  front_.Izz = w.IFxx + w.mH*z[25] + pow(z[16], 2)*z[31] - z[26]*(-w.IHxx*z[26] + w.IHxz*z[1]) - z[27]*(w.IHxz*z[0] - w.IHzz*z[27]);
  front_.Ixz = -w.mH*z[22]*z[23] - z[10]*z[26] - z[11]*z[27] - z[16]*z[18]*z[31];
  front_.J = w.IFyy;
  front_.m = z[3];
  front_.R = w.rF;
  front_.r = w.tF;
  front_.a = -z[16]*z[29];
  front_.b = -z[18]*z[29];
  front_.c = -w.c*z[27] + z[26]*z[4];
  ls_ = -w.w*z[26] + z[27]*(z[4] - z[5]);

  delete [] z;
}

} // namespace bicycle
