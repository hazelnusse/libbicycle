#include "bicycle.h"
#include "whipple.h"

void Bicycle::set_parameters_from_whipple(const Whipple & w) {
  double * z = new double[61];

  z[0] = sin(w.lambda);
  z[1] = cos(w.lambda);
  z[2] = pow(w.mB, 2);
  z[3] = pow(w.mH, 2);
  z[4] = w.mB + w.mR;
  z[5] = w.mF + w.mH;
  z[6] = w.rF + w.tF;
  z[7] = w.rR + w.tR;
  z[8] = w.zH + z[6];
  z[9] = w.zB + z[7];
  z[10] = pow(z[0], 2);
  z[11] = pow(z[1], 2);
  z[12] = pow(z[4], 2);
  z[13] = pow(z[5], 2);
  z[14] = w.w - w.xH;
  z[15] = w.xB*z[1] - z[0]*z[9];
  z[16] = -w.xB*z[0] - z[1]*z[9];
  z[17] = w.mB*w.xB*z[0] + w.mB*z[1]*z[9];
  z[18] = w.mB*w.xB*z[1] - w.mB*z[0]*z[9];
  z[19] = pow(z[17], 2);
  z[20] = z[0]*z[14] - z[1]*z[8];
  z[21] = pow(z[18], 2);
  z[22] = -z[0]*z[8] - z[1]*z[14];
  z[23] = -w.mH*z[0]*z[14] + w.mH*z[1]*z[8];
  z[24] = w.mH*z[22];
  z[25] = -w.mH*z[0]*z[8] - w.mH*z[1]*z[14];
  z[26] = -w.mH*z[0]*z[8] - w.mH*z[1]*z[14];
  z[27] = -w.mH*z[20];
  z[28] = z[23];
  z[29] = pow(z[25], 2);
  z[30] = pow(z[28], 2);
  z[31] = -z[15] + z[18]/z[4];
  z[32] = -pow(w.xB, 2) + w.xB*(z[0]*z[17]/z[4] + z[1]*z[18]/z[4]) + z[15]*z[18]/z[4] - z[16]*z[17]/z[4] - pow(z[9], 2) + z[9]*(-z[0]*z[18]/z[4] + z[1]*z[17]/z[4]);
  z[33] = pow(z[14], 2) - z[14]*(-z[0]*z[23]/z[5] - z[1]*z[26]/z[5]) + z[20]*z[28]/z[5] - z[22]*z[25]/z[5] + pow(z[8], 2) + z[8]*(z[0]*z[26]/z[5] - z[1]*z[23]/z[5]);
  z[34] = -w.IBxz;
  z[35] = -w.IHxz;
  z[36] = -w.mH;
  z[37] = w.xB*z[0];
  z[38] = w.xB*z[1];
  z[39] = -z[0];
  z[40] = -z[1]*z[39];
  z[41] = w.mB*z[9];
  z[42] = -z[36]*z[8];
  z[43] = w.mB/z[4];
  z[44] = w.mF/z[13];
  z[45] = -z[14]*z[36];
  z[46] = -1/z[5];
  z[47] = w.mR*z[2]/z[12];
  z[48] = z[17]/z[4];
  z[49] = -z[16]*z[2]/z[4];
  z[50] = z[19]/z[12];
  z[51] = z[15]*z[2]/z[4];
  z[52] = z[21]/z[12];
  z[53] = -z[24]*z[46];
  z[54] = -z[26]*z[46];
  z[55] = -z[27]*z[46];
  z[56] = z[20]*z[3]*z[46];
  z[57] = z[29]/z[13];
  z[58] = z[30]/z[13];
  z[59] = w.mB*z[31];
  z[60] = -z[20]*z[22];

  rear_.Ixx = w.IBxx*z[11] + w.IBzz*z[10] + w.IRxx + w.mB*(-z[32] + z[50]) + pow(z[16], 2)*z[47] + z[31]*z[39]*z[41] + 2*z[34]*z[40] + z[38]*z[51] + z[38]*z[59] + z[39]*z[51]*z[9];
  rear_.Iyy = w.IByy + w.IRyy + w.mB*(-z[32] + z[50] + z[52]) + w.mR*(z[19] + z[21])/z[12];
  rear_.Izz = w.IBxx*z[10] + w.IBzz*z[11] + w.IRxx - w.mB*z[37]*(-z[16] - z[48]) + w.mB*(-z[32] + z[52]) - z[1]*z[41]*(-z[16] - z[48]) + z[1]*z[49]*z[9] + pow(z[15], 2)*z[47] - 2*z[34]*z[40] + z[37]*z[49];
  rear_.Ixz = w.IBxx*z[40] - w.IBzz*z[40] + pow(w.mB, 3)*z[15]*z[16]/z[12] + z[1]*z[31]*z[41] + z[10]*z[34] - z[11]*z[34] + z[15]*z[16]*z[47] + z[37]*z[59] + z[38]*z[49] + z[39]*z[49]*z[9];
  rear_.J = w.IRyy;
  rear_.m = z[4];
  rear_.R = w.rR;
  rear_.r = w.tR;
  rear_.a = z[15]*z[43];
  rear_.b = -z[16]*z[43];
  rear_.c = z[1]*(w.c + w.w) + z[39]*z[7];
  front_.Ixx = w.IFxx + w.IHxx*z[11] + w.IHzz*z[10] - z[1]*z[45]*z[53] - z[1]*z[45]*(-z[22] + z[54]) + pow(z[27], 2)*z[44] + 2*z[35]*z[40] - z[36]*(z[33] + z[58]) + z[39]*z[42]*z[53] + z[39]*z[42]*(-z[22] - z[25]*z[46]);
  front_.Iyy = w.IFyy + w.IHyy - z[36]*(z[33] + z[57] + z[58]) + z[44]*(z[29] + z[30]);
  front_.Izz = w.IFxx + w.IHxx*z[10] + w.IHzz*z[11] + z[1]*z[42]*z[55] + z[1]*z[42]*(z[20] - z[28]*z[46]) + pow(z[24], 2)*z[44] - 2*z[35]*z[40] - z[36]*(z[33] + z[57]) + z[39]*z[45]*z[55] + z[39]*z[45]*(z[20] - z[23]*z[46]);
  front_.Ixz = w.IHxx*z[40] - w.IHzz*z[40] - z[1]*z[14]*z[56] + z[1]*z[42]*(-z[22] + z[54]) + z[10]*z[35] - z[11]*z[35] - z[3]*z[44]*z[60] - z[39]*z[45]*(z[22] - z[54]) + z[39]*z[56]*z[8] + pow(z[36], 3)*z[60]/z[13];
  front_.J = w.IFyy;
  front_.m = z[5];
  front_.R = w.rF;
  front_.r = w.tF;
  front_.a = z[53];
  front_.b = z[55];
  front_.c = w.c*z[1] + z[39]*z[6];
  ls_ = -w.w*z[39] - z[1]*(z[6] - z[7]);

  delete [] z;
}

