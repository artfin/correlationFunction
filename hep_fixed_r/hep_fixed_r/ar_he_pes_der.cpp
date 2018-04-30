#include "math.h"
#include "ar_he_pes_der.hpp"

double ar_he_pot_der(double R)
{
    double t1, t3, t5, t7, t9, t12, t18, t19, t20, t27, t32, t37, t38, t70, t75;
    t1 = R * R;
    t3 = t1 * t1;
    t5 = t3 * t1;
    t7 = t3 * t3;
    t9 = t7 * t1;
    t12 = exp(-0.16274792500000e1 * R);
    t18 = exp(-0.40000000000000e-7 * R * (0.40762761e8 + 0.1170004e7 * R));
    t19 = t7 * t7;
    t20 = t19 * R;
    t27 = t1 * R;
    t32 = t3 * R;
    t37 = -0.65276120000000e21 * t1 - 0.13982159988000e20 * t3 - 0.37971579600000e18 * t5 - 0.13241841440000e17 * t7 - 0.56450130780000e15 * t9 + 0.37897760000000e23 * t12 + 0.37551611466691e15 * t18 * t20 + 0.61677818021480e23 * t12 * R + 0.50842445707618e23 * t12 * t1 + 0.28289912008269e23 * t12 * t27 + 0.11956533659492e23 * t12 * t3 + 0.40975965274158e22 * t12 * t32 + 0.11877872891918e22 * t12 * t5;
    t38 = t3 * t27;
    t70 = -0.37897760000000e23 + 0.30017223672567e21 * t12 * t38 + 0.67711176330861e20 * t12 * t7 + 0.13906982420029e20 * t12 * t7 * R + 0.26482299316268e19 * t12 * t9 + 0.47597985110546e18 * t12 * t7 * t27 + 0.82255007548410e17 * t12 * t7 * t3 + 0.13919887346146e17 * t12 * t7 * t32 + 0.23407022926629e16 * t12 * t7 * t5 + 0.39334588708829e15 * t12 * t7 * t38 + 0.65688923130369e14 * t12 * t19 + 0.12069498412857e14 * t12 * t20 + 0.21556702512116e14 * t18 * t19 * t1;
    t75 = -0.10000000000000e-12 * (t37 + t70) / t20;
    return t75;
}
