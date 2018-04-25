#include "ar_he_pes.hpp"

double ar_he_pot( const double & R )
{
    double t2, t5, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t25, t26, t27, t28,
        t37, t38, t39, t40, t48, t49, t50, t51, t59, t60, t61, t62, t71, t73, t80;

    t2 = R * R;
    t5 = exp(-0.163051044e1 * R - 0.4680016e-1 * t2);
    t7 = 0.16274792500000e1 * R;
    t8 = exp(-t7);
    t9 = 0.13243443545903e1 * t2;
    t10 = t2 * R;
    t11 = 0.71844765231677e0 * t10;
    t12 = t2 * t2;
    t13 = 0.29231466158919e0 * t12;
    t14 = t12 * R;
    t15 = 0.95147209241436e-1 * t14;
    t16 = t12 * t2;
    t17 = 0.25808351455974e-1 * t16;
    t25 = t12 * t10;
    t26 = 0.60003652101865e-2 * t25;
    t27 = t12 * t12;
    t28 = 0.12206837340001e-2 * t27;
    t37 = 0.22073749421084e-3 * t27 * R;
    t38 = t27 * t2;
    t39 = 0.35924569152514e-4 * t38;
    t40 = 0.1e1 + t7 + t9 + t11 + t13 + t15 + t17 + t26 + t28 + t37 + t39;
    t48 = 0.53151355328098e-5 * t27 * t10;
    t49 = t27 * t12;
    t50 = 0.72085606588213e-6 * t49;
    t51 = 0.1e1 + t7 + t9 + t11 + t13 + t15 + t17 + t26 + t28 + t37 + t39 + t48 + t50;
    t59 = 0.90244483804600e-7 * t27 * t14;
    t60 = t27 * t16;
    t61 = 0.10490787487068e-7 * t60;
    t62 = 0.1e1 + t7 + t9 + t11 + t13 + t15 + t17 + t26 + t28 + t37 + t39 + t48 + t50 + t59 + t61;
    t71 = t27 * t27;
    t73 = 0.1e1 + t7 + t9 + t11 + t13 + t15 + t17 + t26 + t28 + t37 + t39 + t48 + t50 + t59 + t61 + 0.11382359300908e-8 * t27 * t25 + 0.11577845986420e-9 * t71;
    t80 = 0.2303058634e2 * t5 - 0.940835513e1 * (0.1e1 - 0.1e1 * t8 * (0.1e1 + t7 + t9 + t11 + t13 + t15 + t17)) / t16 - 0.165523018e3 * (0.1e1 - 0.1e1 * t8 * (0.1e1 + t7 + t9 + t11 + t13 + t15 + t17 + t26 + t28)) / t27 - 0.379715796e4 * (0.1e1 - 0.1e1 * t8 * t40) / t38 - 0.1165179999e6 * (0.1e1 - 0.1e1 * t8 * t51) / t49 - 0.466258000e7 * (0.1e1 - 0.1e1 * t8 * t62) / t60 - 0.236861000e9 * (0.1e1 - 0.1e1 * t8 * t73) / t71;

    return t80;
}


