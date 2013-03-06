#ifndef BENCHMARK_STEADY_TURNS_H
#define BENCHMARK_STEADY_TURNS_H
/** Table 1 of Basu-Mandal 2007.
 *
 * Shows input coordinates and coordinate first time derivatives and the resulting
 * output second time derivatives as calculated by Lagrange and Newton-Euler
 * methods.
 */

double q[] = {0.0,             /*< \f$x\f$ */
              0.0,             /*< \f$y\f$ */
              0.2440472102925, /*< \f$z\f$ */
              0.0,             /*< \f$\theta\f$ */
              0.9501292851472, /*< \f$\psi\f$ */
              3.1257073014894, /*< \f$\phi\f$ */
              0.2311385135743, /*< \f$\psi_f\f$ */
              0.0,             /*< \f$\beta_r\f$ */
              0.0};            /*< \f$\beta_f\f$ */

double q_dot[] = {-2.8069345714545, /*< \f$\dot{x}\f$ */
                  -0.1480982396001, /*< \f$\dot{y}\f$ */
                   0.1058778746261, /*< \f$\dot{z}\f$ */
                   0.7830033527065, /*< \f$\dot{\theta}\f$ */
                   0.6068425835418, /*< \f$\dot{\psi}\f$ */
                  -0.0119185528069, /*< \f$\dot{\phi}\f$ */
                   0.4859824687093, /*< \f$\dot{\psi}_f\f$ */
                   8.9129896614890, /*< \f$\dot{\beta}_r\f$ */
                   8.0133620584155};/*< \f$\dot{\beta}_f\f$ */

double q_dot_dot[] = {-0.5041626315047, /*< \f$\ddot{x}\f$ */
                      -0.3449706619454, /*< \f$\ddot{y}\f$ */
                      -1.4604528332980, /*< \f$\ddot{z}\f$ */
                       0.8353281706379, /*< \f$\ddot{\theta}\f$ */
                      -7.8555281128244, /*< \f$\ddot{\psi}\f$ */
                       0.1205543897884, /*< \f$\ddot{\phi}\f$ */
                      -4.6198904039403, /*< \f$\ddot{\psi}_f\f$ */
                       1.8472554144217, /*< \f$\ddot{\beta}_r\f$ */
                       2.4548072904550};/*< \f$\ddot{\beta}_f\f$ */
#endif

