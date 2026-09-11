package BSimReservoirStage11;

import bsim.ode.BSimOdeSystem;

/**
 * AC Internal ODE System -- 2-state leaky biochemical processor.
 * Ported from reservoir_new (AC_OurReservoir_Scenario.ACInternalDynamics).
 *
 *   dx/dt = k_u * u(t) - gamma_x * x
 *   ds/dt = k_x * x    - gamma_s * s
 */
public class ACInternalDynamics implements BSimOdeSystem {

    public double externalInput = 0.0;

    public static double K_U     = 2.0;
    public static double GAMMA_X = 1.0;
    public static double K_X     = 1.5;
    public static double GAMMA_S = 0.5;

    @Override
    public double[] derivativeSystem(double t, double[] y) {
        double[] dy = new double[2];
        dy[0] = K_U * externalInput - GAMMA_X * y[0];
        dy[1] = K_X * y[0]          - GAMMA_S * y[1];
        return dy;
    }

    @Override
    public int getNumEq() { return 2; }

    @Override
    public double[] getICs() { return new double[]{0.0, 0.0}; }
}
