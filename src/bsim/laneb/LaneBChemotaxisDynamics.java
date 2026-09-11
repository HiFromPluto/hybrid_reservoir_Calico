package bsim.laneb;

import bsim.BSimRandom;
import bsim.ode.BSimOdeSystem;

/**
 * Barkai–Leibler 3-state chemotaxis ODE ported from paper-1 ChemotaxisDynamics.
 * Ligand L = attractant − repellent. AHL does not enter. B0_MOTILITY_MEMORY.
 */
public final class LaneBChemotaxisDynamics implements BSimOdeSystem {

    public double ligandConcentration = 0.0;

    private final BSimRandom rng;

    public LaneBChemotaxisDynamics(BSimRandom rng) {
        if (rng == null) throw new IllegalArgumentException("rng is required");
        this.rng = rng;
    }

    @Override
    public double[] derivativeSystem(double t, double[] y) {
        double m = y[0];
        double a = y[1];
        double cheyP = y[2];
        double L = ligandConcentration;
        double aInf = 1.0 / (1.0 + Math.exp(LaneBIdentity.ALPHA * (L - m)));
        return new double[] {
                LaneBIdentity.K_R * (1.0 - a) - LaneBIdentity.K_B * a,
                (aInf - a) / LaneBIdentity.TAU_A,
                LaneBIdentity.K_A * a - LaneBIdentity.GAMMA_Y * cheyP
        };
    }

    @Override
    public int getNumEq() {
        return 3;
    }

    @Override
    public double[] getICs() {
        double dm = 0.10 * (rng.nextDouble() - 0.5);
        double da = 0.10 * (rng.nextDouble() - 0.5);
        return new double[] {
                Math.max(0.0, 0.0 + dm),
                Math.max(0.0, Math.min(1.0, 0.5 + da)),
                LaneBIdentity.Y0
        };
    }
}
