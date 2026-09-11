package bsim.laned;

/**
 * Well-mixed scalar glucose pool for LANE_D_CLOSED_BATH.
 * Not a NutrientField mode. No PDE. Sink uses rho_cell, not rho_0.
 */
public final class LaneDClosedBath {

    private double cMm;
    private double uMmUm3;
    private boolean sinkOn;

    public LaneDClosedBath(double c0Mm, boolean sinkOn) {
        if (c0Mm < 0.0) {
            throw new IllegalArgumentException("C(0) must be non-negative");
        }
        this.cMm = c0Mm;
        this.uMmUm3 = 0.0;
        this.sinkOn = sinkOn;
    }

    public double concentrationMm() {
        return cMm;
    }

    public double cumulativeSinkMmUm3() {
        return uMmUm3;
    }

    public double deltaCvMmUm3() {
        return (LaneDD0Identity.C_S_MM - cMm) * LaneDD0Identity.V_BATH_UM3;
    }

    public double massResidualRel() {
        double denom = LaneDD0Identity.C_S_MM * LaneDD0Identity.V_BATH_UM3;
        return Math.abs(deltaCvMmUm3() - uMmUm3) / denom;
    }

    /**
     * Death-on while C &gt; C_cut is invalid on domain. This extra never
     * requests it; the throw is the executable refuse.
     */
    public void refuseDeath(boolean deathOn) {
        if (deathOn && cMm > LaneDD0Identity.C_CUT_MM) {
            throw new IllegalStateException(
                    "LANE_D_CLOSED_BATH refuses death-on while C > C_cut");
        }
    }

    /**
     * One Euler sink step. Elongation is not this class.
     * Same (-dC/dt) updates C and U so the bath ledger can telescope.
     */
    public void stepSink(double dtS, double sumVolumeUm3) {
        if (!sinkOn || dtS <= 0.0 || sumVolumeUm3 <= 0.0 || cMm <= 0.0) {
            return;
        }
        double monod = LaneDD0Identity.monodFactor(cMm);
        double lambdaSPerS = LaneDD0Identity.LAMBDA_S_PER_S;
        double massSink = (lambdaSPerS / LaneDD0Identity.YIELD_GCDW_PER_G)
                * LaneDD0Identity.RHO_CELL_GCDW_PER_UM3
                * sumVolumeUm3
                * monod;
        double dCdt = -massSink
                / (LaneDD0Identity.V_BATH_UM3 * LaneDD0Identity.MM_TO_G_PER_UM3);
        uMmUm3 += -dCdt * LaneDD0Identity.V_BATH_UM3 * dtS;
        cMm += dCdt * dtS;
        if (cMm < 0.0) {
            cMm = 0.0;
        }
    }
}
