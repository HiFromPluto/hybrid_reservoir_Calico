package bsim.laned;

/**
 * Frozen D5_GATES identity. Composition gates, not a new dish.
 * Does not edit D0/D4 feast. Does not import StarvationViability.
 */
public final class LaneDD5Identity {

    public static final String OBJECT = "LANE_D_LIFE_CYCLE";
    public static final String GATE = "D5_GATES";
    public static final String DEVICE = LaneDD0Identity.DEVICE;

    public static final double[] SWEEP_MU = {
        0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.89
    };

    public static final double[] EV1_MU = {0.10, 0.30, 0.50, 0.70};
    public static final double[] EV1_GAMMA = {0.24, 0.29, 0.33, 0.40};
    public static final double[] EV1_SD = {0.02, 0.03, 0.03, 0.04};

    public static final double T_PLUS_H = 3.0;
    public static final double T_MINUS_D = 6.0;
    public static final double OLD_T_MINUS_D = 3.0;
    public static final double OPT_GRID = 1.0e-4;
    public static final double MU_STAR_LOCKED = 0.8675;
    public static final double GAMMA_STAR_LOCKED = 0.5000;

    public static final int N0_FAMINE = 64;

    private LaneDD5Identity() { }

    public static double phiRbFromMu(double muPerHour) {
        return LaneDD4Identity.PHI_RB0
                + muPerHour / LaneDD4Identity.GAMMA_TR_PER_H;
    }

    public static double optimumMu(double feastHours, double starvationDays) {
        double ratio = feastHours
                / (LaneDD4Identity.DEATH_SLOPE_H
                * LaneDD4Identity.GAMMA0_PER_DAY
                * starvationDays);
        return Math.log(ratio) / LaneDD4Identity.DEATH_SLOPE_H;
    }

    public static double numericalOptimum(double feastHours,
            double starvationDays) {
        double bestMu = 0.0;
        double bestF = Double.NEGATIVE_INFINITY;
        for (double mu = 0.0; mu <= 3.0 + 0.5 * OPT_GRID; mu += OPT_GRID) {
            double g = LaneDD4Identity.gammaPerDay(mu);
            double f = mu * feastHours - g * starvationDays;
            if (f > bestF) {
                bestF = f;
                bestMu = mu;
            }
        }
        return bestMu;
    }
}
