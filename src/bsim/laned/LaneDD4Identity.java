package bsim.laned;

/**
 * Frozen D4_THREE_PHASE identity. Sequential feast / handoff / famine.
 * Does not edit D0. Does not import StarvationViability.
 */
public final class LaneDD4Identity {

    public static final String OBJECT = "LANE_D_LIFE_CYCLE";
    public static final String GATE = "D4_THREE_PHASE";
    public static final String DEVICE = LaneDD0Identity.DEVICE;

    public static final double PHI_RB0 = 0.049;
    public static final double GAMMA_TR_PER_H = 11.02;
    public static final double LAMBDA_C_PER_H = 1.17;

    public static final double GAMMA0_PER_DAY = 0.21;
    public static final double DEATH_SLOPE_H = 1.0;
    public static final double T_FAMINE_D = 4.0;
    public static final double FAMINE_DT_D = 0.01;
    public static final double FAMINE_CLOCK_DELTA = 1.0e-6;

    public static final double D0_T_END_S = 10455.0;
    public static final double D0_C_END_MM = 9.977e-4;
    public static final int D0_N_END = 64;
    public static final double D0_C_END_TOL = 5.0e-7;

    /** D3.2 standalone sigma(T) at lambda_cut=1e-4. Comparison only. */
    public static final double D32_SIGMA_T_COMPARISON = 0.8510519648186876;

    public static final double[] FAMINE_OBS_D = {0.0, 1.0, 2.0, 3.0, 4.0};

    private LaneDD4Identity() { }

    public static double lambdaMonodPerHour(double cMm) {
        return LaneDD0Identity.LAMBDA_S_PER_HOUR
                * LaneDD0Identity.monodFactor(cMm);
    }

    public static double muFcr(double lambdaMonodPerHour) {
        double lam = Math.max(0.0, lambdaMonodPerHour);
        double denom = 1.0 - lam / LAMBDA_C_PER_H;
        if (denom <= 1.0e-15) {
            throw new IllegalStateException("FCR mu_f undefined at lambda -> lambda_C");
        }
        return lam / denom;
    }

    public static double phiRbStar(double lamPerHour) {
        return PHI_RB0 + lamPerHour / GAMMA_TR_PER_H;
    }

    public static double sigmaStar(double lamPerHour) {
        return lamPerHour / phiRbStar(lamPerHour);
    }

    public static double muEffFromPhiRb(double phiRb) {
        return GAMMA_TR_PER_H * (phiRb - PHI_RB0);
    }

    public static double gammaPerDay(double muEffPerHour) {
        return GAMMA0_PER_DAY * Math.exp(DEATH_SLOPE_H * muEffPerHour);
    }
}
