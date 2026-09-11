package bsim.laned;

/**
 * Frozen D0_CLOSED_BATH identity. Phase-1 closed bath only.
 * Not LC0_OPEN_BATH. Not NutrientField. Not Erickson. Not famine.
 * Numbers cited from ChassisParameters; this class does not import it.
 */
public final class LaneDD0Identity {

    public static final String OBJECT = "LANE_D_LIFE_CYCLE";
    public static final String GATE = "D0_CLOSED_BATH";
    public static final String DEVICE = "LANE_D_CLOSED_BATH";

    public static final int N0 = 8;
    public static final long RNG_SEED = 101L;

    public static final double LAMBDA_S_PER_HOUR = 1.0;
    public static final double LAMBDA_S_PER_S = LAMBDA_S_PER_HOUR / 3600.0;
    public static final double C_S_MM = 0.5;
    public static final double K_S_MM = 0.02;
    public static final double YIELD_GCDW_PER_G = 0.5;
    public static final double RHO_CELL_GCDW_PER_UM3 = 0.137e-12;
    public static final double MM_TO_G_PER_UM3 = 180.0 / 1.0e15 * 1.0e-3;

    public static final double ELL0_UM = 1.0;
    public static final double ELL_DIV_UM = 3.0;
    public static final double RADIUS_UM = 0.5;
    public static final double SIGMA_PER_S = LAMBDA_S_PER_S
            * (Math.log(ELL_DIV_UM / ELL0_UM) / Math.log(2.0));

    public static final double DT_S = 1.0;
    public static final double LOG_DT_S = 10.0;
    public static final double T_CAP_S = 36000.0;
    public static final double C_CUT_MM = 1.0e-3;
    public static final double DELTA_MASS = 0.01;
    public static final int N_EVER_COMPUTE_CAP = 256;

    public static final double BOX_X_UM = 200.0;
    public static final double BOX_Y_UM = 200.0;
    public static final double BOX_Z_UM = 10.0;
    public static final double V_BATH_UM3 = BOX_X_UM * BOX_Y_UM * BOX_Z_UM;
    public static final int GRID_X = 4;
    public static final int GRID_Y = 2;

    /** Spoken F3 note (protocol). Not a mass-gate input. */
    public static final double F3_V_DIV_UM3 = 2.8798;
    public static final double F3_TWO_DAUGHTERS_UM3 = 2.6180;
    public static final double F3_DROP_UM3 = 0.2618;
    public static final double F3_DROP_FRAC = -0.091;

    private LaneDD0Identity() { }

    public static double monodFactor(double nutrientMm) {
        double n = Math.max(0.0, nutrientMm);
        return n / (n + K_S_MM);
    }

    public static double lambdaPerHour(double nutrientMm) {
        return LAMBDA_S_PER_HOUR * monodFactor(nutrientMm);
    }

    /** Capsule volume. Warren rod: π r² ℓ + (4/3) π r³. */
    public static double capsuleVolumeUm3(double lengthUm) {
        double r = RADIUS_UM;
        return Math.PI * r * r * lengthUm + (4.0 / 3.0) * Math.PI * r * r * r;
    }

    public static double founderX(int i) {
        int col = i % GRID_X;
        return (col + 0.5) * (BOX_X_UM / GRID_X);
    }

    public static double founderY(int i) {
        int row = i / GRID_X;
        return (row + 0.5) * (BOX_Y_UM / GRID_Y);
    }
}
