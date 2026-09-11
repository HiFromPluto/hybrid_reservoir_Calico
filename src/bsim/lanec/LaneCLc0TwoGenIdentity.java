package bsim.lanec;

/**
 * Frozen LC0_TWO_GENERATION identity. Replication only. Two generations.
 * Does not rewrite LC0_REPLICATION_VISIBLE. Not Lane A, not Lane B,
 * not P0 packing, not NARMA, not death.
 */
public final class LaneCLc0TwoGenIdentity {

    public static final String OBJECT = "LANE_C_LIVING_CLOCKS";
    public static final String GATE = "LC0_TWO_GENERATION";
    public static final String DEVICE = "LC0_OPEN_BATH";
    public static final String PARENT_GATE = "LC0_REPLICATION_VISIBLE";

    public static final int N0 = 64;
    public static final long RNG_SEED = 101L;

    public static final double LAMBDA_S_PER_S = 1.0 / 3600.0;
    public static final double C_S_MM = 0.5;
    public static final double K_S_MM = 0.02;
    public static final double MONOD = C_S_MM / (C_S_MM + K_S_MM);
    /** ln 2 / (lambda_S * Monod). Spoken 2595.1 s. Same as LC0. */
    public static final double T_DIV_S = Math.log(2.0) / (LAMBDA_S_PER_S * MONOD);
    /** Two generations. Spoken 5190.3 s. The only horizon change from LC0. */
    public static final double T_HORIZON_S = 2.0 * T_DIV_S;

    public static final double ELL0_UM = 1.0;
    public static final double ELL_DIV_UM = 3.0;
    public static final double ALPHA_PER_S = Math.log(ELL_DIV_UM / ELL0_UM) / T_DIV_S;

    public static final double WINDOW_S = 300.0;
    /** Full windows on [0, 5100). Stub [5100, T] is logged, not scored. */
    public static final int N_FULL_WINDOWS = 17;
    public static final double DT_S = 1.0;

    public static final double CLOCK_DELTA = 0.20;
    public static final int VISIBILITY_MIN_BIRTHS = 5;
    public static final int N_EVER_COMPUTE_CAP = 512;

    public static final double BOX_X_UM = 200.0;
    public static final double BOX_Y_UM = 200.0;
    public static final double BOX_Z_UM = 10.0;
    public static final int GRID = 8;

    private LaneCLc0TwoGenIdentity() { }

    public static double nPred(double t) {
        return N0 * Math.pow(2.0, t / T_DIV_S);
    }

    public static double founderX(int i) {
        int col = i % GRID;
        return (col + 0.5) * (BOX_X_UM / GRID);
    }

    public static double founderY(int i) {
        int row = i / GRID;
        return (row + 0.5) * (BOX_Y_UM / GRID);
    }

    /**
     * Inverse CDF of stable age on [0, T_div):
     * p(a)=(2 ln 2 / T) 2^{-a/T}, F(a)=2(1-2^{-a/T}), a=-T log2(1-u/2).
     * Same as LC0. Do not use a=-T log2(1-u).
     */
    public static double stableAge(double u) {
        double uu = Math.min(Math.max(u, 0.0), 1.0 - 1.0e-15);
        return -T_DIV_S * (Math.log(1.0 - 0.5 * uu) / Math.log(2.0));
    }
}
