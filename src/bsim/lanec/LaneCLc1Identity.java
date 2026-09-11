package bsim.lanec;

/**
 * Frozen LC1_DEATH_VISIBLE identity. Death only. Schink WT gamma.
 * Not LC0 with death on. Not Job 6 RK4. Not Paper 1 rates. Not NARMA.
 */
public final class LaneCLc1Identity {

    public static final String OBJECT = "LANE_C_LIVING_CLOCKS";
    public static final String GATE = "LC1_DEATH_VISIBLE";
    public static final String DEVICE = "LC0_OPEN_BATH";
    public static final String PARENT_GATE = "LC0_REPLICATION_VISIBLE";

    public static final int N0 = 64;
    public static final long RNG_SEED = 101L;

    /** Schink 2019 WT after transfer to zero carbon. Not 0.005 /s. */
    public static final double GAMMA_PER_D = 0.43;
    public static final double SECONDS_PER_DAY = 86400.0;
    public static final double GAMMA_PER_S = GAMMA_PER_D / SECONDS_PER_DAY;

    public static final double T_HORIZON_D = 4.0;
    public static final double T_HORIZON_S = T_HORIZON_D * SECONDS_PER_DAY;
    public static final double T_HALF_D = Math.log(2.0) / GAMMA_PER_D;

    public static final double WINDOW_D = 1.0;
    public static final int N_FULL_WINDOWS = 4;

    public static final double CLOCK_DELTA = 0.20;
    public static final int VISIBILITY_MIN_DEATHS = 5;

    public static final double BOX_X_UM = 200.0;
    public static final double BOX_Y_UM = 200.0;
    public static final double BOX_Z_UM = 10.0;
    public static final int GRID = 8;

    private LaneCLc1Identity() { }

    public static double nPred(double tDays) {
        return N0 * Math.exp(-GAMMA_PER_D * tDays);
    }

    public static double expectedDeaths(double t0Days, double dtDays) {
        return nPred(t0Days) * (1.0 - Math.exp(-GAMMA_PER_D * dtDays));
    }

    /**
     * Inverse CDF of Exp(gamma): t = -ln(1-u)/gamma days.
     */
    public static double deathTimeDays(double u) {
        double uu = Math.min(Math.max(u, 0.0), 1.0 - 1.0e-15);
        return -Math.log(1.0 - uu) / GAMMA_PER_D;
    }

    public static double founderX(int i) {
        int col = i % GRID;
        return (col + 0.5) * (BOX_X_UM / GRID);
    }

    public static double founderY(int i) {
        int row = i / GRID;
        return (row + 0.5) * (BOX_Y_UM / GRID);
    }
}
