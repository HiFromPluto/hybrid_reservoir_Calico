package bsim.laneb;

/**
 * Frozen B2 identity for LANE_B_CHEMOTACTIC_SPATIAL / B2_STRIPE_OCCUPY.
 * Parent FAILs B0_MOTILITY_MEMORY and B1_PAINT_HOLD are kept.
 * Not Fig. 4b. Not Lane A. Not NARMA. Not a B1 J_max hunt.
 * L_slab, t_occ, T_hold are the Python field-scout values, frozen before motile traces.
 */
public final class LaneB2Identity {

    public static final String OBJECT = "LANE_B_CHEMOTACTIC_SPATIAL";
    public static final String GATE = "B2_STRIPE_OCCUPY";
    public static final String DEVICE = "LANE_B_MM_DISH_FLOW0";
    public static final String PARENT_FAIL_B0 = "B0_MOTILITY_MEMORY";
    public static final String PARENT_FAIL_B1 = "B1_PAINT_HOLD";
    public static final String NOT_FIG4B = "NOT_FIG4B";

    public static final double LAMBDA_UM = 200.0;
    public static final double SLAB_UM = 100.0;
    public static final double ELL_RUN_UM = 17.2;
    public static final double L_SLAB = 1.3891295989117673;
    public static final double DELTA_L_T1 = 0.3;
    public static final double T_OCC = 1.1;
    public static final double T_HOLD = 3.1;
    public static final double T_CAP = 40.0;
    public static final double KAPPA_FIELD_OCC_MIN = 0.40;
    public static final double KAPPA_FIELD_HOLD_MAX = 0.15;
    public static final double KAPPA_MOTILE_MIN = 0.12;
    public static final double KAPPA_MINUS_FIELD = 0.08;
    public static final double KAPPA_MINUS_OFF = 0.08;
    public static final double KAPPA_SILENT_MAX = 0.08;
    public static final double LEDGER_REL = 1e-3;
    public static final long RNG_SEED = 0L;

    private LaneB2Identity() { }

    public static boolean evenAttSlab(double x) {
        return Math.floor(x / SLAB_UM) % 2.0 == 0.0;
    }

    public static double paintDuration() {
        double t = T_HOLD;
        if (t > T_CAP) t = T_CAP;
        if (t < T_OCC) t = T_OCC;
        return LaneBIdentity.DT * Math.round(t / LaneBIdentity.DT);
    }
}
