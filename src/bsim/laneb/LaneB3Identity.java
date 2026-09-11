package bsim.laneb;

/**
 * Frozen B3 identity for LANE_B_CHEMOTACTIC_SPATIAL / B3_SUSTAINED_STRIPE.
 * Parent FAILs B0, B1, B2 are kept. Not Fig. 4b. Not Lane A. Not NARMA.
 * Not a B2 L_slab raise. J_slab and T_hold are Python field-scout values,
 * frozen before motile traces.
 */
public final class LaneB3Identity {

    public static final String OBJECT = "LANE_B_CHEMOTACTIC_SPATIAL";
    public static final String GATE = "B3_SUSTAINED_STRIPE";
    public static final String DEVICE = "LANE_B_MM_DISH_FLOW0";
    public static final String PARENT_FAIL_B0 = "B0_MOTILITY_MEMORY";
    public static final String PARENT_FAIL_B1 = "B1_PAINT_HOLD";
    public static final String PARENT_FAIL_B2 = "B2_STRIPE_OCCUPY";
    public static final String NOT_FIG4B = "NOT_FIG4B";

    public static final double LAMBDA_UM = 200.0;
    public static final double SLAB_UM = 100.0;
    public static final double ELL_RUN_UM = 17.2;
    public static final double J_SLAB = 0.5022121815796962;
    public static final double DELTA_L_T25 = 0.3;
    public static final double T_LIVE = 25.0;
    public static final double T_OCC = 25.0;
    public static final double T_HOLD = 25.1;
    public static final double T_CAP = 40.0;
    public static final double KAPPA_FIELD_HOLD_MAX = 0.15;
    public static final double KAPPA_MOTILE_MIN = 0.12;
    public static final double KAPPA_MINUS_FIELD = 0.08;
    public static final double KAPPA_MINUS_OFF = 0.08;
    public static final double KAPPA_SILENT_MAX = 0.08;
    public static final double LEDGER_REL = 1e-3;
    public static final long RNG_SEED = 0L;

    private LaneB3Identity() { }

    public static double paintDuration() {
        double t = T_HOLD;
        if (t > T_CAP) t = T_CAP;
        if (t < T_OCC) t = T_OCC;
        return LaneBIdentity.DT * Math.round(t / LaneBIdentity.DT);
    }
}
