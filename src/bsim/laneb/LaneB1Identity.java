package bsim.laneb;

/**
 * Frozen B1 identity for LANE_B_CHEMOTACTIC_SPATIAL / B1_PAINT_HOLD.
 * Parent B0_MOTILITY_MEMORY FAIL is kept. Not Fig. 4b. Not Lane A. Not NARMA.
 * T_hold is the Python field-scout value, frozen before motile traces.
 */
public final class LaneB1Identity {

    public static final String OBJECT = "LANE_B_CHEMOTACTIC_SPATIAL";
    public static final String GATE = "B1_PAINT_HOLD";
    public static final String DEVICE = "LANE_B_MM_DISH_FLOW0";
    public static final String PARENT_FAIL = "B0_MOTILITY_MEMORY";
    public static final String NOT_FIG4B = "NOT_FIG4B";

    public static final double T_OFF = 60.0;
    public static final double T_OCC = 60.0;
    public static final double T_MIX = 180.0;
    public static final double T_HOLD = 275.3;
    public static final double T_CAP = 400.0;
    public static final double STRIPE_WAVELENGTH_UM = 200.0;
    public static final double STRIPE_HIGH_UM = 100.0;
    public static final double KAPPA_FIELD_HOLD_MAX = 0.15;
    public static final double KAPPA_MOTILE_MIN = 0.12;
    public static final double KAPPA_MINUS_FIELD = 0.08;
    public static final double KAPPA_MINUS_OFF = 0.08;
    public static final double KAPPA_SILENT_MAX = 0.08;
    public static final double LEDGER_REL = 1e-3;
    public static final long RNG_SEED = 0L;

    private LaneB1Identity() { }

    public static double paintDuration() {
        double t = T_HOLD;
        if (t > T_CAP) t = T_CAP;
        if (t < T_OFF) t = T_OFF;
        return LaneBIdentity.DT * Math.round(t / LaneBIdentity.DT);
    }
}
