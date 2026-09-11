package bsim.laneb;

/**
 * Frozen B4 single-cell CheY probe. Not a colony-kappa gate.
 * Parents B0–B3 occupation FAILs are kept.
 */
public final class LaneBCheyIdentity {

    public static final String OBJECT = "LANE_B_CHEMOTACTIC_SPATIAL";
    public static final String GATE = "B4_CHEY_PROBE";
    public static final String DEVICE = "LANE_B_CHEY_PROBE";
    public static final double DL = 0.3;
    public static final double T_ADAPT = 5.0;
    public static final double T_TOTAL = 15.0;
    public static final double HOLD_LO = 0.85;
    public static final double HOLD_HI = 1.15;
    public static final double STEP_UP_MAX = 0.88;
    public static final double STEP_DOWN_MIN = 1.12;
    public static final double WIN_LO = 5.4;
    public static final double WIN_HI = 7.0;
    public static final long RNG_SEED = 0L;

    private LaneBCheyIdentity() { }
}
