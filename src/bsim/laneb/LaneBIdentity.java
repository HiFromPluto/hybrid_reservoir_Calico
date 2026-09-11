package bsim.laneb;

/**
 * Frozen B0 identity for LANE_B_CHEMOTACTIC_SPATIAL / B0_MOTILITY_MEMORY.
 * Not Fig. 4b. Not Lane A. Not C1c. LANES_REVIEW §4 is not this freeze.
 */
public final class LaneBIdentity {

    public static final String OBJECT = "LANE_B_CHEMOTACTIC_SPATIAL";
    public static final String GATE = "B0_MOTILITY_MEMORY";
    public static final String DEVICE = "LANE_B_MM_DISH_FLOW0";
    public static final String DEVICE_FLOW8 = "LANE_B_MM_DISH_FLOW8";
    public static final String NOT_FIG4B = "NOT_FIG4B";

    public static final double LX = 1000.0;
    public static final double LY = 500.0;
    public static final double LZ = 10.0;
    public static final double FLOW_X = 0.0;
    public static final double FLOW_Y = 0.0;
    public static final double FLOW_Z = 0.0;

    public static final int NX = 100;
    public static final int NY = 50;
    public static final int NZ = 1;
    public static final int RX = 20;
    public static final int RY = 10;
    public static final double LAMBDA_UM = 50.0;
    public static final double DX = LX / NX;
    public static final double DY = LY / NY;
    public static final double DZ = LZ / NZ;

    public static final double T_K = 305.0;
    public static final double ETA = 2.7e-3;
    public static final double DT = 0.02;
    public static final double SAMPLE_DT = 0.10;

    public static final double D_ATT = 800.0;
    public static final double D_REP = 800.0;
    public static final double K_ATT = 0.0;
    public static final double K_REP = 0.0;

    public static final double J_MAX = 20000.0;
    public static final double T_OFF = 10.0;
    public static final int AC_COUNT = 4;
    public static final double[][] AC_XYZ = {
            {200.0, 150.0, 5.0},
            {200.0, 350.0, 5.0},
            {800.0, 150.0, 5.0},
            {800.0, 350.0, 5.0}
    };
    public static final LaneBA0Source.Payload[] AC_PAYLOAD = {
            LaneBA0Source.Payload.ATTRACTANT,
            LaneBA0Source.Payload.ATTRACTANT,
            LaneBA0Source.Payload.REPELLENT,
            LaneBA0Source.Payload.REPELLENT
    };

    public static final double ALPHA = 2.0;
    public static final double K_R = 0.05;
    public static final double K_B = 0.05;
    public static final double K_A = 1.0;
    public static final double GAMMA_Y = 0.5;
    public static final double TAU_A = 0.5;
    public static final double Y0 = K_A * 0.5 / GAMMA_Y;
    public static final double P_END_RUN_ELSE = 1.0 / 0.86;
    public static final double P_END_TUMBLE = 1.0 / 0.14;
    public static final double FORCE_PN = 1.0;
    public static final double RADIUS = 1.0;

    public static final int N_MSD = 400;
    public static final int N_DENSITY = 1200;
    public static final double T_MSD = 12.0;
    public static final double T_MIX = 16.0;
    public static final double T_PAINT_BASE = 22.0;
    public static final double MSD_FIT_T0 = 2.0;
    public static final double MSD_FIT_T1 = 10.0;
    public static final double BLOB_SIGMA = 40.0;
    public static final double SPAWN_MARGIN = 80.0;

    public static final double W_OVER_TAU = 0.3;
    public static final double W_BAND_LO = 0.1;
    public static final double W_BAND_HI = 1.0;
    public static final double KAPPA_MOTILE_MIN = 0.12;
    public static final double KAPPA_MINUS_FIELD = 0.08;
    public static final double KAPPA_MINUS_OFF = 0.08;
    public static final double KAPPA_SILENT_MAX = 0.08;
    public static final double MU_FAIL_BELOW = 1.0;
    public static final double MU_FAIL_ABOVE = 5000.0;
    public static final double LEDGER_REL = 1e-3;

    public static final double MU0_MIDDLEBROOKS = 120.0;
    public static final double MU0_ZHAO = 130.0;
    public static final double MU0_ZHAO_SD = 21.0;
    public static final double D_SE_BERG = 0.083;
    public static final double D_SE_WATER37 = 0.33;

    public static final long RNG_SEED = 0L;

    public static final double KB = 1.38e-23;

    private LaneBIdentity() { }

    public static double stokesEinsteinUm2s() {
        double rM = RADIUS * 1e-6;
        double dM2s = KB * T_K / (6.0 * Math.PI * ETA * rM);
        return dM2s * 1e12;
    }

    public static int pdeI(double x) {
        int i = (int) (x / DX);
        if (i < 0) return 0;
        if (i >= NX) return NX - 1;
        return i;
    }

    public static int pdeJ(double y) {
        int j = (int) (y / DY);
        if (j < 0) return 0;
        if (j >= NY) return NY - 1;
        return j;
    }

    public static int readI(double x) {
        int i = (int) (x / (LX / RX));
        if (i < 0) return 0;
        if (i >= RX) return RX - 1;
        return i;
    }

    public static int readJ(double y) {
        int j = (int) (y / (LY / RY));
        if (j < 0) return 0;
        if (j >= RY) return RY - 1;
        return j;
    }
}
