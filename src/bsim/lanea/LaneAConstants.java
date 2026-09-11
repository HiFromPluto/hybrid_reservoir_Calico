package bsim.lanea;

/**
 * Frozen Lane A millimetre identity. Do not retune after traces.
 * LANE_A_OCCUPIED_MILLIMETRE. Not paper 1. Not C1c. Not Fig. 4b. Not NARMA.
 */
public final class LaneAConstants {

    public static final String LABEL = "LANE_A_OCCUPIED_MILLIMETRE";
    public static final String GATE = "LaneA_N0_OCCUPANCY";

    public static final double BOUND_X = 1000.0;
    public static final double BOUND_Y = 500.0;
    public static final double BOUND_Z = 10.0;
    public static final double DT = 0.05;
    public static final int GRID_X = 50;
    public static final int GRID_Y = 25;
    public static final int GRID_Z = 1;
    public static final int READOUT_X = 20;
    public static final int READOUT_Y = 10;
    public static final int READOUT_Z = 1;

    public static final double AC_X = 500.0;
    public static final double AC_Y = 250.0;
    public static final double AC_Z = 5.0;
    public static final int I_AC = 25;
    public static final int J_AC = 12;
    public static final int K_AC = 0;

    public static final double FLOW = 0.0;
    public static final double D_AHL = 159.0;
    public static final double K_AHL = 0.0033;
    public static final double J_MAX = 128000000.0;
    public static final double U_ON = 0.5;
    public static final double MOLECULES_PER_UM3_PER_UM = 602.2;

    public static final double HILL_K_UM = 1.6;
    public static final double HILL_N = 2.0;
    public static final double TAU_R_S = 15.0;
    public static final double TAU_L_S = 1500.0;

    public static final double WINDOW_S = 300.0;
    public static final double PULSE_S = 75.0;
    public static final int NUM_WINDOWS = 8;
    public static final double T_END_DRIVEN = NUM_WINDOWS * WINDOW_S;
    public static final double SAMPLE_INTERVAL_S = 20.0;
    public static final int SAMPLES_PER_WINDOW = 16;
    public static final int OCCUPANCY_WINDOW_LO = 4;
    public static final int OCCUPANCY_WINDOW_HI = 7;
    public static final double OCCUPANCY_ALIVE_MIN_R = 0.05;

    public static final double T_END_CLOSED = 10.0;
    public static final double CLOSED_IMPULSE_S = 1.0;
    public static final double CLOSED_DELTA_M = J_MAX * CLOSED_IMPULSE_S;

    public static final int N_CELLS = 1800;
    public static final long RNG_SEED = 101L;
    public static final double SEED_X0 = 300.0;
    public static final double SEED_XSPAN = 400.0;
    public static final double SEED_Y0 = 150.0;
    public static final double SEED_YSPAN = 200.0;
    public static final double SEED_Z = 5.0;

    public static final double LEDGER_CMD_REL = 1e-3;
    public static final double N0_MASS_REL = 1e-11;
    public static final double EXPECTED_DRIVEN_MASS = J_MAX * U_ON * NUM_WINDOWS * PULSE_S;

    public static final String CARRIER_LABEL = "LANE_A_CARRIER";
    public static final String CARRIER_GATE = "LaneA_CARRIER";
    public static final int CARRIER_DISCARD_LO = 0;
    public static final int CARRIER_DISCARD_HI = 1;
    public static final int CARRIER_TRAIN_LO = 2;
    public static final int CARRIER_TRAIN_HI = 5;
    public static final int CARRIER_TEST_LO = 6;
    public static final int CARRIER_TEST_HI = 7;
    public static final double CARRIER_MARGIN = 0.02;
    public static final int READOUT_BINS = READOUT_X * READOUT_Y * READOUT_Z;

    public static final String NARMA10_LABEL = "LANE_A_NARMA10";
    public static final String NARMA10_GATE = "LaneA_NARMA10";
    public static final int NARMA10_WINDOWS = 200;
    public static final double NARMA10_T_END = NARMA10_WINDOWS * WINDOW_S;
    public static final int NARMA10_WASHOUT_HI = 39;
    public static final int NARMA10_TRAIN_LO = 40;
    public static final int NARMA10_TRAIN_HI = 149;
    public static final int NARMA10_TEST_LO = 150;
    public static final int NARMA10_TEST_HI = 199;
    public static final int NARMA10_OCC_LO = 150;
    public static final int NARMA10_OCC_HI = 199;
    public static final double NARMA10_MARGIN = 0.02;
    public static final String NARMA10_U_SHA256 =
            "2fb692bb50bc8b2a9b8e28d2faf6ee63add7e6c980335c78f4f2c0ef3f07f06b";
    public static final String NARMA10B_U_SHA256_FORBIDDEN =
            "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e";

    public static final String MG_LABEL = "LANE_A_MG";
    public static final String MG_GATE = "LaneA_MG";
    public static final String MG_U_SHA256 =
            "e06810ea70f61bbe4baf61939f93b1e1cae9d92d9eee8d757865c4935dbac780";

    public static final String LORENZ_LABEL = "LANE_A_LORENZ";
    public static final String LORENZ_GATE = "LaneA_LORENZ";
    public static final String LORENZ_U_SHA256 =
            "69b861f256fb39527c037b994885d119ccc8ab2873a3d2bca59f7a8802c44aff";

    public static final String WAVEFORM_LABEL = "LANE_A_WAVEFORM";
    public static final String WAVEFORM_GATE = "LaneA_WAVEFORM";
    public static final String WAVEFORM_U_SHA256 =
            "0979090c051c513dce14632e83ae9bf025c4b666e497dc879ae1db4bba30d465";
    public static final int WAVEFORM_WINDOWS = 400;
    public static final double WAVEFORM_T_END = WAVEFORM_WINDOWS * WINDOW_S;
    public static final int WAVEFORM_OCC_LO = 272;
    public static final int WAVEFORM_OCC_HI = 399;
    public static final double WAVEFORM_AUC_MARGIN = 0.05;

    public static final String KR_GR_LABEL = "LANE_A_KR_GR";
    public static final String KR_GR_GATE = "LaneA_KR_GR";
    public static final int KR_GR_WINDOWS = 840;
    public static final double KR_GR_T_END = KR_GR_WINDOWS * WINDOW_S;
    public static final int KR_GR_WASHOUT_HI = 39;
    public static final int KR_GR_OCC_LO = 600;
    public static final int KR_GR_OCC_HI = 839;
    public static final String KR_IID_U_SHA256 =
            "1134e3ac2031cf33f699341090b3d186f9a776a81cfe072b18e9d0c9553ba1b5";
    public static final String GR_CONST_U_SHA256 =
            "b7503301c6a2fabb9c4324c44d7da57ebea8632442d8681e862886643c4a8f72";

    public static final String FIELD_NULL_LABEL = "LANE_A_FIELD_NULL_TASK";
    public static final String FIELD_NULL_GATE = "LaneA_FIELD_NULL";
    public static final String FIELD_NULL_U_SHA256 =
            "3025394aa72732228f4db2a23ba33c84b4a2b2b6b3a215d098d7d3b6f8072d53";
    public static final int FIELD_NULL_WINDOWS = WAVEFORM_WINDOWS;
    public static final double FIELD_NULL_T_END = FIELD_NULL_WINDOWS * WINDOW_S;
    public static final int FIELD_NULL_OCC_LO = WAVEFORM_OCC_LO;
    public static final int FIELD_NULL_OCC_HI = WAVEFORM_OCC_HI;
    public static final double FIELD_NULL_AUC_MARGIN = 0.05;

    public static int readoutIndex(int x, int y, int z) {
        return x * READOUT_Y * READOUT_Z + y * READOUT_Z + z;
    }

    public static double commandU(int window, double t) {
        if (window < 0 || window >= NUM_WINDOWS) {
            return 0.0;
        }
        double tIn = t - window * WINDOW_S;
        if (tIn < 0.0) {
            tIn = 0.0;
        }
        return tIn < PULSE_S ? U_ON : 0.0;
    }

    private LaneAConstants() { }
}
