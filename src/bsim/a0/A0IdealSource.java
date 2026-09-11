package bsim.a0;

import bsim.c1c.DaninoSIC1cFilledPocket;
import bsim.transport.BSimTransportField;

/**
 * A0_IDEAL_SOURCE: named immobilized one-way releasing module.
 * No internal AC ODE. Payload = extracellular AHL. NOT_FIG4B.
 *
 * J(t) = J_max * u(t). Deposit in the scheduler flux-deposition phase
 * into the containing C1c pocket voxel only.
 */
public final class A0IdealSource {

    public static final String STATUS = "A0_IDEAL_SOURCE";
    public static final String NOT_FIG4B = "NOT_FIG4B";
    public static final String PAYLOAD = "extracellular_AHL";
    public static final String DEVICE = "C1c_FILLED_POCKET";

    public static final double X_UM = 50.0;
    public static final double Y_UM = 50.0;
    public static final double Z_UM = 0.825;
    public static final int I_AC = 25;
    public static final int J_AC = 25;
    public static final int K_AC = 0;
    public static final int I_FAR = 0;
    public static final int J_FAR = 0;
    public static final int K_FAR = 0;

    public static final double HE_KICK = 0.05;
    public static final double DELTA_M_KICK = HE_KICK * DaninoSIC1cFilledPocket.V_E_POCKET;
    public static final double T_ON_PULSE = 0.0;
    public static final double T_OFF_PULSE = 1.0;
    public static final double T_ON_STEP = 0.0;
    public static final double T_OFF_STEP = 8.0;
    public static final double J_MAX = DELTA_M_KICK / (T_OFF_PULSE - T_ON_PULSE);

    public static final double T_END_FIELD = 10.0;
    public static final double SAMPLE_DT_FIELD = 0.1;
    public static final double T_END_BACTERIA = 1000.0;
    public static final double SAMPLE_DT_BACTERIA = 0.5;
    public static final double RK_DT = 0.001;

    public enum Command {
        U_PULSE(T_ON_PULSE, T_OFF_PULSE),
        U_STEP(T_ON_STEP, T_OFF_STEP),
        U_ZERO(Double.NaN, Double.NaN);

        public final double tOn;
        public final double tOff;

        Command(double tOn, double tOff) {
            this.tOn = tOn;
            this.tOff = tOff;
        }
    }

    private final Command command;
    private final double jMax;
    private double commandedMass;

    public A0IdealSource(Command command) {
        if (Math.abs(DELTA_M_KICK - 825.0) > 1e-12) {
            throw new IllegalStateException(
                    "A0 kick mass must be 825 A0_IDEAL_SOURCE NOT_FIG4B, got " + DELTA_M_KICK);
        }
        if (Math.abs(J_MAX - 825.0) > 1e-12) {
            throw new IllegalStateException(
                    "A0 J_max must be 825 A0_IDEAL_SOURCE NOT_FIG4B, got " + J_MAX);
        }
        if (I_AC == 0 || J_AC == 0 || I_AC == DaninoSIC1cFilledPocket.NX - 1
                || J_AC == DaninoSIC1cFilledPocket.NY - 1) {
            throw new IllegalStateException("A0 AC voxel is on a wall/mouth A0_IDEAL_SOURCE NOT_FIG4B");
        }
        this.command = command;
        this.jMax = command == Command.U_ZERO ? 0.0 : J_MAX;
        this.commandedMass = 0.0;
    }

    public Command command() {
        return command;
    }

    public double commandedMass() {
        return commandedMass;
    }

    public double jMax() {
        return jMax;
    }

    /** Exact overlap of [tStart, tStart+dt) with the frozen on-window. */
    public static double overlap(double tStart, double dt, double tOn, double tOff) {
        if (!Double.isFinite(tOn) || !Double.isFinite(tOff) || dt <= 0.0) {
            return 0.0;
        }
        double lo = Math.max(tStart, tOn);
        double hi = Math.min(tStart + dt, tOff);
        return Math.max(0.0, hi - lo);
    }

    public double u(double tStart, double dt) {
        if (command == Command.U_ZERO || jMax == 0.0) {
            return 0.0;
        }
        double ov = overlap(tStart, dt, command.tOn, command.tOff);
        return ov > 0.0 ? 1.0 : 0.0;
    }

    public double j(double tStart, double dt) {
        if (command == Command.U_ZERO || jMax == 0.0) {
            return 0.0;
        }
        return jMax * overlap(tStart, dt, command.tOn, command.tOff) / dt;
    }

    /**
     * Deposit commanded mass for this scheduler step into the AC voxel.
     * Must be called from {@code depositFluxes}, not from a bacterium action.
     */
    public double deposit(BSimTransportField field, double tStart, double dt) {
        if (command == Command.U_ZERO || jMax == 0.0) {
            return 0.0;
        }
        double dM = jMax * overlap(tStart, dt, command.tOn, command.tOff);
        if (dM > 0.0) {
            field.addQuantity(I_AC, J_AC, K_AC, dM);
            commandedMass += dM;
        }
        return dM;
    }

    /** Expected closed-form commanded mass for a run of duration tEnd. */
    public double expectedCommandedMass(double tEnd) {
        if (command == Command.U_ZERO || jMax == 0.0) {
            return 0.0;
        }
        return jMax * overlap(0.0, tEnd, command.tOn, command.tOff);
    }
}
