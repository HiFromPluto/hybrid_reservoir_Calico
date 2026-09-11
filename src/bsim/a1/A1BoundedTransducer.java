package bsim.a1;

import bsim.a0.A0IdealSource;
import bsim.c1c.DaninoSIC1cFilledPocket;
import bsim.transport.BSimTransportField;

/**
 * A1_BOUNDED_TRANSDUCER: delayed, leaky, saturating, payload-limited
 * one-way AC envelope. HYPOTHETICAL_DESIGN_ENVELOPE. NOT_FIG4B.
 * Not calibrated. Not Lentini TX–TL. Payload = extracellular AHL.
 */
public final class A1BoundedTransducer {

    public static final String STATUS = "A1_BOUNDED_TRANSDUCER";
    public static final String PROVENANCE = "HYPOTHETICAL_DESIGN_ENVELOPE";
    public static final String NOT_FIG4B = "NOT_FIG4B";
    public static final String PAYLOAD = "extracellular_AHL";
    public static final String DEVICE = "C1c_FILLED_POCKET";

    public static final double TAU_AC = 2.0;
    public static final double N_HILL = 2.0;
    public static final double K = 0.5;
    public static final double J_MAX = A0IdealSource.J_MAX;
    public static final double J_LEAK = 0.1;
    public static final double M0 = 5000.0;
    public static final double M_HALF = 500.0;
    public static final double J_SYN = 0.0;
    public static final double X1_0 = 0.0;

    public static final double T_ON_PULSE = A0IdealSource.T_ON_PULSE;
    public static final double T_OFF_PULSE = A0IdealSource.T_OFF_PULSE;
    public static final double T_ON_STEP = A0IdealSource.T_ON_STEP;
    public static final double T_OFF_STEP = A0IdealSource.T_OFF_STEP;
    public static final double T_ON_EXHAUST = 0.0;
    public static final double T_OFF_EXHAUST = 50.0;

    public static final double T_END_FIELD = 10.0;
    public static final double T_END_EXHAUST = 50.0;
    public static final double SAMPLE_DT_FIELD = 0.1;
    public static final double T_END_BACTERIA = 1000.0;
    public static final double SAMPLE_DT_BACTERIA = 0.5;
    public static final double RK_DT = 0.001;

    public enum Command {
        U_PULSE(T_ON_PULSE, T_OFF_PULSE),
        U_STEP(T_ON_STEP, T_OFF_STEP),
        U_ZERO(Double.NaN, Double.NaN),
        U_EXHAUST(T_ON_EXHAUST, T_OFF_EXHAUST);

        public final double tOn;
        public final double tOff;

        Command(double tOn, double tOff) {
            this.tOn = tOn;
            this.tOff = tOff;
        }
    }

    private final Command command;
    private final boolean acOff;
    private double x1;
    private double payload;
    private double released;
    private double lastJ;
    private double lastU;

    public A1BoundedTransducer(Command command, boolean acOff) {
        freezeCheck();
        this.command = command;
        this.acOff = acOff;
        this.x1 = X1_0;
        this.payload = M0;
        this.released = 0.0;
        this.lastJ = 0.0;
        this.lastU = 0.0;
    }

    public static void freezeCheck() {
        if (Math.abs(J_MAX - 825.0) > 1e-12 || Math.abs(J_LEAK - 0.1) > 1e-15) {
            throw new IllegalStateException("A1 J_max/J_leak freeze broken A1_BOUNDED_TRANSDUCER NOT_FIG4B");
        }
        if (Math.abs(TAU_AC - 2.0) > 0.0 || Math.abs(N_HILL - 2.0) > 0.0 || Math.abs(K - 0.5) > 0.0) {
            throw new IllegalStateException("A1 tau/n/K freeze broken A1_BOUNDED_TRANSDUCER NOT_FIG4B");
        }
        if (Math.abs(M0 - 5000.0) > 0.0 || Math.abs(M_HALF - 500.0) > 0.0 || J_SYN != 0.0) {
            throw new IllegalStateException("A1 payload freeze broken A1_BOUNDED_TRANSDUCER NOT_FIG4B");
        }
        if (J_LEAK >= 0.01 * J_MAX) {
            throw new IllegalStateException("A1 J_leak must be << J_max A1_BOUNDED_TRANSDUCER NOT_FIG4B");
        }
        if (A0IdealSource.I_AC != 25 || A0IdealSource.J_AC != 25) {
            throw new IllegalStateException("A1 AC voxel must match A0 A1_BOUNDED_TRANSDUCER NOT_FIG4B");
        }
    }

    public Command command() {
        return command;
    }

    public boolean acOff() {
        return acOff;
    }

    public double x1() {
        return x1;
    }

    public double payload() {
        return payload;
    }

    public double released() {
        return released;
    }

    public double lastJ() {
        return lastJ;
    }

    public double lastU() {
        return lastU;
    }

    public double payloadDrop() {
        return M0 - payload;
    }

    public double uOn(double tStart, double dt) {
        if (acOff || command == Command.U_ZERO) {
            return 0.0;
        }
        return A0IdealSource.overlap(tStart, dt, command.tOn, command.tOff) > 0.0 ? 1.0 : 0.0;
    }

    public static double a0J(Command command, double t) {
        if (command == Command.U_ZERO || command == Command.U_EXHAUST) {
            if (command == Command.U_EXHAUST && t >= command.tOn && t < command.tOff) {
                return J_MAX;
            }
            return 0.0;
        }
        if (t > 0.0 && t < command.tOff) {
            return J_MAX;
        }
        return 0.0;
    }

    /**
     * Advance the AC ODE over [tStart, tStart+dt) and deposit J_S dt.
     * Must be called from {@code depositFluxes}.
     */
    public double deposit(BSimTransportField field, double tStart, double dt) {
        if (acOff) {
            lastJ = 0.0;
            lastU = 0.0;
            return 0.0;
        }
        double u = uOn(tStart, dt);
        lastU = u;
        x1 = u + (x1 - u) * Math.exp(-dt / TAU_AC);
        double hillNum = Math.pow(x1, N_HILL);
        double hill = hillNum / (Math.pow(K, N_HILL) + hillNum);
        double pay = payload / (payload + M_HALF);
        double j = J_LEAK + J_MAX * hill * pay;
        double dM = j * dt;
        if (dM > payload) {
            dM = payload;
            j = dt > 0.0 ? dM / dt : 0.0;
        }
        payload -= dM;
        if (payload < 0.0) {
            payload = 0.0;
        }
        lastJ = j;
        if (dM > 0.0) {
            field.addQuantity(A0IdealSource.I_AC, A0IdealSource.J_AC, A0IdealSource.K_AC, dM);
            released += dM;
        }
        return dM;
    }
}
