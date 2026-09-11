package bsim.lanea;

import bsim.transport.BSimTransportField;

/**
 * A0-class ideal current on the Lane A millimetre dish: J = J_max u.
 * Not C1c J_max=825. No internal AC ODE. LANE_A_OCCUPIED_MILLIMETRE.
 */
public final class LaneAIdealSource {

    public enum Command {
        U_PULSE_TRAIN,
        U_ZERO,
        U_CLOSED_IMPULSE,
        U_NARMA10
    }

    private final Command command;
    private final double[] narmaU;
    private double commandedMass;

    public LaneAIdealSource(Command command) {
        this(command, null);
    }

    public LaneAIdealSource(double[] narmaU) {
        this(Command.U_NARMA10, narmaU);
    }

    private LaneAIdealSource(Command command, double[] narmaU) {
        if (Math.abs(LaneAConstants.J_MAX - 128000000.0) > 1e-6) {
            throw new IllegalStateException("Lane A J_max must stay 1.28e8 " + LaneAConstants.LABEL);
        }
        if (command == Command.U_NARMA10) {
            if (narmaU == null
                    || (narmaU.length != LaneAConstants.NARMA10_WINDOWS
                    && narmaU.length != LaneAConstants.WAVEFORM_WINDOWS
                    && narmaU.length != LaneAConstants.KR_GR_WINDOWS)) {
                throw new IllegalStateException(
                        "per-window u must have "
                                + LaneAConstants.NARMA10_WINDOWS + " or "
                                + LaneAConstants.WAVEFORM_WINDOWS + " or "
                                + LaneAConstants.KR_GR_WINDOWS + " windows "
                                + LaneAConstants.LABEL);
            }
            for (double u : narmaU) {
                if (u < 0.0 || u > 0.5) {
                    throw new IllegalStateException(
                            "per-window u escaped [0, 0.5] " + LaneAConstants.LABEL);
                }
            }
        } else if (narmaU != null) {
            throw new IllegalStateException("per-window u only belongs on U_NARMA10");
        }
        this.command = command;
        this.narmaU = narmaU;
        this.commandedMass = 0.0;
    }

    public Command command() {
        return command;
    }

    public double commandedMass() {
        return commandedMass;
    }

    public static double overlap(double tStart, double dt, double tOn, double tOff) {
        if (!Double.isFinite(tOn) || !Double.isFinite(tOff) || dt <= 0.0) {
            return 0.0;
        }
        double lo = Math.max(tStart, tOn);
        double hi = Math.min(tStart + dt, tOff);
        return Math.max(0.0, hi - lo);
    }

    public double onOverlap(double tStart, double dt) {
        if (command == Command.U_ZERO) {
            return 0.0;
        }
        if (command == Command.U_CLOSED_IMPULSE) {
            return overlap(tStart, dt, 0.0, LaneAConstants.CLOSED_IMPULSE_S);
        }
        int nWin = command == Command.U_NARMA10
                ? narmaU.length : LaneAConstants.NUM_WINDOWS;
        double tEnd = nWin * LaneAConstants.WINDOW_S;
        if (tStart >= tEnd || tStart + dt <= 0.0) {
            return 0.0;
        }
        int w = (int) Math.floor(tStart / LaneAConstants.WINDOW_S);
        if (w < 0 || w >= nWin) {
            return 0.0;
        }
        double on = w * LaneAConstants.WINDOW_S;
        double off = on + LaneAConstants.PULSE_S;
        return overlap(tStart, dt, on, off);
    }

    public double uOn() {
        if (command == Command.U_ZERO) {
            return 0.0;
        }
        if (command == Command.U_CLOSED_IMPULSE) {
            return 1.0;
        }
        if (command == Command.U_NARMA10) {
            throw new IllegalStateException("U_NARMA10 amplitude is per-window; use deposit");
        }
        return LaneAConstants.U_ON;
    }

    public double commandU(int window, double t) {
        if (command == Command.U_ZERO) {
            return 0.0;
        }
        if (command == Command.U_CLOSED_IMPULSE) {
            return t < LaneAConstants.CLOSED_IMPULSE_S ? 1.0 : 0.0;
        }
        if (command == Command.U_NARMA10) {
            if (window < 0 || window >= narmaU.length) {
                return 0.0;
            }
            double tIn = t - window * LaneAConstants.WINDOW_S;
            if (tIn < 0.0) {
                tIn = 0.0;
            }
            return tIn < LaneAConstants.PULSE_S ? narmaU[window] : 0.0;
        }
        return LaneAConstants.commandU(window, t);
    }

    public double deposit(BSimTransportField field, double tStart, double dt) {
        if (command == Command.U_ZERO) {
            return 0.0;
        }
        if (command == Command.U_NARMA10) {
            int w = (int) Math.floor(tStart / LaneAConstants.WINDOW_S);
            if (w < 0 || w >= narmaU.length) {
                return 0.0;
            }
            double on = w * LaneAConstants.WINDOW_S;
            double off = on + LaneAConstants.PULSE_S;
            double ov = overlap(tStart, dt, on, off);
            double dM = LaneAConstants.J_MAX * narmaU[w] * ov;
            if (dM > 0.0) {
                field.addQuantity(LaneAConstants.I_AC, LaneAConstants.J_AC, LaneAConstants.K_AC, dM);
                commandedMass += dM;
            }
            return dM;
        }
        double ov = onOverlap(tStart, dt);
        double dM = LaneAConstants.J_MAX * uOn() * ov;
        if (dM > 0.0) {
            field.addQuantity(LaneAConstants.I_AC, LaneAConstants.J_AC, LaneAConstants.K_AC, dM);
            commandedMass += dM;
        }
        return dM;
    }

    public double expectedCommandedMass(double tEnd) {
        if (command == Command.U_ZERO) {
            return 0.0;
        }
        if (command == Command.U_CLOSED_IMPULSE) {
            return LaneAConstants.J_MAX * overlap(0.0, tEnd, 0.0, LaneAConstants.CLOSED_IMPULSE_S);
        }
        if (command == Command.U_NARMA10) {
            double mass = 0.0;
            for (int w = 0; w < narmaU.length; w++) {
                double on = w * LaneAConstants.WINDOW_S;
                double off = on + LaneAConstants.PULSE_S;
                mass += LaneAConstants.J_MAX * narmaU[w] * overlap(0.0, tEnd, on, off);
            }
            return mass;
        }
        double mass = 0.0;
        for (int w = 0; w < LaneAConstants.NUM_WINDOWS; w++) {
            double on = w * LaneAConstants.WINDOW_S;
            double off = on + LaneAConstants.PULSE_S;
            mass += LaneAConstants.J_MAX * LaneAConstants.U_ON * overlap(0.0, tEnd, on, off);
        }
        return mass;
    }
}
