package bsim.laneb;

import bsim.transport.BSimTransportField;

/**
 * A0-class immobilized current into attractant or repellent.
 * No internal AC ODE. B0_MOTILITY_MEMORY. Not Lentini TX-TL.
 */
public final class LaneBA0Source {

    public enum Payload { ATTRACTANT, REPELLENT }

    public final double x;
    public final double y;
    public final double z;
    public final int i;
    public final int j;
    public final int k;
    public final Payload payload;
    public final double jMax;
    public final double tOn;
    public final double tOff;
    private double commandedMass;

    public LaneBA0Source(double x, double y, double z, Payload payload,
                         double jMax, double tOn, double tOff) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.i = LaneBIdentity.pdeI(x);
        this.j = LaneBIdentity.pdeJ(y);
        this.k = 0;
        this.payload = payload;
        this.jMax = jMax;
        this.tOn = tOn;
        this.tOff = tOff;
        this.commandedMass = 0.0;
    }

    public double commandedMass() {
        return commandedMass;
    }

    public double deposit(BSimTransportField field, double tStart, double dt) {
        if (jMax == 0.0 || field == null) {
            return 0.0;
        }
        double dM = jMax * overlap(tStart, dt, tOn, tOff);
        if (dM > 0.0) {
            field.addQuantity(i, j, k, dM);
            commandedMass += dM;
        }
        return dM;
    }

    static double overlap(double tStart, double dt, double tOn, double tOff) {
        if (!Double.isFinite(tOn) || !Double.isFinite(tOff) || dt <= 0.0) {
            return 0.0;
        }
        double lo = Math.max(tStart, tOn);
        double hi = Math.min(tStart + dt, tOff);
        return Math.max(0.0, hi - lo);
    }
}
