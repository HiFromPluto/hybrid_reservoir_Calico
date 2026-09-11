package bsim.laneb;

import bsim.transport.BSimTransportField;

/**
 * Distributed slab current into attractant (even 100 µm slabs) or
 * repellent (odd slabs). Uniform in y. No point ACs. B3_SUSTAINED_STRIPE.
 */
public final class LaneBSlabSource {

    public final double jSlab;
    public final double tOn;
    public final double tOff;
    private double commandedAtt;
    private double commandedRep;

    public LaneBSlabSource(double jSlab, double tOn, double tOff) {
        this.jSlab = jSlab;
        this.tOn = tOn;
        this.tOff = tOff;
        this.commandedAtt = 0.0;
        this.commandedRep = 0.0;
    }

    public double commandedAtt() {
        return commandedAtt;
    }

    public double commandedRep() {
        return commandedRep;
    }

    public void deposit(BSimTransportField att, BSimTransportField rep, double tStart, double dt) {
        if (jSlab == 0.0 || att == null || rep == null) {
            return;
        }
        double live = LaneBA0Source.overlap(tStart, dt, tOn, tOff);
        if (live <= 0.0) {
            return;
        }
        double dM = jSlab * LaneBIdentity.DX * LaneBIdentity.DY * LaneBIdentity.DZ * live;
        for (int i = 0; i < LaneBIdentity.NX; i++) {
            double xc = (i + 0.5) * LaneBIdentity.DX;
            boolean even = LaneB2Identity.evenAttSlab(xc);
            for (int j = 0; j < LaneBIdentity.NY; j++) {
                if (even) {
                    att.addQuantity(i, j, 0, dM);
                    commandedAtt += dM;
                } else {
                    rep.addQuantity(i, j, 0, dM);
                    commandedRep += dM;
                }
            }
        }
    }
}
