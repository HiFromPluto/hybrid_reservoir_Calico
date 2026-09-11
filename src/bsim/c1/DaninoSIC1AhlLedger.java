package bsim.c1;

import bsim.transport.BSimTransportLedger;

/**
 * Closed AHL mass ledger for Gate C1. NOT_FIG4B.
 *
 * Residual: {@code initial + sources - mu_loss - gammaH_loss - outlet
 * - boundary - remaining}. Membrane into cells and into voxels is
 * equal-and-opposite bookkeeping and is omitted from the residual.
 * Extracellular {@code mu} decay, outlet, and boundary come from the
 * N0 {@link BSimTransportLedger} of the spatial field.
 */
public final class DaninoSIC1AhlLedger {

    private double initialMass;
    private double synthSource;
    private double gammaHLoss;
    private double membraneIntoCells;
    private double membraneIntoExtra;

    public void reset(double initialMass) {
        this.initialMass = initialMass;
        this.synthSource = 0.0;
        this.gammaHLoss = 0.0;
        this.membraneIntoCells = 0.0;
        this.membraneIntoExtra = 0.0;
    }

    public void addWeighted(double synth, double gammaH, double memIntoCells) {
        synthSource += synth;
        gammaHLoss += gammaH;
        membraneIntoCells += memIntoCells;
        membraneIntoExtra -= memIntoCells;
    }

    public double initialMass() {
        return initialMass;
    }

    public double synthSource() {
        return synthSource;
    }

    public double gammaHLoss() {
        return gammaHLoss;
    }

    public double membraneIntoCells() {
        return membraneIntoCells;
    }

    public double membraneIntoExtra() {
        return membraneIntoExtra;
    }

    public double membraneCancel() {
        return membraneIntoCells + membraneIntoExtra;
    }

    public double residual(double remainingMass, BSimTransportLedger extra) {
        return initialMass + synthSource
                - extra.getDecayLoss()
                - gammaHLoss
                - extra.getOutletLoss()
                - extra.getBoundaryLoss()
                - remainingMass;
    }

    public double relativeResidual(double remainingMass, BSimTransportLedger extra) {
        double r = residual(remainingMass, extra);
        double mStar = characteristicMass(initialMass, remainingMass, synthSource);
        return Math.abs(r) / mStar;
    }

    public static double characteristicMass(double initial, double remaining, double sources) {
        return Math.max(Math.max(Math.abs(initial), Math.abs(remaining)),
                Math.max(Math.abs(sources), 1e-15));
    }
}
