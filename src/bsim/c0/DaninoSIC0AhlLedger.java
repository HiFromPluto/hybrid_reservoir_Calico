package bsim.c0;

import bsim.transport.BSimTransportLedger;

/**
 * Closed AHL mass ledger for Gate C0. NOT_FIG4B.
 *
 * Residual: {@code initial + sources - mu_loss - gammaH_loss - remaining}.
 * Membrane into cells and into extra is equal-and-opposite bookkeeping
 * and is omitted from the residual. {@link BSimTransportLedger} records
 * extracellular {@code mu} decay of the single He compartment.
 *
 * Enzymatic {@code gammaH} is intracellular. It is not an extra
 * extracellular sink.
 */
public final class DaninoSIC0AhlLedger {

    private final BSimTransportLedger extra = new BSimTransportLedger();
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
        extra.reset();
    }

    /**
     * Accumulate one RK4-weighted increment. {@code memIntoCells} is AHL
     * mass into intracellular volumes from the membrane term; extra gets
     * the exact negative.
     */
    public void addWeighted(double synth, double gammaH, double muLoss, double memIntoCells) {
        synthSource += synth;
        gammaHLoss += gammaH;
        extra.addDecayLoss(muLoss);
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

    public double muLoss() {
        return extra.getDecayLoss();
    }

    public double membraneIntoCells() {
        return membraneIntoCells;
    }

    public double membraneIntoExtra() {
        return membraneIntoExtra;
    }

    /** Should be ~0; membrane is equal-and-opposite by construction. */
    public double membraneCancel() {
        return membraneIntoCells + membraneIntoExtra;
    }

    public BSimTransportLedger extraCompartmentLedger() {
        return extra;
    }

    public double residual(double remainingMass) {
        return initialMass + synthSource - extra.getDecayLoss() - gammaHLoss - remainingMass;
    }

    public static double characteristicMass(double initial, double remaining, double sources) {
        return Math.max(Math.max(Math.abs(initial), Math.abs(remaining)), Math.max(Math.abs(sources), 1e-15));
    }
}
