package bsim.transport;

/** Cumulative molecule accounting for one chemical field. */
public final class BSimTransportLedger {
    private double sourceAdded;
    private double decayLoss;
    private double outletLoss;
    private double boundaryLoss;

    public double getSourceAdded() {
        return sourceAdded;
    }

    public double getDecayLoss() {
        return decayLoss;
    }

    public double getOutletLoss() {
        return outletLoss;
    }

    public double getBoundaryLoss() {
        return boundaryLoss;
    }

    public double residual(double initialQuantity, double currentQuantity) {
        return initialQuantity + sourceAdded
                - decayLoss - outletLoss - boundaryLoss - currentQuantity;
    }

    public void reset() {
        sourceAdded = 0.0;
        decayLoss = 0.0;
        outletLoss = 0.0;
        boundaryLoss = 0.0;
    }

    public void addSource(double quantity) {
        sourceAdded += quantity;
    }

    public void addDecayLoss(double quantity) {
        decayLoss += quantity;
    }

    public void addOutletLoss(double quantity) {
        outletLoss += quantity;
    }

    public void addBoundaryLoss(double quantity) {
        boundaryLoss += quantity;
    }
}
