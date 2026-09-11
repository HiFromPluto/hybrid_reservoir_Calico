package bsim.laned;

/**
 * Scalar sizer body for LANE_D_CLOSED_BATH.
 * Incremental elongation only. Does not call elongateCited (F4).
 */
public final class LaneDClosedBathCell {

    final double x;
    final double y;
    final double z;
    private double lengthUm;

    public LaneDClosedBathCell(double x, double y, double z, double lengthUm) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.lengthUm = lengthUm;
    }

    public double lengthUm() {
        return lengthUm;
    }

    public double volumeUm3() {
        return LaneDD0Identity.capsuleVolumeUm3(lengthUm);
    }

    /**
     * One incremental step: dℓ/dt = σ Monod(C) ℓ at the current C.
     * Not analyticLengthUm(t_since_birth, C). Never shrinks for C ≥ 0.
     *
     * @return true if this step decreased length (F4 violation)
     */
    public boolean elongateIncremental(double dtS, double nutrientMm) {
        if (dtS < 0.0) {
            throw new IllegalArgumentException("dt must be non-negative");
        }
        double before = lengthUm;
        double monod = LaneDD0Identity.monodFactor(nutrientMm);
        lengthUm = lengthUm * Math.exp(LaneDD0Identity.SIGMA_PER_S * monod * dtS);
        return lengthUm < before - 1.0e-15;
    }

    public boolean shouldDivide() {
        return lengthUm + 1.0e-15 >= LaneDD0Identity.ELL_DIV_UM;
    }

    /** Parent resets to ℓ₀. Daughter is constructed by the job at this position. */
    public void resetAfterFission() {
        lengthUm = LaneDD0Identity.ELL0_UM;
    }
}
