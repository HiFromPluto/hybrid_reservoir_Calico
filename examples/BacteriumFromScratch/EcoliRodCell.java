package BacteriumFromScratch;

import javax.vecmath.Vector3d;

import bsim.BSim;
import bsim.capsule.BSimCapsuleBacterium;

/**
 * Isolated E. coli sphero-cylinder. Shape and poles use BSim capsule
 * kinematics only. Growth is Valdez (4)–(5), not BSimCapsuleBacterium.k_growth.
 */
public class EcoliRodCell extends BSimCapsuleBacterium {

    /**
     * How daughters are placed at division.
     *
     * <p>{@code COLLINEAR} — Jobs 2/3/3b. Daughters inherit the mother's axis
     * exactly and both take length ℓ₀. Deterministic and bit-reproducible, but
     * in an OPEN chamber nothing ever breaks the founder's axis, so a colony
     * grows as a 1-D chain (measured: R ~ N, d_centres pinned at 2·ℓ₀, lateral
     * Hertzian contact never engages). Job 3's CLUSTER_BOX hid this because
     * walls forced packing after ~4 cells.
     *
     * <p>{@code SYMMETRY_BROKEN} — Job 3c onward. Small daughter length
     * asymmetry and small axis misalignment at each division, from a seeded
     * RNG. Mechanism per Melke et al. 2010 Methods / Cho et al. 2007;
     * magnitudes are ENGINEERING (see ChassisParameters).
     */
    public enum DivisionMode { COLLINEAR, SYMMETRY_BROKEN }

    /** Placement rule at division. Inherited by daughters. */
    public DivisionMode divisionMode = DivisionMode.COLLINEAR;

    /** Seeded source for SYMMETRY_BROKEN. Shared with daughters. */
    public java.util.Random divisionRng;

    /**
     * Time integrated by {@link #elongateCited} since birth (s).
     * Independent of BSim's tick-then-log clock.
     */
    public double grownSinceBirthS;

    /** Translational velocity. μm/h. Native Valdez/RodCellVSPhage time. */
    public final Vector3d velocity = new Vector3d();

    /** Angular velocity. rad/h. */
    public final Vector3d angularVelocity = new Vector3d();

    /** Accumulator for Valdez net force. Paper native units. */
    public final Vector3d netForce = new Vector3d();

    /** Accumulator for Valdez net torque. */
    public final Vector3d netTorque = new Vector3d();

    public EcoliRodCell(BSim sim, Vector3d pole1, Vector3d pole2) {
        super(sim, pole1, pole2);
        this.radius = ChassisParameters.RADIUS_UM;
        this.L_initial = ChassisParameters.L0_UM;
        this.L = ChassisParameters.L0_UM;
        this.L_th = ChassisParameters.LDIV_UM;
        this.L_max = ChassisParameters.LDIV_UM;
        this.k_growth = 0.0;
        this.grownSinceBirthS = 0.0;
        stretchPolesToLength();
    }

    /**
     * Exact exponential step of Valdez eq. (4) at constant local nutrient.
     * Does not call {@link BSimCapsuleBacterium#grow()}.
     */
    public void elongateCited(double dtS, double nutrientMm) {
        if (dtS < 0.0) {
            throw new IllegalArgumentException("dt must be non-negative");
        }
        double n = Math.max(0.0, nutrientMm);
        this.grownSinceBirthS += dtS;
        this.L = ChassisParameters.analyticLengthUm(this.grownSinceBirthS, n);
        if (this.L < 0.0) {
            this.L = 0.0;
        }
        stretchPolesToLength();
    }

    public boolean shouldDivide() {
        return this.L >= ChassisParameters.LDIV_UM;
    }

    /**
     * Warren Fig. 12B: mother hemisphere centres become daughter centres;
     * each daughter cylindrical length is ℓ₀. Zero length noise (ENGINEERING).
     */
    public EcoliRodCell divideCited() {
        Vector3d axis = new Vector3d();
        axis.sub(this.x2, this.x1);
        double actual = axis.length();
        if (actual < 1e-12) {
            axis.set(1.0, 0.0, 0.0);
            actual = 1.0;
        } else {
            axis.scale(1.0 / actual);
        }

        double lMother = ChassisParameters.L0_UM;
        double lChild = ChassisParameters.L0_UM;
        Vector3d axisMother = axis;
        Vector3d axisChild = axis;

        boolean broken = divisionMode == DivisionMode.SYMMETRY_BROKEN
                && divisionRng != null;
        if (broken) {
            double e = ChassisParameters.DIVISION_LENGTH_ASYM * divisionRng.nextGaussian();
            e = Math.max(-0.5, Math.min(0.5, e));
            lMother = ChassisParameters.L0_UM * (1.0 + e);
            lChild = ChassisParameters.L0_UM * (1.0 - e);
            axisMother = tiltAxis(axis, divisionRng);
            axisChild = tiltAxis(axis, divisionRng);
        }

        Vector3d motherX1 = new Vector3d(this.x1);
        Vector3d motherX2Old = new Vector3d(this.x2);

        Vector3d motherX2 = new Vector3d();
        motherX2.scaleAdd(lMother, axisMother, motherX1);

        Vector3d childX1 = new Vector3d();
        childX1.scaleAdd(-lChild, axisChild, motherX2Old);

        this.initialise(lMother, motherX1, motherX2);
        this.radius = ChassisParameters.RADIUS_UM;
        this.k_growth = 0.0;
        this.grownSinceBirthS = 0.0;
        this.L = lMother;
        stretchPolesToLength();

        EcoliRodCell child = new EcoliRodCell(sim, childX1, new Vector3d(motherX2Old));
        child.grownSinceBirthS = 0.0;
        child.divisionMode = this.divisionMode;
        child.divisionRng = this.divisionRng;
        child.L = lChild;
        child.stretchPolesToLength();
        return child;
    }

    /**
     * Tilt a unit axis by a small Gaussian angle about a uniformly random
     * perpendicular direction. Full 3-D: the chamber geometry, not this
     * routine, decides whether the colony stays a monolayer.
     */
    private static Vector3d tiltAxis(Vector3d axis, java.util.Random rng) {
        Vector3d probe = new Vector3d(rng.nextGaussian(), rng.nextGaussian(),
                rng.nextGaussian());
        Vector3d perp = new Vector3d();
        perp.cross(axis, probe);
        if (perp.length() < 1e-9) {
            probe.set(axis.z, axis.x, axis.y);
            perp.cross(axis, probe);
            if (perp.length() < 1e-9) {
                return new Vector3d(axis);
            }
        }
        perp.normalize();
        double theta = rng.nextGaussian() * ChassisParameters.DIVISION_ANGLE_SD_RAD;
        Vector3d out = new Vector3d(axis);
        out.scale(Math.cos(theta));
        perp.scale(Math.sin(theta));
        out.add(perp);
        double n = out.length();
        if (n < 1e-12) {
            return new Vector3d(axis);
        }
        out.scale(1.0 / n);
        return out;
    }

    public double analyticLengthUm(double nutrientMm) {
        return ChassisParameters.analyticLengthUm(this.grownSinceBirthS, nutrientMm);
    }

    public boolean hasNegativeState(double nutrientMm) {
        return this.L < 0.0 || this.radius < 0.0 || nutrientMm < 0.0;
    }

    public Vector3d centre() {
        Vector3d c = new Vector3d();
        c.interpolate(this.x1, this.x2, 0.5);
        this.position.set(c);
        return c;
    }

    void stretchPolesToLength() {
        Vector3d axis = new Vector3d();
        axis.sub(this.x2, this.x1);
        double actual = axis.length();
        if (actual < 1e-12) {
            axis.set(1.0, 0.0, 0.0);
        } else {
            axis.scale(1.0 / actual);
        }
        Vector3d centre = new Vector3d();
        centre.interpolate(this.x1, this.x2, 0.5);
        this.x1.scaleAdd(-0.5 * this.L, axis, centre);
        this.x2.scaleAdd(0.5 * this.L, axis, centre);
        this.position.set(centre);
    }
}
