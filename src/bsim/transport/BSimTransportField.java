package bsim.transport;

import bsim.BSim;
import bsim.BSimChemicalField;

import javax.vecmath.Vector3d;

/**
 * Geometry-masked conservative advection-diffusion-decay field.
 * Solid mask cells exchange no mass with fluid cells.
 */
public class BSimTransportField extends BSimChemicalField {
    private final boolean[][][] fluid;
    private BSimTransportBoundary boundary;
    private BSimFaceVelocity velocity;
    private double[][][] quantityScratch;
    private double[][][] fluxX;
    private double[][][] fluxY;
    private double[][][] fluxZ;

    public BSimTransportField(BSim sim, int[] boxes, double diffusivity,
                              double decayRate, boolean[][][] fluid) {
        this(sim, boxes, diffusivity, decayRate, fluid,
                BSimTransportBoundary.fromSimulation(sim));
    }

    public BSimTransportField(BSim sim, int[] boxes, double diffusivity,
                              double decayRate, boolean[][][] fluid,
                              BSimTransportBoundary boundary) {
        super(sim, boxes, diffusivity, decayRate);
        if (fluid == null) throw new IllegalArgumentException("fluid mask is required");
        this.fluid = BSimConservativeTransport.copyMask(fluid);
        validateMaskDimensions();
        this.boundary = boundary;
        this.boundary.validatePeriodicPairs();
        this.velocity = new BSimFaceVelocity(boxes);
    }

    public void setBoundary(BSimTransportBoundary boundary) {
        if (boundary == null) throw new IllegalArgumentException("boundary is required");
        boundary.validatePeriodicPairs();
        this.boundary = boundary;
    }

    public BSimTransportBoundary getBoundary() {
        return boundary;
    }

    public void setVelocity(BSimFaceVelocity velocity) {
        if (velocity == null) throw new IllegalArgumentException("velocity is required");
        if (velocity.getNx() != boxes[0] || velocity.getNy() != boxes[1]
                || velocity.getNz() != boxes[2])
            throw new IllegalArgumentException("velocity dimensions must match field boxes");
        this.velocity = velocity;
    }

    public BSimFaceVelocity getVelocity() {
        return velocity;
    }

    public boolean isFluid(int i, int j, int k) {
        return fluid[i][j][k];
    }

    /**
     * In the N0 field, positive addQuantity calls are runtime sources and are
     * therefore ledgered. Call ledger.reset() after constructing initial state.
     */
    @Override
    public void addQuantity(int i, int j, int k, double molecules) {
        if (!fluid[i][j][k])
            throw new IllegalArgumentException("cannot add quantity to a solid voxel");
        requireNonnegativeFinite(molecules, "added quantity");
        quantity[i][j][k] += molecules;
        transportLedger.addSource(molecules);
    }

    @Override
    public void addQuantity(Vector3d position, double molecules) {
        int[] b = boxCoords(position);
        addQuantity(b[0], b[1], b[2], molecules);
    }

    /** Add and account for a nonnegative external source in molecules. */
    @Override
    public void addSourceQuantity(int i, int j, int k, double molecules) {
        addQuantity(i, j, k, molecules);
    }

    /** Add and account for a nonnegative external source at a position. */
    @Override
    public void addSourceQuantity(Vector3d position, double molecules) {
        addQuantity(position, molecules);
    }

    /**
     * Conservative signed mass transfer. Not an external source and not a
     * boundary loss. Used for membrane exchange: equal-and-opposite intracellular
     * mass must not enter {@link BSimTransportLedger#addSource}.
     */
    public void transferQuantity(int i, int j, int k, double molecules) {
        if (!fluid[i][j][k])
            throw new IllegalArgumentException("cannot transfer quantity to a solid voxel");
        if (!Double.isFinite(molecules))
            throw new IllegalArgumentException("transferred quantity must be finite: " + molecules);
        double next = quantity[i][j][k] + molecules;
        if (next < -1e-12)
            throw new IllegalArgumentException(
                    "transfer would make voxel quantity negative: " + next);
        quantity[i][j][k] = next < 0.0 ? 0.0 : next;
    }

    /** Conservative signed mass transfer at a position. */
    public void transferQuantity(Vector3d position, double molecules) {
        int[] b = boxCoords(position);
        transferQuantity(b[0], b[1], b[2], molecules);
    }

    /**
     * Explicitly imposed concentration changes are accounted as external
     * source or boundary removal. This method cannot write solid voxels.
     */
    @Override
    public void setConc(int i, int j, int k, double concentration) {
        if (!fluid[i][j][k])
            throw new IllegalArgumentException("cannot set concentration in a solid voxel");
        requireNonnegativeFinite(concentration, "concentration");
        double replacement = concentration * boxVolume;
        double change = replacement - quantity[i][j][k];
        quantity[i][j][k] = replacement;
        if (change >= 0.0) transportLedger.addSource(change);
        else transportLedger.addBoundaryLoss(-change);
    }

    @Override
    public void setConc(Vector3d position, double concentration) {
        int[] b = boxCoords(position);
        setConc(b[0], b[1], b[2], concentration);
    }

    @Override
    public void setConc(double concentration) {
        requireNonnegativeFinite(concentration, "concentration");
        for (int i = 0; i < boxes[0]; i++)
            for (int j = 0; j < boxes[1]; j++)
                for (int k = 0; k < boxes[2]; k++)
                    if (fluid[i][j][k]) setConc(i, j, k, concentration);
    }

    /** Apply the legacy linear-gradient convention to fluid voxels only. */
    @Override
    public void linearGradient(int axis, double startConc, double endConc) {
        if (axis < 0 || axis > 2)
            throw new IllegalArgumentException("gradient axis must be in [0,2]");
        requireNonnegativeFinite(startConc, "gradient start concentration");
        requireNonnegativeFinite(endConc, "gradient end concentration");
        double gradient = (endConc - startConc) / boxes[axis];
        int[] index = {0, 0, 0};
        for (index[0] = 0; index[0] < boxes[0]; index[0]++)
            for (index[1] = 0; index[1] < boxes[1]; index[1]++)
                for (index[2] = 0; index[2] < boxes[2]; index[2]++)
                    if (fluid[index[0]][index[1]][index[2]])
                        setConc(index[0], index[1], index[2],
                                startConc + index[axis] * gradient);
    }

    @Override
    public void diffuse(double duration) {
        double beforeLoss = transportLedger.getBoundaryLoss();
        ensureScratch();
        double[][][] src = quantity;
        double[][][] dest = quantityScratch;
        double[][][] result = BSimConservativeTransport.diffuseReused(
                src, dest, fluxX, fluxY, fluxZ,
                fluid, box, diffusivity, duration, boundary, transportLedger);
        if (result != src) {
            quantityScratch = src;
            quantity = result;
        }
        lastBoundaryLoss = transportLedger.getBoundaryLoss() - beforeLoss;
    }

    private void ensureScratch() {
        if (quantityScratch != null) {
            return;
        }
        int nx = boxes[0], ny = boxes[1], nz = boxes[2];
        quantityScratch = new double[nx][ny][nz];
        fluxX = new double[nx + 1][ny][nz];
        fluxY = new double[nx][ny + 1][nz];
        fluxZ = new double[nx][ny][nz + 1];
    }

    public void advect() {
        advect(sim.getDt());
    }

    public void advect(double duration) {
        quantity = BSimConservativeTransport.advect(
                quantity, fluid, box, duration, velocity, boundary, transportLedger);
    }

    @Override
    public void decay(double duration) {
        requireNonnegativeFinite(duration, "decay duration");
        requireNonnegativeFinite(decayRate, "decay rate");
        double before = totalFluidQuantity();
        double keep = Math.exp(-decayRate * duration);
        for (int i = 0; i < boxes[0]; i++)
            for (int j = 0; j < boxes[1]; j++)
                for (int k = 0; k < boxes[2]; k++)
                    if (fluid[i][j][k]) quantity[i][j][k] *= keep;
        lastDecayLoss = before - totalFluidQuantity();
        transportLedger.addDecayLoss(lastDecayLoss);
    }

    /**
     * Advance all transport terms for an explicit duration. This is the
     * operation intended for each scheduler half-transport phase.
     */
    public void advance(double duration) {
        diffuse(duration);
        advect(duration);
        decay(duration);
    }

    @Override
    public void update() {
        advance(sim.getDt());
    }

    public double totalFluidQuantity() {
        double total = 0.0;
        for (int i = 0; i < boxes[0]; i++)
            for (int j = 0; j < boxes[1]; j++)
                for (int k = 0; k < boxes[2]; k++)
                    if (fluid[i][j][k]) total += quantity[i][j][k];
        return total;
    }

    public double minimumFluidQuantity() {
        double minimum = Double.POSITIVE_INFINITY;
        for (int i = 0; i < boxes[0]; i++)
            for (int j = 0; j < boxes[1]; j++)
                for (int k = 0; k < boxes[2]; k++)
                    if (fluid[i][j][k]) minimum = Math.min(minimum, quantity[i][j][k]);
        return minimum;
    }

    private void validateMaskDimensions() {
        if (fluid.length != boxes[0]) throw new IllegalArgumentException("mask x size mismatch");
        for (int i = 0; i < boxes[0]; i++) {
            if (fluid[i].length != boxes[1])
                throw new IllegalArgumentException("mask y size mismatch at x=" + i);
            for (int j = 0; j < boxes[1]; j++)
                if (fluid[i][j].length != boxes[2])
                    throw new IllegalArgumentException(
                            "mask z size mismatch at x=" + i + ", y=" + j);
        }
    }

    private static void requireNonnegativeFinite(double value, String name) {
        if (!Double.isFinite(value) || value < 0.0)
            throw new IllegalArgumentException(name + " must be finite and nonnegative: " + value);
    }
}
