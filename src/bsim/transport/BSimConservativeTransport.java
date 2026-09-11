package bsim.transport;

/**
 * Conservative explicit finite-volume operators on a regular Cartesian grid.
 * Quantities are molecules per voxel; concentrations are molecules/um^3.
 */
public final class BSimConservativeTransport {
    private static final double ROUND_OFF_TOLERANCE = 1e-12;

    private BSimConservativeTransport() { }

    public static double[][][] diffuse(double[][][] initial, boolean[][][] fluid,
                                       double[] box, double diffusivity, double duration,
                                       BSimTransportBoundary boundary,
                                       BSimTransportLedger ledger) {
        validateInputs(initial, fluid, box, duration);
        requireNonnegativeFinite(diffusivity, "diffusivity");
        boundary.validatePeriodicPairs();
        if (duration == 0.0 || diffusivity == 0.0 && !hasLeak(boundary)) {
            return copy(initial);
        }

        int subcycles = diffusionSubcycles(initial, fluid, box, diffusivity,
                duration, boundary);
        double h = duration / subcycles;
        double[][][] state = copy(initial);
        for (int step = 0; step < subcycles; step++) {
            state = diffuseStableStep(state, fluid, box, diffusivity, h, boundary, ledger);
        }
        return state;
    }

    /**
     * Same FTCS operator as {@link #diffuse}, writing into caller-owned buffers
     * so a long C1 run does not allocate a new field every substep.
     * Returns the buffer that holds the final state ({@code src} or {@code dest}).
     */
    public static double[][][] diffuseReused(
            double[][][] src, double[][][] dest,
            double[][][] fx, double[][][] fy, double[][][] fz,
            boolean[][][] fluid, double[] box, double diffusivity, double duration,
            BSimTransportBoundary boundary, BSimTransportLedger ledger) {
        validateInputs(src, fluid, box, duration);
        requireNonnegativeFinite(diffusivity, "diffusivity");
        boundary.validatePeriodicPairs();
        if (dest == src)
            throw new IllegalArgumentException("diffuseReused dest must be a second buffer");
        if (duration == 0.0 || diffusivity == 0.0 && !hasLeak(boundary)) {
            return src;
        }
        int subcycles = diffusionSubcycles(src, fluid, box, diffusivity,
                duration, boundary);
        double h = duration / subcycles;
        double[][][] a = src;
        double[][][] b = dest;
        for (int step = 0; step < subcycles; step++) {
            diffuseStableStepInto(a, b, fx, fy, fz, fluid, box, diffusivity, h, boundary, ledger);
            double[][][] tmp = a;
            a = b;
            b = tmp;
        }
        return a;
    }

    public static double[][][] advect(double[][][] initial, boolean[][][] fluid,
                                      double[] box, double duration,
                                      BSimFaceVelocity velocity,
                                      BSimTransportBoundary boundary,
                                      BSimTransportLedger ledger) {
        validateInputs(initial, fluid, box, duration);
        boundary.validatePeriodicPairs();
        if (duration == 0.0) return copy(initial);

        int subcycles = advectionSubcycles(initial, fluid, box, duration, velocity, boundary);
        double h = duration / subcycles;
        double[][][] state = copy(initial);
        for (int step = 0; step < subcycles; step++) {
            state = advectStableStep(state, fluid, box, h, velocity, boundary, ledger);
        }
        return state;
    }

    public static int diffusionSubcycles(double[][][] quantity, boolean[][][] fluid,
                                         double[] box, double diffusivity, double duration,
                                         BSimTransportBoundary boundary) {
        validateInputs(quantity, fluid, box, duration);
        requireNonnegativeFinite(diffusivity, "diffusivity");
        boundary.validatePeriodicPairs();
        double inverseSquareSum = 1.0 / (box[0] * box[0])
                + 1.0 / (box[1] * box[1])
                + 1.0 / (box[2] * box[2]);
        double prescribedFtcsLoad = 2.0 * diffusivity * duration * inverseSquareSum;
        double localLoad = duration * maximumDiffusiveOutgoingRate(
                quantity, fluid, box, diffusivity, boundary);
        return requiredSubcycles(Math.max(prescribedFtcsLoad, localLoad));
    }

    public static int advectionSubcycles(double[][][] quantity, boolean[][][] fluid,
                                         double[] box, double duration,
                                         BSimFaceVelocity velocity,
                                         BSimTransportBoundary boundary) {
        validateInputs(quantity, fluid, box, duration);
        boundary.validatePeriodicPairs();
        double load = duration * maximumAdvectiveOutgoingRate(
                quantity, fluid, box, velocity, boundary);
        return requiredSubcycles(load);
    }

    public static boolean[][][] allFluid(int[] boxes) {
        boolean[][][] fluid = new boolean[boxes[0]][boxes[1]][boxes[2]];
        for (int i = 0; i < boxes[0]; i++)
            for (int j = 0; j < boxes[1]; j++)
                for (int k = 0; k < boxes[2]; k++)
                    fluid[i][j][k] = true;
        return fluid;
    }

    public static boolean[][][] copyMask(boolean[][][] source) {
        boolean[][][] out = new boolean[source.length][][];
        for (int i = 0; i < source.length; i++) {
            out[i] = new boolean[source[i].length][];
            for (int j = 0; j < source[i].length; j++) {
                out[i][j] = source[i][j].clone();
            }
        }
        return out;
    }

    private static double[][][] diffuseStableStep(
            double[][][] old, boolean[][][] fluid, double[] box, double diffusivity,
            double h, BSimTransportBoundary boundary, BSimTransportLedger ledger) {
        int nx = old.length, ny = old[0].length, nz = old[0][0].length;
        double[][][] fx = new double[nx + 1][ny][nz];
        double[][][] fy = new double[nx][ny + 1][nz];
        double[][][] fz = new double[nx][ny][nz + 1];
        double[][][] next = copy(old);
        diffuseStableStepInto(old, next, fx, fy, fz, fluid, box, diffusivity, h, boundary, ledger);
        return next;
    }

    private static void diffuseStableStepInto(
            double[][][] old, double[][][] next,
            double[][][] fx, double[][][] fy, double[][][] fz,
            boolean[][][] fluid, double[] box, double diffusivity, double h,
            BSimTransportBoundary boundary, BSimTransportLedger ledger) {
        int nx = old.length, ny = old[0].length, nz = old[0][0].length;
        zero(fx);
        zero(fy);
        zero(fz);
        double volume = box[0] * box[1] * box[2];
        double cx = diffusivity * h / (box[0] * box[0]);
        double cy = diffusivity * h / (box[1] * box[1]);
        double cz = diffusivity * h / (box[2] * box[2]);

        for (int i = 1; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++)
                    if (fluid[i - 1][j][k] && fluid[i][j][k])
                        fx[i][j][k] = cx * (old[i - 1][j][k] - old[i][j][k]);
        for (int i = 0; i < nx; i++)
            for (int j = 1; j < ny; j++)
                for (int k = 0; k < nz; k++)
                    if (fluid[i][j - 1][k] && fluid[i][j][k])
                        fy[i][j][k] = cy * (old[i][j - 1][k] - old[i][j][k]);
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 1; k < nz; k++)
                    if (fluid[i][j][k - 1] && fluid[i][j][k])
                        fz[i][j][k] = cz * (old[i][j][k - 1] - old[i][j][k]);

        diffusionXBoundaries(old, fluid, fx, cx, h, box[0], volume, boundary, ledger);
        diffusionYBoundaries(old, fluid, fy, cy, h, box[1], volume, boundary, ledger);
        diffusionZBoundaries(old, fluid, fz, cz, h, box[2], volume, boundary, ledger);

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++) {
                    if (!fluid[i][j][k]) {
                        next[i][j][k] = old[i][j][k];
                        continue;
                    }
                    double value = old[i][j][k]
                            + fx[i][j][k] - fx[i + 1][j][k]
                            + fy[i][j][k] - fy[i][j + 1][k]
                            + fz[i][j][k] - fz[i][j][k + 1];
                    next[i][j][k] = nonnegative(value, old[i][j][k]);
                }
    }

    private static void diffusionXBoundaries(
            double[][][] old, boolean[][][] fluid, double[][][] flux,
            double coefficient, double h, double spacing, double volume,
            BSimTransportBoundary boundary, BSimTransportLedger ledger) {
        int nx = old.length, ny = old[0].length, nz = old[0][0].length;
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                if (periodic(boundary, BSimTransportBoundary.X_LOW) && nx > 1
                        && fluid[0][j][k] && fluid[nx - 1][j][k]) {
                    double f = coefficient * (old[nx - 1][j][k] - old[0][j][k]);
                    flux[0][j][k] = f;
                    flux[nx][j][k] = f;
                } else {
                    if (fluid[0][j][k])
                        flux[0][j][k] = lowBoundaryFlux(old[0][j][k], coefficient,
                                h, spacing, volume, boundary,
                                BSimTransportBoundary.X_LOW, ledger);
                    if (fluid[nx - 1][j][k])
                        flux[nx][j][k] = highBoundaryFlux(old[nx - 1][j][k],
                                coefficient, h, spacing, volume, boundary,
                                BSimTransportBoundary.X_HIGH, ledger);
                }
            }
    }

    private static void diffusionYBoundaries(
            double[][][] old, boolean[][][] fluid, double[][][] flux,
            double coefficient, double h, double spacing, double volume,
            BSimTransportBoundary boundary, BSimTransportLedger ledger) {
        int nx = old.length, ny = old[0].length, nz = old[0][0].length;
        for (int i = 0; i < nx; i++)
            for (int k = 0; k < nz; k++) {
                if (periodic(boundary, BSimTransportBoundary.Y_LOW) && ny > 1
                        && fluid[i][0][k] && fluid[i][ny - 1][k]) {
                    double f = coefficient * (old[i][ny - 1][k] - old[i][0][k]);
                    flux[i][0][k] = f;
                    flux[i][ny][k] = f;
                } else {
                    if (fluid[i][0][k])
                        flux[i][0][k] = lowBoundaryFlux(old[i][0][k], coefficient,
                                h, spacing, volume, boundary,
                                BSimTransportBoundary.Y_LOW, ledger);
                    if (fluid[i][ny - 1][k])
                        flux[i][ny][k] = highBoundaryFlux(old[i][ny - 1][k],
                                coefficient, h, spacing, volume, boundary,
                                BSimTransportBoundary.Y_HIGH, ledger);
                }
            }
    }

    private static void diffusionZBoundaries(
            double[][][] old, boolean[][][] fluid, double[][][] flux,
            double coefficient, double h, double spacing, double volume,
            BSimTransportBoundary boundary, BSimTransportLedger ledger) {
        int nx = old.length, ny = old[0].length, nz = old[0][0].length;
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++) {
                if (periodic(boundary, BSimTransportBoundary.Z_LOW) && nz > 1
                        && fluid[i][j][0] && fluid[i][j][nz - 1]) {
                    double f = coefficient * (old[i][j][nz - 1] - old[i][j][0]);
                    flux[i][j][0] = f;
                    flux[i][j][nz] = f;
                } else {
                    if (fluid[i][j][0])
                        flux[i][j][0] = lowBoundaryFlux(old[i][j][0], coefficient,
                                h, spacing, volume, boundary,
                                BSimTransportBoundary.Z_LOW, ledger);
                    if (fluid[i][j][nz - 1])
                        flux[i][j][nz] = highBoundaryFlux(old[i][j][nz - 1],
                                coefficient, h, spacing, volume, boundary,
                                BSimTransportBoundary.Z_HIGH, ledger);
                }
            }
    }

    private static double lowBoundaryFlux(
            double inside, double diffusionCoefficient, double h, double spacing,
            double volume, BSimTransportBoundary boundary, int face,
            BSimTransportLedger ledger) {
        BSimTransportBoundary.Type type = boundary.getType(face);
        double flux = 0.0;
        if (type == BSimTransportBoundary.Type.OPEN) {
            double outside = boundary.getExteriorConcentration(face) * volume;
            flux = diffusionCoefficient * (outside - inside);
            recordBoundaryExchange(flux, ledger);
        } else if (type == BSimTransportBoundary.Type.LEAKY) {
            flux = -boundary.getLeakRate(face) * h / (spacing * spacing) * inside;
            ledger.addBoundaryLoss(-flux);
        }
        return flux;
    }

    private static double highBoundaryFlux(
            double inside, double diffusionCoefficient, double h, double spacing,
            double volume, BSimTransportBoundary boundary, int face,
            BSimTransportLedger ledger) {
        BSimTransportBoundary.Type type = boundary.getType(face);
        double flux = 0.0;
        if (type == BSimTransportBoundary.Type.OPEN) {
            double outside = boundary.getExteriorConcentration(face) * volume;
            flux = diffusionCoefficient * (inside - outside);
            recordBoundaryExchange(-flux, ledger);
        } else if (type == BSimTransportBoundary.Type.LEAKY) {
            flux = boundary.getLeakRate(face) * h / (spacing * spacing) * inside;
            ledger.addBoundaryLoss(flux);
        }
        return flux;
    }

    private static void recordBoundaryExchange(double domainChange,
                                               BSimTransportLedger ledger) {
        if (domainChange >= 0.0) ledger.addSource(domainChange);
        else ledger.addBoundaryLoss(-domainChange);
    }

    private static double[][][] advectStableStep(
            double[][][] old, boolean[][][] fluid, double[] box, double h,
            BSimFaceVelocity velocity, BSimTransportBoundary boundary,
            BSimTransportLedger ledger) {
        int nx = old.length, ny = old[0].length, nz = old[0][0].length;
        double[][][] fx = new double[nx + 1][ny][nz];
        double[][][] fy = new double[nx][ny + 1][nz];
        double[][][] fz = new double[nx][ny][nz + 1];
        double volume = box[0] * box[1] * box[2];

        for (int i = 1; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++)
                    if (fluid[i - 1][j][k] && fluid[i][j][k]) {
                        double v = velocity.getX(i, j, k);
                        fx[i][j][k] = h * v / box[0]
                                * (v >= 0.0 ? old[i - 1][j][k] : old[i][j][k]);
                    }
        for (int i = 0; i < nx; i++)
            for (int j = 1; j < ny; j++)
                for (int k = 0; k < nz; k++)
                    if (fluid[i][j - 1][k] && fluid[i][j][k]) {
                        double v = velocity.getY(i, j, k);
                        fy[i][j][k] = h * v / box[1]
                                * (v >= 0.0 ? old[i][j - 1][k] : old[i][j][k]);
                    }
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 1; k < nz; k++)
                    if (fluid[i][j][k - 1] && fluid[i][j][k]) {
                        double v = velocity.getZ(i, j, k);
                        fz[i][j][k] = h * v / box[2]
                                * (v >= 0.0 ? old[i][j][k - 1] : old[i][j][k]);
                    }

        advectionXBoundaries(old, fluid, fx, h, box[0], volume, velocity,
                boundary, ledger);
        advectionYBoundaries(old, fluid, fy, h, box[1], volume, velocity,
                boundary, ledger);
        advectionZBoundaries(old, fluid, fz, h, box[2], volume, velocity,
                boundary, ledger);

        double[][][] next = copy(old);
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++) {
                    if (!fluid[i][j][k]) continue;
                    double value = old[i][j][k]
                            + fx[i][j][k] - fx[i + 1][j][k]
                            + fy[i][j][k] - fy[i][j + 1][k]
                            + fz[i][j][k] - fz[i][j][k + 1];
                    next[i][j][k] = nonnegative(value, old[i][j][k]);
                }
        return next;
    }

    private static void advectionXBoundaries(
            double[][][] old, boolean[][][] fluid, double[][][] flux,
            double h, double spacing, double volume, BSimFaceVelocity velocity,
            BSimTransportBoundary boundary, BSimTransportLedger ledger) {
        int nx = old.length, ny = old[0].length, nz = old[0][0].length;
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                if (periodic(boundary, BSimTransportBoundary.X_LOW) && nx > 1
                        && fluid[0][j][k] && fluid[nx - 1][j][k]) {
                    double v = periodicVelocity(velocity.getX(0, j, k),
                            velocity.getX(nx, j, k));
                    double f = h * v / spacing
                            * (v >= 0.0 ? old[nx - 1][j][k] : old[0][j][k]);
                    flux[0][j][k] = f;
                    flux[nx][j][k] = f;
                } else {
                    if (fluid[0][j][k])
                        flux[0][j][k] = lowAdvectionFlux(old[0][j][k],
                                velocity.getX(0, j, k), h, spacing, volume,
                                boundary, BSimTransportBoundary.X_LOW, ledger);
                    if (fluid[nx - 1][j][k])
                        flux[nx][j][k] = highAdvectionFlux(old[nx - 1][j][k],
                                velocity.getX(nx, j, k), h, spacing, volume,
                                boundary, BSimTransportBoundary.X_HIGH, ledger);
                }
            }
    }

    private static void advectionYBoundaries(
            double[][][] old, boolean[][][] fluid, double[][][] flux,
            double h, double spacing, double volume, BSimFaceVelocity velocity,
            BSimTransportBoundary boundary, BSimTransportLedger ledger) {
        int nx = old.length, ny = old[0].length, nz = old[0][0].length;
        for (int i = 0; i < nx; i++)
            for (int k = 0; k < nz; k++) {
                if (periodic(boundary, BSimTransportBoundary.Y_LOW) && ny > 1
                        && fluid[i][0][k] && fluid[i][ny - 1][k]) {
                    double v = periodicVelocity(velocity.getY(i, 0, k),
                            velocity.getY(i, ny, k));
                    double f = h * v / spacing
                            * (v >= 0.0 ? old[i][ny - 1][k] : old[i][0][k]);
                    flux[i][0][k] = f;
                    flux[i][ny][k] = f;
                } else {
                    if (fluid[i][0][k])
                        flux[i][0][k] = lowAdvectionFlux(old[i][0][k],
                                velocity.getY(i, 0, k), h, spacing, volume,
                                boundary, BSimTransportBoundary.Y_LOW, ledger);
                    if (fluid[i][ny - 1][k])
                        flux[i][ny][k] = highAdvectionFlux(old[i][ny - 1][k],
                                velocity.getY(i, ny, k), h, spacing, volume,
                                boundary, BSimTransportBoundary.Y_HIGH, ledger);
                }
            }
    }

    private static void advectionZBoundaries(
            double[][][] old, boolean[][][] fluid, double[][][] flux,
            double h, double spacing, double volume, BSimFaceVelocity velocity,
            BSimTransportBoundary boundary, BSimTransportLedger ledger) {
        int nx = old.length, ny = old[0].length, nz = old[0][0].length;
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++) {
                if (periodic(boundary, BSimTransportBoundary.Z_LOW) && nz > 1
                        && fluid[i][j][0] && fluid[i][j][nz - 1]) {
                    double v = periodicVelocity(velocity.getZ(i, j, 0),
                            velocity.getZ(i, j, nz));
                    double f = h * v / spacing
                            * (v >= 0.0 ? old[i][j][nz - 1] : old[i][j][0]);
                    flux[i][j][0] = f;
                    flux[i][j][nz] = f;
                } else {
                    if (fluid[i][j][0])
                        flux[i][j][0] = lowAdvectionFlux(old[i][j][0],
                                velocity.getZ(i, j, 0), h, spacing, volume,
                                boundary, BSimTransportBoundary.Z_LOW, ledger);
                    if (fluid[i][j][nz - 1])
                        flux[i][j][nz] = highAdvectionFlux(old[i][j][nz - 1],
                                velocity.getZ(i, j, nz), h, spacing, volume,
                                boundary, BSimTransportBoundary.Z_HIGH, ledger);
                }
            }
    }

    private static double lowAdvectionFlux(
            double inside, double velocity, double h, double spacing, double volume,
            BSimTransportBoundary boundary, int face, BSimTransportLedger ledger) {
        if (boundary.getType(face) != BSimTransportBoundary.Type.OPEN) return 0.0;
        double upwind = velocity >= 0.0
                ? boundary.getExteriorConcentration(face) * volume : inside;
        double flux = h * velocity / spacing * upwind;
        if (flux >= 0.0) ledger.addSource(flux);
        else ledger.addOutletLoss(-flux);
        return flux;
    }

    private static double highAdvectionFlux(
            double inside, double velocity, double h, double spacing, double volume,
            BSimTransportBoundary boundary, int face, BSimTransportLedger ledger) {
        if (boundary.getType(face) != BSimTransportBoundary.Type.OPEN) return 0.0;
        double upwind = velocity >= 0.0
                ? inside : boundary.getExteriorConcentration(face) * volume;
        double flux = h * velocity / spacing * upwind;
        if (flux >= 0.0) ledger.addOutletLoss(flux);
        else ledger.addSource(-flux);
        return flux;
    }

    private static double maximumDiffusiveOutgoingRate(
            double[][][] quantity, boolean[][][] fluid, double[] box,
            double diffusivity, BSimTransportBoundary boundary) {
        int nx = quantity.length, ny = quantity[0].length, nz = quantity[0][0].length;
        double maximum = 0.0;
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++) {
                    if (!fluid[i][j][k]) continue;
                    double rate = diffusionSideRate(i > 0 && fluid[i - 1][j][k],
                            i == 0, nx, fluid[nx - 1][j][k], diffusivity, box[0],
                            boundary, BSimTransportBoundary.X_LOW);
                    rate += diffusionSideRate(i + 1 < nx && fluid[i + 1][j][k],
                            i == nx - 1, nx, fluid[0][j][k], diffusivity, box[0],
                            boundary, BSimTransportBoundary.X_HIGH);
                    rate += diffusionSideRate(j > 0 && fluid[i][j - 1][k],
                            j == 0, ny, fluid[i][ny - 1][k], diffusivity, box[1],
                            boundary, BSimTransportBoundary.Y_LOW);
                    rate += diffusionSideRate(j + 1 < ny && fluid[i][j + 1][k],
                            j == ny - 1, ny, fluid[i][0][k], diffusivity, box[1],
                            boundary, BSimTransportBoundary.Y_HIGH);
                    rate += diffusionSideRate(k > 0 && fluid[i][j][k - 1],
                            k == 0, nz, fluid[i][j][nz - 1], diffusivity, box[2],
                            boundary, BSimTransportBoundary.Z_LOW);
                    rate += diffusionSideRate(k + 1 < nz && fluid[i][j][k + 1],
                            k == nz - 1, nz, fluid[i][j][0], diffusivity, box[2],
                            boundary, BSimTransportBoundary.Z_HIGH);
                    maximum = Math.max(maximum, rate);
                }
        return maximum;
    }

    private static double diffusionSideRate(
            boolean internalFluid, boolean atBoundary, int extent,
            boolean periodicFluid, double diffusivity, double spacing,
            BSimTransportBoundary boundary, int face) {
        if (internalFluid) return diffusivity / (spacing * spacing);
        if (!atBoundary) return 0.0;
        BSimTransportBoundary.Type type = boundary.getType(face);
        if (type == BSimTransportBoundary.Type.PERIODIC)
            return extent > 1 && periodicFluid ? diffusivity / (spacing * spacing) : 0.0;
        if (type == BSimTransportBoundary.Type.OPEN)
            return diffusivity / (spacing * spacing);
        if (type == BSimTransportBoundary.Type.LEAKY)
            return boundary.getLeakRate(face) / (spacing * spacing);
        return 0.0;
    }

    private static double maximumAdvectiveOutgoingRate(
            double[][][] quantity, boolean[][][] fluid, double[] box,
            BSimFaceVelocity velocity, BSimTransportBoundary boundary) {
        int nx = quantity.length, ny = quantity[0].length, nz = quantity[0][0].length;
        double maximum = 0.0;
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++) {
                    if (!fluid[i][j][k]) continue;
                    double rate = 0.0;
                    if (advectionFaceOpen(i > 0 && fluid[i - 1][j][k], i == 0,
                            nx, fluid[nx - 1][j][k], boundary,
                            BSimTransportBoundary.X_LOW))
                        rate += Math.max(0.0, -velocity.getX(i, j, k)) / box[0];
                    if (advectionFaceOpen(i + 1 < nx && fluid[i + 1][j][k], i == nx - 1,
                            nx, fluid[0][j][k], boundary,
                            BSimTransportBoundary.X_HIGH))
                        rate += Math.max(0.0, velocity.getX(i + 1, j, k)) / box[0];
                    if (advectionFaceOpen(j > 0 && fluid[i][j - 1][k], j == 0,
                            ny, fluid[i][ny - 1][k], boundary,
                            BSimTransportBoundary.Y_LOW))
                        rate += Math.max(0.0, -velocity.getY(i, j, k)) / box[1];
                    if (advectionFaceOpen(j + 1 < ny && fluid[i][j + 1][k], j == ny - 1,
                            ny, fluid[i][0][k], boundary,
                            BSimTransportBoundary.Y_HIGH))
                        rate += Math.max(0.0, velocity.getY(i, j + 1, k)) / box[1];
                    if (advectionFaceOpen(k > 0 && fluid[i][j][k - 1], k == 0,
                            nz, fluid[i][j][nz - 1], boundary,
                            BSimTransportBoundary.Z_LOW))
                        rate += Math.max(0.0, -velocity.getZ(i, j, k)) / box[2];
                    if (advectionFaceOpen(k + 1 < nz && fluid[i][j][k + 1], k == nz - 1,
                            nz, fluid[i][j][0], boundary,
                            BSimTransportBoundary.Z_HIGH))
                        rate += Math.max(0.0, velocity.getZ(i, j, k + 1)) / box[2];
                    maximum = Math.max(maximum, rate);
                }
        return maximum;
    }

    private static boolean advectionFaceOpen(
            boolean internalFluid, boolean atBoundary, int extent,
            boolean periodicFluid, BSimTransportBoundary boundary, int face) {
        if (internalFluid) return true;
        if (!atBoundary) return false;
        BSimTransportBoundary.Type type = boundary.getType(face);
        if (type == BSimTransportBoundary.Type.PERIODIC)
            return extent > 1 && periodicFluid;
        return type == BSimTransportBoundary.Type.OPEN;
    }

    private static boolean periodic(BSimTransportBoundary boundary, int face) {
        return boundary.getType(face) == BSimTransportBoundary.Type.PERIODIC;
    }

    private static double periodicVelocity(double low, double high) {
        double scale = Math.max(1.0, Math.max(Math.abs(low), Math.abs(high)));
        if (Math.abs(low - high) > 1e-12 * scale) {
            throw new IllegalArgumentException(
                    "periodic face velocities must match: " + low + " versus " + high);
        }
        return 0.5 * (low + high);
    }

    private static int requiredSubcycles(double load) {
        if (!Double.isFinite(load) || load < 0.0) {
            throw new IllegalArgumentException("stability load must be finite and nonnegative: " + load);
        }
        if (load <= 1.0) return 1;
        if (load > Integer.MAX_VALUE - 1.0) {
            throw new IllegalArgumentException("requested step needs too many stability subcycles: " + load);
        }
        return (int) Math.ceil(load);
    }

    private static boolean hasLeak(BSimTransportBoundary boundary) {
        for (int face = 0; face < 6; face++)
            if (boundary.getType(face) == BSimTransportBoundary.Type.LEAKY
                    && boundary.getLeakRate(face) > 0.0)
                return true;
        return false;
    }

    private static double nonnegative(double value, double reference) {
        if (value >= 0.0) return value;
        double tolerance = ROUND_OFF_TOLERANCE * Math.max(1.0, Math.abs(reference));
        if (value >= -tolerance) return 0.0;
        throw new IllegalStateException(
                "stable transport produced a negative quantity " + value);
    }

    private static void validateInputs(double[][][] quantity, boolean[][][] fluid,
                                       double[] box, double duration) {
        if (quantity == null || quantity.length == 0
                || quantity[0].length == 0 || quantity[0][0].length == 0) {
            throw new IllegalArgumentException("quantity grid must be non-empty");
        }
        if (box == null || box.length != 3) {
            throw new IllegalArgumentException("box must contain three spacings");
        }
        for (double spacing : box)
            if (!Double.isFinite(spacing) || spacing <= 0.0)
                throw new IllegalArgumentException("box spacing must be finite and positive");
        requireNonnegativeFinite(duration, "duration");
        if (fluid.length != quantity.length
                || fluid[0].length != quantity[0].length
                || fluid[0][0].length != quantity[0][0].length) {
            throw new IllegalArgumentException("fluid mask dimensions must match quantity grid");
        }
        for (int i = 0; i < quantity.length; i++)
            for (int j = 0; j < quantity[i].length; j++)
                for (int k = 0; k < quantity[i][j].length; k++)
                    requireNonnegativeFinite(quantity[i][j][k], "quantity");
    }

    private static void requireNonnegativeFinite(double value, String name) {
        if (!Double.isFinite(value) || value < 0.0) {
            throw new IllegalArgumentException(name + " must be finite and nonnegative: " + value);
        }
    }

    private static double[][][] copy(double[][][] source) {
        double[][][] out = new double[source.length][][];
        for (int i = 0; i < source.length; i++) {
            out[i] = new double[source[i].length][];
            for (int j = 0; j < source[i].length; j++) {
                out[i][j] = source[i][j].clone();
            }
        }
        return out;
    }

    private static void zero(double[][][] a) {
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a[i].length; j++)
                java.util.Arrays.fill(a[i][j], 0.0);
    }
}
