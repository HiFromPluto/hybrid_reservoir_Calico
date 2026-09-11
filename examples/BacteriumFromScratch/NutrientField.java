package BacteriumFromScratch;

import java.util.Arrays;
import java.util.List;

import javax.vecmath.Vector3d;

/**
 * Valdez (16)–(17) / Warren A1.2.2 two-domain nutrient field on a 3D grid.
 * Colony Ω_C: D_C diffusion + Monod consumption (λ_S/Y)ρ₀ N/(N+K_S).
 * Agar Ω_A: D_A diffusion only. Interface: flux-matched faces; N continuous.
 * ENGINEERING BC (declared in ledger): low-index agar faces i=j=k=0 are
 * Dirichlet N = C_s; high faces Neumann. Agar is a corner shell (not Warren's
 * full z&lt;0 slab). Quasi-steady Gauss–Seidel each tick (Warren A1.2.7).
 */
public final class NutrientField {

    /** Colony voxels use D_C; agar padding uses D_A. */
    public enum VoxelKind { AGAR, COLONY }

    /**
     * Boundary-condition geometry.
     *
     * <p>{@code CORNER_SHELL} — Jobs 2/3/3b. Agar padding on the three LOW
     * faces only; Dirichlet C_s at i=0, j=0, k=0; the three HIGH faces
     * flux-free. An asymmetric corner shell, not Warren's z&lt;0 slab. Kept
     * exactly as locked by 3b.3.
     *
     * <p>{@code DISH_LATERAL} — Job 3c. Agar rim on BOTH lateral faces in x
     * and y, none in z; Dirichlet C_s at i=0, i=nx-1, j=0, j=ny-1; the z
     * floor and ceiling flux-free. A microfluidic chamber fed from its rim.
     * Required for a RADIAL gradient: an agar-below BC supplies every cell
     * through the same short vertical path, so no radial gradient can form.
     */
    public enum BcMode { CORNER_SHELL, DISH_LATERAL }

    /**
     * How a capsule is stamped onto the grid as consuming biomass.
     *
     * <p>{@code LEGACY_VOXEL} — Jobs 2/3/3b. The capsule is dilated by
     * {@code +0.5*dx} and every voxel whose CENTRE lands inside gets the full
     * {@code rho_0}. The dilation is a fixed fraction of the GRID, not of the
     * cell, so the consumed volume scales with dx: measured at 48.9x the true
     * capsule volume at dx=2.0, 6.11x at dx=1.0, 1.53x at dx=0.5, 1.43x at
     * dx=0.25. The source term therefore CHANGES UNDER REFINEMENT and no grid
     * study on this mode can converge. Retained only to keep the locked 3b
     * CSVs reproducible.
     *
     * <p>{@code CONSERVATIVE_MASS} — Job 3c onward. No dilation; each voxel
     * receives {@code rho_0 x (occupied volume fraction)}, estimated by
     * sub-sampling, accumulated across cells and capped at {@code rho_0}.
     * Total consumed mass equals the true capsule volume to sub-sampling
     * accuracy, independent of dx.
     */
    public enum OccupancyMode { LEGACY_VOXEL, CONSERVATIVE_MASS }

    /** Sub-samples per axis per voxel for CONSERVATIVE_MASS (8^3 = 512). */
    private static final int OCC_SUBSAMPLE = 8;

    private final double boundX;
    private final double boundY;
    private final double boundZ;
    private final double dx;
    private final BcMode bcMode;
    private OccupancyMode occMode = OccupancyMode.LEGACY_VOXEL;
    private final int loX;
    private final int loY;
    private final int loZ;
    private final int hiX;
    private final int hiY;
    private final int hiZ;
    private final int nx;
    private final int ny;
    private final int nz;
    private final double[] n;
    private final double[] rho;
    private final VoxelKind[] kind;

    /** When true, skip PDE and return C_s everywhere (ρ=0 regression). */
    private final boolean uniformBath;

    /** When true, occupancy is tracked but sink is zero (ρ=0). */
    private final boolean consumptionEnabled;

    private int lastSolveIter;
    private boolean lastSolveConverged;
    private int clampEventCount;

    /**
     * Successive over-relaxation factor. omega = 1.0 is EXACTLY Gauss-Seidel
     * and is the frozen default for Jobs 2/3/3b -- the sweep short-circuits to
     * the plain GS assignment at omega == 1.0, so the locked path is
     * bit-identical by construction, not merely by test.
     *
     * <p>omega > 1 accelerates convergence: plain GS needs O(1/dx^2) sweeps on
     * a 3D Laplace problem, which is why the dx = 0.5 acceptance sweep could
     * not converge within NUTRIENT_MAX_ITER. Optimal omega turns that into
     * O(1/dx). ENGINEERING: the value must be MEASURED per geometry, not taken
     * from the textbook estimate -- this operator carries a Monod sink and two
     * diffusivities.
     */
    private double omega = 1.0;

    /** ENGINEERING study overrides; <=0 means use frozen ChassisParameters. */
    private double convTolOverride = -1.0;
    private double fluxTolOverride = -1.0;

    /**
     * ENGINEERING only: set the over-relaxation factor. 1.0 (default) is plain
     * Gauss-Seidel. Values must satisfy 0 < omega < 2 for SOR to be stable.
     */
    public void setOverRelaxation(double omegaValue) {
        if (!(omegaValue > 0.0) || omegaValue >= 2.0) {
            throw new IllegalArgumentException("omega must lie in (0, 2): " + omegaValue);
        }
        this.omega = omegaValue;
    }

    public double overRelaxation() {
        return omega;
    }

    /** ENGINEERING only: tighten quasi-steady acceptance for a solver study. */
    public void setAcceptanceOverride(double convTolMm, double fluxTolRel) {
        this.convTolOverride = convTolMm;
        this.fluxTolOverride = fluxTolRel;
    }

    public NutrientField(double boundX, double boundY, double boundZ,
                         boolean uniformBath, boolean consumptionEnabled) {
        this(boundX, boundY, boundZ, uniformBath, consumptionEnabled,
             ChassisParameters.NUTRIENT_DX_UM, ChassisParameters.NUTRIENT_AGAR_PAD);
    }

    public NutrientField(double boundX, double boundY, double boundZ,
                         boolean uniformBath, boolean consumptionEnabled,
                         double dxUm, int agarPadVoxels) {
        this(boundX, boundY, boundZ, uniformBath, consumptionEnabled,
             dxUm, agarPadVoxels, BcMode.CORNER_SHELL);
    }

    /**
     * Grid-study constructor. ENGINEERING only. NUTRIENT_AGAR_PAD is a VOXEL
     * count, so refining dx silently thins the physical agar shell (nagar*dx)
     * and moves the Dirichlet plane. Pass both explicitly to hold the physical
     * agar thickness fixed across a refinement study. Job 3b gate runs use the
     * frozen ChassisParameters values via the 5-arg constructor.
     */
    public NutrientField(double boundX, double boundY, double boundZ,
                         boolean uniformBath, boolean consumptionEnabled,
                         double dxUm, int agarPadVoxels, BcMode mode) {
        this.boundX = boundX;
        this.boundY = boundY;
        this.boundZ = boundZ;
        this.uniformBath = uniformBath;
        this.consumptionEnabled = consumptionEnabled;
        this.dx = dxUm;
        this.bcMode = mode;
        if (mode == BcMode.DISH_LATERAL) {
            this.loX = agarPadVoxels; this.hiX = agarPadVoxels;
            this.loY = agarPadVoxels; this.hiY = agarPadVoxels;
            this.loZ = 0;             this.hiZ = 0;
        } else {
            this.loX = agarPadVoxels; this.hiX = 0;
            this.loY = agarPadVoxels; this.hiY = 0;
            this.loZ = agarPadVoxels; this.hiZ = 0;
        }
        this.nx = loX + Math.max(1, (int) Math.ceil(boundX / dx)) + hiX;
        this.ny = loY + Math.max(1, (int) Math.ceil(boundY / dx)) + hiY;
        this.nz = loZ + Math.max(1, (int) Math.ceil(boundZ / dx)) + hiZ;
        int nVox = nx * ny * nz;
        this.n = new double[nVox];
        this.rho = new double[nVox];
        this.kind = new VoxelKind[nVox];
        classifyVoxels();
        fillUniform(ChassisParameters.NUTRIENT_BATH_MM);
    }

    private void classifyVoxels() {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int idx = index(i, j, k);
                    if (i < loX || i >= nx - hiX
                            || j < loY || j >= ny - hiY
                            || k < loZ || k >= nz - hiZ) {
                        kind[idx] = VoxelKind.AGAR;
                    } else {
                        kind[idx] = VoxelKind.COLONY;
                    }
                }
            }
        }
    }

    public void fillUniform(double csMm) {
        Arrays.fill(n, csMm);
    }

    public void clearOccupancy() {
        Arrays.fill(rho, 0.0);
    }

    /** Mark colony voxels overlapping each capsule with constant ρ₀ (Warren v1). */
    public void markOccupancy(List<EcoliRodCell> cells) {
        clearOccupancy();
        if (uniformBath || !consumptionEnabled) {
            return;
        }
        for (EcoliRodCell cell : cells) {
            markCapsule(cell);
        }
    }

    /** Voxels whose centres fall inside the sphero-cylinder get ρ₀. */
    private void markCapsule(EcoliRodCell cell) {
        Vector3d axis = new Vector3d();
        axis.sub(cell.x2, cell.x1);
        double axisLen = axis.length();
        if (axisLen < 1e-12) {
            axis.set(1.0, 0.0, 0.0);
            axisLen = 1.0;
        } else {
            axis.scale(1.0 / axisLen);
        }
        Vector3d c1 = cell.x1;
        boolean legacy = (occMode == OccupancyMode.LEGACY_VOXEL);
        // LEGACY dilates by half a voxel (grid dependent -- see OccupancyMode).
        double r = legacy ? cell.radius + 0.5 * dx : cell.radius;
        double halfL = 0.5 * cell.L + 0.5 * dx;
        Vector3d centre = cell.centre();

        int iMin = Math.max(loX, worldToI(centre.x - r - halfL));
        int iMax = Math.min(nx - 1, worldToI(centre.x + r + halfL));
        int jMin = Math.max(loY, worldToJ(centre.y - r - halfL));
        int jMax = Math.min(ny - 1, worldToJ(centre.y + r + halfL));
        int kMin = Math.max(loZ, worldToK(centre.z - r - halfL));
        int kMax = Math.min(nz - 1, worldToK(centre.z + r + halfL));

        for (int i = iMin; i <= iMax; i++) {
            for (int j = jMin; j <= jMax; j++) {
                for (int k = kMin; k <= kMax; k++) {
                    if (kind[index(i, j, k)] != VoxelKind.COLONY) {
                        continue;
                    }
                    double px = (i - loX + 0.5) * dx;
                    double py = (j - loY + 0.5) * dx;
                    double pz = (k - loZ + 0.5) * dx;
                    if (legacy) {
                        if (insideCapsule(px, py, pz, c1, axis, axisLen, r)) {
                            rho[index(i, j, k)] = ChassisParameters.RHO_0_GCDW_PER_UM3;
                        }
                        continue;
                    }
                    double frac = occupiedFraction(px, py, pz, c1, axis, axisLen, r);
                    if (frac > 0.0) {
                        int idx = index(i, j, k);
                        double add = ChassisParameters.RHO_0_GCDW_PER_UM3 * frac;
                        rho[idx] = Math.min(ChassisParameters.RHO_0_GCDW_PER_UM3,
                                rho[idx] + add);
                    }
                }
            }
        }
    }

    /**
     * Fraction of the voxel centred at (px,py,pz) lying inside the capsule,
     * by regular sub-sampling. Grid independent by construction: the estimate
     * converges to the true overlap volume as OCC_SUBSAMPLE grows, and the sum
     * over voxels converges to the true capsule volume for any dx.
     */
    private double occupiedFraction(double px, double py, double pz,
            Vector3d c1, Vector3d axis, double axisLen, double r) {
        int sN = OCC_SUBSAMPLE;
        double step = dx / sN;
        double x0 = px - 0.5 * dx + 0.5 * step;
        double y0 = py - 0.5 * dx + 0.5 * step;
        double z0 = pz - 0.5 * dx + 0.5 * step;
        int hits = 0;
        for (int a = 0; a < sN; a++) {
            for (int b = 0; b < sN; b++) {
                for (int c = 0; c < sN; c++) {
                    if (insideCapsule(x0 + a * step, y0 + b * step, z0 + c * step,
                            c1, axis, axisLen, r)) {
                        hits++;
                    }
                }
            }
        }
        return hits / (double) (sN * sN * sN);
    }

    private static boolean insideCapsule(double px, double py, double pz,
            Vector3d c1, Vector3d axis, double axisLen, double r) {
        double vx = px - c1.x;
        double vy = py - c1.y;
        double vz = pz - c1.z;
        double t = vx * axis.x + vy * axis.y + vz * axis.z;
        if (t < 0.0) {
            return vx * vx + vy * vy + vz * vz <= r * r;
        }
        if (t > axisLen) {
            double wx = px - (c1.x + axis.x * axisLen);
            double wy = py - (c1.y + axis.y * axisLen);
            double wz = pz - (c1.z + axis.z * axisLen);
            return wx * wx + wy * wy + wz * wz <= r * r;
        }
        double cx = vx - t * axis.x;
        double cy = vy - t * axis.y;
        double cz = vz - t * axis.z;
        return cx * cx + cy * cy + cz * cz <= r * r;
    }

    /**
     * Warren Appendix A1.2.7: solve ∂N/∂t = 0 each mechanical step.
     * Returns iteration count (0 if uniform / no consumption).
     */
    public int solveQuasiSteady() {
        if (uniformBath) {
            fillUniform(ChassisParameters.NUTRIENT_BATH_MM);
            lastSolveIter = 0;
            lastSolveConverged = true;
            return 0;
        }
        if (!consumptionEnabled && maxRho() <= 0.0) {
            fillUniform(ChassisParameters.NUTRIENT_BATH_MM);
            lastSolveIter = 0;
            lastSolveConverged = true;
            return 0;
        }

        applyDirichletBoundaries();
        clampEventCount = 0;
        double tol = convTolOverride > 0.0
                ? convTolOverride : ChassisParameters.NUTRIENT_CONV_TOL_MM;
        double fluxTol = fluxTolOverride > 0.0
                ? fluxTolOverride : ChassisParameters.NUTRIENT_FLUX_TOL_REL;
        int needStable = ChassisParameters.NUTRIENT_GS_STABLE_SWEEPS;
        int maxIter = ChassisParameters.NUTRIENT_MAX_ITER;
        int iter = 0;
        int stableSweeps = 0;
        double maxDelta = Double.POSITIVE_INFINITY;
        for (; iter < maxIter; iter++) {
            maxDelta = oneGaussSeidelSweep();
            applyDirichletBoundaries();
            if (maxDelta < tol) {
                stableSweeps++;
            } else {
                stableSweeps = 0;
            }
            if (stableSweeps >= needStable
                    && fluxBalance().relativeError <= fluxTol) {
                iter++;
                break;
            }
        }
        lastSolveIter = iter;
        lastSolveConverged = stableSweeps >= needStable
                && fluxBalance().relativeError <= fluxTol;
        return iter;
    }

    /**
     * One full sweep; returns max |ΔN| over interior voxels (mM).
     *
     * <p>At omega == 1.0 the Gauss-Seidel value is stored verbatim. Writing
     * {@code old + 1.0 * (gs - old)} would be algebraically identical but NOT
     * bit-identical in IEEE-754, which would break the locked 3b CSVs -- hence
     * the explicit branch.
     */
    private double oneGaussSeidelSweep() {
        double maxDelta = 0.0;
        double bath = ChassisParameters.NUTRIENT_BATH_MM;
        boolean plainGs = (omega == 1.0);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    if (isDirichlet(i, j, k)) {
                        continue;
                    }
                    int idx = index(i, j, k);
                    double old = n[idx];
                    double gs = solveVoxel(i, j, k);
                    double neu = plainGs ? gs : old + omega * (gs - old);
                    if (neu < 0.0 || neu > bath) {
                        clampEventCount++;
                    }
                    neu = Math.max(0.0, Math.min(neu, bath));
                    n[idx] = neu;
                    maxDelta = Math.max(maxDelta, Math.abs(neu - old));
                }
            }
        }
        return maxDelta;
    }

    public int lastSolveIter() {
        return lastSolveIter;
    }

    public boolean lastSolveConverged() {
        return lastSolveConverged;
    }

    /**
     * Steady-state check: Dirichlet inflow must equal colony Monod consumption.
     * Units mM·µm³/s (both sides).
     */
    public FluxBalance fluxBalance() {
        double vol = dx * dx * dx;
        double consume = 0.0;
        for (int i = loX; i < nx - hiX; i++) {
            for (int j = loY; j < ny - hiY; j++) {
                for (int k = loZ; k < nz - hiZ; k++) {
                    int idx = index(i, j, k);
                    if (kind[idx] == VoxelKind.COLONY && consumptionEnabled && rho[idx] > 0.0) {
                        consume += ChassisParameters.consumptionSinkMmPerS(n[idx], rho[idx]) * vol;
                    }
                }
            }
        }

        double inflow = 0.0;
        if (bcMode == BcMode.DISH_LATERAL) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    inflow += dirichletInflowFlux(1, j, k, 0, j, k);
                    inflow += dirichletInflowFlux(nx - 2, j, k, nx - 1, j, k);
                }
            }
            for (int i = 0; i < nx; i++) {
                for (int k = 0; k < nz; k++) {
                    inflow += dirichletInflowFlux(i, 1, k, i, 0, k);
                    inflow += dirichletInflowFlux(i, ny - 2, k, i, ny - 1, k);
                }
            }
        } else {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    inflow += dirichletInflowFlux(1, j, k, 0, j, k);
                }
            }
            for (int i = 0; i < nx; i++) {
                for (int k = 0; k < nz; k++) {
                    inflow += dirichletInflowFlux(i, 1, k, i, 0, k);
                }
            }
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    inflow += dirichletInflowFlux(i, j, 1, i, j, 0);
                }
            }
        }

        double denom = Math.max(consume, 1.0e-30);
        double relErr = Math.abs(inflow - consume) / denom;
        return new FluxBalance(inflow, consume, relErr);
    }

    /** Diffusive flux from Dirichlet face into interior neighbour (mM·µm³/s). */
    private double dirichletInflowFlux(int iIn, int jIn, int kIn,
            int iDir, int jDir, int kDir) {
        if (!inBounds(iIn, jIn, kIn)) {
            return 0.0;
        }
        double dFace = faceDiffusivity(iDir, jDir, kDir, iIn, jIn, kIn);
        double nIn = n[index(iIn, jIn, kIn)];
        return dFace * (ChassisParameters.NUTRIENT_BATH_MM - nIn) * dx;
    }

    private double solveVoxel(int i, int j, int k) {
        double sink = 0.0;
        int idx = index(i, j, k);
        if (kind[idx] == VoxelKind.COLONY && consumptionEnabled && rho[idx] > 0.0) {
            sink = ChassisParameters.consumptionSinkMmPerS(n[idx], rho[idx]);
        }

        double sum = 0.0;
        sum += faceContribution(i - 1, j, k, i, j, k);
        sum += faceContribution(i + 1, j, k, i, j, k);
        sum += faceContribution(i, j - 1, k, i, j, k);
        sum += faceContribution(i, j + 1, k, i, j, k);
        sum += faceContribution(i, j, k - 1, i, j, k);
        sum += faceContribution(i, j, k + 1, i, j, k);

        // Recover denom from the six face loop (each face adds dFace once to denom).
        double denom = faceDiffusivity(i, j, k, i - 1, j, k)
                + faceDiffusivity(i, j, k, i + 1, j, k)
                + faceDiffusivity(i, j, k, i, j - 1, k)
                + faceDiffusivity(i, j, k, i, j + 1, k)
                + faceDiffusivity(i, j, k, i, j, k - 1)
                + faceDiffusivity(i, j, k, i, j, k + 1);

        if (denom <= 0.0) {
            return n[idx];
        }
        return (sum - sink * dx * dx) / denom;
    }

    private double faceContribution(int ni, int nj, int nk, int i, int j, int k) {
        double dFace = faceDiffusivity(i, j, k, ni, nj, nk);
        double nNeighbor = neighborValue(ni, nj, nk, i, j, k);
        return dFace * nNeighbor;
    }

    /** Harmonic mean of D across a face (flux-matched interface). */
    private double faceDiffusivity(int i, int j, int k, int ni, int nj, int nk) {
        if (!inBounds(ni, nj, nk)) {
            if (isNeumannOuter(i, j, k, ni, nj, nk)) {
                return diffusivity(i, j, k);
            }
            return 0.0;
        }
        double d0 = diffusivity(i, j, k);
        double d1 = diffusivity(ni, nj, nk);
        if (d0 <= 0.0 || d1 <= 0.0) {
            return 0.0;
        }
        return 2.0 * d0 * d1 / (d0 + d1);
    }

    private double diffusivity(int i, int j, int k) {
        if (!inBounds(i, j, k)) {
            return ChassisParameters.D_A_UM2_PER_S;
        }
        return kind[index(i, j, k)] == VoxelKind.AGAR
                ? ChassisParameters.D_A_UM2_PER_S
                : ChassisParameters.D_C_UM2_PER_S;
    }

    private double neighborValue(int ni, int nj, int nk, int i, int j, int k) {
        if (isDirichlet(ni, nj, nk)) {
            return ChassisParameters.NUTRIENT_BATH_MM;
        }
        if (!inBounds(ni, nj, nk)) {
            if (isNeumannOuter(i, j, k, ni, nj, nk)) {
                return n[index(i, j, k)];
            }
            return ChassisParameters.NUTRIENT_BATH_MM;
        }
        return n[index(ni, nj, nk)];
    }

    private boolean isNeumannOuter(int i, int j, int k, int ni, int nj, int nk) {
        boolean oob = ni < 0 || ni >= nx || nj < 0 || nj >= ny || nk < 0 || nk >= nz;
        if (!oob) {
            return false;
        }
        if (bcMode == BcMode.DISH_LATERAL) {
            // only the z floor and ceiling are flux-free; x/y outers are Dirichlet
            return nk < 0 || nk >= nz;
        }
        return i == nx - 1 || j == ny - 1 || k == nz - 1;
    }

    private boolean isDirichlet(int i, int j, int k) {
        if (bcMode == BcMode.DISH_LATERAL) {
            return i == 0 || i == nx - 1 || j == 0 || j == ny - 1;
        }
        return i == 0 || j == 0 || k == 0;
    }

    private void applyDirichletBoundaries() {
        double cs = ChassisParameters.NUTRIENT_BATH_MM;
        if (bcMode == BcMode.DISH_LATERAL) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    n[index(0, j, k)] = cs;
                    n[index(nx - 1, j, k)] = cs;
                }
            }
            for (int i = 0; i < nx; i++) {
                for (int k = 0; k < nz; k++) {
                    n[index(i, 0, k)] = cs;
                    n[index(i, ny - 1, k)] = cs;
                }
            }
            return;
        }
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                n[index(0, j, k)] = cs;
            }
        }
        for (int i = 0; i < nx; i++) {
            for (int k = 0; k < nz; k++) {
                n[index(i, 0, k)] = cs;
            }
        }
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                n[index(i, j, 0)] = cs;
            }
        }
    }

    /** Trilinear sample at a physical point (µm). */
    public double sampleMm(Vector3d point) {
        if (uniformBath) {
            return ChassisParameters.NUTRIENT_BATH_MM;
        }
        double fx = point.x / dx + loX - 0.5;
        double fy = point.y / dx + loY - 0.5;
        double fz = point.z / dx + loZ - 0.5;
        int i0 = (int) Math.floor(fx);
        int j0 = (int) Math.floor(fy);
        int k0 = (int) Math.floor(fz);
        double tx = fx - i0;
        double ty = fy - j0;
        double tz = fz - k0;
        double c000 = valueAtGrid(i0, j0, k0);
        double c100 = valueAtGrid(i0 + 1, j0, k0);
        double c010 = valueAtGrid(i0, j0 + 1, k0);
        double c110 = valueAtGrid(i0 + 1, j0 + 1, k0);
        double c001 = valueAtGrid(i0, j0, k0 + 1);
        double c101 = valueAtGrid(i0 + 1, j0, k0 + 1);
        double c011 = valueAtGrid(i0, j0 + 1, k0 + 1);
        double c111 = valueAtGrid(i0 + 1, j0 + 1, k0 + 1);
        double c00 = c000 * (1 - tx) + c100 * tx;
        double c10 = c010 * (1 - tx) + c110 * tx;
        double c01 = c001 * (1 - tx) + c101 * tx;
        double c11 = c011 * (1 - tx) + c111 * tx;
        double c0 = c00 * (1 - ty) + c10 * ty;
        double c1 = c01 * (1 - ty) + c11 * ty;
        return c0 * (1 - tz) + c1 * tz;
    }

    private double valueAtGrid(int i, int j, int k) {
        if (!inBounds(i, j, k)) {
            return ChassisParameters.NUTRIENT_BATH_MM;
        }
        return n[index(i, j, k)];
    }

    public double dxUm() {
        return dx;
    }

    /** Physical agar shell thickness (µm) = pad voxels * dx. */
    public double agarThicknessUm() {
        return loX * dx;
    }

    public BcMode bcMode() {
        return bcMode;
    }

    public OccupancyMode occupancyMode() {
        return occMode;
    }

    /**
     * ENGINEERING. Default LEGACY_VOXEL reproduces the locked 3b behaviour
     * exactly; CONSERVATIVE_MASS makes total consumption dx-independent.
     */
    public void setOccupancyMode(OccupancyMode mode) {
        this.occMode = mode;
    }

    /** Effective consuming volume actually stamped on the grid (um^3). */
    public double markedMassUm3() {
        double sum = 0.0;
        for (double v : rho) {
            sum += v;
        }
        return sum / ChassisParameters.RHO_0_GCDW_PER_UM3 * dx * dx * dx;
    }

    /** Count of voxels currently carrying rho_0 (ENGINEERING diagnostic). */
    public int markedVoxelCount() {
        int c = 0;
        for (double v : rho) {
            if (v > 0.0) {
                c++;
            }
        }
        return c;
    }

    public int voxelCount() {
        return nx * ny * nz;
    }

    public FieldStats stats() {
        double min = Double.POSITIVE_INFINITY;
        double sum = 0.0;
        int count = 0;
        for (int i = loX; i < nx - hiX; i++) {
            for (int j = loY; j < ny - hiY; j++) {
                for (int k = loZ; k < nz - hiZ; k++) {
                    double v = n[index(i, j, k)];
                    min = Math.min(min, v);
                    sum += v;
                    count++;
                }
            }
        }
        if (count == 0) {
            return new FieldStats(ChassisParameters.NUTRIENT_BATH_MM,
                    ChassisParameters.NUTRIENT_BATH_MM);
        }
        return new FieldStats(min, sum / count);
    }

    public int countClampEvents() {
        return clampEventCount;
    }

    public int countNonFinite() {
        int bad = 0;
        for (double v : n) {
            if (!Double.isFinite(v) || v < 0.0) {
                bad++;
            }
        }
        return bad;
    }

    public double maxRho() {
        double m = 0.0;
        for (double r : rho) {
            m = Math.max(m, r);
        }
        return m;
    }

    public boolean isUniformBath() {
        return uniformBath;
    }

    private int worldToI(double x) {
        return (int) Math.floor(x / dx) + loX;
    }

    private int worldToJ(double y) {
        return (int) Math.floor(y / dx) + loY;
    }

    private int worldToK(double z) {
        return (int) Math.floor(z / dx) + loZ;
    }

    private boolean inBounds(int i, int j, int k) {
        return i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz;
    }

    private int index(int i, int j, int k) {
        return i + nx * (j + ny * k);
    }

    public static final class FieldStats {
        public final double minMm;
        public final double meanMm;

        public FieldStats(double minMm, double meanMm) {
            this.minMm = minMm;
            this.meanMm = meanMm;
        }
    }

    public static final class FluxBalance {
        public final double inflowMm3PerS;
        public final double consumeMm3PerS;
        public final double relativeError;

        public FluxBalance(double inflowMm3PerS, double consumeMm3PerS, double relativeError) {
            this.inflowMm3PerS = inflowMm3PerS;
            this.consumeMm3PerS = consumeMm3PerS;
            this.relativeError = relativeError;
        }
    }
}
