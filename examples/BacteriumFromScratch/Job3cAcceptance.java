package BacteriumFromScratch;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import javax.vecmath.Vector3d;

import bsim.BSim;
import bsim.BSimUtils;

/**
 * Job 3c step 3 — acceptance study for the DISH_LATERAL chamber.
 *
 * <p>ENGINEERING. Not a gate. Derives {@code conv_tol} / {@code flux_tol} at
 * {@code dx = 1.0} and runs the mandated refinement sweep.
 *
 * <p>Design note: {@code flux_tol} is set to a NON-BINDING 1.0 throughout, so
 * acceptance is decided purely by the two-consecutive-sweep {@code conv_tol}
 * rule and the resulting flux residual is an INDEPENDENT measurement. Letting
 * flux_tol gate acceptance while also reporting it is the circularity that
 * weakened the 3b.2/3b.3 gate; it is deliberately avoided here.
 *
 * <p>The colony sweep includes very small N because the probe showed the
 * residual pinning at its ceiling when consumption is small — the early ticks,
 * not the converged field, are the hard case.
 */
public final class Job3cAcceptance {

    private Job3cAcceptance() {}

    static final double BX = 60.0, BY = 60.0, BZ = 4.0;
    static final String DIR = "results/job3c_acceptance/";
    static final double FLUX_NON_BINDING = 1.0;

    /** Frozen from the Job 3c omega bracket. ENGINEERING; see JOB3C_STANDING. */
    static final double OMEGA = 1.96;

    private static BSim sim() {
        BSim s = new BSim();
        s.setDt(ChassisParameters.DT_S);
        s.setSimulationTime(100.0);
        s.setBound(BX, BY, BZ);
        s.setSolid(true, true, true);
        return s;
    }

    /** Deterministic disc of rods, radius R, centred in the chamber. */
    private static List<EcoliRodCell> disc(BSim s, double R) {
        List<EcoliRodCell> cells = new ArrayList<EcoliRodCell>();
        double cx = BX / 2.0, cy = BY / 2.0, cz = BZ / 2.0;
        double halfL = ChassisParameters.L0_UM / 2.0;
        if (R <= 0.0) {
            cells.add(new EcoliRodCell(s,
                    new Vector3d(cx - halfL, cy, cz),
                    new Vector3d(cx + halfL, cy, cz)));
            return cells;
        }
        for (double y = cy - R; y <= cy + R + 1e-9; y += 1.2) {
            for (double x = cx - R; x <= cx + R + 1e-9; x += 2.4) {
                double ddx = x - cx, ddy = y - cy;
                if (ddx * ddx + ddy * ddy <= R * R) {
                    cells.add(new EcoliRodCell(s,
                            new Vector3d(x - halfL, y, cz),
                            new Vector3d(x + halfL, y, cz)));
                }
            }
        }
        return cells;
    }

    public static void main(String[] args) {
        BSimUtils.generateDirectoryPath(DIR);
        // Tiered: the full ladder is affordable only on the working grid.
        // dx = 0.25 is ~1M voxels and exists only to confirm that N_centre
        // has converged, so it runs one representative point.
        double[] dxs = {1.0, 0.5, 0.25};
        double[][] radiiFor = {
            {0.0, 4.0, 12.0, 16.0},
            {0.0, 16.0},
            {16.0},
        };
        double[][] convFor = {
            {1e-4, 1e-6, 1e-8, 1e-10, 1e-12},
            {1e-6, 1e-8, 1e-10},
            {1e-8},
        };

        PrintStream csv;
        try {
            csv = new PrintStream(DIR + "acceptance_sweep.csv", "UTF-8");
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        csv.println("dx_um;pad;rim_um;voxels;N_cells;R_um;conv_tol;omega;occupancy;gs_iter;"
                + "gs_converged;flux_rel_err;solve_ms;marked_um3;N_centre_mM;N_edge_mM");

        System.out.println("Job 3c acceptance study   chamber "
                + BX + " x " + BY + " x " + BZ);
        System.out.println("  DISH_LATERAL   CONSERVATIVE_MASS   omega = " + OMEGA
                + "   flux_tol NON-BINDING (1.0)");

        for (int di = 0; di < dxs.length; di++) {
            double dx = dxs[di];
            double[] radii = radiiFor[di];
            double[] convTols = convFor[di];
            int pad = (int) Math.round(1.0 / dx);
            System.out.println();
            System.out.printf(Locale.US, "dx = %.2f um   pad = %d   rim = %.2f um%n",
                    dx, pad, pad * dx);
            System.out.printf(Locale.US, "  %6s %8s %10s %9s %14s %11s %14s%n",
                    "R_um", "N", "conv_tol", "gs_iter", "flux_rel_err", "solve_ms", "N_centre");
            for (double R : radii) {
                BSim s = sim();
                List<EcoliRodCell> cells = disc(s, R);
                for (double ct : convTols) {
                    NutrientField f = new NutrientField(BX, BY, BZ, false, true, dx, pad,
                            NutrientField.BcMode.DISH_LATERAL);
                    f.setOccupancyMode(NutrientField.OccupancyMode.CONSERVATIVE_MASS);
                    f.setOverRelaxation(OMEGA);
                    f.setAcceptanceOverride(ct, FLUX_NON_BINDING);
                    f.markOccupancy(cells);
                    long t0 = System.nanoTime();
                    f.solveQuasiSteady();
                    double ms = (System.nanoTime() - t0) / 1.0e6;
                    double flux = f.fluxBalance().relativeError;
                    double nc = f.sampleMm(new Vector3d(BX / 2.0, BY / 2.0, BZ / 2.0));
                    double best = -1.0, ne = Double.NaN;
                    for (EcoliRodCell c : cells) {
                        Vector3d p = c.centre();
                        double d = Math.hypot(p.x - BX / 2.0, p.y - BY / 2.0);
                        if (d > best) { best = d; ne = f.sampleMm(p); }
                    }
                    System.out.printf(Locale.US,
                            "  %6.1f %8d %10.0e %9d %14.3e %11.1f %14.9f%n",
                            R, cells.size(), ct, f.lastSolveIter(), flux, ms, nc);
                    csv.printf(Locale.US,
                            "%.2f;%d;%.2f;%d;%d;%.1f;%.0e;%.3f;%s;%d;%d;%.6e;%.1f;%.4f;%.9f;%.9f%n",
                            dx, pad, pad * dx, f.voxelCount(), cells.size(), R, ct,
                            OMEGA, f.occupancyMode(),
                            f.lastSolveIter(), f.lastSolveConverged() ? 1 : 0,
                            flux, ms, f.markedMassUm3(), nc, ne);
                }
            }
        }
        csv.close();
        System.out.println();
        System.out.println("wrote " + DIR + "acceptance_sweep.csv");
    }
}
