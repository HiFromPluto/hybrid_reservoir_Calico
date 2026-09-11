package BacteriumFromScratch;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import javax.vecmath.Vector3d;

import bsim.BSim;

/**
 * Job 3c step 3a — measure the over-relaxation factor.
 *
 * <p>ENGINEERING. Textbook optimal omega for a Laplace box of width n is
 * {@code 2/(1 + sin(pi/n))} ~ 1.90 at n = 62, but this operator carries a
 * Monod sink and two diffusivities, so omega is measured, not assumed.
 *
 * <p>Cost is reported at a FIXED convergence target so the comparison is like
 * for like; flux_tol is left non-binding so the flux residual stays an
 * independent measurement.
 */
public final class Job3cOmega {

    private Job3cOmega() {}

    static final double BX = 60.0, BY = 60.0, BZ = 4.0;
    static final double CONV_TARGET = 1.0e-9;
    static final double FLUX_NON_BINDING = 1.0;

    private static BSim sim() {
        BSim s = new BSim();
        s.setDt(ChassisParameters.DT_S);
        s.setSimulationTime(100.0);
        s.setBound(BX, BY, BZ);
        s.setSolid(true, true, true);
        return s;
    }

    private static List<EcoliRodCell> disc(BSim s, double R) {
        List<EcoliRodCell> cells = new ArrayList<EcoliRodCell>();
        double cx = BX / 2.0, cy = BY / 2.0, cz = BZ / 2.0;
        double halfL = ChassisParameters.L0_UM / 2.0;
        if (R <= 0.0) {
            cells.add(new EcoliRodCell(s, new Vector3d(cx - halfL, cy, cz),
                    new Vector3d(cx + halfL, cy, cz)));
            return cells;
        }
        for (double y = cy - R; y <= cy + R + 1e-9; y += 1.2) {
            for (double x = cx - R; x <= cx + R + 1e-9; x += 2.4) {
                double ddx = x - cx, ddy = y - cy;
                if (ddx * ddx + ddy * ddy <= R * R) {
                    cells.add(new EcoliRodCell(s, new Vector3d(x - halfL, y, cz),
                            new Vector3d(x + halfL, y, cz)));
                }
            }
        }
        return cells;
    }

    public static void main(String[] args) {
        double[] omegas = {1.90, 1.94, 1.96, 1.97, 1.98, 1.985, 1.99, 1.995};
        double[] radii = {0.0, 16.0};
        double[] dxs = {1.0, 0.5};

        System.out.printf(Locale.US,
                "Job 3c omega bracket   chamber %.0f x %.0f x %.0f   conv target %.0e   CONSERVATIVE_MASS   flux_tol non-binding%n",
                BX, BY, BZ, CONV_TARGET);
        for (double dx : dxs) {
            int pad = (int) Math.round(1.0 / dx);
            int n = (int) Math.round(BX / dx);
            double theory = 2.0 / (1.0 + Math.sin(Math.PI / n));
            System.out.println();
            System.out.printf(Locale.US,
                    "dx = %.2f  pad = %d  rim = %.2f um   textbook omega_opt(n=%d) = %.3f%n",
                    dx, pad, pad * dx, n, theory);
            for (double R : radii) {
                BSim s = sim();
                List<EcoliRodCell> cells = disc(s, R);
                System.out.printf(Locale.US, "  colony R = %.0f um, N = %d%n", R, cells.size());
                System.out.printf(Locale.US, "    %7s %9s %9s %14s %14s %9s%n",
                        "omega", "iter", "conv", "flux_rel_err", "N_centre", "ms");
                double bestMs = Double.MAX_VALUE;
                double bestOmega = Double.NaN;
                for (double w : omegas) {
                    NutrientField f = new NutrientField(BX, BY, BZ, false, true, dx, pad,
                            NutrientField.BcMode.DISH_LATERAL);
                    f.setOverRelaxation(w);
                    f.setOccupancyMode(NutrientField.OccupancyMode.CONSERVATIVE_MASS);
                    f.setAcceptanceOverride(CONV_TARGET, FLUX_NON_BINDING);
                    f.markOccupancy(cells);
                    long t0 = System.nanoTime();
                    f.solveQuasiSteady();
                    double ms = (System.nanoTime() - t0) / 1.0e6;
                    double flux = f.fluxBalance().relativeError;
                    double nc = f.sampleMm(new Vector3d(BX / 2.0, BY / 2.0, BZ / 2.0));
                    boolean conv = f.lastSolveConverged();
                    if (conv && ms < bestMs) { bestMs = ms; bestOmega = w; }
                    System.out.printf(Locale.US, "    %7.3f %9d %9s %14.3e %14.9f %9.1f%n",
                            w, f.lastSolveIter(), conv ? "yes" : "NO", flux, nc, ms);
                }
                System.out.printf(Locale.US, "    -> fastest converged: omega = %.3f at %.1f ms%n",
                        bestOmega, bestMs);
            }
        }
    }
}
