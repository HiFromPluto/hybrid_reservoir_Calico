package BacteriumFromScratch;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import javax.vecmath.Vector3d;

import bsim.BSim;
import bsim.BSimTicker;

/**
 * Job 3c step 4 — morphology smoke test.
 *
 * <p>Does a colony seeded from ONE founder in an OPEN chamber grow as a 2-D
 * disc or a 1-D filament? Nutrient is a uniform bath here on purpose: this
 * isolates division geometry + Hertzian packing from the PDE.
 *
 * <p>Discriminator: a disc gives R ~ sqrt(N); a chain gives R ~ N. Also
 * reported is the minimum centre distance -- lateral contact means Hertzian is
 * actually engaging, which never happened under COLLINEAR.
 */
public final class Job3cMorph {

    private Job3cMorph() {}

    static final double BX = 60.0, BY = 60.0;
    static double BZ = 4.0;
    static final int N_TARGET = 256;

    private static Vector3d com(List<EcoliRodCell> cells) {
        Vector3d c = new Vector3d();
        for (EcoliRodCell cell : cells) {
            c.add(cell.centre());
        }
        if (!cells.isEmpty()) {
            c.scale(1.0 / cells.size());
        }
        return c;
    }

    private static double radius(List<EcoliRodCell> cells, Vector3d c) {
        double r = 0.0;
        for (EcoliRodCell cell : cells) {
            r = Math.max(r, Math.hypot(cell.x1.x - c.x, cell.x1.y - c.y));
            r = Math.max(r, Math.hypot(cell.x2.x - c.x, cell.x2.y - c.y));
        }
        return r;
    }

    /** Fraction of cell centres lying outside the mid-plane +/- one radius. */
    private static double outOfPlane(List<EcoliRodCell> cells) {
        int out = 0;
        for (EcoliRodCell cell : cells) {
            if (Math.abs(cell.centre().z - BZ / 2.0) > ChassisParameters.RADIUS_UM) {
                out++;
            }
        }
        return cells.isEmpty() ? 0.0 : out / (double) cells.size();
    }

    private static void run(EcoliRodCell.DivisionMode mode, double bz) {
        BZ = bz;
        BSim sim = new BSim();
        sim.setDt(ChassisParameters.DT_S);
        sim.setSimulationTime(40000.0);
        sim.setBound(BX, BY, BZ);
        sim.setSolid(true, true, true);

        double halfL = ChassisParameters.L0_UM / 2.0;
        final List<EcoliRodCell> cells = new ArrayList<EcoliRodCell>();
        EcoliRodCell founder = new EcoliRodCell(sim,
                new Vector3d(BX / 2.0 - halfL, BY / 2.0, BZ / 2.0),
                new Vector3d(BX / 2.0 + halfL, BY / 2.0, BZ / 2.0));
        founder.divisionMode = mode;
        founder.divisionRng = new java.util.Random(ChassisParameters.RNG_SEED);
        cells.add(founder);

        final Vector3d bound = sim.getBound();
        final double bath = ChassisParameters.NUTRIENT_BATH_MM;

        System.out.println();
        System.out.println("mode = " + mode
                + "   chamber " + BX + "x" + BY + "x" + BZ + " um, uniform bath");
        System.out.printf(Locale.US, "  %6s %9s %11s %11s %8s %10s %10s %9s%n",
                "N", "R_um", "R~sqrt(N)", "R~N chain", "R/Rdisc", "d_min_um", "dcc_max", "offplane");

        int nextReport = 2;
        double t = 0.0;
        while (cells.size() < N_TARGET && t < 40000.0) {
            for (EcoliRodCell cell : cells) {
                cell.elongateCited(ChassisParameters.DT_S, bath);
            }
            ValdezHertzian.relaxContacts(cells, bound);
            List<EcoliRodCell> newborns = new ArrayList<EcoliRodCell>();
            for (EcoliRodCell cell : cells) {
                if (cell.shouldDivide()) {
                    newborns.add(cell.divideCited());
                }
            }
            cells.addAll(newborns);
            if (!newborns.isEmpty()) {
                ValdezHertzian.relaxContacts(cells, bound);
            }
            t += ChassisParameters.DT_S;

            if (cells.size() >= nextReport) {
                Vector3d c = com(cells);
                double R = radius(cells, c);
                int n = cells.size();
                double rDisc = Math.sqrt(n * 3.0 / 0.85 / Math.PI);
                double rChain = n * 3.0 / 4.0;
                ValdezHertzian.PackingStats ps = ValdezHertzian.stats(cells, bound);
                System.out.printf(Locale.US,
                        "  %6d %9.2f %11.2f %11.2f %8.2f %10.3f %10.4f %9.2f%n",
                        n, R, rDisc, rChain, R / rDisc, ps.minCentreDist,
                        ps.maxDeltaCc, outOfPlane(cells));
                nextReport *= 2;
            }
        }
        System.out.printf(Locale.US, "  reached N = %d at t = %.0f s%n", cells.size(), t);
    }

    public static void main(String[] args) {
        System.out.println("Job 3c morphology smoke test  (uniform bath, no PDE)");
        System.out.println("  disc  -> R ~ sqrt(N), R/Rdisc ~ 1, lateral contact engages");
        System.out.println("  chain -> R ~ N,       R/Rdisc grows, d_min pinned at 2*L0");
        run(EcoliRodCell.DivisionMode.SYMMETRY_BROKEN, 4.0);
        run(EcoliRodCell.DivisionMode.SYMMETRY_BROKEN, 1.0);
    }
}
