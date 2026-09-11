package ChassisPocket;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import javax.vecmath.Vector3d;

import BacteriumFromScratch.ChassisParameters;
import BacteriumFromScratch.EcoliRodCell;
import BacteriumFromScratch.ValdezHertzian;

import bsim.BSim;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.export.BSimLogger;

/**
 * Job I0 — packed monolayer fill in a closed PocketDish-A-class garage.
 *
 * <p>Every parameter is frozen in PROTOCOL.md before this file existed.
 * Chamber 100x100x1 um, neck closed, uniform Warren bath (no PDE),
 * SYMMETRY_BROKEN division seed 101, Job 3 Hertzian walls and contacts.
 *
 * <p>OFF: AHL, Hill R/L, LuxI, NARMA, motility, death, nutrient PDE.
 * This is not HybridDish, not a chemostat, not an exact Danino blueprint.
 */
public final class ChassisPocketI0 {

    private ChassisPocketI0() {}

    static final String EXPORT_DIR = "results/i0_seed101/";
    static final String CSV_NAME = "colony_timeseries.csv";

    static final double BX = 100.0;
    static final double BY = 100.0;
    static final double BZ = 1.0;
    static final double T_END = 23400.0;
    static final int SMOKE_N = 32;
    static final double SPILL_TOL_UM = 1.0e-3;
    static final double FILAMENT_R_RDISC = 2.0;

    private static Vector3d centreOfMass(List<EcoliRodCell> cells) {
        Vector3d c = new Vector3d();
        for (EcoliRodCell cell : cells) {
            c.add(cell.centre());
        }
        if (!cells.isEmpty()) {
            c.scale(1.0 / cells.size());
        }
        return c;
    }

    /**
     * Fraction of centres more than half a cell radius off the mid-plane.
     * Same definition as Job 3c. At b_z = 2r this is a regression guard.
     */
    private static double outOfPlane(List<EcoliRodCell> cells) {
        int out = 0;
        for (EcoliRodCell cell : cells) {
            if (Math.abs(cell.centre().z - BZ / 2.0) > 0.5 * ChassisParameters.RADIUS_UM) {
                out++;
            }
        }
        return cells.isEmpty() ? 0.0 : out / (double) cells.size();
    }

    private static double discRadius(int n) {
        double footprint = ChassisParameters.W0_UM
                * (0.5 * (ChassisParameters.L0_UM + ChassisParameters.LDIV_UM)
                   + 2.0 * ChassisParameters.RADIUS_UM);
        return Math.sqrt(n * footprint / 0.85 / Math.PI);
    }

    private static double colonyRadius(List<EcoliRodCell> cells, Vector3d com) {
        double r = 0.0;
        for (EcoliRodCell cell : cells) {
            r = Math.max(r, Math.hypot(cell.x1.x - com.x, cell.x1.y - com.y));
            r = Math.max(r, Math.hypot(cell.x2.x - com.x, cell.x2.y - com.y));
        }
        return r;
    }

    /** Count cells with any pole outside the closed garage. */
    static int spillCount(List<EcoliRodCell> cells) {
        int n = 0;
        for (EcoliRodCell cell : cells) {
            if (poleOutside(cell.x1) || poleOutside(cell.x2)) {
                n++;
            }
        }
        return n;
    }

    private static boolean poleOutside(Vector3d p) {
        return p.x < -SPILL_TOL_UM || p.x > BX + SPILL_TOL_UM
                || p.y < -SPILL_TOL_UM || p.y > BY + SPILL_TOL_UM
                || p.z < -SPILL_TOL_UM || p.z > BZ + SPILL_TOL_UM;
    }

    private static void reportMorphology(PrintStream out, List<EcoliRodCell> cells,
            Vector3d bound, double t) {
        Vector3d com = centreOfMass(cells);
        double R = colonyRadius(cells, com);
        int n = cells.size();
        double rDisc = discRadius(n);
        ValdezHertzian.PackingStats pst = ValdezHertzian.stats(cells, bound);
        out.printf(Locale.US,
                "  t=%7.0f  N=%4d  R=%6.2f  R/R_disc=%5.2f  offplane=%.3f  d_min=%.3f  dcc=%.4f  spill=%d%n",
                t, n, R, rDisc > 0.0 ? R / rDisc : 0.0, outOfPlane(cells),
                pst.minCentreDist, pst.maxDeltaCc, spillCount(cells));
    }

    static void run(boolean smoke) {
        final PrintStream progress = System.out;
        BSim sim = new BSim();
        sim.setDt(ChassisParameters.DT_S);
        sim.setSimulationTime(T_END);
        sim.setTimeFormat("0.00");
        sim.setBound(BX, BY, BZ);
        sim.setSolid(true, true, true);

        double halfL = ChassisParameters.L0_UM / 2.0;
        final List<EcoliRodCell> cells = new ArrayList<EcoliRodCell>();
        EcoliRodCell founder = new EcoliRodCell(sim,
                new Vector3d(BX / 2.0 - halfL, BY / 2.0, BZ / 2.0),
                new Vector3d(BX / 2.0 + halfL, BY / 2.0, BZ / 2.0));
        founder.divisionMode = EcoliRodCell.DivisionMode.SYMMETRY_BROKEN;
        founder.divisionRng = new java.util.Random(ChassisParameters.RNG_SEED);
        cells.add(founder);

        final Vector3d bound = sim.getBound();
        final double bathMm = ChassisParameters.NUTRIENT_BATH_MM;
        final int[] nextReport = {2};
        final boolean[] filament = {false};

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                for (EcoliRodCell cell : cells) {
                    cell.elongateCited(sim.getDt(), bathMm);
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
                if (cells.size() >= nextReport[0]) {
                    reportMorphology(progress, cells, bound, sim.getTime() + sim.getDt());
                    double rDisc = discRadius(cells.size());
                    double R = colonyRadius(cells, centreOfMass(cells));
                    double rr = rDisc > 0.0 ? R / rDisc : 0.0;
                    if (rr > FILAMENT_R_RDISC) {
                        filament[0] = true;
                        progress.println("FILAMENT: R/R_disc=" + rr
                                + "  check SYMMETRY_BROKEN + seed 101. Stopping.");
                        sim.setSimulationTime(Math.max(sim.getDt(), sim.getTime()));
                    }
                    nextReport[0] *= 2;
                }
                if (smoke && cells.size() >= SMOKE_N) {
                    sim.setSimulationTime(Math.max(sim.getDt(), sim.getTime()));
                }
            }
        });

        BSimUtils.generateDirectoryPath(EXPORT_DIR);
        BSimLogger logger = new BSimLogger(sim, EXPORT_DIR + CSV_NAME) {
            @Override
            public void before() {
                super.before();
                write("t_s;seed;N;R_um;R_disc_um;R_Rdisc;offplane;"
                        + "delta_cc_max_um;d_centers_min_um;spill");
                logRow(0.0);
            }

            @Override
            public void during() {
                logRow(sim.getTime());
            }

            private void logRow(double t) {
                ValdezHertzian.PackingStats pst = ValdezHertzian.stats(cells, bound);
                Vector3d com = centreOfMass(cells);
                double R = colonyRadius(cells, com);
                double rDisc = discRadius(cells.size());
                write(String.format(Locale.US,
                        "%.1f;%d;%d;%.6f;%.6f;%.6f;%.6f;%.9e;%.9f;%d",
                        t, ChassisParameters.RNG_SEED, cells.size(), R,
                        rDisc, rDisc > 0.0 ? R / rDisc : 0.0, outOfPlane(cells),
                        pst.maxDeltaCc, pst.minCentreDist, spillCount(cells)));
            }
        };
        logger.setDt(ChassisParameters.LOG_DT_S);
        sim.addExporter(logger);

        System.out.println("ChassisPocket I0 -> " + EXPORT_DIR + CSV_NAME);
        System.out.println("  garage " + BX + "x" + BY + "x" + BZ
                + " um  neck=CLOSED  bath=" + bathMm + " mM  (no PDE)");
        System.out.println("  division=SYMMETRY_BROKEN  seed=" + ChassisParameters.RNG_SEED
                + "  Hertzian Job 3  motility=OFF  death=OFF  AHL/Hill/NARMA=OFF");
        System.out.println("  mode=" + (smoke ? "SMOKE N=" + SMOKE_N : "FULL t_end=" + (int) T_END + " s"));

        PrintStream originalOut = System.out;
        try {
            System.setOut(new PrintStream(OutputStream.nullOutputStream()));
            sim.export();
        } finally {
            System.setOut(originalOut);
        }

        progress.println("I0 finished  t=" + (int) sim.getTime() + "  N=" + cells.size()
                + "  spill=" + spillCount(cells));
        reportMorphology(progress, cells, bound, sim.getTime());
        if (filament[0]) {
            throw new IllegalStateException("I0 filament at N=" + cells.size());
        }
        if (smoke) {
            Vector3d com = centreOfMass(cells);
            double R = colonyRadius(cells, com);
            double rr = R / discRadius(cells.size());
            double off = outOfPlane(cells);
            if (off > 0.05 || rr > FILAMENT_R_RDISC) {
                throw new IllegalStateException("I0 smoke FAIL  offplane=" + off
                        + "  R/R_disc=" + rr);
            }
            progress.printf(Locale.US,
                    "I0 smoke OK  N=%d  offplane=%.3f  R/R_disc=%.2f  (expect ~1.4 falling toward 1)%n",
                    cells.size(), off, rr);
        }
    }

    public static void main(String[] args) {
        boolean smoke = args.length > 0 && "smoke".equalsIgnoreCase(args[0]);
        ChassisParameters.printLedger();
        run(smoke);
    }
}
