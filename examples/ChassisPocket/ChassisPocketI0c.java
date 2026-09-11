package ChassisPocket;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
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
 * Job I0c — DoorWeir at T4/T5 width W=20 µm, one extra doubling.
 *
 * <p>Copy of ChassisPocketI0b. The only allowed change is t_end = 25950 s.
 * W, founder (50, 80), and WallSpec.openNeckYPlus(20, 50) are unchanged.
 * Frozen in PROTOCOL.md before this file existed.
 *
 * <p>OFF: AHL, Hill R/L, LuxI, NARMA, motility, death, nutrient PDE, bus.
 * This is not HybridDish, not a chemostat, not an exact Danino blueprint.
 */
public final class ChassisPocketI0c {

    private ChassisPocketI0c() {}

    static final String EXPORT_DIR = "results/i0c_seed101/";
    static final String CSV_NAME = "colony_timeseries.csv";

    static final double BX = 100.0;
    static final double BY = 100.0;
    static final double BZ = 1.0;
    static final double T_END = 25950.0;
    static final int SMOKE_N = 32;
    static final double SPILL_TOL_UM = 1.0e-3;
    static final double FILAMENT_R_RDISC = 2.0;
    static final double W_NECK_UM = 20.0;
    static final double NECK_CENTER_X = BX / 2.0;
    static final double FOUNDER_Y = 80.0;

    private static final ValdezHertzian.WallSpec WALLS =
            ValdezHertzian.WallSpec.openNeckYPlus(W_NECK_UM, NECK_CENTER_X);

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

    private static double maxPoleY(List<EcoliRodCell> cells) {
        double y = Double.NEGATIVE_INFINITY;
        for (EcoliRodCell cell : cells) {
            y = Math.max(y, cell.x1.y);
            y = Math.max(y, cell.x2.y);
        }
        return cells.isEmpty() ? 0.0 : y;
    }

    /** Remaining cell has a pole on or through the +y face (door test fired). */
    static int doorContact(List<EcoliRodCell> cells) {
        double thresh = BY - ChassisParameters.RADIUS_UM;
        for (EcoliRodCell cell : cells) {
            if (cell.x1.y >= thresh || cell.x2.y >= thresh) {
                return 1;
            }
        }
        return 0;
    }

    static boolean centreLeftThroughYPlus(EcoliRodCell cell) {
        return cell.centre().y > BY + SPILL_TOL_UM;
    }

    /** Centre left through a face that is not the open +y mouth. */
    static boolean centreWallLeak(EcoliRodCell cell) {
        Vector3d c = cell.centre();
        return c.x < -SPILL_TOL_UM || c.x > BX + SPILL_TOL_UM
                || c.y < -SPILL_TOL_UM
                || c.z < -SPILL_TOL_UM || c.z > BZ + SPILL_TOL_UM;
    }

    /**
     * Remove cells whose centre left through +y. Returns spill this call.
     * wallLeak[0] accumulates closed-face leaks (must stay 0).
     */
    static int removeSpilled(List<EcoliRodCell> cells, int[] wallLeak) {
        int spilled = 0;
        Iterator<EcoliRodCell> it = cells.iterator();
        while (it.hasNext()) {
            EcoliRodCell cell = it.next();
            if (centreWallLeak(cell)) {
                wallLeak[0]++;
                it.remove();
            } else if (centreLeftThroughYPlus(cell)) {
                spilled++;
                it.remove();
            }
        }
        return spilled;
    }

    private static void reportMorphology(PrintStream out, List<EcoliRodCell> cells,
            Vector3d bound, double t, int spillCum, int wallLeak) {
        Vector3d com = centreOfMass(cells);
        double R = colonyRadius(cells, com);
        int n = cells.size();
        double rDisc = discRadius(n);
        ValdezHertzian.PackingStats pst = ValdezHertzian.stats(cells, bound, WALLS);
        out.printf(Locale.US,
                "  t=%7.0f  N=%4d  R=%6.2f  R/R_disc=%5.2f  offplane=%.3f  d_min=%.3f  "
                        + "dcc=%.4f  spill_cum=%d  door=%d  y_max=%.2f  wall_leak=%d%n",
                t, n, R, rDisc > 0.0 ? R / rDisc : 0.0, outOfPlane(cells),
                pst.minCentreDist, pst.maxDeltaCc, spillCum, doorContact(cells),
                maxPoleY(cells), wallLeak);
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
                new Vector3d(BX / 2.0 - halfL, FOUNDER_Y, BZ / 2.0),
                new Vector3d(BX / 2.0 + halfL, FOUNDER_Y, BZ / 2.0));
        founder.divisionMode = EcoliRodCell.DivisionMode.SYMMETRY_BROKEN;
        founder.divisionRng = new java.util.Random(ChassisParameters.RNG_SEED);
        cells.add(founder);

        final Vector3d bound = sim.getBound();
        final double bathMm = ChassisParameters.NUTRIENT_BATH_MM;
        final int[] nextReport = {2};
        final boolean[] filament = {false};
        final int[] spillCum = {0};
        final int[] spillSinceLog = {0};
        final int[] wallLeak = {0};

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                for (EcoliRodCell cell : cells) {
                    cell.elongateCited(sim.getDt(), bathMm);
                }
                ValdezHertzian.relaxContacts(cells, bound, WALLS);
                spillSinceLog[0] += removeSpilled(cells, wallLeak);
                List<EcoliRodCell> newborns = new ArrayList<EcoliRodCell>();
                for (EcoliRodCell cell : cells) {
                    if (cell.shouldDivide()) {
                        newborns.add(cell.divideCited());
                    }
                }
                cells.addAll(newborns);
                if (!newborns.isEmpty()) {
                    ValdezHertzian.relaxContacts(cells, bound, WALLS);
                    spillSinceLog[0] += removeSpilled(cells, wallLeak);
                }
                if (cells.size() >= nextReport[0]) {
                    reportMorphology(progress, cells, bound, sim.getTime() + sim.getDt(),
                            spillCum[0] + spillSinceLog[0], wallLeak[0]);
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
                        + "delta_cc_max_um;d_centers_min_um;spill_tick;spill_cum;"
                        + "N_ever;door_contact;y_max_um;wall_leak");
                logRow(0.0);
            }

            @Override
            public void during() {
                logRow(sim.getTime());
            }

            private void logRow(double t) {
                spillCum[0] += spillSinceLog[0];
                int tickSpill = spillSinceLog[0];
                spillSinceLog[0] = 0;
                ValdezHertzian.PackingStats pst = ValdezHertzian.stats(cells, bound, WALLS);
                Vector3d com = centreOfMass(cells);
                double R = colonyRadius(cells, com);
                double rDisc = discRadius(cells.size());
                int n = cells.size();
                write(String.format(Locale.US,
                        "%.1f;%d;%d;%.6f;%.6f;%.6f;%.6f;%.9e;%.9f;%d;%d;%d;%d;%.6f;%d",
                        t, ChassisParameters.RNG_SEED, n, R,
                        rDisc, rDisc > 0.0 ? R / rDisc : 0.0, outOfPlane(cells),
                        pst.maxDeltaCc, pst.minCentreDist, tickSpill, spillCum[0],
                        n + spillCum[0], doorContact(cells), maxPoleY(cells),
                        wallLeak[0]));
            }
        };
        logger.setDt(ChassisParameters.LOG_DT_S);
        sim.addExporter(logger);

        System.out.println("ChassisPocket I0c -> " + EXPORT_DIR + CSV_NAME);
        System.out.println("  garage " + BX + "x" + BY + "x" + BZ
                + " um  W=" + (int) W_NECK_UM + " on +y  founder_y=" + (int) FOUNDER_Y
                + "  bath=" + bathMm + " mM  (no PDE)");
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

        spillCum[0] += spillSinceLog[0];
        progress.println("I0c finished  t=" + (int) sim.getTime() + "  N=" + cells.size()
                + "  spill_cum=" + spillCum[0] + "  door=" + doorContact(cells)
                + "  wall_leak=" + wallLeak[0]);
        reportMorphology(progress, cells, bound, sim.getTime(), spillCum[0], wallLeak[0]);
        if (filament[0]) {
            throw new IllegalStateException("I0c filament at N=" + cells.size());
        }
        if (smoke) {
            Vector3d com = centreOfMass(cells);
            double R = colonyRadius(cells, com);
            double rr = R / discRadius(cells.size());
            double off = outOfPlane(cells);
            if (off > 0.05 || rr > FILAMENT_R_RDISC) {
                throw new IllegalStateException("I0c smoke FAIL  offplane=" + off
                        + "  R/R_disc=" + rr);
            }
            progress.printf(Locale.US,
                    "I0c smoke OK  N=%d  offplane=%.3f  R/R_disc=%.2f  spill_cum=%d  y_max=%.2f%n",
                    cells.size(), off, rr, spillCum[0], maxPoleY(cells));
        }
    }

    public static void main(String[] args) {
        boolean smoke = args.length > 0 && "smoke".equalsIgnoreCase(args[0]);
        ChassisParameters.printLedger();
        run(smoke);
    }
}
