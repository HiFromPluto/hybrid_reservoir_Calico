package BacteriumFromScratch;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import javax.vecmath.Vector3d;

import bsim.BSim;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.export.BSimLogger;

/**
 * Job 3c — colony dish, radial nutrient limitation.
 *
 * <p>Every parameter here is frozen in PROTOCOL.md. Chamber 60x60x1 um
 * (b_z = one cell diameter, so the monolayer is enforced by the wall Hertzian),
 * DISH_LATERAL rim at C_s, CONSERVATIVE_MASS occupancy, SYMMETRY_BROKEN
 * division, dx = 1.0 um, pad = 1, omega = 1.96, conv_tol = 1e-10,
 * flux_tol = 1e-3. Seed 101. One founder at chamber centre.
 *
 * <p>COLLINEAR division and b_z = 4 um were the first attempt; they produced a
 * 1-D filament and a 3-layer stack respectively. See
 * results/job3c_archive_collinear/.
 *
 * <p>OFF: quorum sensing, phage, chemotaxis, metabolism ODE, Brownian, F_s,
 * stochastic P(l), cell death, second species, any reservoir readout.
 */
public final class Job3cSims {

    private Job3cSims() {}

    static final String EXPORT_DIR = "results/job3c_seed101/";
    static final String CSV_NAME = "colony_timeseries.csv";

    // ---- frozen chamber + solver stack (PROTOCOL.md Job 3c) ----
    static final double BX = 60.0, BY = 60.0, BZ = 1.0;   // b_z = 2r: monolayer enforced
    static final double DX = 1.0;
    static final int PAD = 1;
    static final double OMEGA = 1.96;
    static final double CONV_TOL = 1.0e-10;
    static final double FLUX_TOL = 1.0e-3;
    static final double T_END = 23400.0;

    /** Per-bin elongation statistics for the Gate 4 radial gradient. */
    private static final class BinStats {
        int n;
        double mean;
        double sem;
    }

    /**
     * Instantaneous fractional elongation rate of a cell, 1/s.
     * Valdez (4)-(5): dL/dt = sigma * monod(N_local) * L, so the rate per unit
     * length is sigma * monod -- independent of the cell's current length and
     * therefore comparable across cells at any point in their cycle. This is
     * the "matched time-since-birth" comparison the protocol asks for, done
     * exactly rather than by binning ages.
     */
    private static double fractionalRate(double nLocal) {
        return ChassisParameters.sigmaPerSecond() * ChassisParameters.monodFactor(nLocal);
    }

    private static BinStats stats(List<Double> vals) {
        BinStats b = new BinStats();
        b.n = vals.size();
        if (b.n == 0) {
            b.mean = Double.NaN;
            b.sem = Double.NaN;
            return b;
        }
        double sum = 0.0;
        for (double v : vals) {
            sum += v;
        }
        b.mean = sum / b.n;
        if (b.n < 2) {
            b.sem = 0.0;
            return b;
        }
        double ss = 0.0;
        for (double v : vals) {
            ss += (v - b.mean) * (v - b.mean);
        }
        b.sem = Math.sqrt(ss / (b.n - 1)) / Math.sqrt(b.n);
        return b;
    }

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
     * Fraction of cell centres more than half a cell radius off the mid-plane.
     *
     * <p>At b_z = 2r this is enforced to zero by the wall Hertzian, so it is a
     * REGRESSION GUARD rather than a discriminating test: it fires only if the
     * chamber height is changed. The load-bearing parts of the morphology gate
     * are R/R_disc and d_centres_min, which the filament run failed badly
     * (4.94 and 1.98 um respectively).
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

    /** Monolayer disc radius for N cells at packing fraction 0.85. */
    private static double discRadius(int n) {
        double footprint = ChassisParameters.W0_UM
                * (0.5 * (ChassisParameters.L0_UM + ChassisParameters.LDIV_UM)
                   + 2.0 * ChassisParameters.RADIUS_UM);
        return Math.sqrt(n * footprint / 0.85 / Math.PI);
    }

    /** Frozen definition: max pole distance from the colony centre of mass. */
    private static double colonyRadius(List<EcoliRodCell> cells, Vector3d com) {
        double r = 0.0;
        for (EcoliRodCell cell : cells) {
            r = Math.max(r, Math.hypot(cell.x1.x - com.x, cell.x1.y - com.y));
            r = Math.max(r, Math.hypot(cell.x2.x - com.x, cell.x2.y - com.y));
        }
        return r;
    }

    static void run() {
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
        final NutrientField field = new NutrientField(BX, BY, BZ, false, true, DX, PAD,
                NutrientField.BcMode.DISH_LATERAL);
        field.setOccupancyMode(NutrientField.OccupancyMode.CONSERVATIVE_MASS);
        field.setOverRelaxation(OMEGA);
        field.setAcceptanceOverride(CONV_TOL, FLUX_TOL);

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                field.markOccupancy(cells);
                field.solveQuasiSteady();
                for (EcoliRodCell cell : cells) {
                    cell.elongateCited(sim.getDt(), field.sampleMm(cell.centre()));
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
            }
        });

        BSimUtils.generateDirectoryPath(EXPORT_DIR);
        BSimLogger logger = new BSimLogger(sim, EXPORT_DIR + CSV_NAME) {
            @Override
            public void before() {
                super.before();
                write("t_s;seed;N;R_um;R_disc_um;R_Rdisc;offplane;N_centre_mM;N_edge_mM;"
                        + "field_min_mM;field_mean_mM;"
                        + "delta_cc_max_um;d_centers_min_um;flux_rel_err;gs_iter;gs_converged;"
                        + "n_clamp;marked_um3;inner_n;inner_rate;inner_sem;outer_n;outer_rate;outer_sem");
                logRow(0.0);
            }

            @Override
            public void during() {
                logRow(sim.getTime());
            }

            private void logRow(double t) {
                ValdezHertzian.PackingStats pst = ValdezHertzian.stats(cells, bound);
                NutrientField.FieldStats fst = field.stats();
                NutrientField.FluxBalance fb = field.fluxBalance();
                Vector3d com = centreOfMass(cells);
                double R = colonyRadius(cells, com);

                List<Double> inner = new ArrayList<Double>();
                List<Double> outer = new ArrayList<Double>();
                List<Double> edgeN = new ArrayList<Double>();
                for (EcoliRodCell cell : cells) {
                    Vector3d p = cell.centre();
                    double d = Math.hypot(p.x - com.x, p.y - com.y);
                    double nLocal = field.sampleMm(p);
                    double rate = fractionalRate(nLocal);
                    if (R > 0.0 && d < R / 3.0) {
                        inner.add(rate);
                    } else if (R > 0.0 && d > 2.0 * R / 3.0) {
                        outer.add(rate);
                    }
                    if (R > 0.0 && d > 0.75 * R) {
                        edgeN.add(nLocal);
                    }
                }
                BinStats bi = stats(inner);
                BinStats bo = stats(outer);
                double nCentre = field.sampleMm(new Vector3d(com.x, com.y, BZ / 2.0));
                double nEdge = Double.NaN;
                if (!edgeN.isEmpty()) {
                    double s = 0.0;
                    for (double v : edgeN) {
                        s += v;
                    }
                    nEdge = s / edgeN.size();
                }

                double rDisc = discRadius(cells.size());
                write(String.format(Locale.US,
                        "%.1f;%d;%d;%.6f;%.6f;%.6f;%.6f;%.9f;%.9f;%.9f;%.9f;%.9e;%.9f;%.6e;%d;%d;%d;%.4f;"
                        + "%d;%.9e;%.9e;%d;%.9e;%.9e",
                        t, ChassisParameters.RNG_SEED, cells.size(), R,
                        rDisc, rDisc > 0.0 ? R / rDisc : 0.0, outOfPlane(cells),
                        nCentre, nEdge, fst.minMm, fst.meanMm,
                        pst.maxDeltaCc, pst.minCentreDist,
                        fb.relativeError, field.lastSolveIter(),
                        field.lastSolveConverged() ? 1 : 0, field.countClampEvents(),
                        field.markedMassUm3(),
                        bi.n, bi.mean, bi.sem, bo.n, bo.mean, bo.sem));
            }
        };
        logger.setDt(ChassisParameters.LOG_DT_S);
        sim.addExporter(logger);
        System.out.println("Job 3c colony dish -> " + EXPORT_DIR + CSV_NAME);
        System.out.println("  chamber " + BX + "x" + BY + "x" + BZ
                + "  dx=" + DX + " pad=" + PAD + "  DISH_LATERAL + CONSERVATIVE_MASS"
                + " + SYMMETRY_BROKEN"
                + "  omega=" + OMEGA + "  conv_tol=" + CONV_TOL + "  flux_tol=" + FLUX_TOL);

        PrintStream originalOut = System.out;
        try {
            System.setOut(new PrintStream(OutputStream.nullOutputStream()));
            sim.export();
        } finally {
            System.setOut(originalOut);
        }
        System.out.println("Job 3c finished  t=" + (int) T_END + "  N=" + cells.size());
    }

    public static void main(String[] args) {
        ChassisParameters.printLedger();
        run();
    }
}
