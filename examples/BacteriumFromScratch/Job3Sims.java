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
 * Job 3 runners: two-body Hertzian, growing cluster in a named closed box,
 * isolated packing-force regression. Growth still Valdez (4)–(5).
 */
public final class Job3Sims {

    static final String TWO_BODY_DIR = "results/job3_seed101/";
    static final String TWO_BODY_CSV = "twobody_timeseries.csv";
    static final String CLUSTER_CSV = "cluster_timeseries.csv";
    static final String ISOLATED_CSV = "isolated_packing_timeseries.csv";

    private Job3Sims() {}

    static void runAll() {
        runTwoBody();
        runCluster();
        runIsolatedPacking();
    }

    static void runTwoBody() {
        BSim sim = newSim(
                ChassisParameters.TWO_BODY_BOUND_X,
                ChassisParameters.TWO_BODY_BOUND_Y,
                ChassisParameters.TWO_BODY_BOUND_Z,
                ChassisParameters.TWO_BODY_TIME_S);

        double midX = ChassisParameters.TWO_BODY_BOUND_X / 2.0;
        double midY = ChassisParameters.TWO_BODY_BOUND_Y / 2.0;
        double midZ = ChassisParameters.TWO_BODY_BOUND_Z / 2.0;
        double halfL = ChassisParameters.L0_UM / 2.0;
        double halfSep = (ChassisParameters.W0_UM - ChassisParameters.TWO_BODY_DELTA0_UM) / 2.0;

        EcoliRodCell a = new EcoliRodCell(sim,
                new Vector3d(midX - halfL, midY - halfSep, midZ),
                new Vector3d(midX + halfL, midY - halfSep, midZ));
        EcoliRodCell b = new EcoliRodCell(sim,
                new Vector3d(midX - halfL, midY + halfSep, midZ),
                new Vector3d(midX + halfL, midY + halfSep, midZ));

        final List<EcoliRodCell> cells = new ArrayList<EcoliRodCell>();
        cells.add(a);
        cells.add(b);

        final Vector3d cA0 = new Vector3d(a.centre());
        final Vector3d cB0 = new Vector3d(b.centre());
        final Vector3d bound = sim.getBound();

        double d0 = ValdezHertzian.pairDelta(a, b);
        System.out.println("Job 3 two-body seed delta_cc = " + d0
                + " um (target " + ChassisParameters.TWO_BODY_DELTA0_UM + ")");

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                ValdezHertzian.relaxContacts(cells, bound);
            }
        });

        BSimUtils.generateDirectoryPath(TWO_BODY_DIR);
        BSimLogger logger = new BSimLogger(sim, TWO_BODY_DIR + TWO_BODY_CSV) {
            @Override
            public void before() {
                super.before();
                write("t_s;seed;delta_cc_um;gap_um;disp_a_um;disp_b_um;F_cc_mag;y_a;y_b");
                ValdezHertzian.PackingStats st = ValdezHertzian.stats(cells, bound);
                write(String.format(Locale.US,
                        "%.1f;%d;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e",
                        0.0,
                        ChassisParameters.RNG_SEED,
                        ValdezHertzian.pairDelta(a, b),
                        ValdezHertzian.surfaceGap(a, b),
                        0.0,
                        0.0,
                        st.maxForceMag,
                        a.centre().y,
                        b.centre().y));
            }

            @Override
            public void during() {
                ValdezHertzian.PackingStats st = ValdezHertzian.stats(cells, bound);
                double gap = ValdezHertzian.surfaceGap(a, b);
                write(String.format(Locale.US,
                        "%.1f;%d;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e",
                        sim.getTime(),
                        ChassisParameters.RNG_SEED,
                        ValdezHertzian.pairDelta(a, b),
                        gap,
                        ValdezHertzian.dist(a.centre(), cA0),
                        ValdezHertzian.dist(b.centre(), cB0),
                        st.maxForceMag,
                        a.centre().y,
                        b.centre().y));
            }
        };
        logger.setDt(ChassisParameters.LOG_DT_S);
        sim.addExporter(logger);
        System.out.println("Job 3 two-body export " + TWO_BODY_DIR + TWO_BODY_CSV);
        silentExport(sim);
        System.out.println("Job 3 two-body finished t="
                + (int) ChassisParameters.TWO_BODY_TIME_S);
    }

    static void runCluster() {
        BSim sim = newSim(
                ChassisParameters.CLUSTER_BOUND_X,
                ChassisParameters.CLUSTER_BOUND_Y,
                ChassisParameters.CLUSTER_BOUND_Z,
                ChassisParameters.CLUSTER_TIME_S);

        double midX = ChassisParameters.CLUSTER_BOUND_X / 2.0;
        double midY = ChassisParameters.CLUSTER_BOUND_Y / 2.0;
        double midZ = ChassisParameters.CLUSTER_BOUND_Z / 2.0;
        double halfL = ChassisParameters.L0_UM / 2.0;
        double halfSep = (ChassisParameters.W0_UM - ChassisParameters.CLUSTER_SEED_DELTA_UM) / 2.0;

        final List<EcoliRodCell> cells = new ArrayList<EcoliRodCell>();
        cells.add(new EcoliRodCell(sim,
                new Vector3d(midX - halfL, midY - halfSep, midZ),
                new Vector3d(midX + halfL, midY - halfSep, midZ)));
        cells.add(new EcoliRodCell(sim,
                new Vector3d(midX - halfL, midY + halfSep, midZ),
                new Vector3d(midX + halfL, midY + halfSep, midZ)));

        final Vector3d bound = sim.getBound();
        final double bathMm = ChassisParameters.NUTRIENT_BATH_MM;

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
            }
        });

        BSimUtils.generateDirectoryPath(TWO_BODY_DIR);
        BSimLogger logger = new BSimLogger(sim, TWO_BODY_DIR + CLUSTER_CSV) {
            @Override
            public void before() {
                super.before();
                write("t_s;seed;N;delta_cc_max_um;d_centers_min_um;L_min_um;L_max_um;radius_um;nutrient_mM;n_negative;F_pack_max");
                ValdezHertzian.PackingStats st = ValdezHertzian.stats(cells, bound);
                double lMin = Double.POSITIVE_INFINITY;
                double lMax = Double.NEGATIVE_INFINITY;
                int nNeg = 0;
                for (EcoliRodCell cell : cells) {
                    lMin = Math.min(lMin, cell.L);
                    lMax = Math.max(lMax, cell.L);
                    if (cell.hasNegativeState(bathMm)) {
                        nNeg++;
                    }
                }
                write(String.format(Locale.US,
                        "%.1f;%d;%d;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%d;%.12e",
                        0.0,
                        ChassisParameters.RNG_SEED,
                        cells.size(),
                        st.maxDeltaCc,
                        st.minCentreDist,
                        lMin,
                        lMax,
                        ChassisParameters.RADIUS_UM,
                        bathMm,
                        nNeg,
                        st.maxForceMag));
            }

            @Override
            public void during() {
                ValdezHertzian.PackingStats st = ValdezHertzian.stats(cells, bound);
                double lMin = Double.POSITIVE_INFINITY;
                double lMax = Double.NEGATIVE_INFINITY;
                int nNeg = 0;
                for (EcoliRodCell cell : cells) {
                    lMin = Math.min(lMin, cell.L);
                    lMax = Math.max(lMax, cell.L);
                    if (cell.hasNegativeState(bathMm)) {
                        nNeg++;
                    }
                }
                write(String.format(Locale.US,
                        "%.1f;%d;%d;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%d;%.12e",
                        sim.getTime(),
                        ChassisParameters.RNG_SEED,
                        cells.size(),
                        st.maxDeltaCc,
                        st.minCentreDist,
                        lMin,
                        lMax,
                        ChassisParameters.RADIUS_UM,
                        bathMm,
                        nNeg,
                        st.maxForceMag));
            }
        };
        logger.setDt(ChassisParameters.LOG_DT_S);
        sim.addExporter(logger);
        System.out.println("Job 3 cluster export " + TWO_BODY_DIR + CLUSTER_CSV
                + "  CLUSTER_BOX="
                + ChassisParameters.CLUSTER_BOUND_X + "x"
                + ChassisParameters.CLUSTER_BOUND_Y + "x"
                + ChassisParameters.CLUSTER_BOUND_Z + " um closed");
        silentExport(sim);
        System.out.println("Job 3 cluster finished t="
                + (int) ChassisParameters.CLUSTER_TIME_S
                + " N=" + cells.size());
    }

    static void runIsolatedPacking() {
        BSim sim = newSim(100.0, 100.0, 10.0, ChassisParameters.ISOLATED_PACK_TIME_S);
        Vector3d centre = new Vector3d(50.0, 50.0, 5.0);
        Vector3d x1 = new Vector3d(centre);
        x1.x -= ChassisParameters.L0_UM / 2.0;
        Vector3d x2 = new Vector3d(centre);
        x2.x += ChassisParameters.L0_UM / 2.0;

        final List<EcoliRodCell> cells = new ArrayList<EcoliRodCell>();
        final EcoliRodCell founder = new EcoliRodCell(sim, x1, x2);
        cells.add(founder);
        final Vector3d bound = sim.getBound();
        final double bathMm = ChassisParameters.NUTRIENT_BATH_MM;

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                for (EcoliRodCell cell : cells) {
                    cell.elongateCited(sim.getDt(), bathMm);
                }
                ValdezHertzian.relaxContacts(cells, bound);
            }
        });

        BSimUtils.generateDirectoryPath(TWO_BODY_DIR);
        BSimLogger logger = new BSimLogger(sim, TWO_BODY_DIR + ISOLATED_CSV) {
            @Override
            public void before() {
                super.before();
                write("t_s;seed;N;founder_L_um;founder_L_analytic_um;F_pack_max;radius_um;nutrient_mM;n_negative");
            }

            @Override
            public void during() {
                ValdezHertzian.PackingStats st = ValdezHertzian.stats(cells, bound);
                int nNeg = founder.hasNegativeState(bathMm) ? 1 : 0;
                write(String.format(Locale.US,
                        "%.1f;%d;%d;%.12e;%.12e;%.12e;%.12e;%.12e;%d",
                        sim.getTime(),
                        ChassisParameters.RNG_SEED,
                        cells.size(),
                        founder.L,
                        founder.analyticLengthUm(bathMm),
                        st.maxForceMag,
                        founder.radius,
                        bathMm,
                        nNeg));
            }
        };
        logger.setDt(ChassisParameters.LOG_DT_S);
        sim.addExporter(logger);
        System.out.println("Job 3 isolated packing export " + TWO_BODY_DIR + ISOLATED_CSV);
        silentExport(sim);
        System.out.println("Job 3 isolated packing finished t="
                + (int) ChassisParameters.ISOLATED_PACK_TIME_S);
    }

    private static BSim newSim(double bx, double by, double bz, double tEnd) {
        BSim sim = new BSim();
        sim.setDt(ChassisParameters.DT_S);
        sim.setSimulationTime(tEnd);
        sim.setTimeFormat("0.00");
        sim.setBound(bx, by, bz);
        sim.setSolid(true, true, true);
        return sim;
    }

    private static void silentExport(BSim sim) {
        PrintStream originalOut = System.out;
        try {
            System.setOut(new PrintStream(OutputStream.nullOutputStream()));
            sim.export();
        } finally {
            System.setOut(originalOut);
        }
    }
}
