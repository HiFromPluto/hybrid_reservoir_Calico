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
 * Job 3b: Valdez/Warren spatial nutrient field + local Monod elongation.
 * Job 3 Hertzian packing unchanged. Growth samples N(x) at cell centre.
 */
public final class Job3bSims {

    static final String EXPORT_DIR = "results/job3b_seed101/";
    static final String UNIFORM_CSV = "uniform_timeseries.csv";
    static final String PACKING_CSV = "packing_cluster_timeseries.csv";
    static final String DEPLETION_CSV = "depletion_timeseries.csv";

    private Job3bSims() {}

    static void runAll() {
        runUniformRegression();
        runPackingRegression();
        runDepletionCoupled();
    }

    /** Large box, uniform C_s, ρ=0 → Job 2 analytic unchanged. */
    static void runUniformRegression() {
        BSim sim = newSim(
                ChassisParameters.UNIFORM_REG_BOUND_X,
                ChassisParameters.UNIFORM_REG_BOUND_Y,
                ChassisParameters.UNIFORM_REG_BOUND_Z,
                ChassisParameters.SIM_TIME_S);

        Vector3d centre = new Vector3d(
                ChassisParameters.UNIFORM_REG_BOUND_X / 2.0,
                ChassisParameters.UNIFORM_REG_BOUND_Y / 2.0,
                ChassisParameters.UNIFORM_REG_BOUND_Z / 2.0);
        Vector3d x1 = new Vector3d(centre);
        x1.x -= ChassisParameters.L0_UM / 2.0;
        Vector3d x2 = new Vector3d(centre);
        x2.x += ChassisParameters.L0_UM / 2.0;

        final List<EcoliRodCell> cells = new ArrayList<EcoliRodCell>();
        final EcoliRodCell founder = new EcoliRodCell(sim, x1, x2);
        cells.add(founder);

        final NutrientField field = new NutrientField(
                sim.getBound().x, sim.getBound().y, sim.getBound().z,
                true, false);
        final double bathMm = ChassisParameters.NUTRIENT_BATH_MM;

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                field.markOccupancy(cells);
                field.solveQuasiSteady();
                for (EcoliRodCell cell : cells) {
                    double nLocal = field.sampleMm(cell.centre());
                    cell.elongateCited(sim.getDt(), nLocal);
                }
                List<EcoliRodCell> newborns = new ArrayList<EcoliRodCell>();
                for (EcoliRodCell cell : cells) {
                    if (cell.shouldDivide()) {
                        newborns.add(cell.divideCited());
                    }
                }
                cells.addAll(newborns);
            }
        });

        exportUniform(sim, cells, founder, field, bathMm);
        System.out.println("Job 3b uniform regression finished t="
                + (int) ChassisParameters.SIM_TIME_S);
    }

    /** CLUSTER_BOX with uniform field (ρ=0) + Hertzian — mirrors Job 3 cluster. */
    static void runPackingRegression() {
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
        final NutrientField field = new NutrientField(
                bound.x, bound.y, bound.z, true, false);
        final double bathMm = ChassisParameters.NUTRIENT_BATH_MM;

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

        exportPackingCluster(sim, cells, field, bound, bathMm);
        System.out.println("Job 3b packing regression finished t="
                + (int) ChassisParameters.CLUSTER_TIME_S + " N=" + cells.size());
    }

    /** Closed CLUSTER_BOX, Monod consumption, local elongation + Hertzian. */
    static void runDepletionCoupled() {
        BSim sim = newSim(
                ChassisParameters.CLUSTER_BOUND_X,
                ChassisParameters.CLUSTER_BOUND_Y,
                ChassisParameters.CLUSTER_BOUND_Z,
                ChassisParameters.DEPLETION_TIME_S);

        double midX = ChassisParameters.CLUSTER_BOUND_X / 2.0;
        double midY = ChassisParameters.CLUSTER_BOUND_Y / 2.0;
        double midZ = ChassisParameters.CLUSTER_BOUND_Z / 2.0;
        double halfL = ChassisParameters.L0_UM / 2.0;
        double halfSep = (ChassisParameters.W0_UM - ChassisParameters.CLUSTER_SEED_DELTA_UM) / 2.0;

        final List<EcoliRodCell> cells = new ArrayList<EcoliRodCell>();
        final EcoliRodCell founder = new EcoliRodCell(sim,
                new Vector3d(midX - halfL, midY - halfSep, midZ),
                new Vector3d(midX + halfL, midY - halfSep, midZ));
        cells.add(founder);
        cells.add(new EcoliRodCell(sim,
                new Vector3d(midX - halfL, midY + halfSep, midZ),
                new Vector3d(midX + halfL, midY + halfSep, midZ)));

        final Vector3d bound = sim.getBound();
        final NutrientField field = new NutrientField(
                bound.x, bound.y, bound.z, false, true);
        final double bathMm = ChassisParameters.NUTRIENT_BATH_MM;

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

        exportDepletion(sim, cells, founder, field, bound, bathMm);
        System.out.println("Job 3b depletion/coupled finished t="
                + (int) ChassisParameters.DEPLETION_TIME_S + " N=" + cells.size());
    }

    private static void exportUniform(BSim sim, final List<EcoliRodCell> cells,
            final EcoliRodCell founder, final NutrientField field,
            final double bathMm) {
        BSimUtils.generateDirectoryPath(EXPORT_DIR);
        BSimLogger logger = new BSimLogger(sim, EXPORT_DIR + UNIFORM_CSV) {
            @Override
            public void before() {
                super.before();
                write("t_s;seed;N;founder_L_um;founder_L_analytic_um;N_local_mM;field_min_mM;field_mean_mM;monod_local;radius_um;n_negative");
            }

            @Override
            public void during() {
                NutrientField.FieldStats st = field.stats();
                double nLocal = field.sampleMm(founder.centre());
                double monod = ChassisParameters.monodFactor(nLocal);
                int nNeg = countNegative(cells, nLocal, field);
                write(String.format(Locale.US,
                        "%.1f;%d;%d;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%d",
                        sim.getTime(),
                        ChassisParameters.RNG_SEED,
                        cells.size(),
                        founder.L,
                        founder.analyticLengthUm(bathMm),
                        nLocal,
                        st.minMm,
                        st.meanMm,
                        monod,
                        founder.radius,
                        nNeg));
            }
        };
        logger.setDt(ChassisParameters.LOG_DT_S);
        sim.addExporter(logger);
        System.out.println("Job 3b uniform export " + EXPORT_DIR + UNIFORM_CSV);
        silentExport(sim);
    }

    private static void exportPackingCluster(BSim sim, final List<EcoliRodCell> cells,
            final NutrientField field, final Vector3d bound, final double bathMm) {
        BSimUtils.generateDirectoryPath(EXPORT_DIR);
        BSimLogger logger = new BSimLogger(sim, EXPORT_DIR + PACKING_CSV) {
            @Override
            public void before() {
                super.before();
                write("t_s;seed;N;delta_cc_max_um;d_centers_min_um;L_min_um;L_max_um;radius_um;N_local_mM;field_min_mM;field_mean_mM;monod_local;n_negative;F_pack_max");
                logRow(0.0);
            }

            @Override
            public void during() {
                logRow(sim.getTime());
            }

            private void logRow(double t) {
                ValdezHertzian.PackingStats pst = ValdezHertzian.stats(cells, bound);
                NutrientField.FieldStats fst = field.stats();
                double lMin = Double.POSITIVE_INFINITY;
                double lMax = Double.NEGATIVE_INFINITY;
                double nLocalSum = 0.0;
                int nNeg = 0;
                for (EcoliRodCell cell : cells) {
                    lMin = Math.min(lMin, cell.L);
                    lMax = Math.max(lMax, cell.L);
                    double nl = field.sampleMm(cell.centre());
                    nLocalSum += nl;
                    if (cell.hasNegativeState(nl) || !Double.isFinite(nl)) {
                        nNeg++;
                    }
                }
                double nLocalMean = cells.isEmpty() ? bathMm : nLocalSum / cells.size();
                write(String.format(Locale.US,
                        "%.1f;%d;%d;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%d;%.12e",
                        t,
                        ChassisParameters.RNG_SEED,
                        cells.size(),
                        pst.maxDeltaCc,
                        pst.minCentreDist,
                        lMin,
                        lMax,
                        ChassisParameters.RADIUS_UM,
                        nLocalMean,
                        fst.minMm,
                        fst.meanMm,
                        ChassisParameters.monodFactor(nLocalMean),
                        nNeg,
                        pst.maxForceMag));
            }
        };
        logger.setDt(ChassisParameters.LOG_DT_S);
        sim.addExporter(logger);
        System.out.println("Job 3b packing export " + EXPORT_DIR + PACKING_CSV);
        silentExport(sim);
    }

    private static void exportDepletion(BSim sim, final List<EcoliRodCell> cells,
            final EcoliRodCell founder, final NutrientField field,
            final Vector3d bound, final double bathMm) {
        BSimUtils.generateDirectoryPath(EXPORT_DIR);
        BSimLogger logger = new BSimLogger(sim, EXPORT_DIR + DEPLETION_CSV) {
            @Override
            public void before() {
                super.before();
                write("t_s;seed;N;founder_L_um;founder_L_analytic_uniform_mM;N_local_mM;field_min_mM;field_mean_mM;monod_local;delta_cc_max_um;d_centers_min_um;flux_rel_err;gs_iter;gs_converged;n_clamp;F_pack_max");
                logRow(0.0);
            }

            @Override
            public void during() {
                logRow(sim.getTime());
            }

            private void logRow(double t) {
                ValdezHertzian.PackingStats pst = ValdezHertzian.stats(cells, bound);
                NutrientField.FieldStats fst = field.stats();
                NutrientField.FluxBalance flux = field.fluxBalance();
                double nLocal = field.sampleMm(founder.centre());
                double monod = ChassisParameters.monodFactor(nLocal);
                int nClamp = field.countClampEvents();
                write(String.format(Locale.US,
                        "%.1f;%d;%d;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%d;%d;%d;%.12e",
                        t,
                        ChassisParameters.RNG_SEED,
                        cells.size(),
                        founder.L,
                        founder.analyticLengthUm(bathMm),
                        nLocal,
                        fst.minMm,
                        fst.meanMm,
                        monod,
                        pst.maxDeltaCc,
                        pst.minCentreDist,
                        flux.relativeError,
                        field.lastSolveIter(),
                        field.lastSolveConverged() ? 1 : 0,
                        nClamp,
                        pst.maxForceMag));
            }
        };
        logger.setDt(ChassisParameters.LOG_DT_S);
        sim.addExporter(logger);
        System.out.println("Job 3b depletion export " + EXPORT_DIR + DEPLETION_CSV
                + "  CLUSTER_BOX closed");
        silentExport(sim);
    }

    private static int countNegative(List<EcoliRodCell> cells, double nLocal,
            NutrientField field) {
        int nNeg = field.countNonFinite();
        for (EcoliRodCell cell : cells) {
            double nl = field.sampleMm(cell.centre());
            if (cell.hasNegativeState(nl)) {
                nNeg++;
            }
        }
        return nNeg;
    }

    /**
     * ENGINEERING grid-refinement study. NOT a Job 3b gate run and NOT gated by
     * check_job3b.py. Same depletion configuration; only dx and the agar pad
     * change, with the PHYSICAL agar thickness (pad*dx) held fixed by the caller.
     */
    static void runDxStudy(double dxUm, int agarPadVoxels, String tag) {
        runDxStudy(dxUm, agarPadVoxels, tag, -1.0, -1.0);
    }

    static void runDxStudy(double dxUm, int agarPadVoxels, String tag,
                           double convTolMm, double fluxTolRel) {
        BSim sim = newSim(
                ChassisParameters.CLUSTER_BOUND_X,
                ChassisParameters.CLUSTER_BOUND_Y,
                ChassisParameters.CLUSTER_BOUND_Z,
                ChassisParameters.DEPLETION_TIME_S);

        double midX = ChassisParameters.CLUSTER_BOUND_X / 2.0;
        double midY = ChassisParameters.CLUSTER_BOUND_Y / 2.0;
        double midZ = ChassisParameters.CLUSTER_BOUND_Z / 2.0;
        double halfL = ChassisParameters.L0_UM / 2.0;
        double halfSep = (ChassisParameters.W0_UM - ChassisParameters.CLUSTER_SEED_DELTA_UM) / 2.0;

        final List<EcoliRodCell> cells = new ArrayList<EcoliRodCell>();
        final EcoliRodCell founder = new EcoliRodCell(sim,
                new Vector3d(midX - halfL, midY - halfSep, midZ),
                new Vector3d(midX + halfL, midY - halfSep, midZ));
        cells.add(founder);
        cells.add(new EcoliRodCell(sim,
                new Vector3d(midX - halfL, midY + halfSep, midZ),
                new Vector3d(midX + halfL, midY + halfSep, midZ)));

        final Vector3d bound = sim.getBound();
        final NutrientField field = new NutrientField(
                bound.x, bound.y, bound.z, false, true, dxUm, agarPadVoxels);
        field.setAcceptanceOverride(convTolMm, fluxTolRel);

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

        final String dir = "results/job3b_dxstudy/";
        final String csv = "depletion_dx" + tag + ".csv";
        BSimUtils.generateDirectoryPath(dir);
        BSimLogger logger = new BSimLogger(sim, dir + csv) {
            @Override
            public void before() {
                super.before();
                write("t_s;seed;N;dx_um;agar_um;n_voxels;founder_L_um;N_local_mM;field_min_mM;field_mean_mM;flux_rel_err;gs_iter;gs_converged;n_clamp");
            }

            @Override
            public void during() {
                NutrientField.FieldStats fst = field.stats();
                NutrientField.FluxBalance fb = field.fluxBalance();
                write(String.format(Locale.US,
                        "%.1f;%d;%d;%.4f;%.4f;%d;%.12e;%.12e;%.12e;%.12e;%.12e;%d;%d;%d",
                        sim.getTime(),
                        ChassisParameters.RNG_SEED,
                        cells.size(),
                        field.dxUm(),
                        field.agarThicknessUm(),
                        field.voxelCount(),
                        founder.L,
                        field.sampleMm(founder.centre()),
                        fst.minMm,
                        fst.meanMm,
                        fb.relativeError,
                        field.lastSolveIter(),
                        field.lastSolveConverged() ? 1 : 0,
                        field.countClampEvents()));
            }
        };
        logger.setDt(ChassisParameters.LOG_DT_S);
        sim.addExporter(logger);
        System.out.println("dx study: dx=" + dxUm + " um  pad=" + agarPadVoxels
                + " voxels (" + (dxUm * agarPadVoxels) + " um agar)  -> " + dir + csv);
        silentExport(sim);
        System.out.println("dx study dx=" + dxUm + " finished N=" + cells.size());
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
