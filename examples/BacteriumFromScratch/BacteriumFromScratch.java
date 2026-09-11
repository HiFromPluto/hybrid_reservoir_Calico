package BacteriumFromScratch;

import java.awt.Color;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import javax.vecmath.Vector3d;

import processing.core.PGraphics3D;
import bsim.BSim;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;

/**
 * Job 2: one non-motile E. coli rod, Valdez elongation, Warren geometry.
 * Isolated. No phage, QS, chemotaxis ODE, or uncited k_growth.
 */
public class BacteriumFromScratch {

    static final double BOUND_X = 100.0;
    static final double BOUND_Y = 100.0;
    static final double BOUND_Z = 10.0;

    static final String EXPORT_DIR = "results/job2_seed101/";
    static final String CSV_NAME = "job2_timeseries.csv";

    public static void main(String[] args) {
        String mode = args.length == 0 ? "all" : args[0];
        boolean preview = "preview".equalsIgnoreCase(mode);
        ChassisParameters.printLedger();

        if (preview) {
            runJob2(true);
            return;
        }
        if ("job2".equalsIgnoreCase(mode) || "all".equalsIgnoreCase(mode)) {
            runJob2(false);
        }
        if ("job3".equalsIgnoreCase(mode) || "all".equalsIgnoreCase(mode)) {
            Job3Sims.runAll();
        }
        if ("job3b".equalsIgnoreCase(mode) || "all".equalsIgnoreCase(mode)) {
            Job3bSims.runAll();
        }
        if ("dxstudy".equalsIgnoreCase(mode)) {
            // ENGINEERING refinement study; physical agar shell held at 1.0 um.
            Job3bSims.runDxStudy(0.5, 2, "050");
            Job3bSims.runDxStudy(0.25, 4, "025");
        }
        if ("fluxstudy".equalsIgnoreCase(mode)) {
            // Is flux_rel_err=0.05 a discretization FLOOR or an iteration CEILING?
            Job3bSims.runDxStudy(0.25, 4, "025_tight", 1.0e-9, 0.005);
        }
    }

    static void runJob2(boolean preview) {
        BSim sim = new BSim();
        sim.setDt(ChassisParameters.DT_S);
        sim.setSimulationTime(ChassisParameters.SIM_TIME_S);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);

        Vector3d centre = new Vector3d(BOUND_X / 2.0, BOUND_Y / 2.0, BOUND_Z / 2.0);
        Vector3d x1 = new Vector3d(centre);
        x1.x -= ChassisParameters.L0_UM / 2.0;
        Vector3d x2 = new Vector3d(centre);
        x2.x += ChassisParameters.L0_UM / 2.0;

        final List<EcoliRodCell> cells = new ArrayList<EcoliRodCell>();
        final EcoliRodCell founder = new EcoliRodCell(sim, x1, x2);
        cells.add(founder);

        final double bathMm = ChassisParameters.NUTRIENT_BATH_MM;

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                for (EcoliRodCell cell : cells) {
                    cell.elongateCited(sim.getDt(), bathMm);
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

        sim.setDrawer(new BSimP3DDrawer(sim, 800, 600) {
            @Override
            public void scene(PGraphics3D p3d) {
                p3d.ortho(0, (float) BOUND_X, (float) BOUND_Y, 0, -1000, 10000);
                p3d.camera((float) BOUND_X / 2f, (float) BOUND_Y / 2f, (float) BOUND_Y,
                        (float) BOUND_X / 2f, (float) BOUND_Y / 2f, 0f,
                        0f, 1f, 0f);
                p3d.perspective((float) Math.PI / 2f,
                        (float) BOUND_X / (float) BOUND_Y, 0.1f, 10000f);
                for (EcoliRodCell cell : cells) {
                    draw(cell, Color.ORANGE);
                }
            }
        });

        if (preview) {
            System.out.println("Job 2 preview");
            sim.preview();
            return;
        }

        BSimUtils.generateDirectoryPath(EXPORT_DIR);
        BSimLogger logger = new BSimLogger(sim, EXPORT_DIR + CSV_NAME) {
            @Override
            public void before() {
                super.before();
                write("t_s;N;founder_L_um;founder_L_analytic_um;L_min_um;L_max_um;radius_um;nutrient_mM;n_negative");
            }

            @Override
            public void during() {
                double t = sim.getTime();
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
                        "%.1f;%d;%.12e;%.12e;%.12e;%.12e;%.12e;%.12e;%d",
                        t,
                        cells.size(),
                        founder.L,
                        founder.analyticLengthUm(bathMm),
                        lMin,
                        lMax,
                        founder.radius,
                        bathMm,
                        nNeg));
            }
        };
        logger.setDt(ChassisParameters.LOG_DT_S);
        sim.addExporter(logger);

        System.out.println("Job 2 export " + EXPORT_DIR + CSV_NAME);
        PrintStream originalOut = System.out;
        try {
            System.setOut(new PrintStream(OutputStream.nullOutputStream()));
            sim.export();
        } finally {
            System.setOut(originalOut);
        }
        System.out.println("Job 2 finished. Completeness is last CSV row t="
                + (int) ChassisParameters.SIM_TIME_S);
    }
}
