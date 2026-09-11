package BSimMonodM3;

import java.awt.Color;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Locale;
import java.util.Random;
import java.util.Vector;

import javax.vecmath.Vector3d;

import processing.core.PGraphics3D;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.particle.BSimBacterium;

/**
 * Isolated Monod M3 example: M2 growth+uptake plus one constant point
 * source. No replenisher. No vesicle AC. No glucose chemotaxis. Not a
 * reservoir window experiment.
 *
 * Growth / uptake (frozen, unchanged from M2):
 *   µ = µ_max * G / (K_s + G)
 *   K_s       = 0.18 µM
 *   µ_max     = Math.log(2.0) / 1800.0
 *   K_UPTAKE  = 1.2e3 molecules/s
 *
 * Source (frozen engineering identity rate, not from seeing contrast):
 *   N0 * K_UPTAKE = 50 * 1.2e3 = 6e4 molecules/s at saturating G
 *   K_SOURCE = 1.0e5 molecules/s   (~1.7× that)
 *
 * M1 and M2 remain closed PASSes. Do not copy glucose into
 * BSimReservoirPlanStage6.
 */
public class BSimMonodM3 {

    static final double BOUND_X = 200.0;
    static final double BOUND_Y = 200.0;
    static final double BOUND_Z = 10.0;

    static final int GRID_X = 20;
    static final int GRID_Y = 20;
    static final int GRID_Z = 1;

    static final double GLU_DIFFUSIVITY = 100.0;
    static final double GLU_DECAY_RATE = 0.0;

    static final double DT = 0.05;
    static final double SIM_TIME = 7200.0;
    static final int N_CELLS = 50;
    static final double FLOW_SPEED = 0.0;
    static final long RNG_SEED = 101L;

    /** 1 µM = 602 molecules/µm³. Same conversion as D1 AHL / M1 / M2. */
    static final double MOL_PER_UM3_PER_UM = 602.0;

    static final double K_S_UM = 0.18;
    static final double MU_MAX = Math.log(2.0) / 1800.0;
    static final double GROWTH_RATE_SA = 4.0 * Math.PI / 1800.0;
    static final double K_UPTAKE = 1.2e3;

    static final double G0_UM = 0.18;
    static final double K_SOURCE_ON = 1.0e5;

    static final double SOURCE_X = 50.0;
    static final double SOURCE_Y = 100.0;
    static final double SOURCE_Z = 5.0;
    static final double NEAR_X = 70.0;
    static final double FAR_X = 150.0;

    static Random experimentRng;
    static int birthsNear;
    static int birthsFar;

    /**
     * M2 Monod surface-area growth + uptake. No setGoal. Births counted
     * by child x into near (x<70) vs far (x>150).
     */
    static class MonodBacterium extends BSimBacterium {
        final BSimChemicalField glucoseField;

        static void reseed() {
            rng.setSeed(RNG_SEED);
        }

        public MonodBacterium(BSim sim, Vector3d position, BSimChemicalField glucoseField) {
            super(sim, position);
            this.glucoseField = glucoseField;
        }

        @Override
        public void setRadius() {
            setRadiusFromSurfaceArea(surfaceArea(replicationRadius) / 2
                    + experimentRng.nextDouble() * surfaceArea(replicationRadius) / 2);
        }

        @Override
        public void grow() {
            double gUm = glucoseField.getConc(position) / MOL_PER_UM3_PER_UM;
            double monod = gUm / (K_S_UM + gUm);
            double saved = surfaceAreaGrowthRate;
            surfaceAreaGrowthRate = GROWTH_RATE_SA * monod;
            super.grow();
            surfaceAreaGrowthRate = saved;

            gUm = glucoseField.getConc(position) / MOL_PER_UM3_PER_UM;
            double uptake = K_UPTAKE * gUm / (K_S_UM + gUm);
            glucoseField.addQuantity(position, -uptake * sim.getDt());
        }

        @SuppressWarnings("unchecked")
        @Override
        public void replicate() {
            setRadiusFromSurfaceArea(surfaceArea(replicationRadius) / 2);
            Vector3d childPosition = new Vector3d(position);
            childPosition.x += 2.0 * radius * (experimentRng.nextDouble() - 0.5);
            childPosition.y += 2.0 * radius * (experimentRng.nextDouble() - 0.5);
            if (childPosition.x < NEAR_X) birthsNear++;
            else if (childPosition.x > FAR_X) birthsFar++;
            MonodBacterium child = new MonodBacterium(sim, childPosition, glucoseField);
            child.setRadius(radius);
            child.setSurfaceAreaGrowthRate(GROWTH_RATE_SA);
            child.setChildList(childList);
            childList.add(child);
        }
    }

    static double voxelCenterX(int i) {
        return (i + 0.5) * (BOUND_X / GRID_X);
    }

    static double meanFieldUm(BSimChemicalField field) {
        double sum = 0.0;
        int n = 0;
        for (int i = 0; i < GRID_X; i++) {
            for (int j = 0; j < GRID_Y; j++) {
                for (int k = 0; k < GRID_Z; k++) {
                    sum += field.getConc(i, j, k);
                    n++;
                }
            }
        }
        return (sum / n) / MOL_PER_UM3_PER_UM;
    }

    static double meanRegionUm(BSimChemicalField field, boolean near) {
        double sum = 0.0;
        int n = 0;
        for (int i = 0; i < GRID_X; i++) {
            double x = voxelCenterX(i);
            boolean in = near ? (x < NEAR_X) : (x > FAR_X);
            if (!in) continue;
            for (int j = 0; j < GRID_Y; j++) {
                for (int k = 0; k < GRID_Z; k++) {
                    sum += field.getConc(i, j, k);
                    n++;
                }
            }
        }
        return n == 0 ? 0.0 : (sum / n) / MOL_PER_UM3_PER_UM;
    }

    static void clampNonNegative(BSimChemicalField field) {
        for (int i = 0; i < GRID_X; i++) {
            for (int j = 0; j < GRID_Y; j++) {
                for (int k = 0; k < GRID_Z; k++) {
                    if (field.getConc(i, j, k) < 0.0) {
                        field.setConc(i, j, k, 0.0);
                    }
                }
            }
        }
    }

    static void run(String arm, boolean preview) {
        final boolean on = "on".equals(arm);
        if (!on && !"off".equals(arm)) {
            throw new IllegalArgumentException("arm must be on or off: " + arm);
        }
        final double kSource = on ? K_SOURCE_ON : 0.0;
        final String runLabel = on ? "m3on" : "m3off";
        final String exportDir = "results/" + runLabel + "_seed101/";
        final String csvName = "m3_timeseries.csv";
        final Vector3d source = new Vector3d(SOURCE_X, SOURCE_Y, SOURCE_Z);

        experimentRng = new Random(RNG_SEED);
        MonodBacterium.reseed();
        birthsNear = 0;
        birthsFar = 0;

        BSim sim = new BSim();
        sim.setDt(DT);
        sim.setSimulationTime(SIM_TIME);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);

        final BSimChemicalField glucoseField = new BSimChemicalField(
                sim, new int[]{GRID_X, GRID_Y, GRID_Z}, GLU_DIFFUSIVITY, GLU_DECAY_RATE);
        glucoseField.setConc(G0_UM * MOL_PER_UM3_PER_UM);

        final Vector<MonodBacterium> bacteria = new Vector<MonodBacterium>();
        final Vector<MonodBacterium> children = new Vector<MonodBacterium>();
        while (bacteria.size() < N_CELLS) {
            Vector3d p = new Vector3d(
                    10.0 + experimentRng.nextDouble() * (BOUND_X - 20.0),
                    10.0 + experimentRng.nextDouble() * (BOUND_Y - 20.0),
                    BOUND_Z / 2.0);
            MonodBacterium b = new MonodBacterium(sim, p, glucoseField);
            b.setRadius();
            b.setSurfaceAreaGrowthRate(GROWTH_RATE_SA);
            b.setChildList(children);
            if (!b.intersection(bacteria)) bacteria.add(b);
        }

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                glucoseField.update();
                glucoseField.addQuantity(source, kSource * sim.getDt());
                for (MonodBacterium b : bacteria) {
                    b.action();
                    b.updatePosition();
                }
                clampNonNegative(glucoseField);
                bacteria.addAll(children);
                children.clear();
            }
        });

        sim.setDrawer(new BSimP3DDrawer(sim, 800, 600) {
            @Override
            public void scene(PGraphics3D p3d) {
                p3d.ortho(0, (float) BOUND_X,
                        (float) BOUND_Y, 0,
                        -1000, 10000);
                p3d.camera((float) BOUND_X / 2f, (float) BOUND_Y / 2f, (float) BOUND_Y,
                        (float) BOUND_X / 2f, (float) BOUND_Y / 2f, 0f,
                        0f, 1f, 0f);
                p3d.perspective((float) Math.PI / 2f,
                        (float) BOUND_X / (float) BOUND_Y,
                        0.1f, 10000f);

                draw(glucoseField, Color.GREEN,
                        (float) (255.0 / (G0_UM * MOL_PER_UM3_PER_UM)));
                sphere(source, 4.0, Color.RED, 255);
                for (MonodBacterium b : bacteria) {
                    draw(b, Color.GREEN);
                }
            }
        });

        System.out.println("Monod M3 " + runLabel);
        System.out.println("  mu_max = " + MU_MAX + " /s");
        System.out.println("  K_s = " + K_S_UM + " uM");
        System.out.println("  GROWTH_RATE_SA = " + GROWTH_RATE_SA);
        System.out.println("  K_UPTAKE = " + K_UPTAKE + " molecules/s");
        System.out.println("  K_SOURCE = " + kSource
                + " molecules/s (engineering identity rate; on=" + K_SOURCE_ON + ")");
        System.out.println("  G0 = " + G0_UM + " uM (K_s)");
        System.out.println("  source = (" + SOURCE_X + ", " + SOURCE_Y + ", " + SOURCE_Z + ")");
        System.out.println("  near: voxel centre x < " + NEAR_X + " um");
        System.out.println("  far:  voxel centre x > " + FAR_X + " um");
        System.out.println("  D = " + GLU_DIFFUSIVITY + " um^2/s  decay = " + GLU_DECAY_RATE);
        System.out.println("  FLOW_SPEED = " + FLOW_SPEED);
        System.out.println("  no bath reset; no replenisher; no AC; no setGoal");

        if (preview) {
            System.out.println("Monod M3: preview mode (" + runLabel + ")");
            sim.preview();
            return;
        }

        BSimUtils.generateDirectoryPath(exportDir);
        final String armName = on ? "on" : "off";
        BSimLogger logger = new BSimLogger(sim, exportDir + csvName) {
            @Override
            public void before() {
                super.before();
                write("t_s;G_uM_mean;G_uM_near;G_uM_far;N;births_near_cum;births_far_cum;arm");
            }

            @Override
            public void during() {
                write(String.format(Locale.US, "%.0f;%.12f;%.12f;%.12f;%d;%d;%d;%s",
                        (double) Math.round(sim.getTime()),
                        meanFieldUm(glucoseField),
                        meanRegionUm(glucoseField, true),
                        meanRegionUm(glucoseField, false),
                        bacteria.size(),
                        birthsNear,
                        birthsFar,
                        armName));
            }
        };
        logger.setDt(10.0);
        sim.addExporter(logger);

        System.out.println("Monod M3: exporting to " + exportDir + csvName);
        PrintStream originalOut = System.out;
        try {
            System.setOut(new PrintStream(OutputStream.nullOutputStream()));
            sim.export();
        } finally {
            System.setOut(originalOut);
        }
        System.out.println("Monod M3 " + runLabel
                + ": export finished. Completeness is the CSV last row t=7200.");
    }

    public static void main(String[] args) {
        boolean preview = args.length > 0 && "preview".equals(args[0]);
        if (preview) {
            String arm = args.length > 1 ? args[1] : "on";
            run(arm, true);
            return;
        }
        if (args.length == 0) {
            run("on", false);
            run("off", false);
            return;
        }
        run(args[0], false);
    }
}
