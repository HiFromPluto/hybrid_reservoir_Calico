package BSimMonodM1;

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
 * Isolated Monod M1 example: µ = µ_max * G / (K_s + G) on a uniform glucose
 * bath. No uptake. Not a reservoir window experiment.
 *
 * Growth law (frozen):
 *   µ = µ_max * G / (K_s + G)
 *   K_s    = 0.18 µM   (Senn et al. 1994; BNID 111049 range 0.18–0.55 µM)
 *   µ_max  = Math.log(2.0) / 1800.0   // 30 min doubling, THIS rebuild
 *
 * Do not use reservoir_new's unitless mu_max=1, K_G=5.
 * Do not use Senn's 0.92 h^-1 here — that would be a 45 min doubling and
 * would disagree with Stage 10 / Plan GROWTH_RATE.
 *
 * Stage 11's Monod+nutrient-AC merge stays a failed reservoir path. M1
 * does not reopen it. Do not copy glucose into BSimReservoirPlanStage6.
 */
public class BSimMonodM1 {

    static final double BOUND_X = 200.0;
    static final double BOUND_Y = 200.0;
    static final double BOUND_Z = 10.0;

    static final int GRID_X = 20;
    static final int GRID_Y = 20;
    static final int GRID_Z = 1;

    /** Bath is overwritten every tick; D and decay cannot change G. */
    static final double GLU_DIFFUSIVITY = 0.0;
    static final double GLU_DECAY_RATE = 0.0;

    static final double DT = 0.05;
    static final double SIM_TIME = 7200.0;
    static final int N_CELLS = 50;
    static final double FLOW_SPEED = 0.0;
    static final long RNG_SEED = 101L;

    /** 1 µM = 602 molecules/µm³. Same conversion as D1 AHL. */
    static final double MOL_PER_UM3_PER_UM = 602.0;

    static final double K_S_UM = 0.18;
    static final double MU_MAX = Math.log(2.0) / 1800.0;
    static final double GROWTH_RATE_SA = 4.0 * Math.PI / 1800.0;

    static final double G_M1A_UM = 0.02;
    static final double G_M1B_UM = 0.18;
    static final double G_M1C_UM = 1.80;

    static Random experimentRng;

    /**
     * Bacterium whose surface-area growth is scaled by the Monod factor
     * G/(K_s+G) with G in µM from the bath. Open-loop: reads the field,
     * does not subtract glucose. Replication on. No clamp, no vesicles.
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
        }

        @SuppressWarnings("unchecked")
        @Override
        public void replicate() {
            setRadiusFromSurfaceArea(surfaceArea(replicationRadius) / 2);
            Vector3d childPosition = new Vector3d(position);
            childPosition.x += 2.0 * radius * (experimentRng.nextDouble() - 0.5);
            childPosition.y += 2.0 * radius * (experimentRng.nextDouble() - 0.5);
            MonodBacterium child = new MonodBacterium(sim, childPosition, glucoseField);
            child.setRadius(radius);
            child.setSurfaceAreaGrowthRate(GROWTH_RATE_SA);
            child.setChildList(childList);
            childList.add(child);
        }
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

    static double bathUm(String runLabel) {
        if ("m1a".equals(runLabel)) return G_M1A_UM;
        if ("m1b".equals(runLabel)) return G_M1B_UM;
        if ("m1c".equals(runLabel)) return G_M1C_UM;
        throw new IllegalArgumentException("run label must be m1a, m1b, or m1c: " + runLabel);
    }

    static void run(String runLabel, boolean preview) {
        final double commandedUm = bathUm(runLabel);
        final String exportDir = "results/" + runLabel + "_seed101/";
        final String csvName = "m1_timeseries.csv";

        experimentRng = new Random(RNG_SEED);
        MonodBacterium.reseed();

        BSim sim = new BSim();
        sim.setDt(DT);
        sim.setSimulationTime(SIM_TIME);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);

        final BSimChemicalField glucoseField = new BSimChemicalField(
                sim, new int[]{GRID_X, GRID_Y, GRID_Z}, GLU_DIFFUSIVITY, GLU_DECAY_RATE);
        glucoseField.setConc(commandedUm * MOL_PER_UM3_PER_UM);

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
                glucoseField.setConc(commandedUm * MOL_PER_UM3_PER_UM);
                for (MonodBacterium b : bacteria) {
                    b.action();
                    b.updatePosition();
                }
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
                        (float) (255.0 / (commandedUm * MOL_PER_UM3_PER_UM)));
                for (MonodBacterium b : bacteria) {
                    draw(b, Color.GREEN);
                }
            }
        });

        System.out.println("Monod M1 " + runLabel);
        System.out.println("  mu_max = " + MU_MAX + " /s");
        System.out.println("  K_s = " + K_S_UM + " uM");
        System.out.println("  G = " + commandedUm + " uM");
        System.out.println("  expected mu/mu_max = "
                + (commandedUm / (K_S_UM + commandedUm)));
        System.out.println("  FLOW_SPEED = " + FLOW_SPEED);
        System.out.println("  GROWTH_RATE_SA = " + GROWTH_RATE_SA);
        System.out.println("  no uptake; bath reset every tick");

        if (preview) {
            System.out.println("Monod M1: preview mode (" + runLabel + ")");
            sim.preview();
            return;
        }

        BSimUtils.generateDirectoryPath(exportDir);
        BSimLogger logger = new BSimLogger(sim, exportDir + csvName) {
            @Override
            public void before() {
                super.before();
                write("t_s;G_uM;N;run_label");
            }

            @Override
            public void during() {
                write(String.format(Locale.US, "%.0f;%.12f;%d;%s",
                        (double) Math.round(sim.getTime()),
                        meanFieldUm(glucoseField),
                        bacteria.size(),
                        runLabel));
            }
        };
        logger.setDt(10.0);
        sim.addExporter(logger);

        System.out.println("Monod M1: exporting to " + exportDir + csvName);
        PrintStream originalOut = System.out;
        try {
            System.setOut(new PrintStream(OutputStream.nullOutputStream()));
            sim.export();
        } finally {
            System.setOut(originalOut);
        }
        System.out.println("Monod M1 " + runLabel
                + ": export finished. Completeness is the CSV last row t=7200.");
    }

    public static void main(String[] args) {
        boolean preview = args.length > 0 && "preview".equals(args[0]);
        if (preview) {
            String label = args.length > 1 ? args[1] : "m1b";
            run(label, true);
            return;
        }
        if (args.length == 0) {
            run("m1a", false);
            run("m1b", false);
            run("m1c", false);
            return;
        }
        run(args[0], false);
    }
}
