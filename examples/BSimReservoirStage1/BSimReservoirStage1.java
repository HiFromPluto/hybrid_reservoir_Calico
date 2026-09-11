package BSimReservoirStage1;

import java.awt.Color;
import java.util.Vector;

import javax.vecmath.Vector3d;

import processing.core.PGraphics3D;
import bsim.BSim;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.ode.BSimOdeSolver;
import bsim.ode.BSimOdeSystem;
import bsim.particle.BSimBacterium;

/**
 * Stage 1: Single stationary bacterium with inducible gene expression.
 *
 * ODE model (1 variable):
 *   dy/dt = (alpha * I^n) / (K^n + I^n) + alpha0 - gamma * y
 *
 * where I is the external inducer concentration (set as a step input),
 * y is the reporter protein concentration, and all other symbols are
 * kinetic parameters.
 *
 * Validation criterion: step-input response shape — exponential rise
 * with time constant 1/gamma, from basal level alpha0/gamma to
 * saturated level (alpha + alpha0)/gamma.
 */
public class BSimReservoirStage1 {

    // --- Gene expression parameters ---
    static final double alpha  = 1.0;     // max induced production rate (AU/s)
    static final double alpha0 = 0.01;    // basal production rate (AU/s)
    static final double gamma  = 0.002;   // degradation+dilution rate (1/s) => tau ~500 s
    static final double K      = 50.0;    // Hill half-max inducer concentration (AU)
    static final int    n      = 2;       // Hill coefficient

    // --- Step input parameters ---
    static final double stepTime    = 600.0;  // time of step input (s)
    static final double inducerHigh = 200.0;  // inducer concentration after step (AU)

    // Current inducer level (mutable, updated by ticker)
    static double inducerLevel = 0.0;

    // --- Domain (from plan, fixed across all stages) ---
    static final double BOUND_X = 1000.0; // um
    static final double BOUND_Y = 500.0;
    static final double BOUND_Z = 10.0;

    /**
     * A stationary bacterium with an internal gene-expression ODE.
     * No run/tumble, no growth — purely an ODE container for Stage 1.
     */
    static class ReservoirBacterium extends BSimBacterium {

        double[] y;  // ODE state: y[0] = reporter protein concentration
        BSimOdeSystem grn;

        public ReservoirBacterium(BSim sim, Vector3d position) {
            super(sim, position);
            y = new double[]{ alpha0 / gamma }; // start at basal steady state
            grn = new ReporterGRN();
        }

        @Override
        public void action() {
            // No super.action() — skip run/tumble and growth entirely.
            // Solve ODE for one timestep.
            y = BSimOdeSolver.rungeKutta45(grn, sim.getTime(), y, sim.getDt());
        }

        /** Normalised reporter level in [0,1] for visualisation. */
        public double reporterNormalised() {
            double yMax = (alpha + alpha0) / gamma; // theoretical max steady state
            return Math.min(1.0, Math.max(0.0, y[0] / yMax));
        }

        /**
         * Simple inducible gene expression.
         * Reads the current inducer level from the enclosing class's static field.
         */
        class ReporterGRN implements BSimOdeSystem {
            @Override
            public double[] derivativeSystem(double t, double[] y) {
                double I = inducerLevel;
                double hillTerm = (alpha * Math.pow(I, n)) / (Math.pow(K, n) + Math.pow(I, n));
                return new double[]{ hillTerm + alpha0 - gamma * y[0] };
            }

            @Override
            public int getNumEq() { return 1; }

            @Override
            public double[] getICs() { return new double[]{ alpha0 / gamma }; }
        }
    }

    public static void main(String[] args) {

        // --- Run mode ---
        boolean exportData = true;  // true = headless CSV export; false = GUI preview
        if (args.length > 0 && args[0].equals("preview")) {
            exportData = false;
        }

        String timestamp = BSimUtils.timeStamp();
        String exportPath = "./results/" + timestamp + "/";

        // --- Simulation setup (plan-specified domain) ---
        BSim sim = new BSim();
        sim.setDt(0.01);                        // 10 ms — fine enough for RK4
        sim.setSimulationTime(3600);             // 1 hour
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);

        // --- Single bacterium, placed at domain centre, stationary ---
        final ReservoirBacterium bacterium = new ReservoirBacterium(sim,
                new Vector3d(BOUND_X / 2.0, BOUND_Y / 2.0, BOUND_Z / 2.0));

        // --- Ticker ---
        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                // Step input: inducer jumps at stepTime
                if (sim.getTime() >= stepTime) {
                    inducerLevel = inducerHigh;
                } else {
                    inducerLevel = 0.0;
                }
                bacterium.action();
                // No updatePosition() needed — bacterium is stationary.
            }
        });

        // --- Drawer ---
        sim.setDrawer(new BSimP3DDrawer(sim, 800, 600) {
            @Override
            public void scene(PGraphics3D p3d) {
                // Ortho + camera + perspective, in this order — matches the
                // proven-working setup from the previous project's
                // ReservoirSim_2.java (same 1000x500x10 domain). camera()
                // alone is not sufficient here: the eye's pull-back distance
                // along z uses BOUND_Y (500), not BOUND_Z (10) — any offset
                // relative to the 10-unit-thick z axis puts the eye almost
                // inside the domain, which is the degenerate-sliver failure
                // mode both earlier attempts hit. perspective() has to follow
                // camera() or ortho() mispositions the box in the corner.
                p3d.ortho(0, (float) BOUND_X,
                          (float) BOUND_Y, 0,
                          -1000, 10000);
                p3d.camera((float) BOUND_X / 2f, (float) BOUND_Y / 2f, (float) BOUND_Y,
                           (float) BOUND_X / 2f, (float) BOUND_Y / 2f, 0f,
                           0f, 1f, 0f);
                p3d.perspective((float) Math.PI / 2f,
                        (float) BOUND_X / (float) BOUND_Y,
                        0.1f, 10000f);

                // Bacterium colour: green intensity tracks reporter level.
                // Draw as 10 um sphere (actual radius is ~1 um, sub-pixel
                // in this 1000 um domain at 800 px).
                double r = bacterium.reporterNormalised();
                Color col = new Color(30, (int)(55 + 200 * r), 30);
                sphere(bacterium.getPosition(), 10.0, col, 255);
            }
        });

        // --- Exporters ---
        if (exportData) {
            BSimUtils.generateDirectoryPath(exportPath);

            // Main time-series logger: time, inducer, protein, expected_ss
            BSimLogger tsLogger = new BSimLogger(sim, exportPath + "stage1_timeseries.csv") {
                @Override
                public void before() {
                    super.before();
                    write("time_s,inducer,protein,expected_ss");
                }

                @Override
                public void during() {
                    // Expected steady state for current inducer level
                    double I = inducerLevel;
                    double hillTerm = (alpha * Math.pow(I, n)) / (Math.pow(K, n) + Math.pow(I, n));
                    double expectedSS = (hillTerm + alpha0) / gamma;

                    write(sim.getFormattedTime()
                            + "," + inducerLevel
                            + "," + bacterium.y[0]
                            + "," + expectedSS);
                }
            };
            tsLogger.setDt(1.0); // log every 1 second
            sim.addExporter(tsLogger);

            // Parameters logger
            BSimLogger paramLogger = new BSimLogger(sim, exportPath + "stage1_params.csv") {
                @Override
                public void before() {
                    super.before();
                    write("parameter,value");
                    write("alpha," + alpha);
                    write("alpha0," + alpha0);
                    write("gamma," + gamma);
                    write("K," + K);
                    write("n," + n);
                    write("stepTime_s," + stepTime);
                    write("inducerHigh," + inducerHigh);
                    write("dt," + sim.getDt());
                    write("simTime_s," + sim.getSimulationTime());
                    double hillAtStep = (alpha * Math.pow(inducerHigh, n))
                            / (Math.pow(K, n) + Math.pow(inducerHigh, n));
                    write("basal_ss," + (alpha0 / gamma));
                    write("induced_ss_actual," + ((hillAtStep + alpha0) / gamma));
                    write("induced_ss_max," + ((alpha + alpha0) / gamma));
                    write("tau_s," + (1.0 / gamma));
                    write("domain," + BOUND_X + "x" + BOUND_Y + "x" + BOUND_Z);
                }

                @Override
                public void during() { }
            };
            sim.addExporter(paramLogger);

            System.out.println("Stage 1: exporting to " + exportPath);
            double hillAtStep = (alpha * Math.pow(inducerHigh, n))
                    / (Math.pow(K, n) + Math.pow(inducerHigh, n));
            System.out.println("  Basal SS  = " + (alpha0 / gamma));
            System.out.println("  Induced SS (actual at I=" + inducerHigh + ") = "
                    + ((hillAtStep + alpha0) / gamma));
            System.out.println("  Tau (s)    = " + (1.0 / gamma));
            sim.export();

        } else {
            System.out.println("Stage 1: preview mode (GUI)");
            System.out.println("  Step input at t=" + stepTime + "s, inducer -> " + inducerHigh);
            sim.preview();
        }
    }
}
