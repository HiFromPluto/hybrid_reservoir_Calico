package BSimReservoirStage2;

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
 * Stage 2: Single bacterium lineage — growth, division, death.
 *
 * Builds on Stage 1's gene-expression ODE. Adds:
 *   - Growth via setSurfaceAreaGrowthRate() (built-in BSimBacterium mechanism,
 *     demonstrated in BSimReplication)
 *   - Division via replicate() override — children inherit ODE state
 *   - Death via density-dependent stochastic removal (new code; no bundled
 *     BSim example implements death — confirmed by reading BSimReplication
 *     and BSimVesiculation)
 *
 * Still NO motility — bacterium positions are updated only by Brownian
 * motion (via super.action()), so the population stays near its origin
 * and any growth-curve validation is not confounded by chemotaxis.
 *
 * Validation criterion: exponential growth phase → logistic saturation
 * at carrying capacity; division-interval distribution should be
 * consistent with the deterministic generation time set by the
 * surface-area growth rate.
 */
public class BSimReservoirStage2 {

    // --- Gene expression parameters (carried from Stage 1) ---
    static final double alpha  = 1.0;
    static final double alpha0 = 0.01;
    static final double gamma  = 0.002;
    static final double K      = 50.0;
    static final int    n      = 2;

    // Constant inducer (no step input in Stage 2 — always induced)
    static double inducerLevel = 200.0;

    // --- Domain (fixed across all stages) ---
    static final double BOUND_X = 1000.0;
    static final double BOUND_Y = 500.0;
    static final double BOUND_Z = 10.0;

    // --- Growth / death parameters ---
    static final double GROWTH_RATE = 0.05;      // um^2/s surface-area growth rate
    static final int    CARRYING_CAPACITY = 200;  // target steady-state population

    // --- Shared mutable state ---
    static final Vector<ReservoirBacterium> bacteria = new Vector<>();
    static final Vector<ReservoirBacterium> children = new Vector<>();
    static final Vector<ReservoirBacterium> removals = new Vector<>();

    // --- Division event log (populated during simulation) ---
    static final Vector<String> divisionLog = new Vector<>();

    // Counter for unique bacterium IDs
    static int nextId = 0;

    /**
     * Expected deterministic generation time (seconds).
     * T_gen = surfaceArea(replicationRadius) / (2 * growthRate)
     * With replicationRadius = sqrt(2): SA = 4*pi*2 = 8*pi
     * T_gen = 8*pi / (2 * 0.05) = 8*pi / 0.1 ≈ 251.3 s
     */
    static final double EXPECTED_T_GEN = 4.0 * Math.PI * 2.0 / (2.0 * GROWTH_RATE);

    /**
     * Death rate constant: at carrying capacity, per-cell death rate = 1/T_gen
     * (balances division). Below capacity, death is proportionally lower.
     *
     * P_death per tick = dt * (1/T_gen) * (pop / capacity)
     */
    static final double DEATH_RATE_AT_CAPACITY = 1.0 / EXPECTED_T_GEN;

    // ================================================================
    //  Bacterium with gene expression + growth + division + death
    // ================================================================
    static class ReservoirBacterium extends BSimBacterium {

        final int id;
        double lastDivisionTime;  // reset on each division (or birth)
        double[] y;  // ODE state: y[0] = reporter protein
        BSimOdeSystem grn;

        public ReservoirBacterium(BSim sim, Vector3d position, double birthTime) {
            super(sim, position);
            this.id = nextId++;
            this.lastDivisionTime = birthTime;
            this.y = new double[]{ alpha0 / gamma };
            this.grn = new ReporterGRN();
        }

        @Override
        public void action() {
            // Full BSimBacterium physics: Brownian + run/tumble + growth.
            // No chemotaxis goal set, so run/tumble is unbiased random walk.
            // This matches BSimReplication's approach — cells swim and spread.
            super.action();

            // Gene expression ODE
            y = BSimOdeSolver.rungeKutta45(grn, sim.getTime(), y, sim.getDt());

            // Density-dependent stochastic death
            double pDeath = sim.getDt() * DEATH_RATE_AT_CAPACITY
                    * ((double) bacteria.size() / CARRYING_CAPACITY);
            if (Math.random() < pDeath) {
                removals.add(this);
            }
        }

        @SuppressWarnings("unchecked")
        @Override
        public void replicate() {
            double divisionTime = sim.getTime();
            double interval = divisionTime - this.lastDivisionTime;

            // Parent shrinks to half surface area
            setRadiusFromSurfaceArea(surfaceArea(replicationRadius) / 2);

            // Child: offset by one diameter so parent and child don't overlap
            Vector3d childPos = new Vector3d(position);
            childPos.x += 2.0 * radius * (Math.random() - 0.5);
            childPos.y += 2.0 * radius * (Math.random() - 0.5);
            ReservoirBacterium child = new ReservoirBacterium(sim,
                    childPos, divisionTime);
            child.setRadius(radius);
            child.setSurfaceAreaGrowthRate(surfaceAreaGrowthRate);
            child.setChildList(childList);
            child.y = new double[]{ this.y[0] };  // inherit protein concentration

            childList.add(child);

            // Log the division event, then reset parent's timer
            divisionLog.add(String.format("%.2f,%d,%d,%.2f",
                    divisionTime, this.id, child.id, interval));
            this.lastDivisionTime = divisionTime;
        }

        /** Normalised reporter level in [0,1] for visualisation. */
        public double reporterNormalised() {
            double yMax = (alpha + alpha0) / gamma;
            return Math.min(1.0, Math.max(0.0, y[0] / yMax));
        }

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

    // ================================================================
    //  Main
    // ================================================================
    public static void main(String[] args) {

        boolean exportData = true;
        if (args.length > 0 && args[0].equals("preview")) {
            exportData = false;
        }

        String timestamp = BSimUtils.timeStamp();
        String exportPath = "./results/" + timestamp + "/";

        // --- Simulation setup ---
        BSim sim = new BSim();
        sim.setDt(0.5);                          // coarser dt OK for growth dynamics
        sim.setSimulationTime(5000);              // ~20 generations
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);

        // --- Seed with one bacterium at domain centre ---
        ReservoirBacterium ancestor = new ReservoirBacterium(sim,
                new Vector3d(BOUND_X / 2.0, BOUND_Y / 2.0, BOUND_Z / 2.0), 0.0);
        ancestor.setRadius();  // random initial size between birth and division
        ancestor.setSurfaceAreaGrowthRate(GROWTH_RATE);
        ancestor.setChildList(children);
        bacteria.add(ancestor);

        // --- Ticker ---
        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                for (ReservoirBacterium b : bacteria) {
                    b.action();
                    b.updatePosition();
                }

                // Remove dead bacteria
                bacteria.removeAll(removals);
                removals.clear();

                // Add newborn bacteria
                bacteria.addAll(children);
                children.clear();
            }
        });

        // --- Drawer ---
        sim.setDrawer(new BSimP3DDrawer(sim, 800, 600) {
            @Override
            public void scene(PGraphics3D p3d) {
                // Same camera setup as Stage 1 (user-validated)
                p3d.ortho(0, (float) BOUND_X,
                          (float) BOUND_Y, 0,
                          -1000, 10000);
                p3d.camera((float) BOUND_X / 2f, (float) BOUND_Y / 2f, (float) BOUND_Y,
                           (float) BOUND_X / 2f, (float) BOUND_Y / 2f, 0f,
                           0f, 1f, 0f);
                p3d.perspective((float) Math.PI / 2f,
                        (float) BOUND_X / (float) BOUND_Y,
                        0.1f, 10000f);

                for (ReservoirBacterium b : bacteria) {
                    double r = b.reporterNormalised();
                    Color col = new Color(30, (int)(55 + 200 * r), 30);
                    sphere(b.getPosition(), 10.0, col, 255);
                }
            }
        });

        // --- Exporters ---
        if (exportData) {
            BSimUtils.generateDirectoryPath(exportPath);

            // Population time-series (every 5 seconds)
            BSimLogger popLogger = new BSimLogger(sim, exportPath + "stage2_population.csv") {
                @Override
                public void before() {
                    super.before();
                    write("time_s,population,mean_protein");
                }
                @Override
                public void during() {
                    double meanProtein = 0;
                    for (ReservoirBacterium b : bacteria) {
                        meanProtein += b.y[0];
                    }
                    if (!bacteria.isEmpty()) meanProtein /= bacteria.size();
                    write(sim.getFormattedTime() + "," + bacteria.size() + "," + meanProtein);
                }
            };
            popLogger.setDt(5.0);
            sim.addExporter(popLogger);

            // Division events (flushed at end from divisionLog)
            final String divPath = exportPath + "stage2_divisions.csv";
            BSimLogger divLogger = new BSimLogger(sim, divPath) {
                @Override
                public void before() {
                    super.before();
                    write("time_s,parent_id,child_id,interval_s");
                }
                @Override
                public void during() {
                    // Division events are accumulated in divisionLog; flush here
                    for (String entry : divisionLog) {
                        write(entry);
                    }
                    divisionLog.clear();
                }
            };
            divLogger.setDt(5.0);
            sim.addExporter(divLogger);

            // Parameters
            BSimLogger paramLogger = new BSimLogger(sim, exportPath + "stage2_params.csv") {
                @Override
                public void before() {
                    super.before();
                    write("parameter,value");
                    write("growth_rate_um2_per_s," + GROWTH_RATE);
                    write("expected_T_gen_s," + EXPECTED_T_GEN);
                    write("carrying_capacity," + CARRYING_CAPACITY);
                    write("death_rate_at_capacity," + DEATH_RATE_AT_CAPACITY);
                    write("dt," + sim.getDt());
                    write("simTime_s," + sim.getSimulationTime());
                    write("inducer," + inducerLevel);
                    write("domain," + BOUND_X + "x" + BOUND_Y + "x" + BOUND_Z);
                }
                @Override
                public void during() { }
            };
            sim.addExporter(paramLogger);

            System.out.println("Stage 2: exporting to " + exportPath);
            System.out.println("  Expected T_gen = " + String.format("%.1f", EXPECTED_T_GEN) + " s");
            System.out.println("  Carrying capacity = " + CARRYING_CAPACITY);
            sim.export();

        } else {
            System.out.println("Stage 2: preview mode (GUI)");
            System.out.println("  Expected T_gen = " + String.format("%.1f", EXPECTED_T_GEN) + " s");
            System.out.println("  Carrying capacity = " + CARRYING_CAPACITY);
            sim.preview();
        }
    }
}
