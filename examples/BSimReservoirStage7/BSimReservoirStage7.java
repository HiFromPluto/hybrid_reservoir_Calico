package BSimReservoirStage7;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import javax.vecmath.Vector3d;

import processing.core.PGraphics3D;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.ode.BSimOdeSolver;
import bsim.ode.BSimOdeSystem;
import bsim.particle.BSimBacterium;

/**
 * Stage 7: Wire real chemotaxis + fix growth rate.
 *
 * Changes from Stage 6:
 *   1. setGoal(signalField) — activates BSimBacterium's Berg-cited run/tumble
 *      gradient response (7 citations, best-grounded mechanism in the codebase).
 *      Run/tumble parameters left at defaults: pEndRunUp=1/1.07, pEndRunElse=1/0.86,
 *      pEndTumble=1/0.14, memory 1s/3s — all Berg/Segall cited.
 *
 *   2. GROWTH_RATE fixed from 0.05 (vesicle-derived, T_gen=251s, ~5x too fast)
 *      to 0.00698 (30 min doubling = 1800s, standard E. coli).
 *      Formula: T_gen = 4*pi / rate. Rate = 4*pi / 1800 = 0.00698 um^2/s.
 *      Growth form (constant surface-area accumulation) is Reshes et al.-cited.
 *
 * NOT changed (per roadmap):
 *   - sensitivity=1 (uncited, ~1.7nM, harmless while field is molecule-count units;
 *     flagged for rescaling if/when Stage 8 moves AHL to uM)
 *   - Death mechanism (Stage 9's job — currently density-dependent stochastic,
 *     ungrounded, per audit §4)
 *   - Gene expression ODE (Stage 8's job — currently reads signalField as inducer,
 *     conflating chemotaxis target with QS input; Stage 8 will un-conflate by
 *     introducing a separate bacterium-sourced AHL field)
 *   - Run/tumble parameters (all cited, do not touch)
 */
public class BSimReservoirStage7 {

    // =================================================================
    //  Domain (fixed across all stages)
    // =================================================================
    static final double BOUND_X = 1000.0;
    static final double BOUND_Y = 500.0;
    static final double BOUND_Z = 10.0;

    // =================================================================
    //  Grid
    // =================================================================
    static final int GRID_X = 200;
    static final int GRID_Y = 100;
    static final int GRID_Z = 1;

    // =================================================================
    //  Signal field parameters (Stage 3)
    // =================================================================
    static final double DIFFUSIVITY = 100.0;   // um^2/s
    static final double DECAY_RATE  = 0.01;    // 1/s  (tau = 100s)
    static final double PROD_RATE   = 1e6;     // molecules/s per active AC

    // =================================================================
    //  AC parameters (Stage 4)
    // =================================================================
    static final double DOSE_THRESHOLD = 1.0;
    static final double DOSE_HIGH      = 10.0;
    static final double DOSE_RADIUS    = 30.0;  // um
    static final double IMPULSE_DURATION = 1.0;
    static final boolean IMPULSE_MODE   = false;
    static final String SEQ_FILE_PATTERN = "input_sequence%d.txt";

    // =================================================================
    //  Flow parameters (Stage 5)
    // =================================================================
    static final double FLOW_SPEED = 10.0;     // um/s along +x (gentler than Stage 5)

    // =================================================================
    //  Gene expression parameters (Stage 1)
    // =================================================================
    static final double alpha  = 1.0;
    static final double alpha0 = 0.01;
    static final double gamma  = 0.002;
    static final double K      = 50.0;
    static final int    n_hill = 2;

    // =================================================================
    //  Growth / death parameters
    // =================================================================
    // Growth rate: 30 min doubling (standard E. coli).
    // T_gen = 4*pi / rate => rate = 4*pi / 1800 = 0.00698 um^2/s.
    // Growth form (constant surface-area accumulation) is Reshes et al.-cited
    // in BSimBacterium.java:202. The OLD value 0.05 was vesicle-derived
    // (BSimReplication demo), yielding T_gen=251s (~5x too fast) — audit §4.
    static final double GROWTH_RATE       = 4.0 * Math.PI / 1800.0;  // 0.00698 um^2/s, 30 min doubling
    static final int    CARRYING_CAPACITY = 100;
    static final int    INITIAL_POP       = 20;
    static final double EXPECTED_T_GEN    = 4.0 * Math.PI / GROWTH_RATE;  // = 1800s = 30 min
    // Death: density-dependent stochastic. UNGROUNDED (audit §4, flagged for Stage 9).
    static final double DEATH_RATE_AT_CAP = 1.0 / EXPECTED_T_GEN;

    // =================================================================
    //  AC class (from Stage 4 — per-AC independent sequences)
    // =================================================================
    static class AC {
        final int id;
        final Vector3d position;
        boolean activated = false;
        final int[] sequence;
        final double bitDuration;
        final double totalTime;

        AC(int id, Vector3d position, int[] sequence, double bitDuration) {
            this.id = id;
            this.position = position;
            this.sequence = sequence;
            this.bitDuration = bitDuration;
            this.totalTime = sequence.length * bitDuration;
        }

        boolean shouldProduce(double t) {
            int bitIndex = (int)(t / bitDuration);
            if (bitIndex >= sequence.length) bitIndex = sequence.length - 1;
            int bitValue = sequence[bitIndex];
            if (IMPULSE_MODE) {
                double tInBit = t - bitIndex * bitDuration;
                return (bitValue == 1) && (tInBit < IMPULSE_DURATION);
            } else {
                return (bitValue == 1);
            }
        }
    }

    static Object[] readSequenceFile(String path, int acId) {
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            double bitDuration = Double.parseDouble(br.readLine().trim());
            String line = br.readLine().trim();
            String[] tokens = line.split("\\s+");
            int[] seq = new int[tokens.length];
            for (int i = 0; i < tokens.length; i++)
                seq[i] = Integer.parseInt(tokens[i]);
            return new Object[]{bitDuration, seq};
        } catch (IOException e) {
            System.err.println("Could not read " + path + " for AC " + acId + ", using default");
            return new Object[]{5.0, new int[]{1, 0, 1, 0, 1}};
        }
    }

    // =================================================================
    //  Bacterium with gene expression + growth + division + death
    //  (from Stage 1+2, now sensing chemical field as inducer)
    // =================================================================
    static int nextId = 0;
    static final Vector<ReservoirBacterium> bacteria = new Vector<>();
    static final Vector<ReservoirBacterium> children = new Vector<>();
    static final Vector<ReservoirBacterium> removals = new Vector<>();

    static class ReservoirBacterium extends BSimBacterium {
        final int id;
        double lastDivisionTime;
        double[] y;
        BSimOdeSystem grn;
        final BSimChemicalField signalField;

        public ReservoirBacterium(BSim sim, Vector3d position, double birthTime,
                                  BSimChemicalField signalField) {
            super(sim, position);
            this.id = nextId++;
            this.lastDivisionTime = birthTime;
            this.y = new double[]{ alpha0 / gamma };
            this.signalField = signalField;
            this.grn = new ReporterGRN();
            // Stage 7: wire chemotaxis toward AC-driven signal field.
            // Activates BSimBacterium's Berg-cited run/tumble gradient response.
            // NOTE (flagged for Stage 8): signalField currently serves double duty
            // as both chemotaxis target AND gene-expression inducer. Stage 8 will
            // un-conflate these by introducing a separate bacterium-sourced AHL field.
            setGoal(signalField);
        }

        @Override
        public void action() {
            // Full BSimBacterium physics: Brownian + run/tumble + growth
            super.action();

            // Apply flow drag force: F = stokesCoefficient * flowSpeed
            // updatePosition() converts F → v via Stokes law
            double flowForce = stokesCoefficient() * FLOW_SPEED;
            addForce(new Vector3d(flowForce, 0, 0));

            // Gene expression ODE — inducer = local signal concentration
            y = BSimOdeSolver.rungeKutta45(grn, sim.getTime(), y, sim.getDt());

            // Density-dependent stochastic death
            double pDeath = sim.getDt() * DEATH_RATE_AT_CAP
                    * ((double) bacteria.size() / CARRYING_CAPACITY);
            if (Math.random() < pDeath) {
                removals.add(this);
            }
        }

        @SuppressWarnings("unchecked")
        @Override
        public void replicate() {
            double divisionTime = sim.getTime();
            setRadiusFromSurfaceArea(surfaceArea(replicationRadius) / 2);

            Vector3d childPos = new Vector3d(position);
            childPos.x += 2.0 * radius * (Math.random() - 0.5);
            childPos.y += 2.0 * radius * (Math.random() - 0.5);

            ReservoirBacterium child = new ReservoirBacterium(sim,
                    childPos, divisionTime, signalField);
            child.setRadius(radius);
            child.setSurfaceAreaGrowthRate(surfaceAreaGrowthRate);
            child.setChildList(childList);
            child.y = new double[]{ this.y[0] };
            // setGoal already called in child's constructor

            childList.add(child);
            this.lastDivisionTime = divisionTime;
        }

        public double reporterNormalised() {
            double yMax = (alpha + alpha0) / gamma;
            return Math.min(1.0, Math.max(0.0, y[0] / yMax));
        }

        class ReporterGRN implements BSimOdeSystem {
            @Override
            public double[] derivativeSystem(double t, double[] y) {
                // Inducer = local signal concentration from the field
                double I = signalField.getConc(position);
                double hillTerm = (alpha * Math.pow(I, n_hill))
                        / (Math.pow(K, n_hill) + Math.pow(I, n_hill));
                return new double[]{ hillTerm + alpha0 - gamma * y[0] };
            }
            @Override
            public int getNumEq() { return 1; }
            @Override
            public double[] getICs() { return new double[]{ alpha0 / gamma }; }
        }
    }

    // =================================================================
    //  Upwind advection (from Stage 5)
    // =================================================================
    static void advect(BSimChemicalField field, double flowSpeed, double dt,
                       int nx, int ny, int nz) {
        double dx = field.getBox()[0];
        double courant = flowSpeed * dt / dx;
        double[][] before = new double[nx][ny];
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                before[i][j] = field.getConc(i, j, 0);

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++) {
                double cUpwind = (i > 0) ? before[i - 1][j] : 0;
                double newConc = before[i][j] - courant * (before[i][j] - cUpwind);
                if (newConc < 0) newConc = 0;
                field.setConc(i, j, 0, newConc);
            }
    }

    // =================================================================
    //  Main
    // =================================================================
    public static void main(String[] args) {

        boolean exportData = true;
        if (args.length > 0 && args[0].equals("preview")) {
            exportData = false;
        }

        // --- Load per-AC sequences ---
        Vector3d[] acPositions = {
            new Vector3d(250, 250, 5),
            new Vector3d(500, 250, 5),
            new Vector3d(750, 250, 5)
        };
        final List<AC> acs = new ArrayList<>();
        double maxACTime = 0;
        for (int i = 0; i < acPositions.length; i++) {
            String seqFile = String.format(SEQ_FILE_PATTERN, i);
            Object[] result = readSequenceFile(seqFile, i);
            double bitDur = (Double) result[0];
            int[] seq = (int[]) result[1];
            AC ac = new AC(i, acPositions[i], seq, bitDur);
            acs.add(ac);
            if (ac.totalTime > maxACTime) maxACTime = ac.totalTime;
            System.out.print("AC " + i + " (" + bitDur + "s/bit, "
                    + seq.length + " bits, " + ac.totalTime + "s): ");
            for (int b : seq) System.out.print(b + " ");
            System.out.println();
        }

        // Simulation time: AC sequences + enough for growth/chemotaxis validation.
        // With 1800s doubling, need ~3600s to see at least one division cycle.
        final double simTime = Math.max(maxACTime + 30.0, 3600.0);

        String timestamp = BSimUtils.timeStamp();
        String exportPath = "./results/" + timestamp + "/";

        double voxelX = BOUND_X / GRID_X;
        double courant = FLOW_SPEED * 0.05 / voxelX;

        System.out.println("Stage 7: Chemotaxis + corrected growth rate");
        System.out.println("  Flow: " + FLOW_SPEED + " um/s, CFL=" + courant);
        System.out.println("  Growth rate: " + String.format("%.5f", GROWTH_RATE) + " um^2/s");
        System.out.println("  Expected T_gen: " + String.format("%.0f", EXPECTED_T_GEN) + " s = "
                + String.format("%.1f", EXPECTED_T_GEN / 60.0) + " min");
        System.out.println("  Initial pop: " + INITIAL_POP + ", carrying capacity: " + CARRYING_CAPACITY);
        System.out.println("  Sim time: " + simTime + "s");
        System.out.println("  Chemotaxis: setGoal(signalField) active");

        // --- Simulation ---
        BSim sim = new BSim();
        sim.setDt(0.05);
        sim.setSimulationTime(simTime);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);

        // --- Signal field ---
        final BSimChemicalField signalField = new BSimChemicalField(sim,
                new int[]{GRID_X, GRID_Y, GRID_Z}, DIFFUSIVITY, DECAY_RATE);

        // --- Dose field (for AC activation gating) ---
        final BSimChemicalField doseField = new BSimChemicalField(sim,
                new int[]{GRID_X, GRID_Y, GRID_Z}, 0, 0);

        final double voxelXval = BOUND_X / GRID_X;
        final double voxelYval = BOUND_Y / GRID_Y;

        // --- Seed bacteria scattered across domain centre ---
        for (int i = 0; i < INITIAL_POP; i++) {
            double bx = 300 + Math.random() * 400;  // x: 300-700
            double by = 150 + Math.random() * 200;  // y: 150-350
            double bz = BOUND_Z / 2.0;
            ReservoirBacterium b = new ReservoirBacterium(sim,
                    new Vector3d(bx, by, bz), 0.0, signalField);
            b.setRadius();
            b.setSurfaceAreaGrowthRate(GROWTH_RATE);
            b.setChildList(children);
            bacteria.add(b);
        }

        // --- Ticker ---
        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                double t = sim.getTime();

                // ---- AC activation (Stage 4) ----
                // Clear dose field
                for (int i = 0; i < GRID_X; i++)
                    for (int j = 0; j < GRID_Y; j++)
                        doseField.setConc(i, j, 0, 0);

                // Apply dose per-AC based on independent sequences
                for (AC ac : acs) {
                    if (t < ac.totalTime && ac.shouldProduce(t)) {
                        int cxMin = Math.max(0, (int)((ac.position.x - DOSE_RADIUS) / voxelXval));
                        int cxMax = Math.min(GRID_X - 1, (int)((ac.position.x + DOSE_RADIUS) / voxelXval));
                        int cyMin = Math.max(0, (int)((ac.position.y - DOSE_RADIUS) / voxelYval));
                        int cyMax = Math.min(GRID_Y - 1, (int)((ac.position.y + DOSE_RADIUS) / voxelYval));
                        for (int i = cxMin; i <= cxMax; i++)
                            for (int j = cyMin; j <= cyMax; j++)
                                doseField.setConc(i, j, 0, DOSE_HIGH);
                    }
                }

                // Each AC checks dose and produces signal if activated
                for (AC ac : acs) {
                    double dose = doseField.getConc(ac.position);
                    ac.activated = dose > DOSE_THRESHOLD;
                    if (ac.activated) {
                        signalField.addQuantity(ac.position, PROD_RATE * sim.getDt());
                    }
                }

                // ---- Signal field: diffuse + decay + advect (Stages 3+5) ----
                signalField.update();
                advect(signalField, FLOW_SPEED, sim.getDt(), GRID_X, GRID_Y, GRID_Z);

                // ---- Bacteria: action + position update (Stages 1+2+5) ----
                for (ReservoirBacterium b : bacteria) {
                    b.action();
                    b.updatePosition();
                }

                // Remove dead / out-of-bounds bacteria
                for (ReservoirBacterium b : bacteria) {
                    if (b.getPosition().x > BOUND_X || b.getPosition().x < 0) {
                        if (!removals.contains(b)) removals.add(b);
                    }
                }
                bacteria.removeAll(removals);
                removals.clear();

                // Add children
                bacteria.addAll(children);
                children.clear();
            }
        });

        // --- Drawer ---
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

                // Signal field (cyan)
                draw(signalField, Color.CYAN, (float)(255.0 / 540.0));

                // ACs: green if activated, red if not
                for (AC ac : acs) {
                    Color col = ac.activated ? Color.GREEN : Color.RED;
                    sphere(ac.position, 15.0, col, 255);
                }

                // Bacteria: brightness tracks reporter protein level
                for (ReservoirBacterium b : bacteria) {
                    double r = b.reporterNormalised();
                    Color col = new Color(30, (int)(55 + 200 * r), 30);
                    sphere(b.getPosition(), 8.0, col, 255);
                }

                // Flow direction arrow
                p3d.stroke(255, 255, 0);
                p3d.strokeWeight(2);
                p3d.line(50, 30, 0, 200, 30, 0);
                p3d.line(200, 30, 0, 180, 20, 0);
                p3d.line(200, 30, 0, 180, 40, 0);
                p3d.noStroke();
            }
        });

        // --- Exporters ---
        if (exportData) {
            BSimUtils.generateDirectoryPath(exportPath);

            // AC + population state (every 1s)
            BSimLogger stateLogger = new BSimLogger(sim, exportPath + "stage7_state.csv") {
                @Override
                public void before() {
                    super.before();
                    StringBuilder header = new StringBuilder(
                            "time_s,population,mean_reporter");
                    for (AC ac : acs) {
                        header.append(",ac").append(ac.id).append("_activated");
                        header.append(",ac").append(ac.id).append("_signal_conc");
                    }
                    write(header.toString());
                }
                @Override
                public void during() {
                    double t = sim.getTime();
                    double meanReporter = 0;
                    for (ReservoirBacterium b : bacteria)
                        meanReporter += b.y[0];
                    if (!bacteria.isEmpty()) meanReporter /= bacteria.size();

                    StringBuilder line = new StringBuilder();
                    line.append(sim.getFormattedTime());
                    line.append(",").append(bacteria.size());
                    line.append(",").append(String.format("%.4f", meanReporter));

                    for (AC ac : acs) {
                        line.append(",").append(ac.activated ? 1 : 0);
                        line.append(",").append(String.format("%.2f",
                                signalField.getConc(ac.position)));
                    }
                    write(line.toString());
                }
            };
            stateLogger.setDt(1.0);
            sim.addExporter(stateLogger);

            // Per-bacterium snapshot (every 30s — position + reporter + distance to nearest AC)
            BSimLogger bacLogger = new BSimLogger(sim, exportPath + "stage7_bacteria.csv") {
                @Override
                public void before() {
                    super.before();
                    write("time_s,bac_id,x,y,z,reporter,signal_at_pos,dist_nearest_ac");
                }
                @Override
                public void during() {
                    String t = sim.getFormattedTime();
                    for (ReservoirBacterium b : bacteria) {
                        Vector3d p = b.getPosition();
                        double minDist = Double.MAX_VALUE;
                        for (AC ac : acs) {
                            double dx = p.x - ac.position.x;
                            double dy = p.y - ac.position.y;
                            double d = Math.sqrt(dx * dx + dy * dy);
                            if (d < minDist) minDist = d;
                        }
                        write(t + "," + b.id
                                + "," + String.format("%.1f", p.x)
                                + "," + String.format("%.1f", p.y)
                                + "," + String.format("%.1f", p.z)
                                + "," + String.format("%.4f", b.y[0])
                                + "," + String.format("%.4f", signalField.getConc(p))
                                + "," + String.format("%.1f", minDist));
                    }
                }
            };
            bacLogger.setDt(30.0);
            sim.addExporter(bacLogger);

            // Parameters
            BSimLogger paramLogger = new BSimLogger(sim, exportPath + "stage7_params.csv") {
                @Override
                public void before() {
                    super.before();
                    write("parameter,value");
                    write("num_acs," + acs.size());
                    for (AC ac : acs) {
                        StringBuilder seqStr = new StringBuilder();
                        for (int i = 0; i < ac.sequence.length; i++) {
                            if (i > 0) seqStr.append(" ");
                            seqStr.append(ac.sequence[i]);
                        }
                        write("ac" + ac.id + "_sequence," + seqStr);
                        write("ac" + ac.id + "_bit_duration_s," + ac.bitDuration);
                    }
                    write("flow_speed_um_per_s," + FLOW_SPEED);
                    write("diffusivity," + DIFFUSIVITY);
                    write("decay_rate," + DECAY_RATE);
                    write("prod_rate," + PROD_RATE);
                    write("initial_pop," + INITIAL_POP);
                    write("carrying_capacity," + CARRYING_CAPACITY);
                    write("growth_rate," + GROWTH_RATE);
                    write("expected_T_gen," + EXPECTED_T_GEN);
                    write("alpha," + alpha);
                    write("gamma," + gamma);
                    write("K," + K);
                    write("n_hill," + n_hill);
                    write("dt," + sim.getDt());
                    write("sim_time," + sim.getSimulationTime());
                    write("chemotaxis,setGoal(signalField)");
                    write("sensitivity,1 (uncited; ~1.7 nM; flagged for Stage 8)");
                }
                @Override
                public void during() { }
            };
            sim.addExporter(paramLogger);

            System.out.println("Stage 7: exporting to " + exportPath);
            sim.export();

        } else {
            System.out.println("Stage 7: preview mode");
            sim.preview();
        }
    }
}
