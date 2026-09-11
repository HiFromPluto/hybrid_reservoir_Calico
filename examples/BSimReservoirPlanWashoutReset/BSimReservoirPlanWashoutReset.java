package BSimReservoirPlanWashoutReset;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.draw.BSimP3DDrawer;
import bsim.particle.BSimBacterium;
import bsim.particle.BSimParticle;
import processing.core.PGraphics3D;

import javax.vecmath.Vector3d;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Properties;
import java.util.Random;
import java.util.Vector;

/**
 * WashoutReset: occupied quiescent load, then a declared chemical /
 * hydrodynamic reset, then measure leftovers. Not a task scout.
 * Mid-run epoch switch: FLOW, AHL u, and acid u change at declared
 * windows. SweepS5 / E0.3 Stage-6 upwind when FLOW>0; NO_FLUX when
 * FLOW=0. Stokes +x on bacterium. Claim box 1000x500x10, CENTER AHL
 * (500,250,5), rate 1.28e8, K=1.6, tau_L=1500. No ridge, no NARMA file.
 */
public final class BSimReservoirPlanWashoutReset {
    enum Arm { BROWNIAN, SILENT, DRIVEN }

    static double BOUND_X = 1000.0, BOUND_Y = 500.0, BOUND_Z = 10.0;
    static double DT = .05;
    static int GRID_X = 50, GRID_Y = 25, GRID_Z = 1;
    static int STATE_X = 20, STATE_Y = 10, STATE_Z = 1;
    static int COUNT_X = 4, COUNT_Y = 2, COUNT_Z = 1;
    static double ATT_DIFFUSIVITY = 100, REP_DIFFUSIVITY = 100, AHL_DIFFUSIVITY = 159;
    static double ACID_DIFFUSIVITY = 200;
    static double ATT_DECAY = .0067, REP_DECAY = .033, AHL_DECAY = .0033, ACID_DECAY = .0067;
    static double WARMUP = 18000, WINDOW_DURATION = 300, PULSE_DURATION = 75;
    static double SAMPLING_DURATION = 300, SAMPLING_INTERVAL = 20;
    static double WARMUP_AHL_INPUT = 0.5;
    static int NUM_WINDOWS = 40, INITIAL_POP = 1800, CARRYING_CAPACITY = 2000;
    static boolean HEADLESS = true, WRITE_SAMPLES = true, WRITE_VOXELS = true;
    static long RNG_SEED = 101;
    static Arm ARM = Arm.DRIVEN;
    static String OUTPUT_DIR = "results/wash_decay_w5_seed101";
    static String RUN_LABEL = "wash_decay_w5";
    static String STOCHASTIC_REPLICATE = "driven_seed101";
    static String CONDITION_ID = "WASH_DECAY_W5";
    static double AHL_SOURCE_RATE = 128000000;
    static double ACID_SOURCE_RATE = 2.0e11;
    static double CONFIGURED_AHL_SOURCE_RATE = 128000000;
    static double CONFIGURED_ACID_SOURCE_RATE = 2.0e11;
    static double ACID_CELL_PRODUCTION_RATE = 1e6;
    static double K_MAX_ACID = 0.002;
    static double K_MAX_ALK = 0.003;

    static double FLOW_SPEED = 0.0;
    static double WASH_FLOW_SPEED = 0.0;
    static double TRANSIT_S = Double.POSITIVE_INFINITY;
    static double CHEMICAL_COURANT = 0.0;
    static double WASH_CHEMICAL_COURANT = 0.0;
    static double DT_CFL_MAX_FLOW = Double.POSITIVE_INFINITY;
    static String CHEMICAL_BOUNDARY = "NO_FLUX";
    static int LOAD_LAST_WINDOW = 19;
    static int WASH_LAST_WINDOW = 24;
    static double LOAD_AHL_U = 0.5;
    static double LOAD_ACID_U = 0.5;
    static PrintStream TELEMETRY;
    static BufferedWriter EPOCH_TRACE;
    static final double PROD_RATE = 1e6; // Preserved silent AC0/AC2 architecture.
    static final double GROWTH_RATE = 4.0 * Math.PI / 1800.0;
    static final double EXPECTED_T_GEN = 4.0 * Math.PI / GROWTH_RATE;
    static final double T_REMOVAL = EXPECTED_T_GEN;

    /** Exact frozen Stage 3B receiver; deliberately not configurable. */
    static final double RECEIVER_K_UM = 1.6;
    static final double RECEIVER_HILL_N = 2.0;
    static final double RECEIVER_TAU_S = 15.0;
    static final double RECEIVER_T95_S = 44.936;
    /** Lag, not gain: alpha and delta stay equal. Claim is 1/1500. */
    static double DELTA_LUX = 1.0 / 1500.0;
    static double ALPHA_LUX = 1.0 / 1500.0;
    static double LUM_T95_S = -Math.log(0.05) * 1500.0;
    static final double MOLECULES_PER_UM3_PER_UM = 602.2;
    static final double MM_TO_MOLECULES_PER_UM3 = 602.2 * 1000.0;

    /**
     * Literature-anchored pH map. Buffer capacity is an assumed weakly buffered
     * medium parameter, not a Small/Castanie-Cornet constant.
     */
    static final double PH_BASE = 7.1;
    static final double BUFFER_CAPACITY_MM_PER_PH = 2.0;
    static final double PH_GROWTH_LIMIT = 4.5;
    static final double PH_ACID_EC50 = 3.75;
    static final double PH_ACID_MIN = 2.5;
    static final double N_ACID = 2.0;
    static final double PH_ALK_ONSET = 9.0;
    static final double PH_ALK_EC50 = 9.85;
    static final double PH_ALK_MAX = 10.2;
    static final double N_ALK = 2.0;

    static Vector3d ATT0_POSITION = new Vector3d(250, 250, 5);
    static Vector3d ATT2_POSITION = new Vector3d(750, 250, 5);
    static Vector3d ACID_SOURCE_POSITION = new Vector3d(300, 375, 5);
    static Vector3d AHL_SOURCE_POSITION = new Vector3d(500, 250, 5);
    static double AHL_SOURCE_X = 500, AHL_SOURCE_Y = 250, AHL_SOURCE_Z = 5;
    static double ACID_SOURCE_X = 300, ACID_SOURCE_Y = 375, ACID_SOURCE_Z = 5;
    static double ATT0_X = 250, ATT0_Y = 250, ATT0_Z = 5;
    static double ATT2_X = 750, ATT2_Y = 250, ATT2_Z = 5;
    static double FIELD_DX = 20.0, FIELD_DY = 20.0, DT_CFL_MAX = 0.5;
    static boolean SMOKE_NOT_NARMA = false;

    enum DeathCause { NONE, CLAMP, ACID, OOB }

    static final Vector<ReservoirBacterium> bacteria = new Vector<>();
    static final Vector<ReservoirBacterium> children = new Vector<>();
    static final Vector<ReservoirBacterium> removals = new Vector<>();
    static final Vector<BrownianParticle> particles = new Vector<>();
    static VoxelAnalyzer voxelAnalyzer;
    static Random experimentRng;
    static int nextId, cumulativeDeaths;

    /**
     * Passive null: thermal drift and diffusion only. No sensing, growth,
     * division, death, communication, or luminescence.
     */
    static final class BrownianParticle extends BSimParticle {
        BrownianParticle(BSim sim, Vector3d position) {
            super(sim, position, 1.0);
        }

        @Override
        public void action() {
            super.action();
            addForce(new Vector3d(stokesCoefficient() * FLOW_SPEED, 0, 0));
        }

        static void seedMotionRng(long seed) {
            rng.setSeed(seed);
        }
    }

    public static final class ReservoirBacterium extends BSimBacterium {
        final int id;
        final BSimChemicalField attractantField, repellentField, ahlField, acidField;
        double receiver;
        double luminescence;
        DeathCause deathCause = DeathCause.NONE;

        ReservoirBacterium(BSim sim, Vector3d position,
                           BSimChemicalField attractantField,
                           BSimChemicalField repellentField,
                           BSimChemicalField ahlField,
                           BSimChemicalField acidField) {
            super(sim, position);
            id = nextId++;
            this.attractantField = attractantField;
            this.repellentField = repellentField;
            this.ahlField = ahlField;
            this.acidField = acidField;
            this.receiver = 0.0;
            this.luminescence = 0.0;
            setGoal(attractantField);
        }

        double getResponse() { return receiver; }
        double getQ() { return receiver; }
        double getLuminescence() { return luminescence; }

        void markDeath(DeathCause cause) {
            if (deathCause == DeathCause.NONE) {
                deathCause = cause;
                removals.add(this);
            }
        }

        @Override
        public void action() {
            super.action();
            addForce(new Vector3d(stokesCoefficient() * FLOW_SPEED, 0, 0));

            double concentrationUm = ahlField.getConc(position) / MOLECULES_PER_UM3_PER_UM;
            double concentrationN = Math.pow(concentrationUm, RECEIVER_HILL_N);
            double target = concentrationN
                    / (Math.pow(RECEIVER_K_UM, RECEIVER_HILL_N) + concentrationN);
            receiver = target + (receiver - target) * Math.exp(-sim.getDt() / RECEIVER_TAU_S);
            luminescence += (ALPHA_LUX * receiver - DELTA_LUX * luminescence) * sim.getDt();

            acidField.addQuantity(position, ACID_CELL_PRODUCTION_RATE * sim.getDt());

            double pRemoval = (sim.getDt() / T_REMOVAL)
                    * Math.pow(2.0, -(1.0 - (double) bacteria.size() / CARRYING_CAPACITY));
            boolean clampFires = experimentRng.nextDouble() < pRemoval;
            double pH = pHFromAcidMm(acidMmFromConc(acidField.getConc(position)));
            boolean acidFires = experimentRng.nextDouble() < killRate(pH) * sim.getDt();
            // Competing risks: acid takes priority so clamp cannot mask the death channel.
            if (acidFires) markDeath(DeathCause.ACID);
            else if (clampFires) markDeath(DeathCause.CLAMP);
        }

        @Override
        public void updatePosition() {
            super.updatePosition();
            double boundY = sim.getBound().y;
            if (position.y < 0) position.y = -position.y;
            else if (position.y > boundY) position.y = 2 * boundY - position.y;
            position.y = Math.max(0, Math.min(boundY, position.y));
            position.z = Math.max(0, Math.min(sim.getBound().z, position.z));
        }

        @SuppressWarnings("unchecked")
        @Override
        public void replicate() {
            setRadiusFromSurfaceArea(surfaceArea(replicationRadius) / 2);
            Vector3d childPosition = new Vector3d(position);
            childPosition.x += 2 * radius * (experimentRng.nextDouble() - .5);
            childPosition.y += 2 * radius * (experimentRng.nextDouble() - .5);
            ReservoirBacterium child = new ReservoirBacterium(sim, childPosition,
                    attractantField, repellentField, ahlField, acidField);
            child.setRadius(radius);
            child.setSurfaceAreaGrowthRate(surfaceAreaGrowthRate);
            child.setChildList(childList);
            child.receiver = receiver;
            child.luminescence = luminescence;
            childList.add(child);
            if (voxelAnalyzer != null) voxelAnalyzer.recordBirth(childPosition);
        }
    }

    static final class WindowAccumulator {
        final List<Double> extracellularAhl = new ArrayList<>();
        final List<Double> response = new ArrayList<>();
        final List<Double> pH = new ArrayList<>();
        final List<Double> luminescence = new ArrayList<>();
        double qSum, lumSum, populationSum;
        long cellObservations;
        int qAboveHalf, sampleCount;

        void reset() {
            extracellularAhl.clear();
            response.clear();
            pH.clear();
            luminescence.clear();
            qSum = lumSum = populationSum = 0;
            cellObservations = 0;
            qAboveHalf = sampleCount = 0;
        }

        void observe(BSimChemicalField ahlField, BSimChemicalField acidField) {
            for (ReservoirBacterium bacterium : bacteria) {
                extracellularAhl.add(ahlField.getConc(bacterium.getPosition()));
                response.add(bacterium.getResponse());
                double localAcidMm = acidMmFromConc(acidField.getConc(bacterium.getPosition()));
                pH.add(pHFromAcidMm(localAcidMm));
                luminescence.add(bacterium.getLuminescence());
                lumSum += bacterium.getLuminescence();
                qSum += bacterium.getQ();
                if (bacterium.getQ() > .5) qAboveHalf++;
                cellObservations++;
            }
            populationSum += bacteria.size();
            sampleCount++;
        }

        void observeBrownian() {
            populationSum += particles.size();
            sampleCount++;
        }
    }

    static final class OutputOwner implements AutoCloseable {
        final BufferedWriter summary, samples, voxels;
        int summaryRows, sampleRows, voxelRows;
        private boolean closed;

        OutputOwner(File directory) throws IOException {
            if (!directory.exists() && !directory.mkdirs())
                throw new IOException("Cannot create " + directory);
            summary = new BufferedWriter(new FileWriter(new File(directory, "window_summary.csv")));
            summary.write("Run_Label;Stochastic_Replicate;Arm;Window;Epoch;AHL_Input;Acid_Input;Flow_um_s;"
                    + "Extracellular_AHL_uM_Mean;Extracellular_AHL_uM_P10;"
                    + "Extracellular_AHL_uM_P50;Extracellular_AHL_uM_P90;"
                    + "Mean_q;Fraction_q_gt_0_5;Mean_L;L_Min;L_P10;L_P50;L_P90;L_Max;Lum_Sum;"
                    + "pH_Mean;pH_P10;pH_P50;pH_P90;Population;Births;"
                    + "Total_Deaths;Clamp_Deaths;Input_Driven_Deaths;OOB_Deaths\n");
            samples = WRITE_SAMPLES
                    ? new BufferedWriter(new FileWriter(new File(directory, "results.csv"))) : null;
            if (samples != null) samples.write("Window;Sample;TimeInWindow_s;"
                    + "Input_AC0;Input_AC1_AHL;Input_AC2;Input_Acid;Total_Count;"
                    + "Births_This_Window;Deaths_This_Window;Clamp_Deaths_This_Window;"
                    + "Input_Driven_Deaths_This_Window;OOB_Deaths_This_Window;Mean_L;Lum_Sum\n");
            voxels = WRITE_VOXELS
                    ? new BufferedWriter(new FileWriter(new File(directory, "voxels.csv"))) : null;
            if (voxels != null) writeVoxelHeader(voxels);
            writeFeatureContract(directory);
            writeMatlabMeta(directory);
            writeSourcePositions(directory);
        }

        @Override
        public void close() {
            if (closed) return;
            closed = true;
            try {
                summary.close();
                if (samples != null) samples.close();
                if (voxels != null) voxels.close();
            } catch (IOException e) {
                throw new RuntimeException("Cannot close WashoutReset output", e);
            }
        }
    }

    public static void main(String[] args) {
        String configPath = args.length > 0 ? args[0] : "sim_config_washout.properties";
        loadConfig(configPath);
        if (args.length > 1 && "preview".equals(args[1])) HEADLESS = false;
        assertStep5Locks();
        printTheory();
        printGeometryBanner();
        resetStaticState();
        experimentRng = new Random(RNG_SEED);
        BrownianParticle.seedMotionRng(RNG_SEED);

        final BSim sim = new BSim();
        sim.setDt(DT);
        sim.setSimulationTime(WARMUP + NUM_WINDOWS * WINDOW_DURATION);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);
        int[] fieldGrid = {GRID_X, GRID_Y, GRID_Z};
        final BSimChemicalField attractantField = new BSimChemicalField(
                sim, fieldGrid, ATT_DIFFUSIVITY, ATT_DECAY);
        final BSimChemicalField repellentField = new BSimChemicalField(
                sim, fieldGrid, REP_DIFFUSIVITY, REP_DECAY);
        final BSimChemicalField ahlField = new BSimChemicalField(
                sim, fieldGrid, AHL_DIFFUSIVITY, AHL_DECAY);
        final BSimChemicalField acidField = new BSimChemicalField(
                sim, fieldGrid, ACID_DIFFUSIVITY, ACID_DECAY);
        voxelAnalyzer = new VoxelAnalyzer(new int[]{STATE_X, STATE_Y, STATE_Z},
                new int[]{COUNT_X, COUNT_Y, COUNT_Z}, BOUND_X, BOUND_Y, BOUND_Z);

        for (int i = 0; i < INITIAL_POP; i++) {
            Vector3d p = new Vector3d(
                    0.3 * BOUND_X + experimentRng.nextDouble() * 0.4 * BOUND_X,
                    0.3 * BOUND_Y + experimentRng.nextDouble() * 0.4 * BOUND_Y,
                    BOUND_Z / 2);
            if (ARM == Arm.BROWNIAN) {
                particles.add(new BrownianParticle(sim, p));
            } else {
                ReservoirBacterium bacterium = new ReservoirBacterium(
                        sim, p, attractantField, repellentField, ahlField, acidField);
                bacterium.setRadius();
                bacterium.setSurfaceAreaGrowthRate(GROWTH_RATE);
                bacterium.setChildList(children);
                bacteria.add(bacterium);
            }
        }

        final int warmupSteps = (int) Math.round(WARMUP / DT);
        final int stepsInWindow = (int) Math.round(WINDOW_DURATION / DT);
        final int stepsInPulse = (int) Math.round(PULSE_DURATION / DT);
        final int samplingStart = stepsInWindow - (int) Math.round(SAMPLING_DURATION / DT);
        final int samplingInterval = Math.max(1, (int) Math.round(SAMPLING_INTERVAL / DT));
        final int[] births = {0}, deaths = {0}, clampDeaths = {0};
        final int[] acidDeaths = {0}, oobDeaths = {0}, completed = {0};
        final WindowAccumulator windowAccumulator = new WindowAccumulator();
        final OutputOwner output;
        try {
            output = new OutputOwner(new File(OUTPUT_DIR));
            EPOCH_TRACE = new BufferedWriter(new FileWriter(new File(OUTPUT_DIR, "epoch_trace.txt")));
        } catch (IOException e) {
            throw new RuntimeException("Cannot open WashoutReset output", e);
        }
        applyFlow(0.0);

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                int step = (int) sim.getTimestep();
                boolean warmup = step < warmupSteps;
                int effective = step - warmupSteps;
                int window = warmup ? -1 : effective / stepsInWindow;
                int stepInWindow = warmup ? -1 : effective % stepsInWindow;

                if (!warmup && window < NUM_WINDOWS && stepInWindow == 0) {
                    if (window > 0) {
                        finishWindow(output, window - 1, ahlU(window - 1),
                                acidU(window - 1), births[0], deaths[0], clampDeaths[0],
                                acidDeaths[0], oobDeaths[0], windowAccumulator);
                        completed[0]++;
                    }
                    applyEpoch(window);
                    births[0] = deaths[0] = clampDeaths[0] = acidDeaths[0] = oobDeaths[0] = 0;
                    voxelAnalyzer.resetWindowCounters();
                    windowAccumulator.reset();
                }

                double input0 = 0.0, input2 = 0.0;
                if (warmup && step % stepsInWindow < stepsInPulse) {
                    ahlField.addQuantity(AHL_SOURCE_POSITION,
                            AHL_SOURCE_RATE * WARMUP_AHL_INPUT * DT);
                }
                if (!warmup && window < NUM_WINDOWS && stepInWindow < stepsInPulse) {
                    attractantField.addQuantity(ATT0_POSITION, PROD_RATE * input0 * DT);
                    attractantField.addQuantity(ATT2_POSITION, PROD_RATE * input2 * DT);
                    double ahlInput = ahlU(window);
                    double acidInput = acidU(window);
                    if (ahlInput > 0)
                        ahlField.addQuantity(AHL_SOURCE_POSITION,
                                AHL_SOURCE_RATE * ahlInput * DT);
                    if (acidInput > 0)
                        acidField.addQuantity(ACID_SOURCE_POSITION,
                                ACID_SOURCE_RATE * acidInput * DT);
                }

                attractantField.update();
                repellentField.update();
                ahlField.update();
                acidField.update();
                if (FLOW_SPEED > 0) {
                    advect(attractantField, FLOW_SPEED, DT, GRID_X, GRID_Y, GRID_Z);
                    advect(repellentField, FLOW_SPEED, DT, GRID_X, GRID_Y, GRID_Z);
                    advect(ahlField, FLOW_SPEED, DT, GRID_X, GRID_Y, GRID_Z);
                    advect(acidField, FLOW_SPEED, DT, GRID_X, GRID_Y, GRID_Z);
                }

                if (ARM == Arm.BROWNIAN) {
                    for (BrownianParticle particle : particles) {
                        particle.action();
                        particle.updatePosition();
                    }
                } else {
                    for (ReservoirBacterium bacterium : bacteria) {
                        bacterium.action();
                        bacterium.updatePosition();
                    }
                    int born = children.size();
                    bacteria.addAll(children);
                    children.clear();
                    if (!warmup) births[0] += born;
                    for (ReservoirBacterium bacterium : bacteria) {
                        if (bacterium.getPosition().x < 0
                                || bacterium.getPosition().x > BOUND_X)
                            bacterium.markDeath(DeathCause.OOB);
                    }
                    for (ReservoirBacterium bacterium : removals) {
                        if (!warmup) {
                            voxelAnalyzer.recordDeath(bacterium.getPosition(), bacterium.deathCause);
                            deaths[0]++;
                            if (bacterium.deathCause == DeathCause.CLAMP) clampDeaths[0]++;
                            else if (bacterium.deathCause == DeathCause.ACID) acidDeaths[0]++;
                            else if (bacterium.deathCause == DeathCause.OOB) oobDeaths[0]++;
                        }
                    }
                    cumulativeDeaths += removals.size();
                    bacteria.removeAll(removals);
                    removals.clear();
                }

                if (!warmup && window < NUM_WINDOWS && stepInWindow >= samplingStart) {
                    boolean interval = (stepInWindow - samplingStart) % samplingInterval == 0;
                    boolean last = stepInWindow == stepsInWindow - 1;
                    if (interval || last) {
                        if (ARM == Arm.BROWNIAN) windowAccumulator.observeBrownian();
                        else windowAccumulator.observe(ahlField, acidField);
                        if (WRITE_SAMPLES || WRITE_VOXELS)
                            writeSample(output, window, stepInWindow * DT,
                                    ahlU(window), acidU(window),
                                    births[0], deaths[0], clampDeaths[0],
                                    acidDeaths[0], oobDeaths[0], ahlField, acidField);
                    }
                }
            }
        });

        sim.setDrawer(new BSimP3DDrawer(sim, 800, 600) {
            @Override
            public void scene(PGraphics3D p3d) {
                p3d.ortho(0, (float) BOUND_X, (float) BOUND_Y, 0, -1000, 10000);
                draw(ahlField, Color.ORANGE, (float) (255.0 / 540.0));
                draw(acidField, Color.MAGENTA, (float) (255.0 / 1.0e7));
                if (ARM == Arm.BROWNIAN) {
                    for (BrownianParticle particle : particles)
                        sphere(particle.getPosition(), 8, Color.GRAY, 255);
                } else {
                    for (ReservoirBacterium bacterium : bacteria)
                        sphere(bacterium.getPosition(), 8, Color.GREEN, 255);
                }
            }
        });

        PrintStream originalOut = System.out;
        TELEMETRY = originalOut;
        try {
            if (HEADLESS) {
                System.setOut(new PrintStream(OutputStream.nullOutputStream()));
                sim.export();
            } else {
                sim.preview();
                return;
            }
            while (completed[0] < NUM_WINDOWS) {
                finishWindow(output, completed[0], ahlU(completed[0]),
                        acidU(completed[0]), births[0], deaths[0], clampDeaths[0],
                        acidDeaths[0], oobDeaths[0], windowAccumulator);
                completed[0]++;
            }
        } finally {
            System.setOut(originalOut);
            closeEpochTrace();
            output.close();
        }
        writeRunStatus(output);
        System.out.printf(Locale.US,
                "WashoutReset complete condition=%s arm=%s wash_flow=%.6f final_FLOW=%.6f wash_Courant=%.6f chemical_boundary=%s bound=(%.1f,%.1f,%.1f) ahl=(%.1f,%.1f,%.1f) rate=%.0f summary_rows=%d expected=%d sample_rows=%d voxel_rows=%d%n",
                CONDITION_ID, ARM.name().toLowerCase(Locale.ROOT), WASH_FLOW_SPEED, FLOW_SPEED,
                WASH_CHEMICAL_COURANT, CHEMICAL_BOUNDARY, BOUND_X, BOUND_Y, BOUND_Z,
                AHL_SOURCE_POSITION.x, AHL_SOURCE_POSITION.y, AHL_SOURCE_POSITION.z,
                CONFIGURED_AHL_SOURCE_RATE,
                output.summaryRows, NUM_WINDOWS,
                output.sampleRows, output.voxelRows);
    }

    static double acidMmFromConc(double moleculesPerUm3) {
        return moleculesPerUm3 / MM_TO_MOLECULES_PER_UM3;
    }

    static String epochName(int window) {
        if (window < 0) return "WARMUP";
        if (window <= LOAD_LAST_WINDOW) return "LOAD";
        if (window <= WASH_LAST_WINDOW) return "WASH";
        return "POST";
    }

    static double ahlU(int window) {
        return window <= LOAD_LAST_WINDOW ? LOAD_AHL_U : 0.0;
    }

    static double acidU(int window) {
        return window <= LOAD_LAST_WINDOW ? LOAD_ACID_U : 0.0;
    }

    static double flowForWindow(int window) {
        if (window < 0) return 0.0;
        if (window > LOAD_LAST_WINDOW && window <= WASH_LAST_WINDOW) return WASH_FLOW_SPEED;
        return 0.0;
    }

    static String epochMark(int window) {
        if (window == 0) return "LOAD";
        if (window == LOAD_LAST_WINDOW + 1) return "WASH_START";
        if (window == WASH_LAST_WINDOW) return "WASH_END";
        if (window == NUM_WINDOWS - 1) return "LAST";
        return null;
    }

    static void telemetry(String line) {
        if (TELEMETRY != null) TELEMETRY.println(line);
        if (EPOCH_TRACE != null) {
            try {
                EPOCH_TRACE.write(line);
                EPOCH_TRACE.write('\n');
                EPOCH_TRACE.flush();
            } catch (IOException e) {
                throw new RuntimeException("Cannot write epoch_trace", e);
            }
        }
    }

    static void closeEpochTrace() {
        if (EPOCH_TRACE == null) return;
        try {
            EPOCH_TRACE.close();
        } catch (IOException e) {
            throw new RuntimeException("Cannot close epoch_trace", e);
        }
        EPOCH_TRACE = null;
    }

    static void applyFlow(double speed) {
        if (speed < 0) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset FLOW_SPEED must be >= 0, printed %.6f", speed));
        }
        FLOW_SPEED = speed;
        CHEMICAL_COURANT = FLOW_SPEED * DT / FIELD_DX;
        if (CHEMICAL_COURANT >= 1.0) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset chemical Courant=%.6f >= 1 (FLOW_SPEED=%.6f dt=%.6f dx=%.6f). Do not lower dt.",
                    CHEMICAL_COURANT, FLOW_SPEED, DT, FIELD_DX));
        }
        if (FLOW_SPEED > 0) {
            TRANSIT_S = BOUND_X / FLOW_SPEED;
            DT_CFL_MAX_FLOW = FIELD_DX / FLOW_SPEED;
            CHEMICAL_BOUNDARY = "OUTFLOW";
        } else {
            TRANSIT_S = Double.POSITIVE_INFINITY;
            DT_CFL_MAX_FLOW = Double.POSITIVE_INFINITY;
            CHEMICAL_BOUNDARY = "NO_FLUX";
        }
    }

    static void applyEpoch(int window) {
        applyFlow(flowForWindow(window));
        String epoch = epochName(window);
        telemetry(String.format(Locale.US,
                "WashoutReset WINDOW window=%d epoch=%s ahl_u=%.2f acid_u=%.2f FLOW_SPEED=%.6f chemical_boundary=%s chemical_Courant=%.6f",
                window, epoch, ahlU(window), acidU(window), FLOW_SPEED, CHEMICAL_BOUNDARY, CHEMICAL_COURANT));
        String mark = epochMark(window);
        if (mark != null) {
            telemetry(String.format(Locale.US,
                    "WashoutReset EPOCH_MARK mark=%s window=%d epoch=%s FLOW_SPEED=%.6f chemical_boundary=%s chemical_Courant=%.6f ahl=(%.1f,%.1f,%.1f)",
                    mark, window, epoch, FLOW_SPEED, CHEMICAL_BOUNDARY, CHEMICAL_COURANT,
                    AHL_SOURCE_POSITION.x, AHL_SOURCE_POSITION.y, AHL_SOURCE_POSITION.z));
        }
    }

    static void printTheory() {
        double tauAhl = 1.0 / 0.0033;
        int washWindows = WASH_LAST_WINDOW - LOAD_LAST_WINDOW;
        double dtWash = washWindows * WINDOW_DURATION;
        double ahlLeft = Math.exp(-dtWash / tauAhl);
        double rLeft = Math.exp(-dtWash / RECEIVER_TAU_S);
        double lLeft = Math.exp(-dtWash / 1500.0);
        System.out.printf(Locale.US,
                "WashoutReset THEORY tau_AHL=%.4f s tau_R=%.1f s tau_L=1500 s%n",
                tauAhl, RECEIVER_TAU_S);
        System.out.printf(Locale.US,
                "WashoutReset THEORY W=%d dt=%.1f s AHL_leftover=%.5f R_leftover=%.6g L_leftover=%.6f%n",
                washWindows, dtWash, ahlLeft, rLeft, lLeft);
        System.out.printf(Locale.US,
                "WashoutReset THEORY W=5 ideal AHL=%.5f L=%.6f; W=15 ideal AHL~0 L=%.6f%n",
                Math.exp(-1500.0 / tauAhl), Math.exp(-1.0), Math.exp(-3.0));
        System.out.println("WashoutReset THEORY spatially uniform exponential; living residuals differ (walls, clamp, OOB, occupancy). Do not retune tau_L.");
    }

    static double pHFromAcidMm(double acidMm) {
        return Math.max(1.0, Math.min(14.0, PH_BASE - acidMm / BUFFER_CAPACITY_MM_PER_PH));
    }

    static double killRate(double pH) {
        double acidExcess = Math.max(0.0, PH_GROWTH_LIMIT - pH);
        double acidHalf = PH_GROWTH_LIMIT - PH_ACID_EC50;
        double kAcid = K_MAX_ACID * Math.pow(acidExcess, N_ACID)
                / (Math.pow(acidHalf, N_ACID) + Math.pow(acidExcess, N_ACID));
        double alkExcess = Math.max(0.0, pH - PH_ALK_ONSET);
        double alkHalf = PH_ALK_EC50 - PH_ALK_ONSET;
        double kAlk = K_MAX_ALK * Math.pow(alkExcess, N_ALK)
                / (Math.pow(alkHalf, N_ALK) + Math.pow(alkExcess, N_ALK));
        return kAcid + kAlk;
    }

    static void finishWindow(OutputOwner output, int window, double ahlInput, double acidInput,
                             int births, int totalDeaths, int clampDeaths, int acidDeaths,
                             int oobDeaths, WindowAccumulator accumulator) {
        double[] ext = accumulator.extracellularAhl.stream()
                .mapToDouble(Double::doubleValue).toArray();
        double[] lum = accumulator.luminescence.stream()
                .mapToDouble(Double::doubleValue).toArray();
        double[] pH = accumulator.pH.stream().mapToDouble(Double::doubleValue).toArray();
        Arrays.sort(ext);
        Arrays.sort(lum);
        Arrays.sort(pH);
        int n = lum.length;
        double observations = Math.max(1, accumulator.cellObservations);
        double ahlUmMean = mean(ext) / MOLECULES_PER_UM3_PER_UM;
        int population = ARM == Arm.BROWNIAN
                ? particles.size()
                : (int) Math.round(accumulator.populationSum
                / Math.max(1, accumulator.sampleCount));
        try {
            output.summary.write(String.format(Locale.US,
                    "%s;%s;%s;%d;%s;%.2f;%.2f;%.6f;%.9g;%.9g;%.9g;%.9g;%.9g;%.9g;"
                            + "%.9g;%.9g;%.9g;%.9g;%.9g;%.9g;%.9g;"
                            + "%.9g;%.9g;%.9g;%.9g;%d;%d;%d;%d;%d;%d%n",
                    RUN_LABEL, STOCHASTIC_REPLICATE, ARM.name().toLowerCase(Locale.ROOT),
                    window, epochName(window), ahlInput, acidInput, flowForWindow(window),
                    ARM == Arm.BROWNIAN ? 0.0 : ahlUmMean,
                    ARM == Arm.BROWNIAN ? 0.0 : quantile(ext, .1) / MOLECULES_PER_UM3_PER_UM,
                    ARM == Arm.BROWNIAN ? 0.0 : quantile(ext, .5) / MOLECULES_PER_UM3_PER_UM,
                    ARM == Arm.BROWNIAN ? 0.0 : quantile(ext, .9) / MOLECULES_PER_UM3_PER_UM,
                    ARM == Arm.BROWNIAN ? 0.0 : accumulator.qSum / observations,
                    ARM == Arm.BROWNIAN ? 0.0 : accumulator.qAboveHalf / observations,
                    ARM == Arm.BROWNIAN ? 0.0 : accumulator.lumSum / observations,
                    n == 0 ? 0.0 : lum[0],
                    n == 0 ? 0.0 : quantile(lum, .1),
                    n == 0 ? 0.0 : quantile(lum, .5),
                    n == 0 ? 0.0 : quantile(lum, .9),
                    n == 0 ? 0.0 : lum[n - 1],
                    ARM == Arm.BROWNIAN ? 0.0 : accumulator.lumSum,
                    n == 0 || ARM == Arm.BROWNIAN ? 0.0 : mean(pH),
                    n == 0 || ARM == Arm.BROWNIAN ? 0.0 : quantile(pH, .1),
                    n == 0 || ARM == Arm.BROWNIAN ? 0.0 : quantile(pH, .5),
                    n == 0 || ARM == Arm.BROWNIAN ? 0.0 : quantile(pH, .9),
                    population,
                    births, totalDeaths, clampDeaths, acidDeaths, oobDeaths));
            output.summaryRows++;
        } catch (IOException e) {
            throw new RuntimeException("Cannot write window summary", e);
        }
    }

    static void writeSample(OutputOwner output, int window, double offset,
                            double ahlInput, double acidInput, int births, int deaths,
                            int clampDeaths, int acidDeaths, int oobDeaths,
                            BSimChemicalField ahl, BSimChemicalField acid) {
        try {
            int sample = (int) Math.round(offset / SAMPLING_INTERVAL);
            int count = ARM == Arm.BROWNIAN ? particles.size() : bacteria.size();
            double lumSum = 0;
            if (ARM != Arm.BROWNIAN) {
                for (ReservoirBacterium bacterium : bacteria) lumSum += bacterium.getLuminescence();
            }
            double meanL = (ARM == Arm.BROWNIAN || bacteria.isEmpty())
                    ? 0.0 : lumSum / bacteria.size();
            if (output.samples != null) {
                output.samples.write(String.format(Locale.US,
                        "%d;%d;%.2f;0;%.2f;0;%.2f;%d;%d;%d;%d;%d;%d;%.9g;%.9g%n",
                        window, sample, offset, ahlInput, acidInput, count,
                        births, deaths, clampDeaths, acidDeaths, oobDeaths, meanL, lumSum));
                output.sampleRows++;
            }
            if (output.voxels != null) {
                output.voxels.write(String.format(Locale.US,
                        "%d;%d;%.2f;0;%.2f;0;%.2f", window, sample, offset, ahlInput, acidInput));
                if (ARM == Arm.BROWNIAN) {
                    int[] density = voxelAnalyzer.density(particles);
                    for (int v : density) output.voxels.write(";" + v);
                } else {
                    VoxelAnalyzer.VoxelReadout r = voxelAnalyzer.analyze(ahl, acid, bacteria);
                    for (double v : r.ahlUm) output.voxels.write(format(v));
                    for (double v : r.pH) output.voxels.write(format(v));
                    for (int v : r.density) output.voxels.write(";" + v);
                    for (double v : r.meanResponse) output.voxels.write(format(v));
                    for (double v : r.fractionQAboveHalf) output.voxels.write(format(v));
                    for (double v : r.meanLum) output.voxels.write(format(v));
                    for (double v : r.sumLum) output.voxels.write(format(v));
                    for (int v : r.births) output.voxels.write(";" + v);
                    for (int v : r.clampDeaths) output.voxels.write(";" + v);
                    for (int v : r.acidDeaths) output.voxels.write(";" + v);
                    for (int v : r.oobDeaths) output.voxels.write(";" + v);
                }
                output.voxels.write("\n");
                output.voxelRows++;
            }
        } catch (IOException e) {
            throw new RuntimeException("Cannot write sample", e);
        }
    }

    static void writeVoxelHeader(BufferedWriter w) throws IOException {
        StringBuilder h = new StringBuilder(
                "Window;Sample;TimeInWindow_s;Input_AC0;Input_AC1_AHL;Input_AC2;Input_Acid");
        int states = voxelAnalyzer.getStateVoxels();
        if (ARM == Arm.BROWNIAN) {
            for (int i = 0; i < states; i++) h.append(";Den_").append(i);
        } else {
            for (int i = 0; i < states; i++) h.append(";AHL_uM_").append(i);
            for (int i = 0; i < states; i++) h.append(";pH_").append(i);
            for (int i = 0; i < states; i++) h.append(";Den_").append(i);
            for (int i = 0; i < states; i++) h.append(";Receiver_R_").append(i);
            for (int i = 0; i < states; i++) h.append(";Fraction_q_gt_0_5_").append(i);
            for (int i = 0; i < states; i++) h.append(";Lum_Mean_").append(i);
            for (int i = 0; i < states; i++) h.append(";Lum_Sum_").append(i);
            for (int i = 0; i < voxelAnalyzer.getCountBins(); i++) h.append(";Birth_").append(i);
            for (int i = 0; i < voxelAnalyzer.getCountBins(); i++) h.append(";Clamp_Death_").append(i);
            for (int i = 0; i < voxelAnalyzer.getCountBins(); i++) h.append(";Input_Driven_Death_").append(i);
            for (int i = 0; i < voxelAnalyzer.getCountBins(); i++) h.append(";OOB_Death_").append(i);
        }
        w.write(h.append('\n').toString());
    }

    /** Stage 6 upwind +x. Left inlet is 0. No recycling. */
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

    static String formatTransit(double transit) {
        if (Double.isInfinite(transit)) return "inf";
        return String.format(Locale.US, "%.6f", transit);
    }

    static void assertClaimAhlRate() {
        if (Math.abs(CONFIGURED_AHL_SOURCE_RATE - 128000000.0) > 1e-3) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset AHL rate must be 128000000, printed %.6f", CONFIGURED_AHL_SOURCE_RATE));
        }
    }

    static void assertClaimAcidRate() {
        if (Math.abs(CONFIGURED_ACID_SOURCE_RATE - 2.0e11) > 1.0) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset acid rate must be 2e11, printed %.6g", CONFIGURED_ACID_SOURCE_RATE));
        }
    }

    static void assertAlphaEqualsDelta() {
        if (Math.abs(ALPHA_LUX - DELTA_LUX) > 1e-15) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset alpha=%.12g != delta=%.12g", ALPHA_LUX, DELTA_LUX));
        }
        if (!(DELTA_LUX > 0)) {
            throw new IllegalStateException("STOP WashoutReset delta must be positive");
        }
        if (Math.abs(ALPHA_LUX - 1.0 / 1500.0) > 1e-12
                || Math.abs(DELTA_LUX - 1.0 / 1500.0) > 1e-12) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset alpha=delta must be 1/1500, printed alpha=%.12g delta=%.12g",
                    ALPHA_LUX, DELTA_LUX));
        }
    }

    static void assertAhlPositionMatchesConfig() {
        if (Math.abs(AHL_SOURCE_POSITION.x - AHL_SOURCE_X) > 1e-12
                || Math.abs(AHL_SOURCE_POSITION.y - AHL_SOURCE_Y) > 1e-12
                || Math.abs(AHL_SOURCE_POSITION.z - AHL_SOURCE_Z) > 1e-12) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP printed AHL (%.4f,%.4f,%.4f) != config (%.4f,%.4f,%.4f)",
                    AHL_SOURCE_POSITION.x, AHL_SOURCE_POSITION.y, AHL_SOURCE_POSITION.z,
                    AHL_SOURCE_X, AHL_SOURCE_Y, AHL_SOURCE_Z));
        }
    }

    static void assertReadoutFrozen() {
        if (STATE_X != 20 || STATE_Y != 10 || STATE_Z != 1
                || COUNT_X != 4 || COUNT_Y != 2 || COUNT_Z != 1) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset readout must stay 20x10x1 + 4x2x1, printed %dx%dx%d / %dx%dx%d",
                    STATE_X, STATE_Y, STATE_Z, COUNT_X, COUNT_Y, COUNT_Z));
        }
    }

    static void assertBoundZ() {
        if (Math.abs(BOUND_Z - 10.0) > 1e-12) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset bound.z must be 10, printed %.6f", BOUND_Z));
        }
    }

    static void assertDtFrozen() {
        if (Math.abs(DT - 0.05) > 1e-12) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset dt must be 0.05, printed %.6f", DT));
        }
    }

    static void assertCfl() {
        FIELD_DX = BOUND_X / GRID_X;
        FIELD_DY = BOUND_Y / GRID_Y;
        double minD = Math.min(FIELD_DX, FIELD_DY);
        DT_CFL_MAX = (minD * minD) / (4.0 * 200.0);
        if (DT > DT_CFL_MAX) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset dt=%.6f exceeds CFL dt_max=%.6f (dx=%.6f dy=%.6f). Do not lower dt.",
                    DT, DT_CFL_MAX, FIELD_DX, FIELD_DY));
        }
    }

    static void requireInside(String name, double x, double y, double z) {
        if (x < 0 || x > BOUND_X || y < 0 || y > BOUND_Y || z < 0 || z > BOUND_Z) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset %s=(%.4f, %.4f, %.4f) outside [0,bound]=(%.1f, %.1f, %.1f)",
                    name, x, y, z, BOUND_X, BOUND_Y, BOUND_Z));
        }
    }

    static void assertSourcesInsideBox() {
        requireInside("AHL", AHL_SOURCE_POSITION.x, AHL_SOURCE_POSITION.y, AHL_SOURCE_POSITION.z);
        requireInside("acid", ACID_SOURCE_POSITION.x, ACID_SOURCE_POSITION.y, ACID_SOURCE_POSITION.z);
        requireInside("att0", ATT0_POSITION.x, ATT0_POSITION.y, ATT0_POSITION.z);
        requireInside("att2", ATT2_POSITION.x, ATT2_POSITION.y, ATT2_POSITION.z);
    }

    static void assertClaimBox() {
        if (Math.abs(BOUND_X - 1000.0) > 1e-9
                || Math.abs(BOUND_Y - 500.0) > 1e-9
                || Math.abs(BOUND_Z - 10.0) > 1e-9) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset bounds must be 1000x500x10, printed %.4f x %.4f x %.4f",
                    BOUND_X, BOUND_Y, BOUND_Z));
        }
        if (GRID_X != 50 || GRID_Y != 25 || GRID_Z != 1) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset grid must be 50x25x1, printed %dx%dx%d",
                    GRID_X, GRID_Y, GRID_Z));
        }
        if (Math.abs(AHL_SOURCE_POSITION.x - 500.0) > 1e-9
                || Math.abs(AHL_SOURCE_POSITION.y - 250.0) > 1e-9
                || Math.abs(AHL_SOURCE_POSITION.z - 5.0) > 1e-9) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset AHL must be CENTER (500, 250, 5), printed (%.4f, %.4f, %.4f)",
                    AHL_SOURCE_POSITION.x, AHL_SOURCE_POSITION.y, AHL_SOURCE_POSITION.z));
        }
        if (Math.abs(ACID_SOURCE_POSITION.x - 300.0) > 1e-9
                || Math.abs(ACID_SOURCE_POSITION.y - 375.0) > 1e-9
                || Math.abs(ACID_SOURCE_POSITION.z - 5.0) > 1e-9) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset acid must be (300, 375, 5), printed (%.4f, %.4f, %.4f)",
                    ACID_SOURCE_POSITION.x, ACID_SOURCE_POSITION.y, ACID_SOURCE_POSITION.z));
        }
    }

    static void assertFlowCourant() {
        if (WASH_FLOW_SPEED < 0) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset wash FLOW_SPEED must be >= 0, printed %.6f", WASH_FLOW_SPEED));
        }
        applyFlow(0.0);
        WASH_CHEMICAL_COURANT = WASH_FLOW_SPEED * DT / FIELD_DX;
        if (WASH_CHEMICAL_COURANT >= 1.0) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset wash chemical Courant=%.6f >= 1 (wash_flow=%.6f dt=%.6f dx=%.6f). Do not lower dt.",
                    WASH_CHEMICAL_COURANT, WASH_FLOW_SPEED, DT, FIELD_DX));
        }
        boolean flowOk = Math.abs(WASH_FLOW_SPEED) < 1e-9
                || Math.abs(WASH_FLOW_SPEED - 0.25) < 1e-6
                || Math.abs(WASH_FLOW_SPEED - 8.0) < 1e-6;
        if (SMOKE_NOT_NARMA) {
            if (Math.abs(WASH_FLOW_SPEED - 8.0) > 1e-6) {
                throw new IllegalStateException(String.format(Locale.US,
                        "STOP WashoutReset smoke wash flow must be 8, printed %.6f", WASH_FLOW_SPEED));
            }
        } else if (!flowOk) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset wash FLOW_SPEED=%.6f is not 0 / 0.25 / 8", WASH_FLOW_SPEED));
        }
        if (LOAD_LAST_WINDOW < 0 || WASH_LAST_WINDOW < LOAD_LAST_WINDOW
                || WASH_LAST_WINDOW >= NUM_WINDOWS) {
            throw new IllegalStateException(String.format(Locale.US,
                    "STOP WashoutReset epoch windows load_last=%d wash_last=%d num.windows=%d",
                    LOAD_LAST_WINDOW, WASH_LAST_WINDOW, NUM_WINDOWS));
        }
        if (ARM != Arm.DRIVEN) {
            throw new IllegalStateException("STOP WashoutReset production/smoke arm must be driven");
        }
    }

    static void assertStep5Locks() {
        assertAhlPositionMatchesConfig();
        assertClaimAhlRate();
        assertClaimAcidRate();
        assertAlphaEqualsDelta();
        assertReadoutFrozen();
        assertBoundZ();
        assertDtFrozen();
        assertCfl();
        assertSourcesInsideBox();
        assertClaimBox();
        assertFlowCourant();
    }

    static void printGeometryBanner() {
        System.out.printf(Locale.US,
                "WashoutReset bound.x=%.6f bound.y=%.6f bound.z=%.6f%n", BOUND_X, BOUND_Y, BOUND_Z);
        System.out.printf(Locale.US,
                "WashoutReset grid.x=%d grid.y=%d grid.z=%d dx=%.6f dy=%.6f CFL_dt_max=%.6f dt=%.6f%n",
                GRID_X, GRID_Y, GRID_Z, FIELD_DX, FIELD_DY, DT_CFL_MAX, DT);
        System.out.printf(Locale.US,
                "WashoutReset field.ahl.source.x=%.6f field.ahl.source.y=%.6f field.ahl.source.z=%.6f rate=%.0f%n",
                AHL_SOURCE_POSITION.x, AHL_SOURCE_POSITION.y, AHL_SOURCE_POSITION.z, AHL_SOURCE_RATE);
        System.out.printf(Locale.US,
                "WashoutReset field.acid.source.x=%.6f field.acid.source.y=%.6f field.acid.source.z=%.6f rate=%.6g%n",
                ACID_SOURCE_POSITION.x, ACID_SOURCE_POSITION.y, ACID_SOURCE_POSITION.z, ACID_SOURCE_RATE);
        System.out.printf(Locale.US,
                "WashoutReset field.att0.x=%.6f field.att0.y=%.6f field.att0.z=%.6f%n",
                ATT0_POSITION.x, ATT0_POSITION.y, ATT0_POSITION.z);
        System.out.printf(Locale.US,
                "WashoutReset field.att2.x=%.6f field.att2.y=%.6f field.att2.z=%.6f%n",
                ATT2_POSITION.x, ATT2_POSITION.y, ATT2_POSITION.z);
        System.out.printf(Locale.US,
                "WashoutReset initial.pop=%d carrying.capacity=%d luminescence.alpha=%.12g luminescence.delta=%.12g%n",
                INITIAL_POP, CARRYING_CAPACITY, ALPHA_LUX, DELTA_LUX);
        System.out.printf(Locale.US,
                "WashoutReset condition=%s LOAD=0..%d WASH=%d..%d POST=%d..%d wash_flow=%.6f wash_Courant=%.6f%n",
                CONDITION_ID, LOAD_LAST_WINDOW, LOAD_LAST_WINDOW + 1, WASH_LAST_WINDOW,
                WASH_LAST_WINDOW + 1, NUM_WINDOWS - 1, WASH_FLOW_SPEED, WASH_CHEMICAL_COURANT);
        System.out.printf(Locale.US,
                "WashoutReset start FLOW_SPEED=%.6f transit_s=%s chemical_Courant=%.6f CFL_dt_max_flow=%s chemical_boundary=%s%n",
                FLOW_SPEED, formatTransit(TRANSIT_S), CHEMICAL_COURANT, formatTransit(DT_CFL_MAX_FLOW),
                CHEMICAL_BOUNDARY);
        System.out.println("WashoutReset advect=AHL+acid+att+rep when FLOW>0 else off; Stokes=bacterium; OUTFLOW iff FLOW>0");
        System.out.printf(Locale.US,
                "WashoutReset expected Courant at 8 um/s=0.02 (dx=20); at 0.25 um/s=0.000625%n");
    }

    static void writeSourcePositions(File directory) throws IOException {
        assertStep5Locks();
        try (BufferedWriter w = new BufferedWriter(
                new FileWriter(new File(directory, "source_positions.txt")))) {
            w.write(String.format(Locale.US, "bound.x=%.6f%n", BOUND_X));
            w.write(String.format(Locale.US, "bound.y=%.6f%n", BOUND_Y));
            w.write(String.format(Locale.US, "bound.z=%.6f%n", BOUND_Z));
            w.write(String.format(Locale.US, "grid.x=%d%n", GRID_X));
            w.write(String.format(Locale.US, "grid.y=%d%n", GRID_Y));
            w.write(String.format(Locale.US, "grid.z=%d%n", GRID_Z));
            w.write(String.format(Locale.US, "dx=%.6f%n", FIELD_DX));
            w.write(String.format(Locale.US, "dy=%.6f%n", FIELD_DY));
            w.write(String.format(Locale.US, "CFL_dt_max=%.6f%n", DT_CFL_MAX));
            w.write(String.format(Locale.US, "dt=%.6f%n", DT));
            w.write(String.format(Locale.US, "field.ahl.source.x=%.6f%n", AHL_SOURCE_POSITION.x));
            w.write(String.format(Locale.US, "field.ahl.source.y=%.6f%n", AHL_SOURCE_POSITION.y));
            w.write(String.format(Locale.US, "field.ahl.source.z=%.6f%n", AHL_SOURCE_POSITION.z));
            w.write(String.format(Locale.US, "field.ahl.source.rate=%.6f%n", CONFIGURED_AHL_SOURCE_RATE));
            w.write(String.format(Locale.US, "field.acid.source.x=%.6f%n", ACID_SOURCE_POSITION.x));
            w.write(String.format(Locale.US, "field.acid.source.y=%.6f%n", ACID_SOURCE_POSITION.y));
            w.write(String.format(Locale.US, "field.acid.source.z=%.6f%n", ACID_SOURCE_POSITION.z));
            w.write(String.format(Locale.US, "field.acid.source.rate=%.6g%n", CONFIGURED_ACID_SOURCE_RATE));
            w.write(String.format(Locale.US, "field.att0.x=%.6f%n", ATT0_POSITION.x));
            w.write(String.format(Locale.US, "field.att0.y=%.6f%n", ATT0_POSITION.y));
            w.write(String.format(Locale.US, "field.att0.z=%.6f%n", ATT0_POSITION.z));
            w.write(String.format(Locale.US, "field.att2.x=%.6f%n", ATT2_POSITION.x));
            w.write(String.format(Locale.US, "field.att2.y=%.6f%n", ATT2_POSITION.y));
            w.write(String.format(Locale.US, "field.att2.z=%.6f%n", ATT2_POSITION.z));
            w.write(String.format(Locale.US, "initial.pop=%d%n", INITIAL_POP));
            w.write(String.format(Locale.US, "carrying.capacity=%d%n", CARRYING_CAPACITY));
            w.write(String.format(Locale.US, "luminescence.alpha=%.12g%n", ALPHA_LUX));
            w.write(String.format(Locale.US, "luminescence.delta=%.12g%n", DELTA_LUX));
            w.write(String.format(Locale.US, "condition.id=%s%n", CONDITION_ID));
            w.write(String.format(Locale.US, "epoch.load.last.window=%d%n", LOAD_LAST_WINDOW));
            w.write(String.format(Locale.US, "epoch.wash.last.window=%d%n", WASH_LAST_WINDOW));
            w.write(String.format(Locale.US, "load.ahl.u=%.2f%n", LOAD_AHL_U));
            w.write(String.format(Locale.US, "load.acid.u=%.2f%n", LOAD_ACID_U));
            w.write(String.format(Locale.US, "wash.flow.speed.um_s=%.6f%n", WASH_FLOW_SPEED));
            w.write(String.format(Locale.US, "wash_chemical_Courant=%.6f%n", WASH_CHEMICAL_COURANT));
            w.write(String.format(Locale.US, "start.flow.speed.um_s=%.6f%n", FLOW_SPEED));
            w.write(String.format(Locale.US, "transit_s=%s%n", formatTransit(TRANSIT_S)));
            w.write(String.format(Locale.US, "chemical_Courant=%.6f%n", CHEMICAL_COURANT));
            w.write(String.format(Locale.US, "CFL_dt_max_flow=%s%n", formatTransit(DT_CFL_MAX_FLOW)));
            w.write(String.format(Locale.US, "chemical_boundary=%s%n", CHEMICAL_BOUNDARY));
            w.write("advect=ahl,acid,att,rep when FLOW>0 else off\n");
            w.write("stokes=bacterium\n");
            w.write("note=WashoutReset: occupied load then declared wash; claim box 1000x500x10; CENTER AHL (500,250,5) rate 1.28e8; K=1.6; tau_L=1500; SweepS5 Stage 6 upwind when FLOW>0; no recycling; no ridge\n");
        }
        if (SMOKE_NOT_NARMA) {
            try (BufferedWriter w = new BufferedWriter(
                    new FileWriter(new File(directory, "SMOKE_NOT_EVIDENCE.txt")))) {
                w.write("This directory is WashoutReset smoke only. Do not score as residual evidence.\n");
            }
        }
    }

    static void writeMatlabMeta(File directory) throws IOException {
        int nxy = STATE_X * STATE_Y;
        int nCount = COUNT_X * COUNT_Y;
        try (BufferedWriter w = new BufferedWriter(
                new FileWriter(new File(directory, "matlab_meta.txt")))) {
            w.write("delimiter=;\n");
            w.write("samplesPerWindow=16\n");
            w.write("stateGridX=" + STATE_X + "\n");
            w.write("stateGridY=" + STATE_Y + "\n");
            w.write("countGridX=" + COUNT_X + "\n");
            w.write("countGridY=" + COUNT_Y + "\n");
            w.write("nState=" + nxy + "\n");
            w.write("nCount=" + nCount + "\n");
            w.write("prefixCols=Window;Sample;TimeInWindow_s;Input_AC0;Input_AC1_AHL;Input_AC2;Input_Acid\n");
            w.write("prefixCount=7\n");
            if (ARM == Arm.BROWNIAN) {
                w.write("channels=Den\n");
                w.write(String.format(Locale.US,
                        "matlabReshape=Den=reshape(row(8:%d),[%d,%d]);%n",
                        7 + nxy, STATE_X, STATE_Y));
            } else {
                w.write("channels=AHL_uM,pH,Den,Receiver_R,Fraction_q_gt_0_5,Lum_Mean,Lum_Sum,Birth,Clamp_Death,Input_Driven_Death,OOB_Death\n");
                int c = 8;
                StringBuilder reshape = new StringBuilder("matlabReshape=");
                String[] stateNames = {
                        "AHL", "pH", "Den", "R", "q", "L", "Lsum"
                };
                for (String name : stateNames) {
                    int last = c + nxy - 1;
                    reshape.append(String.format(Locale.US,
                            "%s=reshape(row(%d:%d),[%d,%d]); ",
                            name, c, last, STATE_X, STATE_Y));
                    c = last + 1;
                }
                String[] countNames = {"Birth", "Clamp", "Death", "OOB"};
                for (int i = 0; i < countNames.length; i++) {
                    int last = c + nCount - 1;
                    reshape.append(String.format(Locale.US,
                            "%s=reshape(row(%d:%d),[%d,%d])",
                            countNames[i], c, last, COUNT_X, COUNT_Y));
                    if (i < countNames.length - 1) reshape.append("; ");
                    c = last + 1;
                }
                w.write(reshape.append('\n').toString());
                w.write("biologyReadout=last-sample dish-mean AHL/R/L; no ridge\n");
            }
            w.write("note=WashoutReset residual scout. Do not score NARMA.\n");
        }
    }

    static void writeFeatureContract(File directory) throws IOException {
        try (BufferedWriter w = new BufferedWriter(
                new FileWriter(new File(directory, "feature_contract.txt")))) {
            w.write("arm=" + ARM.name().toLowerCase(Locale.ROOT) + "\n");
            w.write("ridge=off\n");
            w.write("condition.id=" + CONDITION_ID + "\n");
            w.write("note=WashoutReset residual scout. Do not score NARMA, MG NRMSE, waveform AUC, or ridge.\n");
        }
    }

    static void writeRunStatus(OutputOwner output) {
        File directory = new File(OUTPUT_DIR);
        try (BufferedWriter w = new BufferedWriter(
                new FileWriter(new File(directory, "run_status.txt")))) {
            w.write(String.format(Locale.US,
                    "condition.id=%s%narm=%s%nseed=%d%nbound.x=%.6f%nbound.y=%.6f%nbound.z=%.6f%ngrid.x=%d%ngrid.y=%d%ngrid.z=%d%ndx=%.6f%ndy=%.6f%nCFL_dt_max=%.6f%nahl.x=%.6f%nahl.y=%.6f%nahl.z=%.6f%nahl.rate=%.6f%nacid.x=%.6f%nacid.y=%.6f%nacid.z=%.6f%natt0.x=%.6f%natt0.y=%.6f%natt0.z=%.6f%natt2.x=%.6f%natt2.y=%.6f%natt2.z=%.6f%ninitial.pop=%d%ncarrying.capacity=%d%nalpha=%.12g%ndelta=%.12g%nepoch.load.last.window=%d%nepoch.wash.last.window=%d%nwash.flow.speed.um_s=%.6f%nwash_chemical_Courant=%.6f%nfinal.flow.speed.um_s=%.6f%ntransit_s=%s%nchemical_Courant=%.6f%nCFL_dt_max_flow=%s%nchemical_boundary=%s%nsummary_rows=%d%nsample_rows=%d%nvoxel_rows=%d%n",
                    CONDITION_ID, ARM.name().toLowerCase(Locale.ROOT), RNG_SEED,
                    BOUND_X, BOUND_Y, BOUND_Z, GRID_X, GRID_Y, GRID_Z,
                    FIELD_DX, FIELD_DY, DT_CFL_MAX,
                    AHL_SOURCE_POSITION.x, AHL_SOURCE_POSITION.y, AHL_SOURCE_POSITION.z,
                    AHL_SOURCE_RATE,
                    ACID_SOURCE_POSITION.x, ACID_SOURCE_POSITION.y, ACID_SOURCE_POSITION.z,
                    ATT0_POSITION.x, ATT0_POSITION.y, ATT0_POSITION.z,
                    ATT2_POSITION.x, ATT2_POSITION.y, ATT2_POSITION.z,
                    INITIAL_POP, CARRYING_CAPACITY, ALPHA_LUX, DELTA_LUX,
                    LOAD_LAST_WINDOW, WASH_LAST_WINDOW, WASH_FLOW_SPEED, WASH_CHEMICAL_COURANT,
                    FLOW_SPEED, formatTransit(TRANSIT_S), CHEMICAL_COURANT,
                    formatTransit(DT_CFL_MAX_FLOW), CHEMICAL_BOUNDARY,
                    output.summaryRows, output.sampleRows, output.voxelRows));
        } catch (IOException e) {
            throw new RuntimeException("Cannot write run status", e);
        }
    }

    static double[] readGradedSequence(String path, String name) {
        List<Double> values = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty() && !line.startsWith("#"))
                    values.add(Double.parseDouble(line));
            }
        } catch (IOException | NumberFormatException e) {
            throw new IllegalArgumentException("Cannot read graded " + name + " input " + path, e);
        }
        double[] result = new double[values.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = values.get(i);
            if (result[i] < 0 || result[i] > 1)
                throw new IllegalArgumentException(name + " inputs must be in [0,1]");
        }
        return result;
    }

    static void loadConfig(String path) {
        Properties p = loadProperties(new File(path));
        DT = getDouble(p, "dt", DT);
        BOUND_X = getDouble(p, "bound.x", BOUND_X);
        BOUND_Y = getDouble(p, "bound.y", BOUND_Y);
        BOUND_Z = getDouble(p, "bound.z", BOUND_Z);
        GRID_X = getInt(p, "grid.x", GRID_X); GRID_Y = getInt(p, "grid.y", GRID_Y);
        GRID_Z = getInt(p, "grid.z", GRID_Z);
        STATE_X = getInt(p, "readout.grid.x", STATE_X);
        STATE_Y = getInt(p, "readout.grid.y", STATE_Y);
        STATE_Z = getInt(p, "readout.grid.z", STATE_Z);
        COUNT_X = getInt(p, "readout.countgrid.x", COUNT_X);
        COUNT_Y = getInt(p, "readout.countgrid.y", COUNT_Y);
        COUNT_Z = getInt(p, "readout.countgrid.z", COUNT_Z);
        ATT_DIFFUSIVITY = getDouble(p, "field.att.diff", ATT_DIFFUSIVITY);
        REP_DIFFUSIVITY = getDouble(p, "field.rep.diff", REP_DIFFUSIVITY);
        AHL_DIFFUSIVITY = getDouble(p, "field.ahl.diff", AHL_DIFFUSIVITY);
        ACID_DIFFUSIVITY = getDouble(p, "field.acid.diff", ACID_DIFFUSIVITY);
        ATT_DECAY = getDouble(p, "field.att.decay", ATT_DECAY);
        REP_DECAY = getDouble(p, "field.rep.decay", REP_DECAY);
        AHL_DECAY = getDouble(p, "field.ahl.decay", AHL_DECAY);
        ACID_DECAY = getDouble(p, "field.acid.decay", ACID_DECAY);
        AHL_SOURCE_RATE = getDouble(p, "field.ahl.source.rate", AHL_SOURCE_RATE);
        AHL_SOURCE_X = getDouble(p, "field.ahl.source.x", AHL_SOURCE_X);
        AHL_SOURCE_Y = getDouble(p, "field.ahl.source.y", AHL_SOURCE_Y);
        AHL_SOURCE_Z = getDouble(p, "field.ahl.source.z", AHL_SOURCE_Z);
        AHL_SOURCE_POSITION = new Vector3d(AHL_SOURCE_X, AHL_SOURCE_Y, AHL_SOURCE_Z);
        ACID_SOURCE_X = getDouble(p, "field.acid.source.x", ACID_SOURCE_X);
        ACID_SOURCE_Y = getDouble(p, "field.acid.source.y", ACID_SOURCE_Y);
        ACID_SOURCE_Z = getDouble(p, "field.acid.source.z", ACID_SOURCE_Z);
        ACID_SOURCE_POSITION = new Vector3d(ACID_SOURCE_X, ACID_SOURCE_Y, ACID_SOURCE_Z);
        ATT0_X = getDouble(p, "field.att0.x", ATT0_X);
        ATT0_Y = getDouble(p, "field.att0.y", ATT0_Y);
        ATT0_Z = getDouble(p, "field.att0.z", ATT0_Z);
        ATT0_POSITION = new Vector3d(ATT0_X, ATT0_Y, ATT0_Z);
        ATT2_X = getDouble(p, "field.att2.x", ATT2_X);
        ATT2_Y = getDouble(p, "field.att2.y", ATT2_Y);
        ATT2_Z = getDouble(p, "field.att2.z", ATT2_Z);
        ATT2_POSITION = new Vector3d(ATT2_X, ATT2_Y, ATT2_Z);
        ACID_SOURCE_RATE = getDouble(p, "field.acid.source.rate", ACID_SOURCE_RATE);
        ACID_CELL_PRODUCTION_RATE = getDouble(p,
                "field.acid.cell.production.rate", ACID_CELL_PRODUCTION_RATE);
        K_MAX_ACID = getDouble(p, "field.acid.kmax", K_MAX_ACID);
        K_MAX_ALK = getDouble(p, "field.acid.kmax.alk", K_MAX_ALK);
        ALPHA_LUX = getDouble(p, "luminescence.alpha", ALPHA_LUX);
        DELTA_LUX = getDouble(p, "luminescence.delta", DELTA_LUX);
        LUM_T95_S = -Math.log(0.05) / DELTA_LUX;
        assertAlphaEqualsDelta();
        WARMUP_AHL_INPUT = getDouble(p, "warmup.ahl.input", WARMUP_AHL_INPUT);
        if (WARMUP_AHL_INPUT < 0 || WARMUP_AHL_INPUT > 1)
            throw new IllegalArgumentException("warmup.ahl.input must be in [0,1]");
        WARMUP = getDouble(p, "warmup.s", WARMUP);
        WINDOW_DURATION = getDouble(p, "window.duration.s", WINDOW_DURATION);
        PULSE_DURATION = Math.min(WINDOW_DURATION, getDouble(p, "pulse.duration.s", PULSE_DURATION));
        SAMPLING_DURATION = Math.min(WINDOW_DURATION,
                getDouble(p, "sampling.duration.s", SAMPLING_DURATION));
        SAMPLING_INTERVAL = getDouble(p, "sampling.interval.s", SAMPLING_INTERVAL);
        NUM_WINDOWS = getInt(p, "num.windows", NUM_WINDOWS);
        INITIAL_POP = getInt(p, "initial.pop", INITIAL_POP);
        CARRYING_CAPACITY = getInt(p, "carrying.capacity", CARRYING_CAPACITY);
        HEADLESS = Boolean.parseBoolean(p.getProperty("headless", String.valueOf(HEADLESS)));
        WRITE_SAMPLES = Boolean.parseBoolean(
                p.getProperty("output.write.samples", String.valueOf(WRITE_SAMPLES)));
        WRITE_VOXELS = Boolean.parseBoolean(
                p.getProperty("output.write.voxels", String.valueOf(WRITE_VOXELS)));
        RNG_SEED = Long.parseLong(p.getProperty("rng.seed", String.valueOf(RNG_SEED)));
        OUTPUT_DIR = p.getProperty("output.dir", OUTPUT_DIR);
        RUN_LABEL = p.getProperty("output.run.label", RUN_LABEL);
        STOCHASTIC_REPLICATE = p.getProperty(
                "output.stochastic.replicate", STOCHASTIC_REPLICATE);
        CONDITION_ID = p.getProperty("condition.id", CONDITION_ID);
        String armName = p.getProperty("arm", "driven").trim().toLowerCase(Locale.ROOT);
        if ("driven".equals(armName)) ARM = Arm.DRIVEN;
        else throw new IllegalArgumentException("WashoutReset arm must be driven");
        SMOKE_NOT_NARMA = Boolean.parseBoolean(p.getProperty("sweep.smoke.not.narma", "false"))
                || Boolean.parseBoolean(p.getProperty("washout.smoke.not.evidence", "false"));
        WASH_FLOW_SPEED = getDouble(p, "wash.flow.speed.um_s", WASH_FLOW_SPEED);
        LOAD_LAST_WINDOW = getInt(p, "epoch.load.last.window", LOAD_LAST_WINDOW);
        WASH_LAST_WINDOW = getInt(p, "epoch.wash.last.window", WASH_LAST_WINDOW);
        LOAD_AHL_U = getDouble(p, "load.ahl.u", LOAD_AHL_U);
        LOAD_ACID_U = getDouble(p, "load.acid.u", LOAD_ACID_U);
        if (LOAD_AHL_U < 0 || LOAD_AHL_U > 1 || LOAD_ACID_U < 0 || LOAD_ACID_U > 1)
            throw new IllegalArgumentException("load AHL/acid u must be in [0,1]");
        CONFIGURED_AHL_SOURCE_RATE = AHL_SOURCE_RATE;
        CONFIGURED_ACID_SOURCE_RATE = ACID_SOURCE_RATE;
        FLOW_SPEED = 0.0;
        CHEMICAL_BOUNDARY = "NO_FLUX";
    }

    static Properties loadProperties(File file) {
        Properties local = new Properties();
        try (FileInputStream in = new FileInputStream(file)) {
            local.load(in);
        } catch (IOException e) {
            throw new IllegalArgumentException("Cannot load config " + file, e);
        }
        String include = local.getProperty("config.include");
        if (include == null) return local;
        File parent = file.getAbsoluteFile().getParentFile();
        Properties merged = loadProperties(new File(parent, include));
        merged.putAll(local);
        return merged;
    }

    static int getInt(Properties p, String k, int d) {
        return Integer.parseInt(p.getProperty(k, String.valueOf(d)));
    }
    static double getDouble(Properties p, String k, double d) {
        return Double.parseDouble(p.getProperty(k, String.valueOf(d)));
    }
    static String format(double v) { return String.format(Locale.US, ";%.4e", v); }
    static double mean(double[] a) {
        if (a.length == 0) return Double.NaN;
        double sum = 0; for (double v : a) sum += v; return sum / a.length;
    }
    static double quantile(double[] sorted, double p) {
        if (sorted.length == 0) return Double.NaN;
        double index = p * (sorted.length - 1);
        int lo = (int) Math.floor(index), hi = (int) Math.ceil(index);
        return sorted[lo] + (index - lo) * (sorted[hi] - sorted[lo]);
    }
    static void resetStaticState() {
        bacteria.clear(); children.clear(); removals.clear(); particles.clear();
        voxelAnalyzer = null; nextId = 0; cumulativeDeaths = 0;
    }
}
