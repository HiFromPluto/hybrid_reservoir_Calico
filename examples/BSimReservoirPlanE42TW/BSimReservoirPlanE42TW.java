package BSimReservoirPlanE42TW;

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
 * E4.2 hypothetical two-way mediator on the frozen Narma10b / A1 dish.
 * MODEL_STATUS = HYPOTHETICAL_DESIGN_ENVELOPE. Not a digital twin.
 * P is an inert reporter, not AHL, not acid, not L. No vesicle store.
 * Kinetics, layout, K_AC, Jmax_A1, and FLOW=0 are frozen.
 */
public final class BSimReservoirPlanE42TW {
    enum Arm { BROWNIAN, SILENT, DRIVEN }
    enum Coupling { CLOSED, REPLAY, SHUFFLE }

    static final double BOUND_X = 1000.0, BOUND_Y = 500.0, BOUND_Z = 10.0;
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
    static int NUM_WINDOWS = 200, INITIAL_POP = 1800, CARRYING_CAPACITY = 2000;
    static boolean HEADLESS = true, WRITE_SAMPLES = true, WRITE_VOXELS = true;
    static long RNG_SEED = 111;
    static Arm ARM = Arm.DRIVEN;
    static String OUTPUT_DIR = "results/e42_closed_seed111";
    static String INPUT_AHL_FILE = "input_ahl_narma200.txt";
    static String INPUT_ACID_FILE = "input_acid_held05_200.txt";
    static String INPUT_P_REPLAY_FILE = "";
    static String RUN_LABEL = "e42_closed";
    static String STOCHASTIC_REPLICATE = "closed_seed111";
    static Coupling COUPLING = Coupling.CLOSED;
    /** A0 HybridDish current; used only for payload-match print. Not the A1 source. */
    static final double JMAX_A0 = 128000000.0;
    static final String MODEL_STATUS = "HYPOTHETICAL_DESIGN_ENVELOPE";
    static final double AC_G0 = 0.0;
    static final double AC_N = 2.0;
    static final double AC_K = 0.25;
    /** Payload-matched A1 Jmax. Frozen before occupancy/NRMSE. Not a wet-lab flux. */
    static final double JMAX_A1_FROZEN = 74677509.75821304;
    /** Inert reporter P. AHL-like transport. Not acid, not a second QS HSL, not L. */
    static final double P_DIFFUSIVITY = 159.0;
    static final double P_DECAY = 0.0033;
    static final double AC_N_P = 2.0;
    static final double AC_K_P = 1.6;
    static final double AC_ALPHA = 0.5;
    /** Well-mixed freeze: 1800 cells at R=0.1915 give P≈K_P. Not an NRMSE fit. */
    static final double K_P_RATE = 46121.49695387294;
    static double AHL_SOURCE_RATE = JMAX_A1_FROZEN;
    static double ACID_SOURCE_RATE = 2.0e11;
    static double ACID_CELL_PRODUCTION_RATE = 1e6;
    static double K_MAX_ACID = 0.002;
    static double K_MAX_ALK = 0.003;

    static final double FLOW_SPEED = 0.0;
    static final double PROD_RATE = 1e6; // Preserved silent AC0/AC2 architecture.
    static final double GROWTH_RATE = 4.0 * Math.PI / 1800.0;
    static final double EXPECTED_T_GEN = 4.0 * Math.PI / GROWTH_RATE;
    static final double T_REMOVAL = EXPECTED_T_GEN;

    /** Exact frozen Stage 3B receiver; deliberately not configurable. */
    static final double RECEIVER_K_UM = 1.6;
    static final double RECEIVER_HILL_N = 2.0;
    static final double RECEIVER_TAU_S = 15.0;
    static final double RECEIVER_T95_S = 44.936;
    static final double DELTA_LUX = 1.0 / 1500.0;
    static final double LUM_T95_S = -Math.log(0.05) * 1500.0;
    static double ALPHA_LUX = 1.0 / 1500.0;
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

    static final Vector3d[] AC_POSITIONS = {
            new Vector3d(250, 250, 5), new Vector3d(500, 250, 5),
            new Vector3d(750, 250, 5)
    };
    /** Dedicated death-input source; not AC0/AC1/AC2. */
    static final Vector3d ACID_SOURCE_POSITION = new Vector3d(300, 375, 5);

    enum DeathCause { NONE, CLAMP, ACID, OOB }

    static final Vector<ReservoirBacterium> bacteria = new Vector<>();
    static final Vector<ReservoirBacterium> children = new Vector<>();
    static final Vector<ReservoirBacterium> removals = new Vector<>();
    static final Vector<BrownianParticle> particles = new Vector<>();
    static VoxelAnalyzer voxelAnalyzer;
    static Random experimentRng;
    static int nextId, cumulativeDeaths;
    static boolean ahlFinite = true;
    static boolean ahlNonNegative = true;
    static boolean pFinite = true;
    static boolean pNonNegative = true;
    static double[] replayTime;
    static double[] replayP;
    static int replayCursor;

    /**
     * Passive null: thermal drift and diffusion only. No sensing, growth,
     * division, death, communication, or luminescence.
     */
    static final class BrownianParticle extends BSimParticle {
        BrownianParticle(BSim sim, Vector3d position) {
            super(sim, position, 1.0);
        }

        static void seedMotionRng(long seed) {
            rng.setSeed(seed);
        }
    }

    public static final class ReservoirBacterium extends BSimBacterium {
        final int id;
        final BSimChemicalField attractantField, repellentField, ahlField, acidField, reporterField;
        double receiver;
        double luminescence;
        DeathCause deathCause = DeathCause.NONE;

        ReservoirBacterium(BSim sim, Vector3d position,
                           BSimChemicalField attractantField,
                           BSimChemicalField repellentField,
                           BSimChemicalField ahlField,
                           BSimChemicalField acidField,
                           BSimChemicalField reporterField) {
            super(sim, position);
            id = nextId++;
            this.attractantField = attractantField;
            this.repellentField = repellentField;
            this.ahlField = ahlField;
            this.acidField = acidField;
            this.reporterField = reporterField;
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
            reporterField.addQuantity(position, K_P_RATE * receiver * sim.getDt());

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
                    attractantField, repellentField, ahlField, acidField, reporterField);
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
        final BufferedWriter summary, samples, voxels, pLog;
        int summaryRows, sampleRows, voxelRows;
        private boolean closed;

        OutputOwner(File directory) throws IOException {
            if (!directory.exists() && !directory.mkdirs())
                throw new IOException("Cannot create " + directory);
            summary = new BufferedWriter(new FileWriter(new File(directory, "window_summary.csv")));
            summary.write("Run_Label;Stochastic_Replicate;Arm;Window;AHL_Input;Acid_Input;"
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
            pLog = new BufferedWriter(new FileWriter(new File(directory, "p_ac_timeseries.txt")));
            pLog.write("# t_s P_ac_uM h P_live_uM\n");
            writeFeatureContract(directory);
        }

        @Override
        public void close() {
            if (closed) return;
            closed = true;
            try {
                summary.close();
                if (samples != null) samples.close();
                if (voxels != null) voxels.close();
                if (pLog != null) pLog.close();
            } catch (IOException e) {
            throw new RuntimeException("Cannot close E42 output", e);
            }
        }
    }

    public static void main(String[] args) {
        String configPath = args.length > 0 ? args[0] : "sim_config_e42.properties";
        loadConfig(configPath);
        if (args.length > 1 && "preview".equals(args[1])) HEADLESS = false;
        resetStaticState();
        if (FLOW_SPEED != 0.0)
            throw new IllegalStateException("E42 claim dish requires FLOW_SPEED=0");
        if (ARM != Arm.SILENT) {
            double relJ = Math.abs(AHL_SOURCE_RATE - JMAX_A1_FROZEN) / JMAX_A1_FROZEN;
            if (relJ > 1e-9)
                throw new IllegalStateException("A1 Jmax is frozen at payload-matched "
                        + JMAX_A1_FROZEN + "; do not retune after scores");
        }
        experimentRng = new Random(RNG_SEED);
        BrownianParticle.seedMotionRng(RNG_SEED);
        final double[] ahlInputs = readGradedSequence(INPUT_AHL_FILE, "AHL");
        final double[] acidInputs = readGradedSequence(INPUT_ACID_FILE, "acid");
        if (ahlInputs.length < NUM_WINDOWS)
            throw new IllegalArgumentException("AHL input has " + ahlInputs.length
                    + " values; num.windows=" + NUM_WINDOWS);
        if (acidInputs.length < NUM_WINDOWS)
            throw new IllegalArgumentException("acid input has " + acidInputs.length
                    + " values; num.windows=" + NUM_WINDOWS);
        double sumU = 60.0 * 0.5;
        double sumG = 60.0 * acGate(0.5);
        int hashed = Math.min(ahlInputs.length, 200);
        for (int i = 0; i < hashed; i++) {
            sumU += ahlInputs[i];
            sumG += acGate(ahlInputs[i]);
        }
        double massA0 = JMAX_A0 * 75.0 * sumU;
        double massA1 = JMAX_A1_FROZEN * 75.0 * sumG;
        double payloadRel = Math.abs(massA1 - massA0) / massA0;
        double h0 = acFeedback(0.0);
        double hK = acFeedback(AC_K_P);
        double jTw0 = JMAX_A1_FROZEN * acGate(0.4) * acFeedback(0.0);
        double jA1 = JMAX_A1_FROZEN * acGate(0.4);
        if (COUPLING != Coupling.CLOSED) loadReplay(INPUT_P_REPLAY_FILE);
        System.out.printf(Locale.US,
                "E42 MODEL_STATUS=%s digital_twin=false dish=%.0fx%.0fx%.0f CENTER "
                        + "FLOW_SPEED=%.1f AHL_source=(500,250,5) coupling=%s%n",
                MODEL_STATUS, BOUND_X, BOUND_Y, BOUND_Z, FLOW_SPEED,
                COUPLING.name().toLowerCase(Locale.ROOT));
        System.out.printf(Locale.US,
                "E42 A1 g0=%.6g n_AC=%.6g K_AC=%.6g Jmax_A1=%.12e g(0)=%.12f g(0.5)=%.12f "
                        + "payload_rel_err=%.3e%n",
                AC_G0, AC_N, AC_K, AHL_SOURCE_RATE, acGate(0.0), acGate(0.5), payloadRel);
        System.out.printf(Locale.US,
                "E42 TW n_P=%.6g K_P=%.6g k_P=%.12e alpha=%.6g h(0)=%.12f h(K_P)=%.12f "
                        + "h(inf)=%.12f P0_identity_rel=%.3e D_P=%.6g decay_P=%.6g%n",
                AC_N_P, AC_K_P, K_P_RATE, AC_ALPHA, h0, hK, acFeedback(1.0e12),
                Math.abs(jTw0 - jA1) / Math.max(1.0, Math.abs(jA1)),
                P_DIFFUSIVITY, P_DECAY);
        System.out.printf(Locale.US,
                "E42 arm=%s seed=%d num.windows=%d pulse=%.0f window=%.0f%n",
                ARM.name().toLowerCase(Locale.ROOT), RNG_SEED, NUM_WINDOWS,
                PULSE_DURATION, WINDOW_DURATION);
        if (payloadRel > 0.01)
            throw new IllegalStateException("A1 payload match exceeded 1%; do not start living run");
        if (Math.abs(h0 - 1.0) > 1e-12)
            throw new IllegalStateException("h(0) must be 1 so P=0 recovers A1");

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
        final BSimChemicalField reporterField = new BSimChemicalField(
                sim, fieldGrid, P_DIFFUSIVITY, P_DECAY);
        voxelAnalyzer = new VoxelAnalyzer(new int[]{STATE_X, STATE_Y, STATE_Z},
                new int[]{COUNT_X, COUNT_Y, COUNT_Z}, BOUND_X, BOUND_Y, BOUND_Z);

        for (int i = 0; i < INITIAL_POP; i++) {
            Vector3d p = new Vector3d(300 + experimentRng.nextDouble() * 400,
                    150 + experimentRng.nextDouble() * 200, BOUND_Z / 2);
            if (ARM == Arm.BROWNIAN) {
                particles.add(new BrownianParticle(sim, p));
            } else {
                ReservoirBacterium bacterium = new ReservoirBacterium(
                        sim, p, attractantField, repellentField, ahlField, acidField, reporterField);
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
        final int pLogInterval = Math.max(1, (int) Math.round(1.0 / DT));
        final int[] births = {0}, deaths = {0}, clampDeaths = {0};
        final int[] acidDeaths = {0}, oobDeaths = {0}, completed = {0};
        final WindowAccumulator windowAccumulator = new WindowAccumulator();
        final OutputOwner output;
        try {
            output = new OutputOwner(new File(OUTPUT_DIR));
        } catch (IOException e) {
            throw new RuntimeException("Cannot open E42 output", e);
        }

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
                        finishWindow(output, window - 1, ahlInputs[window - 1],
                                acidInputs[window - 1], births[0], deaths[0], clampDeaths[0],
                                acidDeaths[0], oobDeaths[0], windowAccumulator);
                        completed[0]++;
                    }
                    births[0] = deaths[0] = clampDeaths[0] = acidDeaths[0] = oobDeaths[0] = 0;
                    voxelAnalyzer.resetWindowCounters();
                    windowAccumulator.reset();
                }

                double pLive = reporterField.getConc(AC_POSITIONS[1]) / MOLECULES_PER_UM3_PER_UM;
                double pSense = (COUPLING == Coupling.CLOSED)
                        ? pLive
                        : replayPAt(step * DT);
                if (pSense < 0.0) pSense = 0.0;
                double h = acFeedback(pSense);
                double input0 = 0.0, input2 = 0.0;
                if (warmup && step % stepsInWindow < stepsInPulse) {
                    ahlField.addQuantity(AC_POSITIONS[1],
                            AHL_SOURCE_RATE * acGate(WARMUP_AHL_INPUT) * h * DT);
                }
                if (!warmup && window < NUM_WINDOWS && stepInWindow < stepsInPulse) {
                    attractantField.addQuantity(AC_POSITIONS[0], PROD_RATE * input0 * DT);
                    attractantField.addQuantity(AC_POSITIONS[2], PROD_RATE * input2 * DT);
                    ahlField.addQuantity(AC_POSITIONS[1],
                            AHL_SOURCE_RATE * acGate(ahlInputs[window]) * h * DT);
                    acidField.addQuantity(ACID_SOURCE_POSITION,
                            ACID_SOURCE_RATE * acidInputs[window] * DT);
                }

                attractantField.update();
                repellentField.update();
                ahlField.update();
                acidField.update();
                reporterField.update();
                double sourceConc = ahlField.getConc(AC_POSITIONS[1]);
                if (!Double.isFinite(sourceConc)) ahlFinite = false;
                if (sourceConc < -1e-12) ahlNonNegative = false;
                double pAfter = reporterField.getConc(AC_POSITIONS[1]) / MOLECULES_PER_UM3_PER_UM;
                if (!Double.isFinite(pAfter) || !Double.isFinite(pLive)) pFinite = false;
                if (pAfter < -1e-12 || pLive < -1e-12) pNonNegative = false;
                if (step % pLogInterval == 0) {
                    try {
                        output.pLog.write(String.format(Locale.US,
                                "%.2f %.12e %.12e %.12e%n",
                                step * DT, pSense, h, pLive));
                    } catch (IOException e) {
                        throw new RuntimeException("Cannot write P_ac timeseries", e);
                    }
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
                                    ahlInputs[window], acidInputs[window],
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
        try {
            if (HEADLESS) {
                System.setOut(new PrintStream(OutputStream.nullOutputStream()));
                sim.export();
            } else {
                sim.preview();
                return;
            }
            while (completed[0] < NUM_WINDOWS) {
                finishWindow(output, completed[0], ahlInputs[completed[0]],
                        acidInputs[completed[0]], births[0], deaths[0], clampDeaths[0],
                        acidDeaths[0], oobDeaths[0], windowAccumulator);
                completed[0]++;
            }
        } finally {
            System.setOut(originalOut);
            output.close();
        }
        writeRunStatus(output);
        System.out.printf(Locale.US,
                "E42 complete arm=%s coupling=%s FLOW_SPEED=%.1f AHL_finite=%s AHL_nonneg=%s "
                        + "P_finite=%s P_nonneg=%s summary_rows=%d expected=%d sample_rows=%d voxel_rows=%d%n",
                ARM.name().toLowerCase(Locale.ROOT),
                COUPLING.name().toLowerCase(Locale.ROOT), FLOW_SPEED,
                ahlFinite ? "true" : "false", ahlNonNegative ? "true" : "false",
                pFinite ? "true" : "false", pNonNegative ? "true" : "false",
                output.summaryRows, NUM_WINDOWS,
                output.sampleRows, output.voxelRows);
    }

    static double acidMmFromConc(double moleculesPerUm3) {
        return moleculesPerUm3 / MM_TO_MOLECULES_PER_UM3;
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
                    "%s;%s;%s;%d;%.2f;%.2f;%.9g;%.9g;%.9g;%.9g;%.9g;%.9g;"
                            + "%.9g;%.9g;%.9g;%.9g;%.9g;%.9g;%.9g;"
                            + "%.9g;%.9g;%.9g;%.9g;%d;%d;%d;%d;%d;%d%n",
                    RUN_LABEL, STOCHASTIC_REPLICATE, ARM.name().toLowerCase(Locale.ROOT),
                    window, ahlInput, acidInput,
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
                    for (double v : r.ahlUm) {
                        if (!Double.isFinite(v) || v < -1e-12) {
                            ahlFinite = false;
                            if (v < -1e-12) ahlNonNegative = false;
                        }
                        output.voxels.write(format(v));
                    }
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

    static void writeFeatureContract(File directory) throws IOException {
        try (BufferedWriter w = new BufferedWriter(
                new FileWriter(new File(directory, "feature_contract.txt")))) {
            w.write("arm=" + ARM.name().toLowerCase(Locale.ROOT) + "\n");
            w.write("coupling=" + COUPLING.name().toLowerCase(Locale.ROOT) + "\n");
            w.write("MODEL_STATUS=HYPOTHETICAL_DESIGN_ENVELOPE\n");
            w.write("digital_twin=false\n");
            w.write("transducer=A1_static_Hill_gate_times_h_P\n");
            w.write("vesicle_store=false\n");
            w.write("feedback_species=inert_reporter_P\n");
            w.write("P_in_AHL_receiver=false\n");
            w.write("stage4_status=FAIL\n");
            w.write("stage5_status=PASS\n");
            if (ARM == Arm.BROWNIAN) {
                w.write("ridge_features=Den (20x10 spatial density)\n");
                w.write("excluded=all chemical-field and biological channels\n");
                w.write("note=chemical fields still update on the dish; they are not readout features\n");
            } else {
                w.write("ridge_features=Receiver_R / Mean_q (20x10); Lum_Mean / Mean_L (20x10); Input_Driven_Death (4x2)\n");
                w.write("excluded=window AHL, occupancy, pH, Births, Total_Deaths, Clamp_Deaths, OOB, Population, Lum_Sum, voxel AHL_uM, voxel pH, voxel Den, voxel Fraction_q, voxel Lum_Sum\n");
            }
        }
    }

    static void writeRunStatus(OutputOwner output) {
        File directory = new File(OUTPUT_DIR);
        try (BufferedWriter w = new BufferedWriter(
                new FileWriter(new File(directory, "run_status.txt")))) {
            w.write(String.format(Locale.US,
                    "arm=%s%ncoupling=%s%nseed=%d%nMODEL_STATUS=%s%nJmax_A1=%.12e%ng0=%.6g%nn_AC=%.6g%nK_AC=%.6g%n"
                            + "n_P=%.6g%nK_P=%.6g%nk_P=%.12e%nalpha=%.6g%n"
                            + "FLOW_SPEED=%.1f%nAHL_finite=%s%nAHL_nonneg=%s%nP_finite=%s%nP_nonneg=%s%n"
                            + "summary_rows=%d%nsample_rows=%d%nvoxel_rows=%d%n",
                    ARM.name().toLowerCase(Locale.ROOT),
                    COUPLING.name().toLowerCase(Locale.ROOT),
                    RNG_SEED, MODEL_STATUS,
                    AHL_SOURCE_RATE, AC_G0, AC_N, AC_K,
                    AC_N_P, AC_K_P, K_P_RATE, AC_ALPHA, FLOW_SPEED,
                    ahlFinite ? "true" : "false", ahlNonNegative ? "true" : "false",
                    pFinite ? "true" : "false", pNonNegative ? "true" : "false",
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
        ACID_SOURCE_RATE = getDouble(p, "field.acid.source.rate", ACID_SOURCE_RATE);
        ACID_CELL_PRODUCTION_RATE = getDouble(p,
                "field.acid.cell.production.rate", ACID_CELL_PRODUCTION_RATE);
        K_MAX_ACID = getDouble(p, "field.acid.kmax", K_MAX_ACID);
        K_MAX_ALK = getDouble(p, "field.acid.kmax.alk", K_MAX_ALK);
        ALPHA_LUX = getDouble(p, "luminescence.alpha", ALPHA_LUX);
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
        INPUT_AHL_FILE = p.getProperty("input.ahl.file", INPUT_AHL_FILE);
        INPUT_ACID_FILE = p.getProperty("input.acid.file", INPUT_ACID_FILE);
        INPUT_P_REPLAY_FILE = p.getProperty("input.p_replay.file", INPUT_P_REPLAY_FILE);
        RUN_LABEL = p.getProperty("output.run.label", RUN_LABEL);
        STOCHASTIC_REPLICATE = p.getProperty(
                "output.stochastic.replicate", STOCHASTIC_REPLICATE);
        String armName = p.getProperty("arm", "driven").trim().toLowerCase(Locale.ROOT);
        if ("brownian".equals(armName)) ARM = Arm.BROWNIAN;
        else if ("silent".equals(armName)) ARM = Arm.SILENT;
        else if ("driven".equals(armName)) ARM = Arm.DRIVEN;
        else throw new IllegalArgumentException("arm must be brownian, silent, or driven");
        String couplingName = p.getProperty("coupling", "closed").trim().toLowerCase(Locale.ROOT);
        if ("closed".equals(couplingName)) COUPLING = Coupling.CLOSED;
        else if ("replay".equals(couplingName)) COUPLING = Coupling.REPLAY;
        else if ("shuffle".equals(couplingName)) COUPLING = Coupling.SHUFFLE;
        else throw new IllegalArgumentException("coupling must be closed, replay, or shuffle");
        if (ARM == Arm.SILENT) {
            AHL_SOURCE_RATE = 0.0;
            ACID_SOURCE_RATE = 0.0;
            WARMUP_AHL_INPUT = 0.0;
        }
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
    static double acGate(double u) {
        if (u < 0.0) u = 0.0;
        double un = Math.pow(u, AC_N);
        double kn = Math.pow(AC_K, AC_N);
        return AC_G0 + (1.0 - AC_G0) * un / (kn + un);
    }

    static double acFeedback(double pUm) {
        if (pUm < 0.0) pUm = 0.0;
        double pn = Math.pow(pUm, AC_N_P);
        double kn = Math.pow(AC_K_P, AC_N_P);
        return 1.0 + AC_ALPHA * pn / (kn + pn);
    }

    static void loadReplay(String path) {
        if (path == null || path.trim().isEmpty())
            throw new IllegalArgumentException("replay/shuffle requires input.p_replay.file");
        List<Double> times = new ArrayList<>();
        List<Double> values = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                String[] parts = line.split("[;,\t ]+");
                times.add(Double.parseDouble(parts[0]));
                values.add(Double.parseDouble(parts[1]));
            }
        } catch (IOException | NumberFormatException e) {
            throw new IllegalArgumentException("Cannot read P replay " + path, e);
        }
        if (times.isEmpty())
            throw new IllegalArgumentException("P replay file is empty: " + path);
        replayTime = new double[times.size()];
        replayP = new double[values.size()];
        for (int i = 0; i < times.size(); i++) {
            replayTime[i] = times.get(i);
            replayP[i] = values.get(i);
        }
        replayCursor = 0;
    }

    static double replayPAt(double t) {
        if (replayP == null || replayP.length == 0) return 0.0;
        if (t <= replayTime[0]) return replayP[0];
        while (replayCursor + 1 < replayTime.length && replayTime[replayCursor + 1] <= t)
            replayCursor++;
        if (replayCursor >= replayTime.length - 1) return replayP[replayP.length - 1];
        double t0 = replayTime[replayCursor];
        double t1 = replayTime[replayCursor + 1];
        double p0 = replayP[replayCursor];
        double p1 = replayP[replayCursor + 1];
        if (t1 <= t0) return p0;
        return p0 + (t - t0) / (t1 - t0) * (p1 - p0);
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
        ahlFinite = true; ahlNonNegative = true;
        pFinite = true; pNonNegative = true;
        replayTime = null; replayP = null; replayCursor = 0;
    }
}
