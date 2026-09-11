package BSimReservoirPlan;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.draw.BSimP3DDrawer;
import bsim.ode.BSimOdeSolver;
import bsim.ode.BSimOdeSystem;
import bsim.particle.BSimBacterium;
import processing.core.PGraphics3D;

import javax.vecmath.Vector3d;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Properties;
import java.util.Random;
import java.util.TreeMap;
import java.util.Vector;

/**
 * Plan Stages 1-2 rebuilt from BSimReservoirStage10.
 *
 * Preserved biology: Stage10 fixed 30-minute surface-area growth, LacOperon
 * homeostatic clamp, Danino QS ODE, direct on/off AC sources, and the archived
 * repellent sign-flip plumbing. No glucose, Monod, nutrient AC, ODE-gated
 * vesicle AC, or attractant-growth boost is present.
 *
 * Added here: independent 50x25 field, 20x10 state readout, 4x2 event-count
 * readout, configurable window protocol, and numeric Stage 1/2 gates.
 */
public final class BSimReservoirPlan {
    static final double BOUND_X = 1000.0;
    static final double BOUND_Y = 500.0;
    static final double BOUND_Z = 10.0;

    static double DT = 0.05;
    static int GRID_X = 50;
    static int GRID_Y = 25;
    static int GRID_Z = 1;
    static int STATE_X = 20;
    static int STATE_Y = 10;
    static int STATE_Z = 1;
    static int COUNT_X = 4;
    static int COUNT_Y = 2;
    static int COUNT_Z = 1;

    static double ATT_DIFFUSIVITY = 100.0;
    static double REP_DIFFUSIVITY = 100.0;
    static double AHL_DIFFUSIVITY = 159.0;
    static double ATT_DECAY = 0.0067;
    static double REP_DECAY = 0.033;
    static double AHL_DECAY = 0.0033;
    static final double PROD_RATE = 1e6;

    static double WARMUP = 18000.0;
    static double WINDOW_DURATION = 300.0;
    static double PULSE_DURATION = 75.0;
    static double SAMPLING_DURATION = 300.0;
    static double SAMPLING_INTERVAL = 20.0;
    static int NUM_WINDOWS = 40;

    static int INITIAL_POP = 1800;
    static int CARRYING_CAPACITY = 2000;
    static boolean HEADLESS = true;
    static boolean METRICS_ONLY = false;
    static long RNG_SEED = 11;
    static String METRICS_TAG = "default";
    static String OUTPUT_DIR = "results";

    static final double FLOW_SPEED = 0.0;
    static final double GROWTH_RATE = 4.0 * Math.PI / 1800.0;
    static final double EXPECTED_T_GEN = 4.0 * Math.PI / GROWTH_RATE;
    static final double T_REMOVAL = EXPECTED_T_GEN;

    static final double TOX_EC50 = 150.0;
    static final double TOX_K_MAX = 0.001;
    static final double TOX_N = 2.0;

    static final double TIME_ADJ = 60.0;
    static final double QS_DELTA1 = 0.8487 / TIME_ADJ;
    static final double QS_DELTA2 = 0.0234 / TIME_ADJ;
    static final double QS_G = 0.0412;
    static final double QS_KP2 = 9.0 / TIME_ADJ;
    static final double QS_KR1OFF = 6e-6 / TIME_ADJ;
    static final double QS_KR1ON = 5.99e-5 / TIME_ADJ;
    static final double QS_KCAT_AIIA = 2631.4 / TIME_ADJ;
    static final double QS_T_A = 0.00276 / TIME_ADJ;
    static final double QS_T_LA = 0.024 / TIME_ADJ;
    static final double QS_A0LI = 7.785e-6 / TIME_ADJ;
    static final double QS_A0AA = 6.183e-6 / TIME_ADJ;
    static final double QS_KPLI = 0.9 / TIME_ADJ;
    static final double QS_KPAA = 0.9 / TIME_ADJ;
    static final double QS_KMLA = 1e-2;
    static final double QS_KMAA = 1200.0;
    static final double QS_LTOT = 15.0;
    static final double QS_N = 2.0;
    static final double CELL_WALL_DIFF = 3.0 / TIME_ADJ;

    // Matches the Stage10 comments: AC0/AC2 attractant, AC1 repellent.
    static final boolean[] AC_IS_REPELLENT = {false, true, false};
    static final Vector3d[] AC_POSITIONS = {
            new Vector3d(250, 250, 5),
            new Vector3d(500, 250, 5),
            new Vector3d(750, 250, 5)
    };

    static final Vector<ReservoirBacterium> bacteria = new Vector<>();
    static final Vector<ReservoirBacterium> children = new Vector<>();
    static final Vector<ReservoirBacterium> removals = new Vector<>();
    static VoxelAnalyzer voxelAnalyzer;
    static Random experimentRng;
    static int nextId;
    static int cumulativeDeaths;

    static final class AC {
        final int id;
        final Vector3d position;
        final int[] sequence;
        final boolean repellent;
        boolean active;

        AC(int id, Vector3d position, int[] sequence, boolean repellent) {
            this.id = id;
            this.position = position;
            this.sequence = sequence;
            this.repellent = repellent;
        }

        int inputForWindow(int window) {
            return sequence[window % sequence.length];
        }
    }

    static int[] readSequenceFile(String path, int acId) {
        try (BufferedReader reader = new BufferedReader(new FileReader(path))) {
            reader.readLine(); // Stage10 bit duration; Plan protocol owns window timing.
            String line = reader.readLine();
            if (line == null) throw new IOException("missing sequence line");
            String[] tokens = line.trim().split("\\s+");
            int[] sequence = new int[tokens.length];
            for (int i = 0; i < tokens.length; i++) sequence[i] = Integer.parseInt(tokens[i]);
            return sequence;
        } catch (IOException | NumberFormatException e) {
            System.err.printf("Could not read %s for AC%d; using 1,0 fallback.%n", path, acId);
            return new int[]{1, 0};
        }
    }

    public static final class ReservoirBacterium extends BSimBacterium {
        final int id;
        final BSimChemicalField attractantField;
        final BSimChemicalField repellentField;
        final BSimChemicalField ahlField;
        final double[] repMemory;
        double[] y = {0.05, 0.05, 0.05, 0.05};
        double diffConc;

        ReservoirBacterium(BSim sim, Vector3d position,
                           BSimChemicalField attractantField,
                           BSimChemicalField repellentField,
                           BSimChemicalField ahlField) {
            super(sim, position);
            id = nextId++;
            this.attractantField = attractantField;
            this.repellentField = repellentField;
            this.ahlField = ahlField;
            setGoal(attractantField);
            int memoryLength = sim.timesteps(shortTermMemoryDuration + longTermMemoryDuration);
            repMemory = new double[memoryLength];
            java.util.Arrays.fill(repMemory, repellentField.getConc(position));
        }

        double getLuxI() {
            return y[0];
        }

        private double repellentGradientDelta() {
            System.arraycopy(repMemory, 0, repMemory, 1, repMemory.length - 1);
            repMemory[0] = repellentField.getConc(position);
            double shortSum = 0;
            double longSum = 0;
            for (int i = 0; i < repMemory.length; i++) {
                if (i < shortTermMemoryLength) shortSum += repMemory[i];
                else longSum += repMemory[i];
            }
            return shortSum / shortTermMemoryLength - longSum / longTermMemoryLength;
        }

        @Override
        public double pEndRun() {
            boolean upAttractant = goal != null && movingUpGradient();
            double repellentDelta = repellentGradientDelta();
            boolean upRepellent = repellentDelta > sensitivity;
            boolean downRepellent = repellentDelta < -sensitivity;
            if (upAttractant && !upRepellent) return pEndRunUp;
            if (upRepellent && !upAttractant) return pEndRunElse;
            if (downRepellent && !upAttractant) return pEndRunUp;
            return pEndRunElse;
        }

        @Override
        public void action() {
            super.action();
            addForce(new Vector3d(stokesCoefficient() * FLOW_SPEED, 0, 0));

            // Preserved verbatim from Stage10. Unit repair belongs to a later QS stage.
            diffConc = y[1] * 1e15 - ahlField.getConc(position);
            y = BSimOdeSolver.rungeKutta45(new QSGRN(), sim.getTime(), y, sim.getDt());
            for (int i = 0; i < y.length; i++) if (y[i] < 0) y[i] = 0;
            double cellVolume = 4.0 * Math.PI * Math.pow(radius, 3) / 3.0;
            ahlField.addQuantity(position, diffConc * CELL_WALL_DIFF * sim.getDt() * cellVolume);

            double pRemoval = (sim.getDt() / T_REMOVAL)
                    * Math.pow(2.0, -(1.0 - (double) bacteria.size() / CARRYING_CAPACITY));
            if (experimentRng.nextDouble() < pRemoval) {
                removals.add(this);
                return;
            }

            double repellent = repellentField.getConc(position);
            if (repellent > 0) {
                double numerator = Math.pow(repellent, TOX_N);
                double killRate = TOX_K_MAX * numerator
                        / (Math.pow(TOX_EC50, TOX_N) + numerator);
                if (experimentRng.nextDouble() < killRate * sim.getDt()) removals.add(this);
            }
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
            childPosition.x += 2 * radius * (experimentRng.nextDouble() - 0.5);
            childPosition.y += 2 * radius * (experimentRng.nextDouble() - 0.5);
            ReservoirBacterium child = new ReservoirBacterium(sim, childPosition,
                    attractantField, repellentField, ahlField);
            child.setRadius(radius);
            child.setSurfaceAreaGrowthRate(surfaceAreaGrowthRate);
            child.setChildList(childList);
            child.y = y.clone();
            childList.add(child);
            if (voxelAnalyzer != null) voxelAnalyzer.recordBirth(childPosition);
        }

        final class QSGRN implements BSimOdeSystem {
            @Override
            public double[] derivativeSystem(double t, double[] state) {
                double exchange = diffConc * 1e-15;
                double[] d = new double[4];
                d[0] = QS_A0LI
                        + QS_KPLI * Math.pow(state[3], QS_N)
                        / (Math.pow(QS_KMLA, QS_N) + Math.pow(state[3], QS_N))
                        - QS_DELTA1 * state[0] / (QS_G * (state[0] + state[2]) + 1);
                d[1] = QS_KP2 * state[0]
                        - QS_KR1ON * (QS_LTOT - state[3]) * state[1]
                        + QS_KR1OFF * state[3]
                        - QS_KCAT_AIIA * state[2] * state[1] / (QS_KMAA + state[1])
                        - QS_T_A * state[1] - CELL_WALL_DIFF * exchange;
                d[2] = QS_A0AA
                        + QS_KPAA * Math.pow(state[3], QS_N)
                        / (Math.pow(QS_KMLA, QS_N) + Math.pow(state[3], QS_N))
                        - QS_DELTA2 * state[2] / (QS_G * (state[0] + state[2]) + 1);
                d[3] = QS_KR1ON * (QS_LTOT - state[3]) * state[1]
                        - QS_KR1OFF * state[3] - QS_T_LA * state[3];
                return d;
            }

            @Override
            public int getNumEq() {
                return 4;
            }

            @Override
            public double[] getICs() {
                return new double[]{0.05, 0.05, 0.05, 0.05};
            }
        }
    }

    static void loadConfig(String path) {
        Properties p = new Properties();
        try (FileInputStream input = new FileInputStream(path)) {
            p.load(input);
        } catch (IOException e) {
            System.out.println("No config file " + path + "; using defaults.");
            return;
        }
        DT = getDouble(p, "dt", DT);
        GRID_X = getInt(p, "grid.x", GRID_X);
        GRID_Y = getInt(p, "grid.y", GRID_Y);
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
        ATT_DECAY = getDouble(p, "field.att.decay", ATT_DECAY);
        REP_DECAY = getDouble(p, "field.rep.decay", REP_DECAY);
        AHL_DECAY = getDouble(p, "field.ahl.decay", AHL_DECAY);
        WARMUP = getDouble(p, "warmup.s", WARMUP);
        WINDOW_DURATION = getDouble(p, "window.duration.s", WINDOW_DURATION);
        PULSE_DURATION = getDouble(p, "pulse.duration.s", PULSE_DURATION);
        SAMPLING_DURATION = getDouble(p, "sampling.duration.s", SAMPLING_DURATION);
        SAMPLING_INTERVAL = getDouble(p, "sampling.interval.s", SAMPLING_INTERVAL);
        NUM_WINDOWS = getInt(p, "num.windows", NUM_WINDOWS);
        INITIAL_POP = getInt(p, "initial.pop", INITIAL_POP);
        CARRYING_CAPACITY = getInt(p, "carrying.capacity", CARRYING_CAPACITY);
        HEADLESS = Boolean.parseBoolean(p.getProperty("headless", String.valueOf(HEADLESS)));
        METRICS_ONLY = Boolean.parseBoolean(p.getProperty("regression.metrics.only",
                String.valueOf(METRICS_ONLY)));
        RNG_SEED = Long.parseLong(p.getProperty("rng.seed", String.valueOf(RNG_SEED)));
        METRICS_TAG = p.getProperty("output.metrics.tag", METRICS_TAG);
        OUTPUT_DIR = p.getProperty("output.dir", OUTPUT_DIR);
        PULSE_DURATION = Math.min(PULSE_DURATION, WINDOW_DURATION);
        SAMPLING_DURATION = Math.min(SAMPLING_DURATION, WINDOW_DURATION);
        System.out.println("Loaded config " + path);
    }

    static double getDouble(Properties p, String key, double fallback) {
        return Double.parseDouble(p.getProperty(key, String.valueOf(fallback)));
    }

    static int getInt(Properties p, String key, int fallback) {
        return Integer.parseInt(p.getProperty(key, String.valueOf(fallback)));
    }

    public static void main(String[] args) {
        String configPath = args.length > 0 ? args[0] : "sim_config.properties";
        loadConfig(configPath);
        if (args.length > 1 && "preview".equals(args[1])) HEADLESS = false;

        resetStaticState();
        experimentRng = new Random(RNG_SEED);
        List<AC> acs = new ArrayList<>();
        for (int i = 0; i < AC_POSITIONS.length; i++) {
            acs.add(new AC(i, AC_POSITIONS[i],
                    readSequenceFile("input_sequence" + i + ".txt", i),
                    AC_IS_REPELLENT[i]));
        }

        double simulationTime = WARMUP + NUM_WINDOWS * WINDOW_DURATION;
        printConfiguration(simulationTime);

        BSim sim = new BSim();
        sim.setDt(DT);
        sim.setSimulationTime(simulationTime);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);

        int[] fieldGrid = {GRID_X, GRID_Y, GRID_Z};
        BSimChemicalField attractantField = new BSimChemicalField(
                sim, fieldGrid, ATT_DIFFUSIVITY, ATT_DECAY);
        BSimChemicalField repellentField = new BSimChemicalField(
                sim, fieldGrid, REP_DIFFUSIVITY, REP_DECAY);
        BSimChemicalField ahlField = new BSimChemicalField(
                sim, fieldGrid, AHL_DIFFUSIVITY, AHL_DECAY);
        voxelAnalyzer = new VoxelAnalyzer(
                new int[]{STATE_X, STATE_Y, STATE_Z},
                new int[]{COUNT_X, COUNT_Y, COUNT_Z},
                BOUND_X, BOUND_Y, BOUND_Z);

        for (int i = 0; i < INITIAL_POP; i++) {
            Vector3d p = new Vector3d(
                    300 + experimentRng.nextDouble() * 400,
                    150 + experimentRng.nextDouble() * 200,
                    BOUND_Z / 2);
            ReservoirBacterium bacterium = new ReservoirBacterium(
                    sim, p, attractantField, repellentField, ahlField);
            bacterium.setRadius();
            bacterium.setSurfaceAreaGrowthRate(GROWTH_RATE);
            bacterium.setChildList(children);
            bacteria.add(bacterium);
        }

        int warmupSteps = (int) Math.round(WARMUP / DT);
        int stepsInWindow = (int) Math.round(WINDOW_DURATION / DT);
        int stepsInPulse = (int) Math.round(PULSE_DURATION / DT);
        int samplingStartStep = stepsInWindow
                - (int) Math.round(SAMPLING_DURATION / DT);
        int samplingIntervalSteps = Math.max(1, (int) Math.round(SAMPLING_INTERVAL / DT));
        int warmupSampleSteps = Math.max(1, (int) Math.round(300.0 / DT));

        File resultsDir = new File(OUTPUT_DIR);
        if (!METRICS_ONLY && !resultsDir.exists() && !resultsDir.mkdirs()) {
            throw new IllegalStateException("Cannot create " + resultsDir);
        }

        final FileWriter resultsWriter;
        final FileWriter voxelsWriter;
        try {
            if (METRICS_ONLY) {
                resultsWriter = null;
                voxelsWriter = null;
            } else {
                resultsWriter = new FileWriter(new File(resultsDir, "results.csv"));
                resultsWriter.write(
                        "Window;Sample;TimeInWindow_s;Input_AC0;Input_AC1;Input_AC2;"
                                + "Total_Count;Births_This_Window;Deaths_This_Window\n");
                voxelsWriter = new FileWriter(new File(resultsDir, "voxels.csv"));
                writeVoxelHeader(voxelsWriter);
            }
        } catch (IOException e) {
            throw new RuntimeException("Cannot open output files", e);
        }

        List<Integer> warmupPop = new ArrayList<>();
        List<Double> warmupTimes = new ArrayList<>();
        List<Integer> windowPop = new ArrayList<>();
        List<Double> windowLuxI = new ArrayList<>();
        List<int[]> birthBins = new ArrayList<>();
        List<int[]> deathBins = new ArrayList<>();
        List<Double> sampleOffsets = new ArrayList<>();
        int[] birthsThisWindow = {0};
        int[] deathsThisWindow = {0};
        int[] completedWindows = {0};

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                int step = (int) sim.getTimestep();
                boolean warmup = step < warmupSteps;
                int effectiveStep = step - warmupSteps;
                int window = warmup ? -1 : effectiveStep / stepsInWindow;
                int stepInWindow = warmup ? -1 : effectiveStep % stepsInWindow;

                if (warmup && (step % warmupSampleSteps == 0 || step == warmupSteps - 1)) {
                    warmupTimes.add(sim.getTime());
                    warmupPop.add(bacteria.size());
                }

                if (!warmup && window < NUM_WINDOWS && stepInWindow == 0) {
                    if (window > 0) finishWindow(
                            birthBins, deathBins, completedWindows,
                            birthsThisWindow, deathsThisWindow);
                    windowPop.add(bacteria.size());
                    windowLuxI.add(meanLuxI());
                    birthsThisWindow[0] = 0;
                    deathsThisWindow[0] = 0;
                    voxelAnalyzer.resetWindowCounters();
                    for (AC ac : acs) ac.active = ac.inputForWindow(window) == 1;
                    if (window % 10 == 0) {
                        System.out.printf("  [W%02d] N=%d inputs=%d,%d,%d%n",
                                window, bacteria.size(),
                                acs.get(0).inputForWindow(window),
                                acs.get(1).inputForWindow(window),
                                acs.get(2).inputForWindow(window));
                    }
                } else if (!warmup && stepInWindow >= stepsInPulse) {
                    for (AC ac : acs) ac.active = false;
                } else if (warmup) {
                    for (AC ac : acs) ac.active = false;
                }

                for (AC ac : acs) {
                    if (ac.active) {
                        BSimChemicalField target = ac.repellent
                                ? repellentField : attractantField;
                        target.addQuantity(ac.position, PROD_RATE * DT);
                    }
                }

                attractantField.update();
                repellentField.update();
                ahlField.update();

                for (ReservoirBacterium bacterium : bacteria) {
                    bacterium.action();
                    bacterium.updatePosition();
                }

                int born = children.size();
                bacteria.addAll(children);
                children.clear();
                if (!warmup) birthsThisWindow[0] += born;

                for (ReservoirBacterium bacterium : bacteria) {
                    if ((bacterium.getPosition().x < 0 || bacterium.getPosition().x > BOUND_X)
                            && !removals.contains(bacterium)) {
                        removals.add(bacterium);
                    }
                }
                if (!removals.isEmpty()) {
                    for (ReservoirBacterium bacterium : removals) {
                        if (!warmup) {
                            voxelAnalyzer.recordDeath(bacterium.getPosition());
                            deathsThisWindow[0]++;
                        }
                    }
                    cumulativeDeaths += removals.size();
                    bacteria.removeAll(removals);
                    removals.clear();
                }

                if (!METRICS_ONLY && !warmup && window < NUM_WINDOWS
                        && stepInWindow >= samplingStartStep) {
                    boolean intervalSample =
                            (stepInWindow - samplingStartStep) % samplingIntervalSteps == 0;
                    boolean finalSample = stepInWindow == stepsInWindow - 1;
                    if (intervalSample || finalSample) {
                        double offset = stepInWindow * DT;
                        sampleOffsets.add(offset);
                        writeSample(resultsWriter, voxelsWriter, window, offset, acs,
                                birthsThisWindow[0], deathsThisWindow[0],
                                attractantField, repellentField, ahlField);
                    }
                }
            }
        });

        sim.setDrawer(new BSimP3DDrawer(sim, 800, 600) {
            @Override
            public void scene(PGraphics3D p3d) {
                p3d.ortho(0, (float) BOUND_X, (float) BOUND_Y, 0, -1000, 10000);
                p3d.camera((float) BOUND_X / 2, (float) BOUND_Y / 2, (float) BOUND_Y,
                        (float) BOUND_X / 2, (float) BOUND_Y / 2, 0,
                        0, 1, 0);
                p3d.perspective((float) Math.PI / 2,
                        (float) BOUND_X / (float) BOUND_Y, 0.1f, 10000f);
                draw(attractantField, Color.CYAN, (float) (255.0 / 540.0));
                draw(repellentField, Color.MAGENTA, (float) (255.0 / 540.0));
                for (AC ac : acs) {
                    sphere(ac.position, 15,
                            ac.active ? (ac.repellent ? Color.RED : Color.GREEN) : Color.DARK_GRAY,
                            255);
                }
                for (ReservoirBacterium bacterium : bacteria) {
                    sphere(bacterium.getPosition(), 8, Color.GREEN, 255);
                }
            }
        });

        long wallStart = System.currentTimeMillis();
        if (HEADLESS) sim.export();
        else {
            sim.preview();
            return;
        }
        long wallMs = System.currentTimeMillis() - wallStart;

        if (completedWindows[0] < NUM_WINDOWS) {
            finishWindow(birthBins, deathBins, completedWindows,
                    birthsThisWindow, deathsThisWindow);
        }

        if (METRICS_ONLY) {
            writeRegressionMetrics(windowPop, windowLuxI, cumulativeDeaths);
        } else {
            close(resultsWriter);
            close(voxelsWriter);
            writeGateSummary(resultsDir, wallMs, warmupTimes, warmupPop,
                    birthBins, deathBins, sampleOffsets);
        }
    }

    static void resetStaticState() {
        bacteria.clear();
        children.clear();
        removals.clear();
        voxelAnalyzer = null;
        nextId = 0;
        cumulativeDeaths = 0;
    }

    static void printConfiguration(double simulationTime) {
        System.out.println("============================================================");
        System.out.println("BSimReservoirPlan — Plan Stages 1-2");
        System.out.println("============================================================");
        System.out.printf(Locale.US, "dt=%.3f field_grid=%dx%dx%d domain=%.0fx%.0fx%.0f%n",
                DT, GRID_X, GRID_Y, GRID_Z, BOUND_X, BOUND_Y, BOUND_Z);
        System.out.printf("state_grid=%dx%dx%d count_grid=%dx%dx%d%n",
                STATE_X, STATE_Y, STATE_Z, COUNT_X, COUNT_Y, COUNT_Z);
        System.out.printf(Locale.US,
                "warmup=%.0fs window=%.0fs pulse=%.0fs sampling=%.0fs/%.0fs windows=%d%n",
                WARMUP, WINDOW_DURATION, PULSE_DURATION,
                SAMPLING_DURATION, SAMPLING_INTERVAL, NUM_WINDOWS);
        System.out.printf(Locale.US, "decay(att,rep,ahl)=%.4f,%.4f,%.4f%n",
                ATT_DECAY, REP_DECAY, AHL_DECAY);
        System.out.printf("initial_pop=%d K=%d seed=%d sim_time=%.0fs%n",
                INITIAL_POP, CARRYING_CAPACITY, RNG_SEED, simulationTime);
        System.out.println("AC types=ATTRACTANT,REPELLENT,ATTRACTANT");
        System.out.println("glucose_fields=0 monod_parameters=0 nutrient_ACs=0");
        System.out.println("============================================================");
    }

    static double meanLuxI() {
        double total = 0;
        for (ReservoirBacterium bacterium : bacteria) total += bacterium.getLuxI();
        return total / Math.max(1, bacteria.size());
    }

    static void finishWindow(List<int[]> birthBins, List<int[]> deathBins,
                             int[] completedWindows,
                             int[] birthsThisWindow, int[] deathsThisWindow) {
        birthBins.add(voxelAnalyzer.snapshotBirths());
        deathBins.add(voxelAnalyzer.snapshotDeaths());
        completedWindows[0]++;
    }

    static void writeVoxelHeader(FileWriter writer) throws IOException {
        StringBuilder header = new StringBuilder(
                "Window;Sample;TimeInWindow_s;Input_AC0;Input_AC1;Input_AC2");
        int stateVoxels = voxelAnalyzer.getStateVoxels();
        int countBins = voxelAnalyzer.getCountBins();
        for (int i = 0; i < stateVoxels; i++) header.append(";Att_").append(i);
        for (int i = 0; i < stateVoxels; i++) header.append(";Rep_").append(i);
        for (int i = 0; i < stateVoxels; i++) header.append(";AHL_").append(i);
        for (int i = 0; i < stateVoxels; i++) header.append(";Den_").append(i);
        for (int i = 0; i < stateVoxels; i++) header.append(";LuxI_").append(i);
        for (int i = 0; i < countBins; i++) header.append(";Birth_").append(i);
        for (int i = 0; i < countBins; i++) header.append(";Death_").append(i);
        writer.write(header.append('\n').toString());
    }

    static void writeSample(FileWriter resultsWriter, FileWriter voxelsWriter,
                            int window, double offset, List<AC> acs,
                            int births, int deaths,
                            BSimChemicalField attractantField,
                            BSimChemicalField repellentField,
                            BSimChemicalField ahlField) {
        try {
            int sample = (int) Math.round(offset / SAMPLING_INTERVAL);
            int u0 = acs.get(0).inputForWindow(window);
            int u1 = acs.get(1).inputForWindow(window);
            int u2 = acs.get(2).inputForWindow(window);
            resultsWriter.write(String.format(Locale.US,
                    "%d;%d;%.2f;%d;%d;%d;%d;%d;%d%n",
                    window, sample, offset, u0, u1, u2,
                    bacteria.size(), births, deaths));

            VoxelAnalyzer.VoxelReadout readout = voxelAnalyzer.analyze(
                    attractantField, repellentField, ahlField, bacteria);
            voxelsWriter.write(String.format(Locale.US, "%d;%d;%.2f;%d;%d;%d",
                    window, sample, offset, u0, u1, u2));
            for (double v : readout.attractant) voxelsWriter.write(format(v));
            for (double v : readout.repellent) voxelsWriter.write(format(v));
            for (double v : readout.ahl) voxelsWriter.write(format(v));
            for (int v : readout.density) voxelsWriter.write(";" + v);
            for (double v : readout.meanLuxI) voxelsWriter.write(format(v));
            for (int v : readout.births) voxelsWriter.write(";" + v);
            for (int v : readout.deaths) voxelsWriter.write(";" + v);
            voxelsWriter.write("\n");
        } catch (IOException e) {
            throw new RuntimeException("Cannot write sample", e);
        }
    }

    static String format(double value) {
        return String.format(Locale.US, ";%.4e", value);
    }

    static void close(FileWriter writer) {
        try {
            if (writer != null) writer.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    static void writeRegressionMetrics(List<Integer> populations,
                                       List<Double> luxI,
                                       int deaths) {
        File directory = new File("results/regression");
        if (!directory.exists() && !directory.mkdirs()) {
            throw new IllegalStateException("Cannot create regression directory");
        }
        Properties metrics = new Properties();
        metrics.setProperty("tag", METRICS_TAG);
        metrics.setProperty("mean_total_count", String.valueOf(meanInts(populations)));
        metrics.setProperty("mean_luxI", String.valueOf(meanDoubles(luxI)));
        metrics.setProperty("total_deaths", String.valueOf(deaths));
        metrics.setProperty("windows", String.valueOf(populations.size()));
        try (java.io.FileOutputStream output = new java.io.FileOutputStream(
                new File(directory, "metrics_" + METRICS_TAG + ".properties"))) {
            metrics.store(output, "Plan Stage 1 regression metrics");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        System.out.printf(Locale.US,
                "METRICS tag=%s Total_Count=%.6f Mean_LuxI=%.9g Total_Deaths=%d%n",
                METRICS_TAG, meanInts(populations), meanDoubles(luxI), deaths);
    }

    static double meanInts(List<Integer> values) {
        double total = 0;
        for (int value : values) total += value;
        return total / Math.max(1, values.size());
    }

    static double meanDoubles(List<Double> values) {
        double total = 0;
        for (double value : values) total += value;
        return total / Math.max(1, values.size());
    }

    static void writeGateSummary(File directory, long wallMs,
                                 List<Double> warmupTimes,
                                 List<Integer> warmupPop,
                                 List<int[]> birthBins,
                                 List<int[]> deathBins,
                                 List<Double> sampleOffsets) {
        File voxels = new File(directory, "voxels.csv");
        double peakWarmupSlope = peakWarmupSlope(warmupTimes, warmupPop);
        // A 3000 s trailing fit suppresses single-interval birth/death shot noise
        // while still measuring the derivative immediately before analysis.
        int plateauFitPoints = Math.min(11, warmupPop.size());
        double analysisStartSlope = trailingSlope(
                warmupTimes, warmupPop, plateauFitPoints);
        double plateauRatio = peakWarmupSlope > 0
                ? Math.abs(analysisStartSlope) / peakWarmupSlope : Double.POSITIVE_INFINITY;
        double meanBirths = meanBins(birthBins);
        double meanDeaths = meanBins(deathBins);
        double minSample = sampleOffsets.isEmpty() ? Double.NaN
                : sampleOffsets.stream().mapToDouble(Double::doubleValue).min().orElse(Double.NaN);
        double maxSample = sampleOffsets.isEmpty() ? Double.NaN
                : sampleOffsets.stream().mapToDouble(Double::doubleValue).max().orElse(Double.NaN);
        boolean samplingPass = minSample <= DT * 1.5
                && maxSample >= WINDOW_DURATION - DT * 1.5;

        StringBuilder report = new StringBuilder();
        report.append("BSimReservoirPlan numeric gates\n");
        report.append("================================\n");
        report.append(String.format(Locale.US,
                "stage1_field_grid=%dx%dx%d target=50x25x1%n", GRID_X, GRID_Y, GRID_Z));
        report.append(String.format("stage1_state_grid=%dx%dx%d target=20x10x1%n",
                STATE_X, STATE_Y, STATE_Z));
        report.append(String.format("stage1_count_grid=%dx%dx%d target=4x2x1%n",
                COUNT_X, COUNT_Y, COUNT_Z));
        report.append(String.format("stage1_windows=%d target=40%n", NUM_WINDOWS));
        report.append(String.format(Locale.US,
                "stage1_wall_s=%.3f threshold_s=1800.000 pass=%s%n",
                wallMs / 1000.0, wallMs < 1_800_000 ? "PASS" : "FAIL"));
        report.append(String.format(Locale.US,
                "stage1_voxels_bytes=%d voxels_mb=%.6f threshold_mb=50.000000 pass=%s%n",
                voxels.length(), voxels.length() / (1024.0 * 1024.0),
                voxels.length() < 50L * 1024 * 1024 ? "PASS" : "FAIL"));
        report.append("stage1_regression=see results/regression/regression_gate_summary.txt\n");
        report.append("stage1_mean_luxI_caveat=may be degenerate until Plan Stage 3; not a blocking gate\n");
        boolean stage2Protocol = WARMUP > 0
                && Math.abs(WINDOW_DURATION - 300.0) < 1e-9
                && Math.abs(PULSE_DURATION - 75.0) < 1e-9
                && Math.abs(SAMPLING_DURATION - 300.0) < 1e-9;
        if (!stage2Protocol) {
            report.append("stage2_gate=NOT_RUN_IN_STAGE1_PROTOCOL\n");
            report.append("plan_stage3_started=0\n");
            writeReport(directory, report);
            return;
        }
        report.append(String.format(Locale.US,
                "stage2_decay_att=%.6f target=0.006700%n", ATT_DECAY));
        report.append(String.format(Locale.US,
                "stage2_decay_rep=%.6f target=0.033000%n", REP_DECAY));
        report.append(String.format(Locale.US,
                "stage2_decay_ahl=%.6f target=0.003300%n", AHL_DECAY));
        report.append("stage2_glucose_fields=0\n");
        report.append(String.format(Locale.US,
                "stage2_warmup_peak_abs_dNdt=%.9f cells/s%n", peakWarmupSlope));
        report.append(String.format(Locale.US,
                "stage2_analysis_start_abs_dNdt=%.9f cells/s%n",
                Math.abs(analysisStartSlope)));
        report.append(String.format(
                "stage2_plateau_fit_points=%d warmup_sample_interval_s=300.000%n",
                plateauFitPoints));
        report.append("stage2_warmup_population_samples=");
        for (int i = 0; i < warmupPop.size(); i++) {
            if (i > 0) report.append(' ');
            report.append(String.format(Locale.US, "%.0f:%d",
                    warmupTimes.get(i), warmupPop.get(i)));
        }
        report.append('\n');
        report.append(String.format(Locale.US,
                "stage2_plateau_ratio=%.9f threshold=0.050000000 pass=%s%n",
                plateauRatio, plateauRatio < 0.05 ? "PASS" : "FAIL"));
        report.append(String.format(Locale.US,
                "stage2_mean_births_per_4x2_bin_window=%.6f threshold=5.000000 pass=%s%n",
                meanBirths, meanBirths >= 5 ? "PASS" : "FAIL"));
        report.append(String.format(Locale.US,
                "stage2_mean_deaths_per_4x2_bin_window=%.6f threshold=5.000000 pass=%s%n",
                meanDeaths, meanDeaths >= 5 ? "PASS" : "FAIL"));
        report.append("stage2_birth_histogram=").append(histogram(birthBins)).append('\n');
        report.append("stage2_death_histogram=").append(histogram(deathBins)).append('\n');
        report.append(String.format(Locale.US,
                "stage2_sampling_min_s=%.3f max_s=%.3f window_s=%.3f "
                        + "sampling_duration_s=%.3f pass=%s%n",
                minSample, maxSample, WINDOW_DURATION, SAMPLING_DURATION,
                samplingPass ? "PASS" : "FAIL"));
        report.append("stage2_birth_input_correlation_gate=NOT_APPLICABLE countability_only=1\n");
        boolean pass = plateauRatio < 0.05 && meanBirths >= 5
                && meanDeaths >= 5 && samplingPass;
        report.append("stage2_gate=").append(pass ? "PASS" : "FAIL").append('\n');
        report.append("plan_stage3_started=0\n");

        writeReport(directory, report);
    }

    static void writeReport(File directory, StringBuilder report) {
        try {
            Files.write(new File(directory, "gate_summary.txt").toPath(),
                    report.toString().getBytes(StandardCharsets.UTF_8));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        System.out.print(report);
    }

    static double peakWarmupSlope(List<Double> times, List<Integer> populations) {
        double peak = 0;
        for (int i = 1; i < populations.size(); i++) {
            double dt = times.get(i) - times.get(i - 1);
            if (dt > 0) {
                peak = Math.max(peak,
                        Math.abs((populations.get(i) - populations.get(i - 1)) / dt));
            }
        }
        return peak;
    }

    static double trailingSlope(List<Double> times, List<Integer> populations, int points) {
        int from = Math.max(0, populations.size() - points);
        int n = populations.size() - from;
        if (n < 2) return Double.NaN;
        double sumT = 0;
        double sumN = 0;
        double sumTT = 0;
        double sumTN = 0;
        for (int i = from; i < populations.size(); i++) {
            double t = times.get(i);
            double population = populations.get(i);
            sumT += t;
            sumN += population;
            sumTT += t * t;
            sumTN += t * population;
        }
        double denominator = n * sumTT - sumT * sumT;
        return denominator == 0 ? 0 : (n * sumTN - sumT * sumN) / denominator;
    }

    static double meanBins(List<int[]> windows) {
        long total = 0;
        int count = 0;
        for (int[] window : windows) {
            for (int value : window) {
                total += value;
                count++;
            }
        }
        return count == 0 ? 0 : (double) total / count;
    }

    static String histogram(List<int[]> windows) {
        TreeMap<Integer, Integer> histogram = new TreeMap<>();
        for (int[] window : windows) {
            for (int value : window) {
                histogram.put(value, histogram.getOrDefault(value, 0) + 1);
            }
        }
        StringBuilder text = new StringBuilder();
        for (java.util.Map.Entry<Integer, Integer> entry : histogram.entrySet()) {
            if (text.length() > 0) text.append(' ');
            text.append(entry.getKey()).append(':').append(entry.getValue());
        }
        return text.toString();
    }
}
