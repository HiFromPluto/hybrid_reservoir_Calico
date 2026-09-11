package BSimReservoirPlanStage3B;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.draw.BSimP3DDrawer;
import bsim.particle.BSimBacterium;
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Properties;
import java.util.Random;
import java.util.Vector;

/**
 * Spatial Plan Stage 3B. This copies the Stage 3 spatial, population, timing,
 * source, and output architecture, but replaces its Danino ODE with the
 * predeclared scalar fast receiver R.
 */
public final class BSimReservoirPlanStage3B {
    static final double BOUND_X = 1000.0, BOUND_Y = 500.0, BOUND_Z = 10.0;
    static double DT = .05;
    static int GRID_X = 50, GRID_Y = 25, GRID_Z = 1;
    static int STATE_X = 20, STATE_Y = 10, STATE_Z = 1;
    static int COUNT_X = 4, COUNT_Y = 2, COUNT_Z = 1;
    static double ATT_DIFFUSIVITY = 100, REP_DIFFUSIVITY = 100, AHL_DIFFUSIVITY = 159;
    static double ATT_DECAY = .0067, REP_DECAY = .033, AHL_DECAY = .0033;
    static double WARMUP = 18000, WINDOW_DURATION = 300, PULSE_DURATION = 75;
    static double SAMPLING_DURATION = 300, SAMPLING_INTERVAL = 20;
    static double WARMUP_AHL_INPUT = 0.5;
    static int NUM_WINDOWS = 40, INITIAL_POP = 1800, CARRYING_CAPACITY = 2000;
    static boolean HEADLESS = true, WRITE_SAMPLES = true, WRITE_VOXELS = true;
    static long RNG_SEED = 11;
    static String OUTPUT_DIR = "results/stage3b_perm1_rep1";
    static String INPUT_AHL_FILE = "input_ahl_balanced40.txt";
    static String RUN_LABEL = "stage3b_fast_receiver";
    static String STOCHASTIC_REPLICATE = "perm1_rep1";
    static double AHL_SOURCE_RATE = 128000000;

    static final double FLOW_SPEED = 0.0;
    static final double PROD_RATE = 1e6; // Preserved silent AC0/AC2 architecture.
    static final double GROWTH_RATE = 4.0 * Math.PI / 1800.0;
    static final double EXPECTED_T_GEN = 4.0 * Math.PI / GROWTH_RATE;
    static final double T_REMOVAL = EXPECTED_T_GEN;

    /** Exact predeclared receiver values; these are deliberately not configurable. */
    static final double RECEIVER_K_UM = 1.6;
    static final double RECEIVER_HILL_N = 2.0;
    static final double RECEIVER_TAU_S = 15.0;
    static final double RECEIVER_T95_S = 44.936;
    static final double MOLECULES_PER_UM3_PER_UM = 602.2;
    static final Vector3d[] AC_POSITIONS = {
            new Vector3d(250, 250, 5), new Vector3d(500, 250, 5),
            new Vector3d(750, 250, 5)
    };

    static final Vector<ReservoirBacterium> bacteria = new Vector<>();
    static final Vector<ReservoirBacterium> children = new Vector<>();
    static final Vector<ReservoirBacterium> removals = new Vector<>();
    static VoxelAnalyzer voxelAnalyzer;
    static Random experimentRng;
    static int nextId, cumulativeDeaths;

    public static final class ReservoirBacterium extends BSimBacterium {
        final int id;
        final BSimChemicalField attractantField, repellentField, ahlField;
        double receiver;

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
        }

        double getResponse() { return receiver; }
        double getQ() { return receiver; } // Compatibility name for the unchanged gate.

        @Override
        public void action() {
            super.action();
            addForce(new Vector3d(stokesCoefficient() * FLOW_SPEED, 0, 0));

            double concentrationUm = ahlField.getConc(position) / MOLECULES_PER_UM3_PER_UM;
            double concentrationN = Math.pow(concentrationUm, RECEIVER_HILL_N);
            double target = concentrationN
                    / (Math.pow(RECEIVER_K_UM, RECEIVER_HILL_N) + concentrationN);
            receiver = target + (receiver - target) * Math.exp(-sim.getDt() / RECEIVER_TAU_S);

            double pRemoval = (sim.getDt() / T_REMOVAL)
                    * Math.pow(2.0, -(1.0 - (double) bacteria.size() / CARRYING_CAPACITY));
            if (experimentRng.nextDouble() < pRemoval) removals.add(this);
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
                    attractantField, repellentField, ahlField);
            child.setRadius(radius);
            child.setSurfaceAreaGrowthRate(surfaceAreaGrowthRate);
            child.setChildList(childList);
            child.receiver = receiver;
            childList.add(child);
            if (voxelAnalyzer != null) voxelAnalyzer.recordBirth(childPosition);
        }
    }

    static final class WindowAccumulator {
        final List<Double> extracellularAhl = new ArrayList<>();
        final List<Double> response = new ArrayList<>();
        double qSum, populationSum;
        long cellObservations;
        int qAboveHalf, sampleCount;

        void reset() {
            extracellularAhl.clear();
            response.clear();
            qSum = populationSum = 0;
            cellObservations = 0;
            qAboveHalf = sampleCount = 0;
        }

        void observe(BSimChemicalField ahlField) {
            for (ReservoirBacterium bacterium : bacteria) {
                extracellularAhl.add(ahlField.getConc(bacterium.getPosition()));
                response.add(bacterium.getResponse());
                qSum += bacterium.getQ();
                if (bacterium.getQ() > .5) qAboveHalf++;
                cellObservations++;
            }
            populationSum += bacteria.size();
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
            summary.write("Run_Label;Stochastic_Replicate;Window;AHL_Input;"
                    + "Extracellular_AHL_Molecules_um3_Mean;Extracellular_AHL_Molecules_um3_P10;"
                    + "Extracellular_AHL_Molecules_um3_P50;Extracellular_AHL_Molecules_um3_P90;"
                    + "Extracellular_AHL_uM_Mean;Extracellular_AHL_uM_P10;"
                    + "Extracellular_AHL_uM_P50;Extracellular_AHL_uM_P90;"
                    + "Receiver_R_Mean;Receiver_R_Min;Receiver_R_P10;Receiver_R_P50;"
                    + "Receiver_R_P90;Receiver_R_Max;Mean_q;Fraction_q_gt_0_5;Population\n");
            samples = WRITE_SAMPLES
                    ? new BufferedWriter(new FileWriter(new File(directory, "results.csv"))) : null;
            if (samples != null) samples.write("Window;Sample;TimeInWindow_s;"
                    + "Input_AC0;Input_AC1_AHL;Input_AC2;Total_Count;"
                    + "Births_This_Window;Deaths_This_Window\n");
            voxels = WRITE_VOXELS
                    ? new BufferedWriter(new FileWriter(new File(directory, "voxels.csv"))) : null;
            if (voxels != null) writeVoxelHeader(voxels);
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
                throw new RuntimeException("Cannot close Stage 3B output", e);
            }
        }
    }

    public static void main(String[] args) {
        String configPath = args.length > 0 ? args[0] : "sim_config_stage3b.properties";
        loadConfig(configPath);
        if (args.length > 1 && "preview".equals(args[1])) HEADLESS = false;
        resetStaticState();
        experimentRng = new Random(RNG_SEED);
        final double[] ahlInputs = readGradedSequence(INPUT_AHL_FILE);
        if (ahlInputs.length < NUM_WINDOWS)
            throw new IllegalArgumentException("AHL input has " + ahlInputs.length
                    + " values; num.windows=" + NUM_WINDOWS);

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
        voxelAnalyzer = new VoxelAnalyzer(new int[]{STATE_X, STATE_Y, STATE_Z},
                new int[]{COUNT_X, COUNT_Y, COUNT_Z}, BOUND_X, BOUND_Y, BOUND_Z);

        for (int i = 0; i < INITIAL_POP; i++) {
            Vector3d p = new Vector3d(300 + experimentRng.nextDouble() * 400,
                    150 + experimentRng.nextDouble() * 200, BOUND_Z / 2);
            ReservoirBacterium bacterium = new ReservoirBacterium(
                    sim, p, attractantField, repellentField, ahlField);
            bacterium.setRadius();
            bacterium.setSurfaceAreaGrowthRate(GROWTH_RATE);
            bacterium.setChildList(children);
            bacteria.add(bacterium);
        }

        final int warmupSteps = (int) Math.round(WARMUP / DT);
        final int stepsInWindow = (int) Math.round(WINDOW_DURATION / DT);
        final int stepsInPulse = (int) Math.round(PULSE_DURATION / DT);
        final int samplingStart = stepsInWindow - (int) Math.round(SAMPLING_DURATION / DT);
        final int samplingInterval = Math.max(1, (int) Math.round(SAMPLING_INTERVAL / DT));
        final int[] births = {0}, deaths = {0}, completed = {0};
        final WindowAccumulator windowAccumulator = new WindowAccumulator();
        final OutputOwner output;
        try {
            output = new OutputOwner(new File(OUTPUT_DIR));
        } catch (IOException e) {
            throw new RuntimeException("Cannot open Stage 3B output", e);
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
                                windowAccumulator);
                        completed[0]++;
                    }
                    births[0] = deaths[0] = 0;
                    voxelAnalyzer.resetWindowCounters();
                    windowAccumulator.reset();
                }

                double input0 = 0.0, input2 = 0.0;
                if (warmup && step % stepsInWindow < stepsInPulse) {
                    ahlField.addQuantity(AC_POSITIONS[1],
                            AHL_SOURCE_RATE * WARMUP_AHL_INPUT * DT);
                }
                if (!warmup && window < NUM_WINDOWS && stepInWindow < stepsInPulse) {
                    attractantField.addQuantity(AC_POSITIONS[0], PROD_RATE * input0 * DT);
                    attractantField.addQuantity(AC_POSITIONS[2], PROD_RATE * input2 * DT);
                    ahlField.addQuantity(AC_POSITIONS[1],
                            AHL_SOURCE_RATE * ahlInputs[window] * DT);
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
                if (!warmup) births[0] += born;
                for (ReservoirBacterium bacterium : bacteria) {
                    if ((bacterium.getPosition().x < 0
                            || bacterium.getPosition().x > BOUND_X)
                            && !removals.contains(bacterium)) removals.add(bacterium);
                }
                for (ReservoirBacterium bacterium : removals) {
                    if (!warmup) {
                        voxelAnalyzer.recordDeath(bacterium.getPosition());
                        deaths[0]++;
                    }
                }
                cumulativeDeaths += removals.size();
                bacteria.removeAll(removals);
                removals.clear();

                if (!warmup && window < NUM_WINDOWS && stepInWindow >= samplingStart) {
                    boolean interval = (stepInWindow - samplingStart) % samplingInterval == 0;
                    boolean last = stepInWindow == stepsInWindow - 1;
                    if (interval || last) {
                        windowAccumulator.observe(ahlField);
                        if (WRITE_SAMPLES || WRITE_VOXELS)
                            writeSample(output, window, stepInWindow * DT,
                                    ahlInputs[window], births[0], deaths[0],
                                    attractantField, repellentField, ahlField);
                    }
                }
            }
        });

        sim.setDrawer(new BSimP3DDrawer(sim, 800, 600) {
            @Override
            public void scene(PGraphics3D p3d) {
                p3d.ortho(0, (float) BOUND_X, (float) BOUND_Y, 0, -1000, 10000);
                draw(ahlField, Color.ORANGE, (float) (255.0 / 540.0));
                for (ReservoirBacterium bacterium : bacteria)
                    sphere(bacterium.getPosition(), 8, Color.GREEN, 255);
            }
        });

        try {
            if (HEADLESS) sim.export();
            else {
                sim.preview();
                return;
            }
            while (completed[0] < NUM_WINDOWS) {
                finishWindow(output, completed[0], ahlInputs[completed[0]],
                        windowAccumulator);
                completed[0]++;
            }
        } finally {
            output.close();
        }
        System.out.printf(Locale.US,
                "Stage3B complete summary_rows=%d expected=%d sample_rows=%d voxel_rows=%d%n",
                output.summaryRows, NUM_WINDOWS, output.sampleRows, output.voxelRows);
    }

    static void finishWindow(OutputOwner output, int window, double input,
                             WindowAccumulator accumulator) {
        double[] ext = accumulator.extracellularAhl.stream()
                .mapToDouble(Double::doubleValue).toArray();
        double[] response = accumulator.response.stream()
                .mapToDouble(Double::doubleValue).toArray();
        Arrays.sort(ext);
        Arrays.sort(response);
        int n = response.length;
        double observations = Math.max(1, accumulator.cellObservations);
        try {
            output.summary.write(String.format(Locale.US,
                    "%s;%s;%d;%.2f;%.9g;%.9g;%.9g;%.9g;%.9g;%.9g;%.9g;%.9g;"
                            + "%.9g;%.9g;%.9g;%.9g;%.9g;%.9g;%.9g;%.9g;%d%n",
                    RUN_LABEL, STOCHASTIC_REPLICATE, window, input,
                    mean(ext), quantile(ext, .1), quantile(ext, .5), quantile(ext, .9),
                    mean(ext) / MOLECULES_PER_UM3_PER_UM,
                    quantile(ext, .1) / MOLECULES_PER_UM3_PER_UM,
                    quantile(ext, .5) / MOLECULES_PER_UM3_PER_UM,
                    quantile(ext, .9) / MOLECULES_PER_UM3_PER_UM,
                    mean(response), n == 0 ? Double.NaN : response[0],
                    quantile(response, .1), quantile(response, .5),
                    quantile(response, .9), n == 0 ? Double.NaN : response[n - 1],
                    accumulator.qSum / observations,
                    accumulator.qAboveHalf / observations,
                    (int) Math.round(accumulator.populationSum
                            / Math.max(1, accumulator.sampleCount))));
            output.summaryRows++;
        } catch (IOException e) {
            throw new RuntimeException("Cannot write window summary", e);
        }
    }

    static void writeSample(OutputOwner output, int window, double offset, double input,
                            int births, int deaths, BSimChemicalField attractant,
                            BSimChemicalField repellent, BSimChemicalField ahl) {
        try {
            int sample = (int) Math.round(offset / SAMPLING_INTERVAL);
            if (output.samples != null) {
                output.samples.write(String.format(Locale.US,
                        "%d;%d;%.2f;0;%.2f;0;%d;%d;%d%n", window, sample, offset,
                        input, bacteria.size(), births, deaths));
                output.sampleRows++;
            }
            if (output.voxels != null) {
                VoxelAnalyzer.VoxelReadout r = voxelAnalyzer.analyze(
                        attractant, repellent, ahl, bacteria);
                output.voxels.write(String.format(Locale.US,
                        "%d;%d;%.2f;0;%.2f;0", window, sample, offset, input));
                for (double v : r.attractant) output.voxels.write(format(v));
                for (double v : r.repellent) output.voxels.write(format(v));
                for (double v : r.ahl) output.voxels.write(format(v));
                for (int v : r.density) output.voxels.write(";" + v);
                for (double v : r.meanResponse) output.voxels.write(format(v));
                for (double v : r.fractionQAboveHalf) output.voxels.write(format(v));
                for (int v : r.births) output.voxels.write(";" + v);
                for (int v : r.deaths) output.voxels.write(";" + v);
                output.voxels.write("\n");
                output.voxelRows++;
            }
        } catch (IOException e) {
            throw new RuntimeException("Cannot write sample", e);
        }
    }

    static void writeVoxelHeader(BufferedWriter w) throws IOException {
        StringBuilder h = new StringBuilder(
                "Window;Sample;TimeInWindow_s;Input_AC0;Input_AC1_AHL;Input_AC2");
        int states = voxelAnalyzer.getStateVoxels();
        for (int i = 0; i < states; i++) h.append(";Att_").append(i);
        for (int i = 0; i < states; i++) h.append(";Rep_").append(i);
        for (int i = 0; i < states; i++) h.append(";AHL_").append(i);
        for (int i = 0; i < states; i++) h.append(";Den_").append(i);
        for (int i = 0; i < states; i++) h.append(";Receiver_R_").append(i);
        for (int i = 0; i < states; i++) h.append(";Fraction_q_gt_0_5_").append(i);
        for (int i = 0; i < voxelAnalyzer.getCountBins(); i++) h.append(";Birth_").append(i);
        for (int i = 0; i < voxelAnalyzer.getCountBins(); i++) h.append(";Death_").append(i);
        w.write(h.append('\n').toString());
    }

    static double[] readGradedSequence(String path) {
        List<Double> values = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty() && !line.startsWith("#"))
                    values.add(Double.parseDouble(line));
            }
        } catch (IOException | NumberFormatException e) {
            throw new IllegalArgumentException("Cannot read graded AHL input " + path, e);
        }
        double[] result = new double[values.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = values.get(i);
            if (result[i] < 0 || result[i] > 1)
                throw new IllegalArgumentException("AHL inputs must be in [0,1]");
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
        ATT_DECAY = getDouble(p, "field.att.decay", ATT_DECAY);
        REP_DECAY = getDouble(p, "field.rep.decay", REP_DECAY);
        AHL_DECAY = getDouble(p, "field.ahl.decay", AHL_DECAY);
        AHL_SOURCE_RATE = getDouble(p, "field.ahl.source.rate", AHL_SOURCE_RATE);
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
        RUN_LABEL = p.getProperty("output.run.label", RUN_LABEL);
        STOCHASTIC_REPLICATE = p.getProperty(
                "output.stochastic.replicate", STOCHASTIC_REPLICATE);
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
        bacteria.clear(); children.clear(); removals.clear();
        voxelAnalyzer = null; nextId = 0; cumulativeDeaths = 0;
    }
}
