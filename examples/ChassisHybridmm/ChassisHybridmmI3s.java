package ChassisHybridmm;

import BacteriumFromScratch.ChassisParameters;
import BacteriumFromScratch.EcoliRodCell;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;

import javax.vecmath.Vector3d;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Locale;
import java.util.Properties;
import java.util.Random;

/**
 * I3s: full-dish, motility-OFF map occupancy using the frozen I3p dynamics.
 */
public final class ChassisHybridmmI3s {
    static final int PLACEMENT_X = 60;
    static final int PLACEMENT_Y = 30;
    static final double PLACEMENT_JITTER_UM = 0.4;
    static final int STATE_X = 20;
    static final int STATE_Y = 10;
    static final int STATE_BINS = STATE_X * STATE_Y;
    static final double SAMPLE_INTERVAL_S = 20.0;
    static final int SAMPLES_PER_WINDOW = 16;
    static final int FAR_LEFT_MAX_X = 3;
    static final int FAR_RIGHT_MIN_X = 16;

    static long rngSeed = 111;
    static String outputDir = "results/i3s_silent_seed111";
    static String runLabel = "i3s_silent";
    static String stochasticReplicate = "silent_seed111";
    static Random experimentRng;
    static int initialOccupiedBins;

    private ChassisHybridmmI3s() {}

    static final class MapAccumulator {
        final boolean[] occupied = new boolean[STATE_BINS];
        int sampleCount;

        void reset() {
            for (int i = 0; i < occupied.length; i++) occupied[i] = false;
            sampleCount = 0;
        }

        int[] observe() {
            int[] density = density();
            for (int i = 0; i < density.length; i++) {
                if (density[i] > 0) occupied[i] = true;
            }
            sampleCount++;
            return density;
        }

        int occupiedCount() {
            int count = 0;
            for (boolean value : occupied) if (value) count++;
            return count;
        }

        int farOccupiedCount() {
            int count = 0;
            for (int x = 0; x < STATE_X; x++) {
                if (x > FAR_LEFT_MAX_X && x < FAR_RIGHT_MIN_X) continue;
                for (int y = 0; y < STATE_Y; y++) {
                    if (occupied[index(x, y)]) count++;
                }
            }
            return count;
        }
    }

    static final class OutputOwner implements AutoCloseable {
        final File directory;
        final BufferedWriter summary;
        final BufferedWriter mapSamples;
        final BufferedWriter mechanics;
        int summaryRows;
        int mapRows;
        int mechanicsRows;
        boolean closed;

        OutputOwner(File directory) throws IOException {
            this.directory = directory;
            if (!directory.exists() && !directory.mkdirs()) {
                throw new IOException("Cannot create " + directory);
            }
            summary = new BufferedWriter(new FileWriter(
                    new File(directory, "window_summary.csv")));
            summary.write("Run_Label;Stochastic_Replicate;Window;Population_End;"
                    + "Births;Total_Deaths;Clamp_Deaths;Acid_Deaths;OOB_Deaths;"
                    + "Occupied_Bins_Window;Far_Occupied_Bins_Window\n");
            mapSamples = new BufferedWriter(new FileWriter(
                    new File(directory, "map_samples.csv")));
            mapSamples.write("Window;Sample;TimeInWindow_s;Population;Occupied_Bins");
            for (int i = 0; i < STATE_BINS; i++) {
                mapSamples.write(";Den_" + i);
            }
            mapSamples.write("\n");
            mechanics = new BufferedWriter(new FileWriter(
                    new File(directory, "mechanics.csv")));
            mechanics.write("t_s;seed;N;offplane;nematic_order;"
                    + "delta_cc_max_um;d_centers_min_um\n");
        }

        @Override
        public void close() {
            if (closed) return;
            closed = true;
            try {
                summary.close();
                mapSamples.close();
                mechanics.close();
            } catch (IOException e) {
                throw new RuntimeException("Cannot close I3s output", e);
            }
        }
    }

    public static void main(String[] args) {
        String configPath = args.length > 0 ? args[0] : "sim_config_i3s.properties";
        loadConfig(configPath);
        validateFrozenConfig();
        ChassisHybridmmI3p.rods.clear();
        experimentRng = new Random(rngSeed);
        ChassisHybridmmI3p.experimentRng = experimentRng;

        final BSim sim = new BSim();
        sim.setDt(ChassisHybridmmI3p.DT);
        sim.setSimulationTime(ChassisHybridmmI3p.WARMUP_S
                + ChassisHybridmmI3p.NUM_WINDOWS
                * ChassisHybridmmI3p.WINDOW_DURATION_S);
        sim.setTimeFormat("0.00");
        sim.setBound(ChassisHybridmmI3p.BOUND_X,
                ChassisHybridmmI3p.BOUND_Y, ChassisHybridmmI3p.BOUND_Z);
        sim.setSolid(true, true, true);
        final BSimChemicalField acidField = new BSimChemicalField(
                sim,
                new int[]{ChassisHybridmmI3p.GRID_X,
                        ChassisHybridmmI3p.GRID_Y,
                        ChassisHybridmmI3p.GRID_Z},
                ChassisHybridmmI3p.ACID_DIFFUSIVITY,
                ChassisHybridmmI3p.ACID_DECAY);
        placeFullDishPopulation(sim);

        final int warmupSteps = ChassisHybridmmI3p.exactSteps(
                ChassisHybridmmI3p.WARMUP_S,
                ChassisHybridmmI3p.DT, "warmup");
        final int stepsInWindow = ChassisHybridmmI3p.exactSteps(
                ChassisHybridmmI3p.WINDOW_DURATION_S,
                ChassisHybridmmI3p.DT, "window");
        final int growthEvery = ChassisHybridmmI3p.exactSteps(
                ChassisParameters.DT_S,
                ChassisHybridmmI3p.DT, "growth cadence");
        final int mechanicsEvery = ChassisHybridmmI3p.exactSteps(
                ChassisHybridmmI3p.MECHANICS_LOG_INTERVAL_S,
                ChassisHybridmmI3p.DT, "mechanics log");
        final int sampleEvery = ChassisHybridmmI3p.exactSteps(
                SAMPLE_INTERVAL_S, ChassisHybridmmI3p.DT, "sample interval");
        final int[] births = {0};
        final int[] deaths = {0};
        final int[] clampDeaths = {0};
        final int[] acidDeaths = {0};
        final int[] oobDeaths = {0};
        final int[] completed = {0};
        final MapAccumulator map = new MapAccumulator();
        final OutputOwner output;
        try {
            output = new OutputOwner(new File(outputDir));
            writeMechanics(output, 0.0, sim.getBound());
        } catch (IOException e) {
            throw new RuntimeException("Cannot open I3s output", e);
        }

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                int step = (int) sim.getTimestep();
                boolean warmup = step < warmupSteps;
                int effective = step - warmupSteps;
                int window = warmup ? -1 : effective / stepsInWindow;
                int stepInWindow = warmup ? -1 : effective % stepsInWindow;

                if (!warmup && window < ChassisHybridmmI3p.NUM_WINDOWS
                        && stepInWindow == 0) {
                    if (window > 0) {
                        finishWindow(output, map, window - 1,
                                births[0], deaths[0], clampDeaths[0],
                                acidDeaths[0], oobDeaths[0]);
                        completed[0]++;
                    }
                    births[0] = deaths[0] = clampDeaths[0] = 0;
                    acidDeaths[0] = oobDeaths[0] = 0;
                    map.reset();
                }

                acidField.update();
                ChassisHybridmmI3p.updateRodChemistry(
                        acidField, warmup, deaths, clampDeaths,
                        acidDeaths, oobDeaths);
                if ((step + 1) % growthEvery == 0) {
                    ChassisHybridmmI3p.growDivideAndPack(
                            sim.getBound(), warmup, births);
                }
                if ((step + 1) % mechanicsEvery == 0) {
                    writeMechanics(output, sim.getTime()
                            + ChassisHybridmmI3p.DT, sim.getBound());
                }

                if (!warmup && window < ChassisHybridmmI3p.NUM_WINDOWS) {
                    boolean interval = stepInWindow % sampleEvery == 0;
                    boolean last = stepInWindow == stepsInWindow - 1;
                    if (interval || last) {
                        int[] density = map.observe();
                        writeMapSample(output, window,
                                stepInWindow * ChassisHybridmmI3p.DT, density);
                    }
                }
            }
        });

        PrintStream originalOut = System.out;
        try {
            System.out.printf(Locale.US,
                    "I3s start seed=%d windows=%d placement=full-dish-60x30 "
                            + "initial_bins=%d motility=OFF%n",
                    rngSeed, ChassisHybridmmI3p.NUM_WINDOWS,
                    initialOccupiedBins);
            System.setOut(new PrintStream(OutputStream.nullOutputStream()));
            sim.export();
            while (completed[0] < ChassisHybridmmI3p.NUM_WINDOWS) {
                finishWindow(output, map, completed[0],
                        births[0], deaths[0], clampDeaths[0],
                        acidDeaths[0], oobDeaths[0]);
                completed[0]++;
            }
        } finally {
            System.setOut(originalOut);
            output.close();
        }
        writeRunStatus(output);
        System.out.printf(Locale.US,
                "I3s complete seed=%d summary=%d map=%d mechanics=%d N=%d%n",
                rngSeed, output.summaryRows, output.mapRows,
                output.mechanicsRows, ChassisHybridmmI3p.rods.size());
    }

    static void placeFullDishPopulation(BSim sim) {
        double dx = ChassisHybridmmI3p.BOUND_X / PLACEMENT_X;
        double dy = ChassisHybridmmI3p.BOUND_Y / PLACEMENT_Y;
        for (int i = 0; i < ChassisHybridmmI3p.INITIAL_POP; i++) {
            int ix = i % PLACEMENT_X;
            int iy = i / PLACEMENT_X;
            double cx = (ix + 0.5) * dx
                    + (experimentRng.nextDouble() - 0.5)
                    * (2.0 * PLACEMENT_JITTER_UM);
            double cy = (iy + 0.5) * dy
                    + (experimentRng.nextDouble() - 0.5)
                    * (2.0 * PLACEMENT_JITTER_UM);
            double theta = 2.0 * Math.PI * experimentRng.nextDouble();
            double phase = ChassisHybridmmI3p.T_DIV_S
                    * experimentRng.nextDouble();
            double length = ChassisParameters.analyticLengthUm(
                    phase, ChassisParameters.NUTRIENT_BATH_MM);
            Vector3d axis = new Vector3d(
                    Math.cos(theta), Math.sin(theta), 0.0);
            Vector3d centre = new Vector3d(cx, cy, 0.5);
            Vector3d p1 = new Vector3d();
            Vector3d p2 = new Vector3d();
            p1.scaleAdd(-0.5 * length, axis, centre);
            p2.scaleAdd(0.5 * length, axis, centre);
            EcoliRodCell body = new EcoliRodCell(sim, p1, p2);
            body.grownSinceBirthS = phase;
            body.elongateCited(0.0, ChassisParameters.NUTRIENT_BATH_MM);
            body.divisionMode = EcoliRodCell.DivisionMode.SYMMETRY_BROKEN;
            body.divisionRng = experimentRng;
            ChassisHybridmmI3p.rods.add(
                    new ChassisHybridmmI3p.RodAgent(body));
        }
        initialOccupiedBins = occupiedCount(density());
        ChassisHybridmmI3p.relaxMechanics(
                ChassisHybridmmI3p.rods, sim.getBound());
    }

    static int[] density() {
        int[] result = new int[STATE_BINS];
        for (ChassisHybridmmI3p.RodAgent rod : ChassisHybridmmI3p.rods) {
            Vector3d position = rod.position();
            int x = Math.max(0, Math.min(STATE_X - 1,
                    (int) (position.x
                    / (ChassisHybridmmI3p.BOUND_X / STATE_X))));
            int y = Math.max(0, Math.min(STATE_Y - 1,
                    (int) (position.y
                    / (ChassisHybridmmI3p.BOUND_Y / STATE_Y))));
            result[index(x, y)]++;
        }
        return result;
    }

    static int index(int x, int y) {
        return x * STATE_Y + y;
    }

    static int occupiedCount(int[] density) {
        int count = 0;
        for (int value : density) if (value > 0) count++;
        return count;
    }

    static void writeMapSample(
            OutputOwner output, int window, double offset, int[] density) {
        try {
            int sample = (int) Math.round(offset / SAMPLE_INTERVAL_S);
            output.mapSamples.write(String.format(Locale.US,
                    "%d;%d;%.2f;%d;%d",
                    window, sample, offset,
                    ChassisHybridmmI3p.rods.size(),
                    occupiedCount(density)));
            for (int value : density) output.mapSamples.write(";" + value);
            output.mapSamples.write("\n");
            output.mapRows++;
        } catch (IOException e) {
            throw new RuntimeException("Cannot write I3s map sample", e);
        }
    }

    static void finishWindow(OutputOwner output, MapAccumulator map,
            int window, int births, int deaths, int clampDeaths,
            int acidDeaths, int oobDeaths) {
        if (map.sampleCount != SAMPLES_PER_WINDOW) {
            throw new IllegalStateException(
                    "I3s window " + window + " has "
                    + map.sampleCount + " samples, expected "
                    + SAMPLES_PER_WINDOW);
        }
        try {
            output.summary.write(String.format(Locale.US,
                    "%s;%s;%d;%d;%d;%d;%d;%d;%d;%d;%d%n",
                    runLabel, stochasticReplicate, window,
                    ChassisHybridmmI3p.rods.size(), births, deaths,
                    clampDeaths, acidDeaths, oobDeaths,
                    map.occupiedCount(), map.farOccupiedCount()));
            output.summaryRows++;
        } catch (IOException e) {
            throw new RuntimeException("Cannot write I3s summary", e);
        }
    }

    static void writeMechanics(
            OutputOwner output, double time, Vector3d bound) {
        ChassisHybridmmI3p.MechanicsStats stats =
                ChassisHybridmmI3p.mechanicsStats(bound);
        double minCentre = Double.isFinite(stats.minCentre)
                ? stats.minCentre
                : ChassisHybridmmI3p.MECHANICS_COMPONENT_CUTOFF_UM;
        try {
            output.mechanics.write(String.format(Locale.US,
                    "%.2f;%d;%d;%.9g;%.9g;%.9g;%.9g%n",
                    time, rngSeed, ChassisHybridmmI3p.rods.size(),
                    stats.offplane, stats.nematic,
                    stats.maxDelta, minCentre));
            output.mechanicsRows++;
        } catch (IOException e) {
            throw new RuntimeException("Cannot write I3s mechanics", e);
        }
    }

    static void loadConfig(String path) {
        Properties properties = new Properties();
        try (FileInputStream input = new FileInputStream(path)) {
            properties.load(input);
        } catch (IOException e) {
            throw new IllegalArgumentException(
                    "Cannot load I3s config " + path, e);
        }
        rngSeed = Long.parseLong(properties.getProperty(
                "rng.seed", String.valueOf(rngSeed)));
        outputDir = properties.getProperty("output.dir", outputDir);
        runLabel = properties.getProperty("output.run.label", runLabel);
        stochasticReplicate = properties.getProperty(
                "output.stochastic.replicate", stochasticReplicate);
        for (String key : properties.stringPropertyNames()) {
            if (!"rng.seed".equals(key)
                    && !"output.dir".equals(key)
                    && !"output.run.label".equals(key)
                    && !"output.stochastic.replicate".equals(key)) {
                throw new IllegalArgumentException(
                        "I3s config may not change frozen key " + key);
            }
        }
    }

    static void validateFrozenConfig() {
        if (PLACEMENT_X * PLACEMENT_Y
                != ChassisHybridmmI3p.INITIAL_POP
                || STATE_X != 20 || STATE_Y != 10
                || FAR_LEFT_MAX_X != 3 || FAR_RIGHT_MIN_X != 16
                || ChassisHybridmmI3p.T_REMOVAL_S
                != ChassisHybridmmI3p.T_DIV_S) {
            throw new IllegalArgumentException(
                    "I3s frozen placement or I3p clock changed");
        }
    }

    static void writeRunStatus(OutputOwner output) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(
                new File(output.directory, "run_status.txt")))) {
            writer.write(String.format(Locale.US,
                    "job=ChassisHybridmm_I3s_map_occupancy%n"
                            + "arm=silent%nseed=%d%n"
                            + "placement=full-dish-60x30%n"
                            + "placement_jitter_um=%.12g%n"
                            + "initial_occupied_bins=%d%n"
                            + "summary_rows=%d%nmap_rows=%d%n"
                            + "mechanics_rows=%d%nT_div_s=%.12g%n"
                            + "T_removal_s=%.12g%n"
                            + "motility=OFF%nridge=OFF%n"
                            + "NARMA_overall=NONE%n",
                    rngSeed, PLACEMENT_JITTER_UM, initialOccupiedBins,
                    output.summaryRows, output.mapRows,
                    output.mechanicsRows, ChassisHybridmmI3p.T_DIV_S,
                    ChassisHybridmmI3p.T_REMOVAL_S));
        } catch (IOException e) {
            throw new RuntimeException("Cannot write I3s status", e);
        }
    }
}
