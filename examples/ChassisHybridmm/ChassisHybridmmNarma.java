package ChassisHybridmm;

import BacteriumFromScratch.ChassisParameters;
import BacteriumFromScratch.EcoliRodCell;
import BacteriumFromScratch.ValdezHertzian;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.particle.BSimParticle;

import javax.vecmath.Vector3d;
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
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Properties;
import java.util.Random;

/**
 * Job I3n: a new millimetre ChassisHybridmm dish.
 *
 * <p>Hertzian E. coli rods are imported from BacteriumFromScratch. This file
 * adds only the frozen HybridDish AHL -> Hill R -> slow L task layer and I/O.
 * It does not modify or inherit the HybridDish claim Java.
 */
public final class ChassisHybridmmNarma {
    enum Arm { BROWNIAN, SILENT, DRIVEN }
    enum DeathCause { CLAMP, ACID, OOB }

    static final double BOUND_X = 1000.0;
    static final double BOUND_Y = 500.0;
    static final double BOUND_Z = 1.0;
    static final double FLOW_SPEED = 0.0;
    static final double RECEIVER_K_UM = 1.6;
    static final double RECEIVER_HILL_N = 2.0;
    static final double RECEIVER_TAU_S = 15.0;
    static final double ALPHA_L = 1.0 / 1500.0;
    static final double DELTA_L = 1.0 / 1500.0;
    static final double MOLECULES_PER_UM3_PER_UM = 602.2;
    static final double MM_TO_MOLECULES_PER_UM3 = 602200.0;
    static final double T_REMOVAL_S = 1800.0;
    static final double PH_BASE = 7.1;
    static final double BUFFER_CAPACITY_MM_PER_PH = 2.0;
    static final double PH_GROWTH_LIMIT = 4.5;
    static final double PH_ACID_EC50 = 3.75;
    static final double PH_ALK_ONSET = 9.0;
    static final double PH_ALK_EC50 = 9.85;
    static final double N_ACID = 2.0;
    static final double N_ALK = 2.0;
    static final double K_MAX_ALK = 0.003;
    static final double MECHANICS_COMPONENT_CUTOFF_UM = 4.5;
    static final int MECHANICS_REBUILDS = 2;
    static final double MECHANICS_LOG_INTERVAL_S = 10.0;

    static double DT = 0.05;
    static int GRID_X = 50, GRID_Y = 25, GRID_Z = 1;
    static int STATE_X = 20, STATE_Y = 10, STATE_Z = 1;
    static int COUNT_X = 4, COUNT_Y = 2, COUNT_Z = 1;
    static double AHL_DIFFUSIVITY = 159.0, AHL_DECAY = 0.0033;
    static double ACID_DIFFUSIVITY = 200.0, ACID_DECAY = 0.0067;
    static double AHL_SOURCE_RATE = 1.28e7;
    static double ACID_SOURCE_RATE = 2.0e10;
    static double ACID_CELL_PRODUCTION_RATE = 1.0e5;
    static double K_MAX_ACID = 0.002;
    static double WARMUP = 18000.0, WINDOW_DURATION = 300.0, PULSE_DURATION = 75.0;
    static double SAMPLING_DURATION = 300.0, SAMPLING_INTERVAL = 20.0;
    static double WARMUP_AHL_INPUT = 0.5;
    static int NUM_WINDOWS = 200, INITIAL_POP = 1800, CARRYING_CAPACITY = 2000;
    static long RNG_SEED = 111;
    static Arm ARM = Arm.DRIVEN;
    static String INPUT_AHL_FILE = "input_ahl_narma200.txt";
    static String INPUT_ACID_FILE = "input_acid_held05_200.txt";
    static String OUTPUT_DIR = "results/i3n_driven_seed111";
    static String RUN_LABEL = "i3n_driven";
    static String STOCHASTIC_REPLICATE = "driven_seed111";

    static final Vector3d AHL_SOURCE_POSITION = new Vector3d(500, 250, 0.5);
    static final Vector3d ACID_SOURCE_POSITION = new Vector3d(300, 375, 0.5);
    static final List<RodAgent> rods = new ArrayList<RodAgent>();
    static final List<BrownianParticle> particles = new ArrayList<BrownianParticle>();
    static Random experimentRng;
    static VoxelAnalyzer voxelAnalyzer;

    private ChassisHybridmmNarma() {}

    static final class RodAgent {
        final EcoliRodCell body;
        double receiver;
        double luminescence;

        RodAgent(EcoliRodCell body) {
            this.body = body;
        }

        RodAgent daughter() {
            RodAgent child = new RodAgent(body.divideCited());
            child.receiver = receiver;
            child.luminescence = luminescence;
            return child;
        }

        Vector3d position() {
            return body.centre();
        }
    }

    /** Narma10b-style passive density null. It has no biological state. */
    static final class BrownianParticle extends BSimParticle {
        BrownianParticle(BSim sim, Vector3d position) {
            super(sim, position, ChassisParameters.RADIUS_UM);
        }

        static void seedMotionRng(long seed) {
            rng.setSeed(seed);
        }
    }

    static final class WindowAccumulator {
        double rSum, lSum, populationSum;
        long observations;
        int samples;

        void reset() {
            rSum = lSum = populationSum = 0.0;
            observations = 0;
            samples = 0;
        }

        void observe() {
            if (ARM == Arm.BROWNIAN) {
                populationSum += particles.size();
            } else {
                for (RodAgent rod : rods) {
                    rSum += rod.receiver;
                    lSum += rod.luminescence;
                    observations++;
                }
                populationSum += rods.size();
            }
            samples++;
        }
    }

    static final class MechanicsStats {
        double maxDelta;
        double minCentre = Double.POSITIVE_INFINITY;
        double offplane;
        double nematic;
    }

    static final class OutputOwner implements AutoCloseable {
        final BufferedWriter summary;
        final BufferedWriter samples;
        final BufferedWriter voxels;
        final BufferedWriter mechanics;
        int summaryRows, sampleRows, voxelRows, mechanicsRows;
        boolean closed;

        OutputOwner(File directory) throws IOException {
            if (!directory.exists() && !directory.mkdirs()) {
                throw new IOException("Cannot create " + directory);
            }
            summary = new BufferedWriter(new FileWriter(new File(directory, "window_summary.csv")));
            summary.write("Run_Label;Stochastic_Replicate;Arm;Window;AHL_Input;Acid_Input;"
                    + "Mean_R;Mean_L;Population;Births;Total_Deaths;Clamp_Deaths;"
                    + "Input_Driven_Deaths;OOB_Deaths\n");
            samples = new BufferedWriter(new FileWriter(new File(directory, "results.csv")));
            samples.write("Window;Sample;TimeInWindow_s;Input_AC1_AHL;Input_Acid;"
                    + "Total_Count;Births_This_Window;Deaths_This_Window;"
                    + "Clamp_Deaths_This_Window;Input_Driven_Deaths_This_Window;"
                    + "OOB_Deaths_This_Window;Mean_R;Mean_L\n");
            voxels = new BufferedWriter(new FileWriter(new File(directory, "voxels.csv")));
            writeVoxelHeader(voxels);
            mechanics = new BufferedWriter(new FileWriter(new File(directory, "mechanics.csv")));
            mechanics.write("t_s;seed;arm;N;mean_R;mean_L;offplane;nematic_order;"
                    + "delta_cc_max_um;d_centers_min_um\n");
            writeFeatureContract(directory);
        }

        @Override
        public void close() {
            if (closed) return;
            closed = true;
            try {
                summary.close();
                samples.close();
                voxels.close();
                mechanics.close();
            } catch (IOException e) {
                throw new RuntimeException("Cannot close I3n output", e);
            }
        }
    }

    public static void main(String[] args) {
        String configPath = args.length > 0 ? args[0] : "sim_config_i3n.properties";
        loadConfig(configPath);
        resetStaticState();
        experimentRng = new Random(RNG_SEED);
        BrownianParticle.seedMotionRng(RNG_SEED);
        double[] ahlInputs = readSequence(INPUT_AHL_FILE, "AHL");
        double[] acidInputs = readSequence(INPUT_ACID_FILE, "acid");
        if (ahlInputs.length < NUM_WINDOWS || acidInputs.length < NUM_WINDOWS) {
            throw new IllegalArgumentException("Input sidecars shorter than num.windows=" + NUM_WINDOWS);
        }

        final BSim sim = new BSim();
        sim.setDt(DT);
        sim.setSimulationTime(WARMUP + NUM_WINDOWS * WINDOW_DURATION);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);
        int[] fieldGrid = {GRID_X, GRID_Y, GRID_Z};
        final BSimChemicalField ahlField = new BSimChemicalField(
                sim, fieldGrid, AHL_DIFFUSIVITY, AHL_DECAY);
        final BSimChemicalField acidField = new BSimChemicalField(
                sim, fieldGrid, ACID_DIFFUSIVITY, ACID_DECAY);
        voxelAnalyzer = new VoxelAnalyzer(
                new int[]{STATE_X, STATE_Y, STATE_Z},
                new int[]{COUNT_X, COUNT_Y, COUNT_Z});
        placeInitialPopulation(sim);

        final int warmupSteps = exactSteps(WARMUP, DT, "warmup");
        final int stepsInWindow = exactSteps(WINDOW_DURATION, DT, "window");
        final int stepsInPulse = exactSteps(PULSE_DURATION, DT, "pulse");
        final int growthEvery = exactSteps(ChassisParameters.DT_S, DT, "growth cadence");
        final int mechanicsLogEvery = exactSteps(MECHANICS_LOG_INTERVAL_S, DT, "mechanics log");
        final int samplingStart = stepsInWindow - exactSteps(SAMPLING_DURATION, DT, "sampling duration");
        final int sampleEvery = exactSteps(SAMPLING_INTERVAL, DT, "sampling interval");
        final int[] births = {0}, deaths = {0}, clampDeaths = {0};
        final int[] acidDeaths = {0}, oobDeaths = {0}, completed = {0};
        final WindowAccumulator accumulator = new WindowAccumulator();
        final OutputOwner output;
        try {
            output = new OutputOwner(new File(OUTPUT_DIR));
            writeMechanics(output, 0.0, sim.getBound());
        } catch (IOException e) {
            throw new RuntimeException("Cannot open I3n output", e);
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
                                acidInputs[window - 1], births[0], deaths[0],
                                clampDeaths[0], acidDeaths[0], oobDeaths[0], accumulator);
                        completed[0]++;
                    }
                    births[0] = deaths[0] = clampDeaths[0] = acidDeaths[0] = oobDeaths[0] = 0;
                    voxelAnalyzer.resetWindowCounters();
                    accumulator.reset();
                }

                if (warmup && step % stepsInWindow < stepsInPulse) {
                    ahlField.addQuantity(AHL_SOURCE_POSITION,
                            AHL_SOURCE_RATE * WARMUP_AHL_INPUT * DT);
                } else if (!warmup && window < NUM_WINDOWS && stepInWindow < stepsInPulse) {
                    ahlField.addQuantity(AHL_SOURCE_POSITION,
                            AHL_SOURCE_RATE * ahlInputs[window] * DT);
                    acidField.addQuantity(ACID_SOURCE_POSITION,
                            ACID_SOURCE_RATE * acidInputs[window] * DT);
                }
                ahlField.update();
                acidField.update();

                if (ARM == Arm.BROWNIAN) {
                    for (BrownianParticle particle : particles) {
                        particle.action();
                        particle.updatePosition();
                    }
                } else {
                    updateRodChemistry(ahlField, acidField, warmup, deaths,
                            clampDeaths, acidDeaths, oobDeaths);
                    if ((step + 1) % growthEvery == 0) {
                        growDivideAndPack(sim.getBound(), warmup, births);
                    }
                }

                if ((step + 1) % mechanicsLogEvery == 0) {
                    writeMechanics(output, sim.getTime() + DT, sim.getBound());
                }

                if (!warmup && window < NUM_WINDOWS && stepInWindow >= samplingStart) {
                    boolean interval = (stepInWindow - samplingStart) % sampleEvery == 0;
                    boolean last = stepInWindow == stepsInWindow - 1;
                    if (interval || last) {
                        accumulator.observe();
                        writeSample(output, window, stepInWindow * DT,
                                ahlInputs[window], acidInputs[window], births[0], deaths[0],
                                clampDeaths[0], acidDeaths[0], oobDeaths[0], ahlField);
                    }
                }
            }
        });

        PrintStream originalOut = System.out;
        try {
            System.out.printf(Locale.US,
                    "I3n start arm=%s seed=%d windows=%d domain=1000x500x1 "
                            + "J_max=1.28e7 rods=%s%n",
                    ARM.name().toLowerCase(Locale.ROOT), RNG_SEED, NUM_WINDOWS,
                    ARM == Arm.BROWNIAN ? "density-null" : "Hertzian");
            System.setOut(new PrintStream(OutputStream.nullOutputStream()));
            sim.export();
            while (completed[0] < NUM_WINDOWS) {
                finishWindow(output, completed[0], ahlInputs[completed[0]],
                        acidInputs[completed[0]], births[0], deaths[0],
                        clampDeaths[0], acidDeaths[0], oobDeaths[0], accumulator);
                completed[0]++;
            }
        } finally {
            System.setOut(originalOut);
            output.close();
        }
        writeRunStatus(output);
        System.out.printf(Locale.US,
                "I3n complete arm=%s seed=%d summary=%d samples=%d voxels=%d mechanics=%d%n",
                ARM.name().toLowerCase(Locale.ROOT), RNG_SEED, output.summaryRows,
                output.sampleRows, output.voxelRows, output.mechanicsRows);
    }

    static void placeInitialPopulation(BSim sim) {
        final int nx = 60, ny = 30;
        double dx = 400.0 / nx, dy = 200.0 / ny;
        double divisionTime = ChassisParameters.analyticDivisionTimeS(
                ChassisParameters.NUTRIENT_BATH_MM);
        for (int i = 0; i < INITIAL_POP; i++) {
            int ix = i % nx, iy = i / nx;
            double cx = 300.0 + (ix + 0.5) * dx + (experimentRng.nextDouble() - 0.5) * 0.8;
            double cy = 150.0 + (iy + 0.5) * dy + (experimentRng.nextDouble() - 0.5) * 0.8;
            if (ARM == Arm.BROWNIAN) {
                particles.add(new BrownianParticle(sim, new Vector3d(cx, cy, 0.5)));
                continue;
            }
            double theta = 2.0 * Math.PI * experimentRng.nextDouble();
            double phase = divisionTime * experimentRng.nextDouble();
            double length = ChassisParameters.analyticLengthUm(
                    phase, ChassisParameters.NUTRIENT_BATH_MM);
            Vector3d axis = new Vector3d(Math.cos(theta), Math.sin(theta), 0.0);
            Vector3d centre = new Vector3d(cx, cy, 0.5);
            Vector3d p1 = new Vector3d(), p2 = new Vector3d();
            p1.scaleAdd(-0.5 * length, axis, centre);
            p2.scaleAdd(0.5 * length, axis, centre);
            EcoliRodCell body = new EcoliRodCell(sim, p1, p2);
            body.grownSinceBirthS = phase;
            body.elongateCited(0.0, ChassisParameters.NUTRIENT_BATH_MM);
            body.divisionMode = EcoliRodCell.DivisionMode.SYMMETRY_BROKEN;
            body.divisionRng = experimentRng;
            rods.add(new RodAgent(body));
        }
        if (ARM != Arm.BROWNIAN) {
            relaxMechanics(rods, sim.getBound());
        }
    }

    static void updateRodChemistry(BSimChemicalField ahl, BSimChemicalField acid,
            boolean warmup, int[] deaths, int[] clampDeaths, int[] acidDeaths,
            int[] oobDeaths) {
        Iterator<RodAgent> iterator = rods.iterator();
        while (iterator.hasNext()) {
            RodAgent rod = iterator.next();
            Vector3d p = rod.position();
            double ahlUm = ahl.getConc(p) / MOLECULES_PER_UM3_PER_UM;
            double cn = Math.pow(ahlUm, RECEIVER_HILL_N);
            double target = cn / (Math.pow(RECEIVER_K_UM, RECEIVER_HILL_N) + cn);
            rod.receiver = target + (rod.receiver - target) * Math.exp(-DT / RECEIVER_TAU_S);
            rod.luminescence += (ALPHA_L * rod.receiver - DELTA_L * rod.luminescence) * DT;
            acid.addQuantity(p, ACID_CELL_PRODUCTION_RATE * DT);

            double pRemoval = (DT / T_REMOVAL_S)
                    * Math.pow(2.0, -(1.0 - rods.size() / (double) CARRYING_CAPACITY));
            boolean clamp = experimentRng.nextDouble() < pRemoval;
            double pH = pHFromAcidMm(acid.getConc(p) / MM_TO_MOLECULES_PER_UM3);
            boolean acidDeath = experimentRng.nextDouble() < killRate(pH) * DT;
            DeathCause cause = acidDeath ? DeathCause.ACID : (clamp ? DeathCause.CLAMP : null);
            if (p.x < 0 || p.x > BOUND_X || p.y < 0 || p.y > BOUND_Y
                    || p.z < 0 || p.z > BOUND_Z) {
                cause = DeathCause.OOB;
            }
            if (cause != null) {
                if (!warmup) {
                    voxelAnalyzer.recordDeath(p, cause);
                    deaths[0]++;
                    if (cause == DeathCause.CLAMP) clampDeaths[0]++;
                    else if (cause == DeathCause.ACID) acidDeaths[0]++;
                    else oobDeaths[0]++;
                }
                iterator.remove();
            }
        }
    }

    static void growDivideAndPack(Vector3d bound, boolean warmup, int[] births) {
        List<RodAgent> newborns = new ArrayList<RodAgent>();
        for (RodAgent rod : rods) {
            rod.body.elongateCited(ChassisParameters.DT_S,
                    ChassisParameters.NUTRIENT_BATH_MM);
            if (rod.body.shouldDivide()) {
                RodAgent child = rod.daughter();
                newborns.add(child);
                if (!warmup) {
                    voxelAnalyzer.recordBirth(child.position());
                    births[0]++;
                }
            }
        }
        rods.addAll(newborns);
        relaxMechanics(rods, bound);
    }

    static void relaxMechanics(List<RodAgent> population, Vector3d bound) {
        for (int pass = 0; pass < MECHANICS_REBUILDS; pass++) {
            for (List<RodAgent> component : contactComponents(population)) {
                List<EcoliRodCell> bodies = new ArrayList<EcoliRodCell>(component.size());
                for (RodAgent rod : component) bodies.add(rod.body);
                ValdezHertzian.relaxContacts(bodies, bound);
            }
        }
    }

    /** Conservative connected components; exact Valdez forces remain unchanged. */
    static List<List<RodAgent>> contactComponents(List<RodAgent> population) {
        int n = population.size();
        int[] parent = new int[n], rank = new int[n];
        for (int i = 0; i < n; i++) parent[i] = i;
        Map<Long, List<Integer>> bins = new HashMap<Long, List<Integer>>();
        for (int i = 0; i < n; i++) {
            Vector3d c = population.get(i).position();
            int bx = (int) Math.floor(c.x / MECHANICS_COMPONENT_CUTOFF_UM);
            int by = (int) Math.floor(c.y / MECHANICS_COMPONENT_CUTOFF_UM);
            int bz = (int) Math.floor(c.z / MECHANICS_COMPONENT_CUTOFF_UM);
            for (int ox = -1; ox <= 1; ox++) {
                for (int oy = -1; oy <= 1; oy++) {
                    for (int oz = -1; oz <= 1; oz++) {
                        List<Integer> candidates = bins.get(binKey(bx + ox, by + oy, bz + oz));
                        if (candidates == null) continue;
                        for (int j : candidates) {
                            Vector3d d = new Vector3d(c);
                            d.sub(population.get(j).position());
                            if (d.lengthSquared() <= MECHANICS_COMPONENT_CUTOFF_UM
                                    * MECHANICS_COMPONENT_CUTOFF_UM) {
                                union(parent, rank, i, j);
                            }
                        }
                    }
                }
            }
            long key = binKey(bx, by, bz);
            if (!bins.containsKey(key)) bins.put(key, new ArrayList<Integer>());
            bins.get(key).add(i);
        }
        Map<Integer, List<RodAgent>> grouped = new LinkedHashMap<Integer, List<RodAgent>>();
        for (int i = 0; i < n; i++) {
            int root = find(parent, i);
            if (!grouped.containsKey(root)) grouped.put(root, new ArrayList<RodAgent>());
            grouped.get(root).add(population.get(i));
        }
        return new ArrayList<List<RodAgent>>(grouped.values());
    }

    static long binKey(int x, int y, int z) {
        return ((long) (x + 2048) << 32) ^ ((long) (y + 2048) << 16) ^ (z + 2048L);
    }

    static int find(int[] parent, int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    static void union(int[] parent, int[] rank, int a, int b) {
        int ra = find(parent, a), rb = find(parent, b);
        if (ra == rb) return;
        if (rank[ra] < rank[rb]) parent[ra] = rb;
        else if (rank[ra] > rank[rb]) parent[rb] = ra;
        else {
            parent[rb] = ra;
            rank[ra]++;
        }
    }

    static MechanicsStats mechanicsStats(Vector3d bound) {
        MechanicsStats result = new MechanicsStats();
        if (rods.isEmpty()) return result;
        int off = 0;
        double c2 = 0.0, s2 = 0.0;
        for (RodAgent rod : rods) {
            Vector3d centre = rod.position();
            if (Math.abs(centre.z - 0.5) > 0.5 * ChassisParameters.RADIUS_UM) off++;
            Vector3d axis = new Vector3d();
            axis.sub(rod.body.x2, rod.body.x1);
            double theta = Math.atan2(axis.y, axis.x);
            c2 += Math.cos(2.0 * theta);
            s2 += Math.sin(2.0 * theta);
        }
        result.offplane = off / (double) rods.size();
        result.nematic = Math.hypot(c2 / rods.size(), s2 / rods.size());
        for (List<RodAgent> component : contactComponents(rods)) {
            List<EcoliRodCell> bodies = new ArrayList<EcoliRodCell>(component.size());
            for (RodAgent rod : component) bodies.add(rod.body);
            ValdezHertzian.PackingStats stats = ValdezHertzian.stats(bodies, bound);
            result.maxDelta = Math.max(result.maxDelta, stats.maxDeltaCc);
            result.minCentre = Math.min(result.minCentre, stats.minCentreDist);
        }
        return result;
    }

    static void writeMechanics(OutputOwner output, double time, Vector3d bound) {
        try {
            MechanicsStats stats = ARM == Arm.BROWNIAN ? new MechanicsStats() : mechanicsStats(bound);
            double r = 0.0, l = 0.0;
            for (RodAgent rod : rods) {
                r += rod.receiver;
                l += rod.luminescence;
            }
            int n = ARM == Arm.BROWNIAN ? particles.size() : rods.size();
            if (!rods.isEmpty()) {
                r /= rods.size();
                l /= rods.size();
            }
            output.mechanics.write(String.format(Locale.US,
                    "%.2f;%d;%s;%d;%.9g;%.9g;%.9g;%.9g;%.9g;%.9g%n",
                    time, RNG_SEED, ARM.name().toLowerCase(Locale.ROOT), n, r, l,
                    stats.offplane, stats.nematic, stats.maxDelta, stats.minCentre));
            output.mechanicsRows++;
        } catch (IOException e) {
            throw new RuntimeException("Cannot write mechanics", e);
        }
    }

    static void finishWindow(OutputOwner output, int window, double ahlInput,
            double acidInput, int births, int deaths, int clampDeaths,
            int acidDeaths, int oobDeaths, WindowAccumulator accumulator) {
        double observations = Math.max(1.0, accumulator.observations);
        int population = (int) Math.round(accumulator.populationSum
                / Math.max(1, accumulator.samples));
        try {
            output.summary.write(String.format(Locale.US,
                    "%s;%s;%s;%d;%.12g;%.12g;%.9g;%.9g;%d;%d;%d;%d;%d;%d%n",
                    RUN_LABEL, STOCHASTIC_REPLICATE, ARM.name().toLowerCase(Locale.ROOT),
                    window, ahlInput, acidInput,
                    ARM == Arm.BROWNIAN ? 0.0 : accumulator.rSum / observations,
                    ARM == Arm.BROWNIAN ? 0.0 : accumulator.lSum / observations,
                    population, births, deaths, clampDeaths, acidDeaths, oobDeaths));
            output.summaryRows++;
        } catch (IOException e) {
            throw new RuntimeException("Cannot write window summary", e);
        }
    }

    static void writeSample(OutputOwner output, int window, double offset,
            double ahlInput, double acidInput, int births, int deaths,
            int clampDeaths, int acidDeaths, int oobDeaths,
            BSimChemicalField ahlField) {
        try {
            int sample = (int) Math.round(offset / SAMPLING_INTERVAL);
            int count = ARM == Arm.BROWNIAN ? particles.size() : rods.size();
            double meanR = 0.0, meanL = 0.0;
            for (RodAgent rod : rods) {
                meanR += rod.receiver;
                meanL += rod.luminescence;
            }
            if (!rods.isEmpty()) {
                meanR /= rods.size();
                meanL /= rods.size();
            }
            output.samples.write(String.format(Locale.US,
                    "%d;%d;%.2f;%.12g;%.12g;%d;%d;%d;%d;%d;%d;%.9g;%.9g%n",
                    window, sample, offset, ahlInput, acidInput, count, births,
                    deaths, clampDeaths, acidDeaths, oobDeaths, meanR, meanL));
            output.sampleRows++;

            output.voxels.write(String.format(Locale.US,
                    "%d;%d;%.2f;%.12g;%.12g", window, sample, offset, ahlInput, acidInput));
            if (ARM == Arm.BROWNIAN) {
                for (int value : voxelAnalyzer.densityParticles(particles)) {
                    output.voxels.write(";" + value);
                }
            } else {
                VoxelReadout readout = voxelAnalyzer.analyze(ahlField, rods);
                for (double value : readout.ahlUm) output.voxels.write(format(value));
                for (int value : readout.density) output.voxels.write(";" + value);
                for (double value : readout.meanR) output.voxels.write(format(value));
                for (double value : readout.meanL) output.voxels.write(format(value));
                for (int value : readout.acidDeaths) output.voxels.write(";" + value);
            }
            output.voxels.write("\n");
            output.voxelRows++;
        } catch (IOException e) {
            throw new RuntimeException("Cannot write I3n sample", e);
        }
    }

    static void writeVoxelHeader(BufferedWriter writer) throws IOException {
        StringBuilder h = new StringBuilder(
                "Window;Sample;TimeInWindow_s;Input_AC1_AHL;Input_Acid");
        int states = STATE_X * STATE_Y * STATE_Z;
        if (ARM == Arm.BROWNIAN) {
            for (int i = 0; i < states; i++) h.append(";Den_").append(i);
        } else {
            for (int i = 0; i < states; i++) h.append(";AHL_uM_").append(i);
            for (int i = 0; i < states; i++) h.append(";Den_").append(i);
            for (int i = 0; i < states; i++) h.append(";Receiver_R_").append(i);
            for (int i = 0; i < states; i++) h.append(";Lum_Mean_").append(i);
            for (int i = 0; i < COUNT_X * COUNT_Y * COUNT_Z; i++) {
                h.append(";Input_Driven_Death_").append(i);
            }
        }
        writer.write(h.append('\n').toString());
    }

    static void writeFeatureContract(File directory) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(
                new FileWriter(new File(directory, "feature_contract.txt")))) {
            writer.write("dish=ChassisHybridmm I3n (new dish; not HybridDish Narma10b)\n");
            writer.write("arm=" + ARM.name().toLowerCase(Locale.ROOT) + "\n");
            writer.write("body=imported EcoliRodCell; Job 3 ValdezHertzian; b_z=1 um\n");
            writer.write("ahl_clocks=HybridDish D=159 um2/s, k=0.0033 /s\n");
            writer.write("source=ENGINEERING volume match J_max=1.28e7 molecules/s\n");
            if (ARM == Arm.BROWNIAN) {
                writer.write("ridge_features=Den (20x10) only\n");
            } else {
                writer.write("ridge_features=Receiver_R (20x10); Lum_Mean (20x10); "
                        + "Input_Driven_Death (4x2) = 408\n");
                writer.write("excluded=AHL_uM and Den from biology ridge\n");
            }
        }
    }

    static final class VoxelReadout {
        final double[] ahlUm, meanR, meanL;
        final int[] density, acidDeaths;

        VoxelReadout(double[] ahlUm, int[] density, double[] meanR,
                double[] meanL, int[] acidDeaths) {
            this.ahlUm = ahlUm;
            this.density = density;
            this.meanR = meanR;
            this.meanL = meanL;
            this.acidDeaths = acidDeaths;
        }
    }

    static final class VoxelAnalyzer {
        final int[] stateGrid, countGrid;
        final int[] acidDeaths;

        VoxelAnalyzer(int[] stateGrid, int[] countGrid) {
            this.stateGrid = stateGrid.clone();
            this.countGrid = countGrid.clone();
            acidDeaths = new int[countGrid[0] * countGrid[1] * countGrid[2]];
        }

        void resetWindowCounters() {
            Arrays.fill(acidDeaths, 0);
        }

        void recordBirth(Vector3d position) {
            // Births are logged in summaries but are not ridge features.
        }

        void recordDeath(Vector3d position, DeathCause cause) {
            if (cause == DeathCause.ACID) acidDeaths[countIndex(position)]++;
        }

        VoxelReadout analyze(BSimChemicalField ahl, List<RodAgent> population) {
            int total = stateGrid[0] * stateGrid[1] * stateGrid[2];
            double[] ahlUm = new double[total], rSum = new double[total];
            double[] lSum = new double[total], meanR = new double[total];
            double[] meanL = new double[total];
            int[] density = new int[total];
            double vx = BOUND_X / stateGrid[0], vy = BOUND_Y / stateGrid[1];
            double vz = BOUND_Z / stateGrid[2];
            for (int x = 0; x < stateGrid[0]; x++) {
                for (int y = 0; y < stateGrid[1]; y++) {
                    for (int z = 0; z < stateGrid[2]; z++) {
                        int id = index(x, y, z, stateGrid);
                        Vector3d c = new Vector3d((x + 0.5) * vx,
                                (y + 0.5) * vy, (z + 0.5) * vz);
                        ahlUm[id] = ahl.getConc(c) / MOLECULES_PER_UM3_PER_UM;
                    }
                }
            }
            for (RodAgent rod : population) {
                int id = stateIndex(rod.position());
                density[id]++;
                rSum[id] += rod.receiver;
                lSum[id] += rod.luminescence;
            }
            for (int i = 0; i < total; i++) {
                if (density[i] > 0) {
                    meanR[i] = rSum[i] / density[i];
                    meanL[i] = lSum[i] / density[i];
                }
            }
            return new VoxelReadout(ahlUm, density, meanR, meanL, acidDeaths.clone());
        }

        int[] densityParticles(List<BrownianParticle> population) {
            int[] density = new int[stateGrid[0] * stateGrid[1] * stateGrid[2]];
            for (BrownianParticle particle : population) density[stateIndex(particle.getPosition())]++;
            return density;
        }

        int stateIndex(Vector3d p) {
            return index(clamp((int) (p.x / (BOUND_X / stateGrid[0])), 0, stateGrid[0] - 1),
                    clamp((int) (p.y / (BOUND_Y / stateGrid[1])), 0, stateGrid[1] - 1),
                    clamp((int) (p.z / (BOUND_Z / stateGrid[2])), 0, stateGrid[2] - 1),
                    stateGrid);
        }

        int countIndex(Vector3d p) {
            return index(clamp((int) (p.x / (BOUND_X / countGrid[0])), 0, countGrid[0] - 1),
                    clamp((int) (p.y / (BOUND_Y / countGrid[1])), 0, countGrid[1] - 1),
                    clamp((int) (p.z / (BOUND_Z / countGrid[2])), 0, countGrid[2] - 1),
                    countGrid);
        }
    }

    static int index(int x, int y, int z, int[] grid) {
        return x * grid[1] * grid[2] + y * grid[2] + z;
    }

    static int clamp(int value, int low, int high) {
        return Math.max(low, Math.min(high, value));
    }

    static double pHFromAcidMm(double acidMm) {
        return Math.max(1.0, Math.min(14.0,
                PH_BASE - acidMm / BUFFER_CAPACITY_MM_PER_PH));
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

    static int exactSteps(double duration, double dt, String name) {
        int steps = (int) Math.round(duration / dt);
        if (Math.abs(steps * dt - duration) > 1e-9) {
            throw new IllegalArgumentException(name + " is not an integer number of field ticks");
        }
        return steps;
    }

    static double[] readSequence(String path, String name) {
        List<Double> values = new ArrayList<Double>();
        try (BufferedReader reader = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty() && !line.startsWith("#")) {
                    double value = Double.parseDouble(line);
                    if (value < 0.0 || value > 0.5) {
                        throw new IllegalArgumentException(name + " escaped frozen [0,0.5]");
                    }
                    values.add(value);
                }
            }
        } catch (IOException | NumberFormatException e) {
            throw new IllegalArgumentException("Cannot read " + name + " input " + path, e);
        }
        double[] result = new double[values.size()];
        for (int i = 0; i < result.length; i++) result[i] = values.get(i);
        return result;
    }

    static void loadConfig(String path) {
        Properties p = new Properties();
        try (FileInputStream input = new FileInputStream(path)) {
            p.load(input);
        } catch (IOException e) {
            throw new IllegalArgumentException("Cannot load config " + path, e);
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
        AHL_DIFFUSIVITY = getDouble(p, "field.ahl.diff", AHL_DIFFUSIVITY);
        AHL_DECAY = getDouble(p, "field.ahl.decay", AHL_DECAY);
        ACID_DIFFUSIVITY = getDouble(p, "field.acid.diff", ACID_DIFFUSIVITY);
        ACID_DECAY = getDouble(p, "field.acid.decay", ACID_DECAY);
        AHL_SOURCE_RATE = getDouble(p, "field.ahl.source.rate", AHL_SOURCE_RATE);
        ACID_SOURCE_RATE = getDouble(p, "field.acid.source.rate", ACID_SOURCE_RATE);
        ACID_CELL_PRODUCTION_RATE = getDouble(p,
                "field.acid.cell.production.rate", ACID_CELL_PRODUCTION_RATE);
        K_MAX_ACID = getDouble(p, "field.acid.kmax", K_MAX_ACID);
        WARMUP = getDouble(p, "warmup.s", WARMUP);
        WARMUP_AHL_INPUT = getDouble(p, "warmup.ahl.input", WARMUP_AHL_INPUT);
        WINDOW_DURATION = getDouble(p, "window.duration.s", WINDOW_DURATION);
        PULSE_DURATION = getDouble(p, "pulse.duration.s", PULSE_DURATION);
        SAMPLING_DURATION = getDouble(p, "sampling.duration.s", SAMPLING_DURATION);
        SAMPLING_INTERVAL = getDouble(p, "sampling.interval.s", SAMPLING_INTERVAL);
        NUM_WINDOWS = getInt(p, "num.windows", NUM_WINDOWS);
        INITIAL_POP = getInt(p, "initial.pop", INITIAL_POP);
        CARRYING_CAPACITY = getInt(p, "carrying.capacity", CARRYING_CAPACITY);
        RNG_SEED = Long.parseLong(p.getProperty("rng.seed", String.valueOf(RNG_SEED)));
        INPUT_AHL_FILE = p.getProperty("input.ahl.file", INPUT_AHL_FILE);
        INPUT_ACID_FILE = p.getProperty("input.acid.file", INPUT_ACID_FILE);
        OUTPUT_DIR = p.getProperty("output.dir", OUTPUT_DIR);
        RUN_LABEL = p.getProperty("output.run.label", RUN_LABEL);
        STOCHASTIC_REPLICATE = p.getProperty(
                "output.stochastic.replicate", STOCHASTIC_REPLICATE);
        String arm = p.getProperty("arm", "driven").trim().toLowerCase(Locale.ROOT);
        if ("driven".equals(arm)) ARM = Arm.DRIVEN;
        else if ("silent".equals(arm)) ARM = Arm.SILENT;
        else if ("brownian".equals(arm)) ARM = Arm.BROWNIAN;
        else throw new IllegalArgumentException("arm must be driven, brownian, or silent");
        if (ARM == Arm.SILENT) {
            AHL_SOURCE_RATE = 0.0;
            ACID_SOURCE_RATE = 0.0;
            WARMUP_AHL_INPUT = 0.0;
        }
        validateFrozenConfig();
    }

    static void validateFrozenConfig() {
        if (DT != 0.05 || GRID_X != 50 || GRID_Y != 25 || GRID_Z != 1
                || STATE_X != 20 || STATE_Y != 10 || STATE_Z != 1
                || COUNT_X != 4 || COUNT_Y != 2 || COUNT_Z != 1
                || WARMUP != 18000.0 || WINDOW_DURATION != 300.0
                || PULSE_DURATION != 75.0 || SAMPLING_INTERVAL != 20.0
                || INITIAL_POP != 1800 || CARRYING_CAPACITY != 2000
                || AHL_DIFFUSIVITY != 159.0 || AHL_DECAY != 0.0033
                || (ARM != Arm.SILENT && AHL_SOURCE_RATE != 1.28e7)) {
            throw new IllegalArgumentException("I3n frozen physical configuration changed");
        }
    }

    static int getInt(Properties p, String key, int fallback) {
        return Integer.parseInt(p.getProperty(key, String.valueOf(fallback)));
    }

    static double getDouble(Properties p, String key, double fallback) {
        return Double.parseDouble(p.getProperty(key, String.valueOf(fallback)));
    }

    static String format(double value) {
        return String.format(Locale.US, ";%.4e", value);
    }

    static void writeRunStatus(OutputOwner output) {
        File directory = new File(OUTPUT_DIR);
        try (BufferedWriter writer = new BufferedWriter(
                new FileWriter(new File(directory, "run_status.txt")))) {
            writer.write(String.format(Locale.US,
                    "dish=ChassisHybridmm_I3n_new_dish%narm=%s%nseed=%d%n"
                            + "summary_rows=%d%nsample_rows=%d%nvoxel_rows=%d%n"
                            + "mechanics_rows=%d%n",
                    ARM.name().toLowerCase(Locale.ROOT), RNG_SEED, output.summaryRows,
                    output.sampleRows, output.voxelRows, output.mechanicsRows));
        } catch (IOException e) {
            throw new RuntimeException("Cannot write run status", e);
        }
    }

    static void resetStaticState() {
        rods.clear();
        particles.clear();
        voxelAnalyzer = null;
    }
}
