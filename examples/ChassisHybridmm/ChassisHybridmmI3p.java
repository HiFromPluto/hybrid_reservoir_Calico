package ChassisHybridmm;

import BacteriumFromScratch.ChassisParameters;
import BacteriumFromScratch.EcoliRodCell;
import BacteriumFromScratch.ValdezHertzian;
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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Properties;
import java.util.Random;

/**
 * I3p: population-clock join only. No NARMA readout or ridge.
 *
 * <p>This copies the frozen I3n silent rod/mechanics/acid scenario and makes
 * one licensed change: T_removal equals the Warren bath division time.
 */
public final class ChassisHybridmmI3p {
    enum DeathCause { CLAMP, ACID, OOB }

    static final double BOUND_X = 1000.0;
    static final double BOUND_Y = 500.0;
    static final double BOUND_Z = 1.0;
    static final double DT = 0.05;
    static final int GRID_X = 50;
    static final int GRID_Y = 25;
    static final int GRID_Z = 1;
    static final double ACID_DIFFUSIVITY = 200.0;
    static final double ACID_DECAY = 0.0067;
    static final double ACID_CELL_PRODUCTION_RATE = 1.0e5;
    static final double K_MAX_ACID = 0.002;
    static final double MM_TO_MOLECULES_PER_UM3 = 602200.0;
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
    static final double WARMUP_S = 18000.0;
    static final double WINDOW_DURATION_S = 300.0;
    static final int NUM_WINDOWS = 200;
    static final int INITIAL_POP = 1800;
    static final int CARRYING_CAPACITY = 2000;
    static final double T_DIV_S = ChassisParameters.analyticDivisionTimeS(
            ChassisParameters.NUTRIENT_BATH_MM);
    static final double T_REMOVAL_S = T_DIV_S;

    static long rngSeed = 111;
    static String outputDir = "results/i3p_silent_seed111";
    static String runLabel = "i3p_silent";
    static String stochasticReplicate = "silent_seed111";

    static final List<RodAgent> rods = new ArrayList<RodAgent>();
    static Random experimentRng;

    private ChassisHybridmmI3p() {}

    static final class RodAgent {
        final EcoliRodCell body;

        RodAgent(EcoliRodCell body) {
            this.body = body;
        }

        RodAgent daughter() {
            return new RodAgent(body.divideCited());
        }

        Vector3d position() {
            return body.centre();
        }
    }

    static final class MechanicsStats {
        double maxDelta;
        double minCentre = Double.POSITIVE_INFINITY;
        double offplane;
        double nematic;
    }

    static final class OutputOwner implements AutoCloseable {
        final File directory;
        final BufferedWriter summary;
        final BufferedWriter mechanics;
        int summaryRows;
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
                    + "Births;Total_Deaths;Clamp_Deaths;Acid_Deaths;OOB_Deaths\n");
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
                mechanics.close();
            } catch (IOException e) {
                throw new RuntimeException("Cannot close I3p output", e);
            }
        }
    }

    public static void main(String[] args) {
        String configPath = args.length > 0 ? args[0] : "sim_config_i3p.properties";
        loadConfig(configPath);
        validateFrozenConfig();
        rods.clear();
        experimentRng = new Random(rngSeed);

        final BSim sim = new BSim();
        sim.setDt(DT);
        sim.setSimulationTime(WARMUP_S + NUM_WINDOWS * WINDOW_DURATION_S);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);
        final BSimChemicalField acidField = new BSimChemicalField(
                sim, new int[]{GRID_X, GRID_Y, GRID_Z},
                ACID_DIFFUSIVITY, ACID_DECAY);
        placeInitialPopulation(sim);

        final int warmupSteps = exactSteps(WARMUP_S, DT, "warmup");
        final int stepsInWindow = exactSteps(WINDOW_DURATION_S, DT, "window");
        final int growthEvery = exactSteps(ChassisParameters.DT_S, DT, "growth cadence");
        final int mechanicsLogEvery = exactSteps(
                MECHANICS_LOG_INTERVAL_S, DT, "mechanics log");
        final int[] births = {0};
        final int[] deaths = {0};
        final int[] clampDeaths = {0};
        final int[] acidDeaths = {0};
        final int[] oobDeaths = {0};
        final int[] completed = {0};
        final OutputOwner output;
        try {
            output = new OutputOwner(new File(outputDir));
            writeMechanics(output, 0.0, sim.getBound());
        } catch (IOException e) {
            throw new RuntimeException("Cannot open I3p output", e);
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
                        finishWindow(output, window - 1, births[0], deaths[0],
                                clampDeaths[0], acidDeaths[0], oobDeaths[0]);
                        completed[0]++;
                    }
                    births[0] = deaths[0] = clampDeaths[0] = 0;
                    acidDeaths[0] = oobDeaths[0] = 0;
                }

                acidField.update();
                updateRodChemistry(acidField, warmup, deaths, clampDeaths,
                        acidDeaths, oobDeaths);
                if ((step + 1) % growthEvery == 0) {
                    growDivideAndPack(sim.getBound(), warmup, births);
                }
                if ((step + 1) % mechanicsLogEvery == 0) {
                    writeMechanics(output, sim.getTime() + DT, sim.getBound());
                }
            }
        });

        PrintStream originalOut = System.out;
        try {
            System.out.printf(Locale.US,
                    "I3p start seed=%d windows=%d T_div=T_removal=%.9f s%n",
                    rngSeed, NUM_WINDOWS, T_REMOVAL_S);
            System.setOut(new PrintStream(OutputStream.nullOutputStream()));
            sim.export();
            while (completed[0] < NUM_WINDOWS) {
                finishWindow(output, completed[0], births[0], deaths[0],
                        clampDeaths[0], acidDeaths[0], oobDeaths[0]);
                completed[0]++;
            }
        } finally {
            System.setOut(originalOut);
            output.close();
        }
        writeRunStatus(output);
        System.out.printf(Locale.US,
                "I3p complete seed=%d summary=%d mechanics=%d N=%d%n",
                rngSeed, output.summaryRows, output.mechanicsRows, rods.size());
    }

    static void placeInitialPopulation(BSim sim) {
        final int nx = 60;
        final int ny = 30;
        double dx = 400.0 / nx;
        double dy = 200.0 / ny;
        for (int i = 0; i < INITIAL_POP; i++) {
            int ix = i % nx;
            int iy = i / nx;
            double cx = 300.0 + (ix + 0.5) * dx
                    + (experimentRng.nextDouble() - 0.5) * 0.8;
            double cy = 150.0 + (iy + 0.5) * dy
                    + (experimentRng.nextDouble() - 0.5) * 0.8;
            double theta = 2.0 * Math.PI * experimentRng.nextDouble();
            double phase = T_DIV_S * experimentRng.nextDouble();
            double length = ChassisParameters.analyticLengthUm(
                    phase, ChassisParameters.NUTRIENT_BATH_MM);
            Vector3d axis = new Vector3d(Math.cos(theta), Math.sin(theta), 0.0);
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
            rods.add(new RodAgent(body));
        }
        relaxMechanics(rods, sim.getBound());
    }

    static void updateRodChemistry(BSimChemicalField acid, boolean warmup,
            int[] deaths, int[] clampDeaths, int[] acidDeaths, int[] oobDeaths) {
        Iterator<RodAgent> iterator = rods.iterator();
        while (iterator.hasNext()) {
            RodAgent rod = iterator.next();
            Vector3d position = rod.position();
            acid.addQuantity(position, ACID_CELL_PRODUCTION_RATE * DT);

            double pRemoval = (DT / T_REMOVAL_S)
                    * Math.pow(2.0,
                    -(1.0 - rods.size() / (double) CARRYING_CAPACITY));
            boolean clamp = experimentRng.nextDouble() < pRemoval;
            double acidMm = acid.getConc(position) / MM_TO_MOLECULES_PER_UM3;
            boolean acidDeath = experimentRng.nextDouble()
                    < killRate(pHFromAcidMm(acidMm)) * DT;
            DeathCause cause = acidDeath
                    ? DeathCause.ACID : (clamp ? DeathCause.CLAMP : null);
            if (position.x < 0.0 || position.x > BOUND_X
                    || position.y < 0.0 || position.y > BOUND_Y
                    || position.z < 0.0 || position.z > BOUND_Z) {
                cause = DeathCause.OOB;
            }
            if (cause != null) {
                if (!warmup) {
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
                if (!warmup) births[0]++;
            }
        }
        rods.addAll(newborns);
        relaxMechanics(rods, bound);
    }

    static void relaxMechanics(List<RodAgent> population, Vector3d bound) {
        for (int pass = 0; pass < MECHANICS_REBUILDS; pass++) {
            for (List<RodAgent> component : contactComponents(population)) {
                List<EcoliRodCell> bodies =
                        new ArrayList<EcoliRodCell>(component.size());
                for (RodAgent rod : component) bodies.add(rod.body);
                ValdezHertzian.relaxContacts(bodies, bound);
            }
        }
    }

    static List<List<RodAgent>> contactComponents(List<RodAgent> population) {
        int n = population.size();
        int[] parent = new int[n];
        int[] rank = new int[n];
        for (int i = 0; i < n; i++) parent[i] = i;
        Map<Long, List<Integer>> bins = new HashMap<Long, List<Integer>>();
        for (int i = 0; i < n; i++) {
            Vector3d centre = population.get(i).position();
            int bx = (int) Math.floor(centre.x / MECHANICS_COMPONENT_CUTOFF_UM);
            int by = (int) Math.floor(centre.y / MECHANICS_COMPONENT_CUTOFF_UM);
            int bz = (int) Math.floor(centre.z / MECHANICS_COMPONENT_CUTOFF_UM);
            for (int ox = -1; ox <= 1; ox++) {
                for (int oy = -1; oy <= 1; oy++) {
                    for (int oz = -1; oz <= 1; oz++) {
                        List<Integer> candidates =
                                bins.get(binKey(bx + ox, by + oy, bz + oz));
                        if (candidates == null) continue;
                        for (int j : candidates) {
                            Vector3d delta = new Vector3d(centre);
                            delta.sub(population.get(j).position());
                            if (delta.lengthSquared()
                                    <= MECHANICS_COMPONENT_CUTOFF_UM
                                    * MECHANICS_COMPONENT_CUTOFF_UM) {
                                union(parent, rank, i, j);
                            }
                        }
                    }
                }
            }
            long key = binKey(bx, by, bz);
            if (!bins.containsKey(key)) {
                bins.put(key, new ArrayList<Integer>());
            }
            bins.get(key).add(i);
        }
        Map<Integer, List<RodAgent>> grouped =
                new LinkedHashMap<Integer, List<RodAgent>>();
        for (int i = 0; i < n; i++) {
            int root = find(parent, i);
            if (!grouped.containsKey(root)) {
                grouped.put(root, new ArrayList<RodAgent>());
            }
            grouped.get(root).add(population.get(i));
        }
        return new ArrayList<List<RodAgent>>(grouped.values());
    }

    static long binKey(int x, int y, int z) {
        return ((long) (x + 2048) << 32)
                ^ ((long) (y + 2048) << 16) ^ (z + 2048L);
    }

    static int find(int[] parent, int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    static void union(int[] parent, int[] rank, int a, int b) {
        int rootA = find(parent, a);
        int rootB = find(parent, b);
        if (rootA == rootB) return;
        if (rank[rootA] < rank[rootB]) parent[rootA] = rootB;
        else if (rank[rootA] > rank[rootB]) parent[rootB] = rootA;
        else {
            parent[rootB] = rootA;
            rank[rootA]++;
        }
    }

    static MechanicsStats mechanicsStats(Vector3d bound) {
        MechanicsStats result = new MechanicsStats();
        if (rods.isEmpty()) return result;
        int offplane = 0;
        double cosine = 0.0;
        double sine = 0.0;
        for (RodAgent rod : rods) {
            Vector3d centre = rod.position();
            if (Math.abs(centre.z - 0.5)
                    > 0.5 * ChassisParameters.RADIUS_UM) {
                offplane++;
            }
            Vector3d axis = new Vector3d();
            axis.sub(rod.body.x2, rod.body.x1);
            double theta = Math.atan2(axis.y, axis.x);
            cosine += Math.cos(2.0 * theta);
            sine += Math.sin(2.0 * theta);
        }
        result.offplane = offplane / (double) rods.size();
        result.nematic = Math.hypot(
                cosine / rods.size(), sine / rods.size());
        for (List<RodAgent> component : contactComponents(rods)) {
            List<EcoliRodCell> bodies =
                    new ArrayList<EcoliRodCell>(component.size());
            for (RodAgent rod : component) bodies.add(rod.body);
            ValdezHertzian.PackingStats stats =
                    ValdezHertzian.stats(bodies, bound);
            result.maxDelta = Math.max(result.maxDelta, stats.maxDeltaCc);
            result.minCentre = Math.min(
                    result.minCentre, stats.minCentreDist);
        }
        return result;
    }

    static void writeMechanics(
            OutputOwner output, double time, Vector3d bound) {
        MechanicsStats stats = mechanicsStats(bound);
        // No pair inside the conservative broadphase means the true minimum is
        // at least the cutoff. Log that finite censored bound, not Infinity.
        double minCentre = Double.isFinite(stats.minCentre)
                ? stats.minCentre : MECHANICS_COMPONENT_CUTOFF_UM;
        try {
            output.mechanics.write(String.format(Locale.US,
                    "%.2f;%d;%d;%.9g;%.9g;%.9g;%.9g%n",
                    time, rngSeed, rods.size(), stats.offplane,
                    stats.nematic, stats.maxDelta, minCentre));
            output.mechanicsRows++;
        } catch (IOException e) {
            throw new RuntimeException("Cannot write I3p mechanics", e);
        }
    }

    static void finishWindow(OutputOwner output, int window,
            int births, int deaths, int clampDeaths,
            int acidDeaths, int oobDeaths) {
        try {
            output.summary.write(String.format(Locale.US,
                    "%s;%s;%d;%d;%d;%d;%d;%d;%d%n",
                    runLabel, stochasticReplicate, window, rods.size(),
                    births, deaths, clampDeaths, acidDeaths, oobDeaths));
            output.summaryRows++;
        } catch (IOException e) {
            throw new RuntimeException("Cannot write I3p summary", e);
        }
    }

    static double pHFromAcidMm(double acidMm) {
        return Math.max(1.0, Math.min(14.0,
                PH_BASE - acidMm / BUFFER_CAPACITY_MM_PER_PH));
    }

    static double killRate(double pH) {
        double acidExcess = Math.max(0.0, PH_GROWTH_LIMIT - pH);
        double acidHalf = PH_GROWTH_LIMIT - PH_ACID_EC50;
        double kAcid = K_MAX_ACID * Math.pow(acidExcess, N_ACID)
                / (Math.pow(acidHalf, N_ACID)
                + Math.pow(acidExcess, N_ACID));
        double alkExcess = Math.max(0.0, pH - PH_ALK_ONSET);
        double alkHalf = PH_ALK_EC50 - PH_ALK_ONSET;
        double kAlk = K_MAX_ALK * Math.pow(alkExcess, N_ALK)
                / (Math.pow(alkHalf, N_ALK)
                + Math.pow(alkExcess, N_ALK));
        return kAcid + kAlk;
    }

    static int exactSteps(double duration, double dt, String name) {
        int steps = (int) Math.round(duration / dt);
        if (Math.abs(steps * dt - duration) > 1e-9) {
            throw new IllegalArgumentException(
                    name + " is not an integer number of field ticks");
        }
        return steps;
    }

    static void loadConfig(String path) {
        Properties properties = new Properties();
        try (FileInputStream input = new FileInputStream(path)) {
            properties.load(input);
        } catch (IOException e) {
            throw new IllegalArgumentException(
                    "Cannot load I3p config " + path, e);
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
                        "I3p config may not change frozen key " + key);
            }
        }
    }

    static void validateFrozenConfig() {
        double expected = ChassisParameters.analyticDivisionTimeS(
                ChassisParameters.NUTRIENT_BATH_MM);
        if (Double.doubleToLongBits(T_REMOVAL_S)
                != Double.doubleToLongBits(expected)
                || T_REMOVAL_S != T_DIV_S) {
            throw new IllegalArgumentException(
                    "I3p frozen population-clock configuration changed");
        }
    }

    static void writeRunStatus(OutputOwner output) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(
                new File(output.directory, "run_status.txt")))) {
            writer.write(String.format(Locale.US,
                    "job=ChassisHybridmm_I3p_population_clock%n"
                            + "arm=silent%nseed=%d%nsummary_rows=%d%n"
                            + "mechanics_rows=%d%nT_div_s=%.12g%n"
                            + "T_removal_s=%.12g%n"
                            + "d_centers_no_pair_sentinel_um=%.12g%n"
                            + "ridge=OFF%nNARMA_overall=NONE%n",
                    rngSeed, output.summaryRows, output.mechanicsRows,
                    T_DIV_S, T_REMOVAL_S,
                    MECHANICS_COMPONENT_CUTOFF_UM));
        } catch (IOException e) {
            throw new RuntimeException("Cannot write I3p run status", e);
        }
    }
}
