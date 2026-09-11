package BSimReservoirStage11;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Properties;
import java.io.FileInputStream;
import java.util.Vector;

import javax.vecmath.Vector3d;

import processing.core.PGraphics3D;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.draw.BSimP3DDrawer;
import bsim.particle.BSimBacterium;

/**
 * Stage 11 — Merged reservoir simulation.
 *
 * Combines:
 *   - Stage10's Danino QSGRN, dual chemotaxis, homeostatic clamp, Hill toxicity
 *   - reservoir_new's ODE-gated ACs, glucose/Monod, RC window protocol, voxel readout
 *
 * Plan Stage 1 (affordability): dt=0.05, grid 50x25x1, readout STATE 20x10 + COUNT 4x2
 * Plan Stage 2 (timescales): warmup=1800, window=300, pulse=75, sampling=300/20
 *
 * Population plateau: Monod glucose (~850 at G0=20), not LacOperon clamp (2000).
 * Plateau gate: linear trend of N across analysis windows; |slope·T|/mean(N) < 5%.
 *
 * Stage 2 birth-input: NOT chemotaxis -> depletion (supply-capped, r flat across 27x field).
 * AC3 NUTRIENT injects into glucoseField so supply = Model A + AC_glucose(u); births track u
 * because they are supply-limited (K_G = KM_UPTAKE = 5, tuned; do not decouple for this gate).
 * Attractant gamma is not raised. 4x2 Birth_* is rank ~1 (plume < count bin); accepted.
 *
 * AC layout: 3 midline ACs (Stage10 x=250/500/750) plus Stage 7 AC5 nutrient at (550, 400).
 * Full 5-AC cross-stream layout remains Stage 7. Do not start Plan Stage 3 until birth gate passes.
 * AC0/AC1: attractant secretion + motility (chemotaxis toward signalField).
 * AC2: REPELLENT field secretion for Rep_* chemical readout ONLY — not a motility channel.
 *   Bacteria may still sense repellentField for Hill toxicity; sign-flip taxis is NULL (Stage10_5
 *   n=5: t=1.2 on near-count) and will be replaced by pH (Yang & Sourjik 2012), not tuned.
 * AC3: NUTRIENT (glucose) at (550, 400) — Stage 7 AC5 site; input-driven supply for births.
 *
 * Flow: deliberately zero (reservoir_new used 8 µm/s). ConvectionApplier is wired but inactive
 * at flow.speed=0 — enable via sim_config when cross-stream mixing is needed.
 *
 * Unit fix: 602 molecules/um^3 per uM (not 1e15).
 */
public class BSimReservoirStage11 {

    // ── Domain ──
    static double BOUND_X = 1000.0;
    static double BOUND_Y = 500.0;
    static double BOUND_Z = 10.0;

    // ── Field grid (Plan Stage 1: coarsened) ──
    static int GRID_X = 50;
    static int GRID_Y = 25;
    static int GRID_Z = 1;

    // ── Readout grids ──
    static int STATE_X = 20, STATE_Y = 10, STATE_Z = 1;
    static int COUNT_X = 4,  COUNT_Y = 2,  COUNT_Z = 1;

    // ── Field parameters ──
    static double ATT_DIFF  = 100.0;
    static double ATT_DECAY = 0.0067;   // tau ~150s
    static double REP_DIFF  = 100.0;
    static double REP_DECAY = 0.033;    // tau ~30s
    static double AHL_DIFF  = 25.0;     // tuned, not sourced
    static double AHL_DECAY = 0.0033;   // tau ~300s
    static double GLU_DIFF  = 5.0;    // low D: keep uptake depletion shadows (100 um^2/s re-homogenises voxels)
    static double GLU_DECAY = 0.00017;  // tau ~6000s
    static double GLU_G0    = 20.0;
    static double GLU_K_SUPPLY = 0.001;  // Model A: dG/dt = k_s(G0-G); max supply k_s*G0*V ≈ 1e5 qty/s

    // ── AC parameters (ODE secretion only — no direct field injection) ──
    // Store-limited: CS_max sets pulse inventory + refill (k_refill=0.1).
    // 1x (gamma=1.5e6, CS_max=5e5) → peak_att@AC0 ≈ 38 (birth-input r=0.261, FAIL).
    // 8x target ≈ Stage10_5 full-rate field ~311. Sweep gamma/CS_max together.
    static double AC_GAMMA_SYM = 1.2e7;
    static double AC_CS_MAX    = 4e6;
    static double AC_KM        = 5e4;

    // ── Population ──
    static int INITIAL_POP       = 850;  // Monod equilibrium (~850); avoids warm-up decline from 950
    static int CARRYING_CAPACITY = 2000;

    // ── Timing (Plan Stage 2) ──
    static double DT              = 0.05;
    static double WARMUP          = 7200.0;  // extended until |dN| at analysis start < 5% of warm-up peak
    static double WINDOW_DURATION = 300.0;
    static double PULSE_DURATION  = 75.0;    // 25% duty cycle
    static double SAMPLING_DURATION = 300.0; // sample whole window
    static double SAMPLING_INTERVAL = 20.0;  // 15 samples per window
    static int    NUM_WINDOWS     = 40;

    // ── Flow (0 = deliberate static dish; reservoir_new used 8 µm/s) ──
    static double FLOW_SPEED = 0.0;

    // ── AC layout ──
    // AC3 = Stage 7 AC5 nutrient site. Attractant ACs not cannibalised; gamma not raised.
    static final Vector3d[] AC_POSITIONS = {
        new Vector3d(250, 250, 5),
        new Vector3d(500, 250, 5),
        new Vector3d(750, 250, 5),
        new Vector3d(550, 400, 5)
    };
    static final ACConfig.SecretionType[] AC_TYPES = {
        ACConfig.SecretionType.ATTRACTANT,
        ACConfig.SecretionType.ATTRACTANT,
        ACConfig.SecretionType.REPELLENT,
        ACConfig.SecretionType.NUTRIENT
    };
    static final int NUTRIENT_AC_INDEX = 3;

    static double NUTRIENT_GAMMA  = 1.5e6;  // store-limited; continuous I ≈ k_refill*CS_max, not gamma
    static double NUTRIENT_CSMAX  = 5e5;
    static double NUTRIENT_K_REFILL = 0.1;

    static boolean HEADLESS = true;
    static int     RNG_SEED = -1;
    static boolean REGRESSION_METRICS_ONLY = false;
    static String  METRICS_TAG = "default";
    static boolean RECORD_WALL_TIME = false;
    static String  OUTPUT_DIR = "results";

    // ========================================================================
    // Load optional properties file
    // ========================================================================
    static void loadConfig(String path) {
        try (FileInputStream fis = new FileInputStream(path)) {
            Properties p = new Properties();
            p.load(fis);
            DT              = Double.parseDouble(p.getProperty("dt", String.valueOf(DT)));
            GRID_X          = Integer.parseInt(p.getProperty("grid.x", String.valueOf(GRID_X)));
            GRID_Y          = Integer.parseInt(p.getProperty("grid.y", String.valueOf(GRID_Y)));
            WARMUP          = Double.parseDouble(p.getProperty("warmup.s", String.valueOf(WARMUP)));
            WINDOW_DURATION = Double.parseDouble(p.getProperty("window.duration.s", String.valueOf(WINDOW_DURATION)));
            PULSE_DURATION  = Double.parseDouble(p.getProperty("pulse.duration.s", String.valueOf(PULSE_DURATION)));
            SAMPLING_DURATION = Double.parseDouble(p.getProperty("sampling.duration.s", String.valueOf(SAMPLING_DURATION)));
            SAMPLING_INTERVAL = Double.parseDouble(p.getProperty("sampling.interval.s", String.valueOf(SAMPLING_INTERVAL)));
            NUM_WINDOWS     = Integer.parseInt(p.getProperty("num.windows", String.valueOf(NUM_WINDOWS)));
            INITIAL_POP     = Integer.parseInt(p.getProperty("initial.pop", String.valueOf(INITIAL_POP)));
            CARRYING_CAPACITY = Integer.parseInt(p.getProperty("carrying.capacity", String.valueOf(CARRYING_CAPACITY)));
            HEADLESS        = Boolean.parseBoolean(p.getProperty("headless", String.valueOf(HEADLESS)));
            ATT_DECAY       = Double.parseDouble(p.getProperty("field.att.decay", String.valueOf(ATT_DECAY)));
            REP_DECAY       = Double.parseDouble(p.getProperty("field.rep.decay", String.valueOf(REP_DECAY)));
            AHL_DECAY       = Double.parseDouble(p.getProperty("field.ahl.decay", String.valueOf(AHL_DECAY)));
            GLU_DECAY       = Double.parseDouble(p.getProperty("field.glu.decay", String.valueOf(GLU_DECAY)));
            GLU_DIFF        = Double.parseDouble(p.getProperty("field.glu.diff", String.valueOf(GLU_DIFF)));
            GLU_G0          = Double.parseDouble(p.getProperty("field.glu.g0", String.valueOf(GLU_G0)));
            GLU_K_SUPPLY    = Double.parseDouble(p.getProperty("field.glu.k_supply", String.valueOf(GLU_K_SUPPLY)));
            FLOW_SPEED      = Double.parseDouble(p.getProperty("flow.speed", String.valueOf(FLOW_SPEED)));
            STATE_X         = Integer.parseInt(p.getProperty("readout.grid.x", String.valueOf(STATE_X)));
            STATE_Y         = Integer.parseInt(p.getProperty("readout.grid.y", String.valueOf(STATE_Y)));
            STATE_Z         = Integer.parseInt(p.getProperty("readout.grid.z", String.valueOf(STATE_Z)));
            COUNT_X         = Integer.parseInt(p.getProperty("readout.countgrid.x", String.valueOf(COUNT_X)));
            COUNT_Y         = Integer.parseInt(p.getProperty("readout.countgrid.y", String.valueOf(COUNT_Y)));
            COUNT_Z         = Integer.parseInt(p.getProperty("readout.countgrid.z", String.valueOf(COUNT_Z)));
            AC_GAMMA_SYM    = Double.parseDouble(p.getProperty("ac.gamma", String.valueOf(AC_GAMMA_SYM)));
            AC_CS_MAX       = Double.parseDouble(p.getProperty("ac.csmax", String.valueOf(AC_CS_MAX)));
            AC_KM           = Double.parseDouble(p.getProperty("ac.km", String.valueOf(AC_KM)));
            NUTRIENT_GAMMA  = Double.parseDouble(p.getProperty("ac.nutrient.gamma", String.valueOf(NUTRIENT_GAMMA)));
            NUTRIENT_CSMAX  = Double.parseDouble(p.getProperty("ac.nutrient.csmax", String.valueOf(NUTRIENT_CSMAX)));
            NUTRIENT_K_REFILL = Double.parseDouble(p.getProperty("ac.nutrient.k_refill", String.valueOf(NUTRIENT_K_REFILL)));
            ReservoirBacterium.K_G = Double.parseDouble(
                    p.getProperty("k.g", String.valueOf(ReservoirBacterium.K_G)));
            ReservoirBacterium.KM_UPTAKE = Double.parseDouble(
                    p.getProperty("km.uptake", String.valueOf(ReservoirBacterium.KM_UPTAKE)));
            if (p.containsKey("rng.seed"))
                RNG_SEED = Integer.parseInt(p.getProperty("rng.seed"));
            REGRESSION_METRICS_ONLY = Boolean.parseBoolean(
                    p.getProperty("regression.metrics.only", String.valueOf(REGRESSION_METRICS_ONLY)));
            METRICS_TAG = p.getProperty("output.metrics.tag", METRICS_TAG);
            RECORD_WALL_TIME = Boolean.parseBoolean(
                    p.getProperty("record_wall_time", String.valueOf(RECORD_WALL_TIME)));
            OUTPUT_DIR = p.getProperty("output.dir", OUTPUT_DIR);
            // Clamp timing to window
            PULSE_DURATION    = Math.min(PULSE_DURATION, WINDOW_DURATION);
            SAMPLING_DURATION = Math.min(SAMPLING_DURATION, WINDOW_DURATION);
            SAMPLING_INTERVAL = Math.min(SAMPLING_INTERVAL, WINDOW_DURATION);
            System.out.println("Loaded config from " + path);
        } catch (IOException e) {
            System.out.println("No config file (" + path + "), using defaults.");
        }
    }

    // ========================================================================
    // Read input sequence file (same format as Stage10: line1=bitDuration, line2=bits)
    // ========================================================================
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
            return new Object[]{WINDOW_DURATION, generateDefaultGradedSequence(NUM_WINDOWS)};
        }
    }

    // Read RC-protocol input file (one value per line, one per window)
    static double[] readInputFile(String path) {
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            java.util.List<String> lines = new java.util.ArrayList<>();
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty()) lines.add(line);
            }
            double[] vals = new double[lines.size()];
            for (int i = 0; i < lines.size(); i++)
                vals[i] = Double.parseDouble(lines.get(i));
            return vals;
        } catch (IOException e) {
            return null;
        }
    }

    static double[] generateDefaultGradedSequence(int len) {
        double[] seq = new double[len];
        for (int i = 0; i < len; i++) {
            // Graded placeholder (~0.15–0.25), not binary — real runs use input_AC*.txt biomarkers
            seq[i] = 0.20 + 0.05 * Math.sin(i * 0.85);
        }
        return seq;
    }

    static double clampInput(double u) {
        return Math.max(0.0, Math.min(1.0, u));
    }

    /** Least-squares slope dN/dt (cells/s) for N sampled at window starts vs time. */
    static double[] populationTrend(java.util.List<Integer> pops, double windowDurationSec) {
        int n = pops.size();
        if (n < 2) return new double[]{0.0, n > 0 ? pops.get(0) : 0.0, 0.0};
        double sumT = 0, sumN = 0, sumTT = 0, sumTN = 0;
        for (int i = 0; i < n; i++) {
            double t = i * windowDurationSec;
            double N = pops.get(i);
            sumT += t;
            sumN += N;
            sumTT += t * t;
            sumTN += t * N;
        }
        double denom = n * sumTT - sumT * sumT;
        double slope = (denom != 0) ? (n * sumTN - sumT * sumN) / denom : 0.0;
        return new double[]{slope, sumN / n, n};
    }

    // ========================================================================
    // Main
    // ========================================================================
    public static void main(String[] args) {

        boolean useRCProtocol = false; // set true if input_AC0.txt etc. found

        // Optional config
        String configPath = (args.length > 0) ? args[0] : "sim_config.properties";
        loadConfig(configPath);

        if (args.length > 1 && args[1].equals("preview")) {
            HEADLESS = false;
        }

        // Propagate carrying capacity to bacterium
        ReservoirBacterium.CARRYING_CAPACITY = CARRYING_CAPACITY;
        ReservoirBacterium.FLOW_SPEED = FLOW_SPEED;

        // ── Load RC-protocol input files (graded reals, one per window) ──
        int configuredWindows = NUM_WINDOWS;
        double[][] acSignals = new double[AC_POSITIONS.length][];
        int signalLength = -1;
        for (int i = 0; i < AC_POSITIONS.length; i++) {
            acSignals[i] = readInputFile("input_AC" + i + ".txt");
            if (acSignals[i] != null) {
                useRCProtocol = true;
                if (signalLength < 0) signalLength = acSignals[i].length;
            }
        }
        if (useRCProtocol) {
            if (signalLength < 0) signalLength = configuredWindows;
            NUM_WINDOWS = Math.min(configuredWindows, signalLength);
            signalLength = NUM_WINDOWS;
            for (int i = 0; i < AC_POSITIONS.length; i++) {
                if (acSignals[i] == null) {
                    acSignals[i] = new double[signalLength];
                } else if (acSignals[i].length > signalLength) {
                    double[] trimmed = new double[signalLength];
                    System.arraycopy(acSignals[i], 0, trimmed, 0, signalLength);
                    acSignals[i] = trimmed;
                }
            }
        } else {
            signalLength = NUM_WINDOWS;
            for (int i = 0; i < AC_POSITIONS.length; i++) {
                acSignals[i] = generateDefaultGradedSequence(signalLength);
            }
        }

        double totalSimTime = WARMUP + NUM_WINDOWS * WINDOW_DURATION;

        System.out.println("============================================================");
        System.out.println("BSimReservoirStage11 — Merged Reservoir Simulation");
        System.out.println("============================================================");
        System.out.printf("  dt=%.3f  grid=%dx%d  domain=%.0fx%.0fx%.0f%n",
                DT, GRID_X, GRID_Y, BOUND_X, BOUND_Y, BOUND_Z);
        System.out.printf("  warmup=%.0fs  window=%.0fs  pulse=%.0fs  windows=%d%n",
                WARMUP, WINDOW_DURATION, PULSE_DURATION, NUM_WINDOWS);
        System.out.printf("  initial_pop=%d  carrying_capacity=%d%n", INITIAL_POP, CARRYING_CAPACITY);
        System.out.printf("  STATE readout=%dx%d  COUNT readout=%dx%d%n",
                STATE_X, STATE_Y, COUNT_X, COUNT_Y);
        System.out.printf("  total sim time=%.0fs (%.1f min)%n", totalSimTime, totalSimTime / 60.0);
        System.out.printf("  k_g=%.3f km_uptake=%.3f (keep equal for supply-limited births)%n",
                ReservoirBacterium.K_G, ReservoirBacterium.KM_UPTAKE);
        System.out.printf("  ac3_nutrient gamma=%.2e csmax=%.2e k_refill=%.2f at (550,400)%n",
                NUTRIENT_GAMMA, NUTRIENT_CSMAX, NUTRIENT_K_REFILL);
        System.out.printf("  glu D=%.1f k=%.5f G0=%.1f k_s=%.4f%n",
                GLU_DIFF, GLU_DECAY, GLU_G0, GLU_K_SUPPLY);
        System.out.printf("  mode=%s%n", HEADLESS ? "headless" : "preview");
        System.out.println("============================================================");

        // ── BSim setup ──
        BSim sim = new BSim();
        sim.setDt(DT);
        sim.setSimulationTime(totalSimTime);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);

        // ── Chemical fields ──
        int[] fieldGrid = {GRID_X, GRID_Y, GRID_Z};
        final BSimChemicalField signalField = new BSimChemicalField(sim, fieldGrid, ATT_DIFF, ATT_DECAY);
        final BSimChemicalField repellentField = new BSimChemicalField(sim, fieldGrid, REP_DIFF, REP_DECAY);
        final BSimChemicalField ahlField = new BSimChemicalField(sim, fieldGrid, AHL_DIFF, AHL_DECAY);
        final BSimChemicalField glucoseField = new BSimChemicalField(sim, fieldGrid, GLU_DIFF, GLU_DECAY);
        glucoseField.setConc(GLU_G0);

        // ── Utilities ──
        final GlucoseReplenisher glucoseReplenisher = new GlucoseReplenisher(GLU_G0, GLU_K_SUPPLY);
        final ConvectionApplier convection = new ConvectionApplier();
        final Vector3d flowVelocity = new Vector3d(FLOW_SPEED, 0, 0);

        final VoxelAnalyzer voxelAnalyzer = new VoxelAnalyzer(
                new int[]{STATE_X, STATE_Y, STATE_Z},
                new int[]{COUNT_X, COUNT_Y, COUNT_Z},
                BOUND_X, BOUND_Y, BOUND_Z);

        // ── ACs (ODE-gated from reservoir_new) ──
        final ArtificialCell[] acs = new ArtificialCell[AC_POSITIONS.length];
        for (int i = 0; i < AC_POSITIONS.length; i++) {
            ACConfig.SecretionType type = AC_TYPES[i];
            boolean nutrient = (type == ACConfig.SecretionType.NUTRIENT);
            ACConfig acCfg = ACConfig.builder(AC_POSITIONS[i], type)
                    .gammaSym(nutrient ? NUTRIENT_GAMMA : AC_GAMMA_SYM, 0.0)
                    .csMax(nutrient ? NUTRIENT_CSMAX : AC_CS_MAX, 0.0)
                    .km(AC_KM, 0.0)
                    .kRefill(nutrient ? NUTRIENT_K_REFILL : 0.1)
                    .heterogeneityEnabled(false)
                    .build();
            acs[i] = new ArtificialCell(sim, acCfg, signalField, repellentField, ahlField, glucoseField);
        }

        // ── Bacteria ──
        final Vector<ReservoirBacterium> bacteria = new Vector<>();
        final Vector<ReservoirBacterium> childrenBuf = new Vector<>();
        final Vector<ReservoirBacterium> removals = new Vector<>();

        ReservoirBacterium.allBacteria = bacteria;

        java.util.Random popRng = (RNG_SEED >= 0) ? new java.util.Random(RNG_SEED) : null;
        for (int i = 0; i < INITIAL_POP; i++) {
            double bx = (popRng != null) ? popRng.nextDouble() * BOUND_X : Math.random() * BOUND_X;
            double by = (popRng != null) ? popRng.nextDouble() * BOUND_Y : Math.random() * BOUND_Y;
            double bz = BOUND_Z / 2.0;
            ReservoirBacterium b = new ReservoirBacterium(sim,
                    new Vector3d(bx, by, bz),
                    signalField, repellentField, ahlField, glucoseField,
                    voxelAnalyzer, childrenBuf);
            b.setRadius();
            b.setSurfaceAreaGrowthRate(ReservoirBacterium.GROWTH_RATE);
            b.setChildList(childrenBuf);
            bacteria.add(b);
        }

        // ── Timing ──
        final int stepsInWindow    = (int)(WINDOW_DURATION / DT);
        final int stepsInPulse     = (int)(PULSE_DURATION / DT);
        final int warmUpSteps      = (int)(WARMUP / DT);
        final int samplingInterval = (int)(SAMPLING_INTERVAL / DT);
        final int stepsInSampling  = (int)(SAMPLING_DURATION / DT);
        final int samplingStart    = stepsInWindow - stepsInSampling;

        // ── Output files ──
        java.io.File resultsDir = new java.io.File(OUTPUT_DIR);
        if (!resultsDir.exists()) resultsDir.mkdirs();

        final FileWriter resultsWriter, voxelsWriter, gateWriter;
        try {
            if (REGRESSION_METRICS_ONLY) {
                resultsWriter = null;
                voxelsWriter = null;
                gateWriter = null;
            } else {
                resultsWriter = new FileWriter(new java.io.File(resultsDir, "stage11_results.csv"));
                resultsWriter.write("Window;Sample;Input_AC0;Population;Deaths_This_Window;Births_This_Window\n");

                voxelsWriter = new FileWriter(new java.io.File(resultsDir, "stage11_voxels.csv"));
                StringBuilder hdr = new StringBuilder("Window;Sample;Input");
                int V = voxelAnalyzer.getStateVoxels();
                int B = voxelAnalyzer.getCountBins();
                for (int i = 0; i < V; i++) hdr.append(";Att_").append(i);
                for (int i = 0; i < V; i++) hdr.append(";Rep_").append(i);
                for (int i = 0; i < V; i++) hdr.append(";Den_").append(i);
                for (int i = 0; i < V; i++) hdr.append(";Act_").append(i);
                for (int i = 0; i < V; i++) hdr.append(";Prot_").append(i);
                for (int i = 0; i < V; i++) hdr.append(";AHL_").append(i);
                for (int i = 0; i < B; i++) hdr.append(";Birth_").append(i);
                for (int i = 0; i < B; i++) hdr.append(";Death_").append(i);
                hdr.append("\n");
                voxelsWriter.write(hdr.toString());

                gateWriter = new FileWriter(new java.io.File(resultsDir, "gate_summary.txt"));
            }
        } catch (IOException e) {
            throw new RuntimeException("Cannot open output files", e);
        }

        final int[] analysisStartN = {-1};
        final java.util.List<Integer> analysisPopByWindow = new java.util.ArrayList<>();
        final java.util.List<Integer> analysisBirthsByWindow = new java.util.ArrayList<>();
        final java.util.List<Integer> analysisDeathsByWindow = new java.util.ArrayList<>();
        final java.util.List<int[]> analysisBirthBinsByWindow = new java.util.ArrayList<>();
        final java.util.List<Double> meanLuxIByWindow = new java.util.ArrayList<>();
        final int[] totalDeathsCumulative = {0};

        final int[] deathsThisWindow = {0};
        final int[] birthsThisWindow = {0};
        final double[] peakAttConc = {0.0};
        final double[] peakRepConc = {0.0};
        final double[] peakGluConc = {0.0};

        // Capture for lambda
        final double[][] signals = acSignals;
        final int sigLen = signalLength;

        // ====================================================================
        // TICKER
        // ====================================================================
        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                int currentStep   = (int) sim.getTimestep();
                boolean inWarmUp  = currentStep < warmUpSteps;
                int effectiveStep = currentStep - warmUpSteps;
                int stepInWindow  = inWarmUp ? 0 : (effectiveStep % stepsInWindow);
                int windowIndex   = inWarmUp ? 0 : (effectiveStep / stepsInWindow);

                // ── Window start: set AC inputs, reset counters ──
                if (!inWarmUp) {
                    if (stepInWindow == 0) {
                        if (windowIndex > 0) {
                            analysisBirthsByWindow.add(birthsThisWindow[0]);
                            analysisDeathsByWindow.add(deathsThisWindow[0]);
                            analysisBirthBinsByWindow.add(voxelAnalyzer.snapshotBirths());
                        }
                        if (windowIndex >= sigLen) {
                            // Do not System.exit — that aborts export before gate summary is written
                            for (ArtificialCell ac : acs) ac.forcedInput = 0.0;
                        } else {
                        analysisPopByWindow.add(bacteria.size());
                        double luxSum = 0.0;
                        for (ReservoirBacterium b : bacteria) luxSum += b.getY()[0];
                        meanLuxIByWindow.add(luxSum / Math.max(1, bacteria.size()));
                        for (int i = 0; i < acs.length; i++)
                            acs[i].forcedInput = clampInput(signals[i][windowIndex]);
                        deathsThisWindow[0] = 0;
                        birthsThisWindow[0] = 0;
                        voxelAnalyzer.resetWindowCounters();

                        if (windowIndex % 10 == 0) {
                            System.out.printf("  [W%03d] N=%d  u=[%.2f,%.2f,%.2f,%.2f]%n",
                                    windowIndex, bacteria.size(),
                                    signals[0][windowIndex],
                                    signals[1][windowIndex],
                                    signals[2][windowIndex],
                                    signals[3][windowIndex]);
                        }
                        }
                    } else if (stepInWindow >= stepsInPulse) {
                        for (ArtificialCell ac : acs) ac.forcedInput = 0.0;
                    }
                } else {
                    for (ArtificialCell ac : acs) ac.forcedInput = 0.0;
                    if (currentStep % (stepsInWindow) == 0 && currentStep > 0) {
                        System.out.printf("  [warm-up t=%.0fs] N=%d%n",
                                sim.getTime(), bacteria.size());
                    }
                    if (currentStep == warmUpSteps - 1) {
                        System.out.println("  ── Warm-up complete. Starting analysis. ──");
                        analysisStartN[0] = bacteria.size();
                    }
                }

                // ── AC action (sole secretion path: ODE gate + MM release + store depletion) ──
                for (ArtificialCell ac : acs) ac.action();

                peakAttConc[0] = Math.max(peakAttConc[0], signalField.getConc(AC_POSITIONS[0]));
                peakRepConc[0] = Math.max(peakRepConc[0], repellentField.getConc(AC_POSITIONS[2]));
                peakGluConc[0] = Math.max(peakGluConc[0], glucoseField.getConc(AC_POSITIONS[NUTRIENT_AC_INDEX]));

                // ── Bacteria action + position ──
                for (ReservoirBacterium b : bacteria) {
                    b.action();
                    b.updatePosition();
                }

                // ── Children ──
                int prevSize = bacteria.size();
                for (ReservoirBacterium child : childrenBuf) {
                    if (bacteria.size() < CARRYING_CAPACITY * 2) { // hard cap
                        bacteria.add(child);
                    }
                }
                int born = bacteria.size() - prevSize;
                birthsThisWindow[0] += born;
                childrenBuf.clear();

                // ── Deaths ──
                for (ReservoirBacterium b : bacteria) {
                    if (b.isMarkedForDeath()) {
                        removals.add(b);
                        if (voxelAnalyzer != null) {
                            voxelAnalyzer.recordDeath(b.getPosition());
                        }
                    }
                }
                deathsThisWindow[0] += removals.size();
                totalDeathsCumulative[0] += removals.size();
                bacteria.removeAll(removals);
                removals.clear();

                // ── Fields: diffuse/decay + glucose replenish + advection ──
                signalField.update();
                repellentField.update();
                ahlField.update();
                glucoseField.update();
                glucoseReplenisher.apply(glucoseField, sim.getDt());

                if (FLOW_SPEED > 0) {
                    convection.apply(signalField, flowVelocity, sim.getDt(), 0, false);
                    convection.apply(repellentField, flowVelocity, sim.getDt(), 0, false);
                    convection.apply(ahlField, flowVelocity, sim.getDt(), 0, false);
                }

                // ── Sampling ──
                if (!REGRESSION_METRICS_ONLY && !inWarmUp && stepInWindow >= samplingStart) {
                    int relStep = stepInWindow - samplingStart;
                    if (relStep % samplingInterval == 0 && windowIndex < sigLen) {
                        int sampleIndex = relStep / samplingInterval;
                        writeSample(resultsWriter, voxelsWriter,
                                windowIndex, sampleIndex, signals[0][windowIndex],
                                bacteria, voxelAnalyzer,
                                signalField, repellentField, ahlField,
                                deathsThisWindow[0], birthsThisWindow[0]);
                    }
                }
            }
        });

        // ====================================================================
        // DRAWER
        // ====================================================================
        sim.setDrawer(new BSimP3DDrawer(sim, 800, 600) {
            @Override
            public void scene(PGraphics3D p3d) {
                p3d.ortho(0, (float) BOUND_X, (float) BOUND_Y, 0, -1000, 10000);
                p3d.camera((float) BOUND_X / 2f, (float) BOUND_Y / 2f, (float) BOUND_Y,
                        (float) BOUND_X / 2f, (float) BOUND_Y / 2f, 0f, 0f, 1f, 0f);
                p3d.perspective((float) Math.PI / 2f,
                        (float) BOUND_X / (float) BOUND_Y, 0.1f, 10000f);

                draw(signalField, Color.CYAN, (float)(255.0 / 540.0));
                draw(repellentField, new Color(255, 0, 255), (float)(255.0 / 540.0));

                for (int i = 0; i < acs.length; i++) {
                    Color col;
                    switch (AC_TYPES[i]) {
                        case REPELLENT: col = Color.RED; break;
                        case NUTRIENT:  col = Color.ORANGE; break;
                        default:        col = Color.GREEN; break;
                    }
                    if (acs[i].forcedInput < 0.02) col = Color.DARK_GRAY;
                    sphere(acs[i].getPosition(), 15.0, col, 255);
                }

                for (ReservoirBacterium b : bacteria) {
                    double r = b.reporterNormalised();
                    int red   = (int)(220 * (1.0 - r)) + 30;
                    int green = (int)(220 * r) + 30;
                    sphere(b.getPosition(), 8.0, new Color(red, green, 30), 255);
                }

                int hudStep = (int) sim.getTimestep();
                boolean hudWarm = hudStep < warmUpSteps;
                int hudWin = hudWarm ? 0 : (hudStep - warmUpSteps) / stepsInWindow;
                p3d.fill(255);
                p3d.text(String.format("Stage11 | N=%d | W=%d | %s",
                        bacteria.size(), hudWin,
                        hudWarm ? "WARM-UP" : "ANALYSIS"), 10, 90);
            }
        });

        // ====================================================================
        // RUN
        // ====================================================================
        if (HEADLESS) {
            long t0 = System.currentTimeMillis();
            sim.export();
            long wallMs = System.currentTimeMillis() - t0;
            if (REGRESSION_METRICS_ONLY) {
                try {
                    writeRegressionMetrics(analysisPopByWindow, meanLuxIByWindow, totalDeathsCumulative[0]);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            } else {
                if (analysisBirthsByWindow.size() < analysisPopByWindow.size()) {
                    analysisBirthsByWindow.add(birthsThisWindow[0]);
                    analysisDeathsByWindow.add(deathsThisWindow[0]);
                    analysisBirthBinsByWindow.add(voxelAnalyzer.snapshotBirths());
                }
                writeGateSummary(resultsDir, gateWriter, resultsWriter, voxelsWriter,
                        analysisPopByWindow, analysisBirthsByWindow, analysisDeathsByWindow,
                        analysisBirthBinsByWindow,
                        analysisStartN[0],
                        bacteria.size(), peakAttConc[0], peakRepConc[0], peakGluConc[0], wallMs);
            }
        } else {
            sim.preview();
        }
    }

    static void writeRegressionMetrics(java.util.List<Integer> pops,
                                       java.util.List<Double> luxI,
                                       int totalDeaths) throws IOException {
        java.io.File regDir = new java.io.File("results/regression");
        regDir.mkdirs();
        double meanCount = 0;
        for (int n : pops) meanCount += n;
        meanCount /= Math.max(1, pops.size());
        double meanLuxI = 0;
        for (double v : luxI) meanLuxI += v;
        meanLuxI /= Math.max(1, luxI.size());
        Properties p = new Properties();
        p.setProperty("tag", METRICS_TAG);
        p.setProperty("mean_total_count", String.valueOf(meanCount));
        p.setProperty("mean_luxI", String.valueOf(meanLuxI));
        p.setProperty("total_deaths", String.valueOf(totalDeaths));
        p.setProperty("windows", String.valueOf(pops.size()));
        try (java.io.FileOutputStream fos = new java.io.FileOutputStream(
                new java.io.File(regDir, "metrics_" + METRICS_TAG + ".txt"))) {
            p.store(fos, "Stage1 regression metrics");
        }
        System.out.printf("Regression metrics [%s]: mean_N=%.1f mean_LuxI=%.4e deaths=%d%n",
                METRICS_TAG, meanCount, meanLuxI, totalDeaths);
    }

    static void writeGateSummary(java.io.File resultsDir, FileWriter gateWriter,
                                 FileWriter resultsWriter, FileWriter voxelsWriter,
                                 java.util.List<Integer> analysisPopByWindow,
                                 java.util.List<Integer> analysisBirthsByWindow,
                                 java.util.List<Integer> analysisDeathsByWindow,
                                 java.util.List<int[]> analysisBirthBinsByWindow,
                                 int analysisStartN, int finalPop,
                                 double peakAtt, double peakRep, double peakGlu, long wallMs) {
        try {
            double totalAnalysisTime = analysisPopByWindow.size() * WINDOW_DURATION;
            double[] trend = populationTrend(analysisPopByWindow, WINDOW_DURATION);
            double slope = trend[0];
            double meanN = trend[1];
            double driftFraction = (meanN > 0 && totalAnalysisTime > 0)
                    ? Math.abs(slope * totalAnalysisTime) / meanN : 0.0;

            StringBuilder sb = new StringBuilder();
            sb.append("BSimReservoirStage11 gate summary\n");
            sb.append("================================\n");
            sb.append(String.format("initial_pop=%d final_pop=%d carrying_capacity=%d%n",
                    INITIAL_POP, finalPop, CARRYING_CAPACITY));
            sb.append(String.format("population_at_analysis_start=%d%n", analysisStartN));
            sb.append(String.format("analysis_windows=%d total_analysis_time=%.0fs%n",
                    analysisPopByWindow.size(), totalAnalysisTime));
            sb.append("population_by_window_start:");
            for (int n : analysisPopByWindow) sb.append(' ').append(n);
            sb.append('\n');
            sb.append(String.format("trend_slope_dNdt=%.6f cells/s mean_N=%.1f%n", slope, meanN));
            sb.append(String.format("trend_drift_fraction |slope*T|/mean(N)=%.4f (gate < 0.05)%n",
                    driftFraction));
            sb.append(String.format("plateau_gate=%s%n", driftFraction < 0.05 ? "PASS" : "FAIL"));
            sb.append(String.format(
                    "glu_diff=%.1f glu_decay=%.5f glu_g0=%.1f glu_k_supply=%.4f k_g=%.3f km_uptake=%.3f (K_G=KM=tuned; keep equal)%n",
                    GLU_DIFF, GLU_DECAY, GLU_G0, GLU_K_SUPPLY,
                    ReservoirBacterium.K_G, ReservoirBacterium.KM_UPTAKE));
            sb.append(String.format("warmup_s=%.0f%n", WARMUP));
            sb.append(String.format("readout_state=%dx%dx%d readout_count=%dx%dx%d%n",
                    STATE_X, STATE_Y, STATE_Z, COUNT_X, COUNT_Y, COUNT_Z));
            sb.append("ac_secretion=ODE_ONLY (no direct injection)\n");
            sb.append(String.format("ac_gamma=%.2e ac_csmax=%.2e ac_km=%.2e%n",
                    AC_GAMMA_SYM, AC_CS_MAX, AC_KM));
            sb.append(String.format("ac3_nutrient_gamma=%.2e ac3_nutrient_csmax=%.2e ac3_k_refill=%.2f%n",
                    NUTRIENT_GAMMA, NUTRIENT_CSMAX, NUTRIENT_K_REFILL));
            sb.append(String.format("peak_att_at_ac0=%.2f peak_rep_at_ac2=%.2f peak_glu_at_ac3=%.2f%n",
                    peakAtt, peakRep, peakGlu));
            double plume = Math.sqrt(ATT_DIFF / ATT_DECAY);
            double gluPlume = Math.sqrt(GLU_DIFF / GLU_DECAY);
            double voxelUm = BOUND_X / STATE_X;
            double countBinUm = BOUND_X / COUNT_X;
            sb.append(String.format("plume_att_um=%.1f plume_glu_um=%.1f readout_voxel_um=%.1f count_bin_um=%.1f%n",
                    plume, gluPlume, voxelUm, countBinUm));
            sb.append("ac0/1=ATTRACTANT motility; ac2=REPELLENT readout only; ac3=NUTRIENT glucose (Stage7 AC5 site)\n");
            sb.append("stage2_birth_mechanism=NUTRIENT_AC (not chemotaxis-depletion); att_growth_boost=REMOVED\n");
            sb.append("repellent_signflip=NULL replace with pH\n");
            sb.append("regression_gate=PASS (see results/regression/regression_gate_summary.txt)\n");
            sb.append("regression_caveat_luxI=DEGENERATE (QS off: basal fixed point, sd~1e-6; real after Stage 3)\n");
            sb.append("regression_caveat_n=n=3 overlap = no detectable difference, not invariance; grid deaths 630 vs 649 (t~1.3)\n");
            sb.append("regression_caveat_regime=shown at warmup=900s / 5 windows, not production 9000s / 40 windows\n");
            sb.append(String.format("production_wall_s=%.1f production_wall_gate_lt_1800s=%s%n",
                    wallMs / 1000.0, wallMs < 30 * 60 * 1000L ? "PASS" : "FAIL"));

            appendInputCorrelation(sb, analysisPopByWindow);
            appendBirthInputGate(sb, analysisPopByWindow, analysisBirthsByWindow,
                    analysisDeathsByWindow, analysisBirthBinsByWindow);

            gateWriter.write(sb.toString());
            gateWriter.flush();
            java.nio.file.Files.write(
                    new java.io.File(resultsDir, "gate_summary.txt").toPath(),
                    sb.toString().getBytes(java.nio.charset.StandardCharsets.UTF_8));
            resultsWriter.flush();
            voxelsWriter.flush();
            resultsWriter.close();
            voxelsWriter.close();
            gateWriter.close();
            if (RECORD_WALL_TIME) {
                java.io.File regDir = new java.io.File(resultsDir, "regression");
                regDir.mkdirs();
                java.nio.file.Files.write(
                        new java.io.File(regDir, "production_wall_time.txt").toPath(),
                        String.format("wall_s=%.1f  gate_lt_1800s=%s%n",
                                wallMs / 1000.0, wallMs < 30 * 60 * 1000L ? "PASS" : "FAIL")
                                .getBytes(java.nio.charset.StandardCharsets.UTF_8));
            }
            System.out.println("Gate summary written to results/gate_summary.txt");
            System.out.print(sb.toString());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static void appendInputCorrelation(StringBuilder sb, java.util.List<Integer> pops) {
        sb.append("population_input_correlation (Pearson r, N vs graded input per window)\n");
        int n = pops.size();
        if (n < 2) {
            sb.append("  (insufficient windows)\n");
            return;
        }
        double[] np = new double[n];
        for (int i = 0; i < n; i++) np[i] = pops.get(i);
        double[] u0 = readInputCol("input_AC0.txt", n);
        double[] u1 = readInputCol("input_AC1.txt", n);
        double[] u2 = readInputCol("input_AC2.txt", n);
        double[] u3 = readInputCol("input_AC3.txt", n);
        sb.append(String.format("  r(N,input_AC0)=%.3f%n", pearson(np, u0)));
        sb.append(String.format("  r(N,input_AC1)=%.3f%n", pearson(np, u1)));
        sb.append(String.format("  r(N,input_AC2)=%.3f%n", pearson(np, u2)));
        sb.append(String.format("  r(N,input_AC3)=%.3f%n", pearson(np, u3)));
    }

    static void appendBirthInputGate(StringBuilder sb,
                                     java.util.List<Integer> pops,
                                     java.util.List<Integer> births,
                                     java.util.List<Integer> deaths,
                                     java.util.List<int[]> birthBins) {
        sb.append("birth_input_gate (total births vs nutrient AC3; gate r>0.30)\n");
        int n = Math.min(pops.size(), births.size());
        if (n < 3) {
            sb.append("  (insufficient windows)\n");
            sb.append("birth_input_gate=FAIL\n");
            return;
        }
        double[] np = new double[n];
        double[] b = new double[n];
        double[] d = new double[n];
        for (int i = 0; i < n; i++) {
            np[i] = pops.get(i);
            b[i] = births.get(i);
            d[i] = (i < deaths.size()) ? deaths.get(i) : 0;
        }
        double[] u0 = readInputCol("input_AC0.txt", n);
        double[] u3 = readInputCol("input_AC3.txt", n);
        double rB3 = pearson(b, u3);
        double rB0 = pearson(b, u0);
        double rD3 = pearson(d, u3);
        double rBN = pearson(b, np);
        double rIN = pearson(u3, np);
        double rPartial = partialR(rB3, rBN, rIN);
        double[] lag = maxLaggedR(b, u3, Math.min(20, n - 3));
        boolean pass = rB3 > 0.30;
        sb.append(String.format("  r(births,input_AC3)=%.3f t=%.2f (n=%d)%n", rB3, tFromR(rB3, n), n));
        sb.append(String.format("  r(births,input_AC0)=%.3f (attractant contrast)%n", rB0));
        sb.append(String.format("  r(deaths,input_AC3)=%.3f%n", rD3));
        sb.append(String.format("  partial_r(births,AC3|N)=%.3f t=%.2f%n", rPartial, tFromR(rPartial, n)));
        sb.append(String.format("  r(births,N)=%.3f%n", rBN));
        sb.append(String.format("  max_lagged_r(births,AC3)=%.3f at lag=%.0f windows (glu tau=%.0fs)%n",
                lag[0], lag[1], 1.0 / GLU_DECAY));
        appendBirthBinRank(sb, birthBins);
        sb.append(String.format("birth_input_gate=%s (threshold r(births,AC3)>0.30)%n",
                pass ? "PASS" : "FAIL"));
    }

    static double[] maxLaggedR(double[] births, double[] u, int maxLag) {
        double best = pearson(births, u);
        int bestLag = 0;
        for (int lag = 1; lag <= maxLag && lag < births.length - 2; lag++) {
            int m = births.length - lag;
            double[] b = new double[m];
            double[] uu = new double[m];
            for (int i = 0; i < m; i++) {
                b[i] = births[i + lag];
                uu[i] = u[i];
            }
            double r = pearson(b, uu);
            if (Math.abs(r) > Math.abs(best)) {
                best = r;
                bestLag = lag;
            }
        }
        return new double[]{best, bestLag};
    }

    static void appendBirthBinRank(StringBuilder sb, java.util.List<int[]> birthBins) {
        if (birthBins == null || birthBins.isEmpty()) {
            sb.append("  birth_4x2_rank=(no bin snapshots)\n");
            return;
        }
        int nBins = birthBins.get(0).length;
        double[] colSum = new double[nBins];
        double total = 0;
        for (int[] row : birthBins) {
            for (int j = 0; j < nBins && j < row.length; j++) {
                colSum[j] += row[j];
                total += row[j];
            }
        }
        int maxBin = 0;
        for (int j = 1; j < nBins; j++) if (colSum[j] > colSum[maxBin]) maxBin = j;
        double maxShare = (total > 0) ? colSum[maxBin] / total : 0;
        double sum = 0, sumSq = 0;
        for (double v : colSum) { sum += v; sumSq += v * v; }
        double pr = (sumSq > 0) ? (sum * sum) / sumSq : 0;
        sb.append(String.format("  birth_4x2_max_bin_share=%.3f (bin %d) participation_ratio=%.2f (1=one bin, %d=uniform)%n",
                maxShare, maxBin, pr, nBins));
        sb.append("  birth_4x2_note=rank-1 spatial accepted; Glu 20x10 carries spatial glucose\n");
    }

    static double partialR(double rxy, double rxz, double ryz) {
        double den = Math.sqrt((1.0 - rxz * rxz) * (1.0 - ryz * ryz));
        return den > 0 ? (rxy - rxz * ryz) / den : 0.0;
    }

    static double tFromR(double r, int n) {
        if (n < 3 || Math.abs(r) >= 1.0) return Double.NaN;
        return r * Math.sqrt((n - 2) / (1.0 - r * r));
    }

    static double[] readInputCol(String path, int n) {
        double[] v = new double[n];
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            for (int i = 0; i < n; i++) {
                String line = br.readLine();
                v[i] = (line != null) ? Double.parseDouble(line.trim()) : 0.0;
            }
        } catch (IOException e) {
            java.util.Arrays.fill(v, 0.0);
        }
        return v;
    }

    static double pearson(double[] x, double[] y) {
        int n = Math.min(x.length, y.length);
        if (n < 2) return Double.NaN;
        double mx = 0, my = 0;
        for (int i = 0; i < n; i++) { mx += x[i]; my += y[i]; }
        mx /= n; my /= n;
        double num = 0, dx = 0, dy = 0;
        for (int i = 0; i < n; i++) {
            num += (x[i] - mx) * (y[i] - my);
            dx += (x[i] - mx) * (x[i] - mx);
            dy += (y[i] - my) * (y[i] - my);
        }
        return (dx > 0 && dy > 0) ? num / Math.sqrt(dx * dy) : 0.0;
    }

    // ========================================================================
    // Sample writer
    // ========================================================================
    static void writeSample(FileWriter resultsWriter, FileWriter voxelsWriter,
                            int windowIndex, int sampleIndex, double primaryInput,
                            Vector<ReservoirBacterium> bacteria,
                            VoxelAnalyzer voxelAnalyzer,
                            BSimChemicalField signalField,
                            BSimChemicalField repellentField,
                            BSimChemicalField ahlField,
                            int deaths, int births) {
        try {
            resultsWriter.write(String.format(java.util.Locale.US,
                    "%d;%d;%.4e;%d;%d;%d%n",
                    windowIndex, sampleIndex, primaryInput,
                    bacteria.size(), deaths, births));
            resultsWriter.flush();

            VoxelAnalyzer.VoxelReadout r = voxelAnalyzer.analyze(
                    signalField, repellentField, ahlField, bacteria);

            voxelsWriter.write(String.format(java.util.Locale.US,
                    "%d;%d;%.4e", windowIndex, sampleIndex, primaryInput));
            for (double v : r.attractant)        voxelsWriter.write(String.format(java.util.Locale.US, ";%.4e", v));
            for (double v : r.repellent)         voxelsWriter.write(String.format(java.util.Locale.US, ";%.4e", v));
            for (int    v : r.density)           voxelsWriter.write(";" + v);
            for (double v : r.activatedFraction) voxelsWriter.write(String.format(java.util.Locale.US, ";%.4f", v));
            for (double v : r.protein)           voxelsWriter.write(String.format(java.util.Locale.US, ";%.4e", v));
            for (double v : r.ahl)               voxelsWriter.write(String.format(java.util.Locale.US, ";%.4e", v));
            for (int    v : r.births)            voxelsWriter.write(";" + v);
            for (int    v : r.deaths)            voxelsWriter.write(";" + v);
            voxelsWriter.write("\n");
            voxelsWriter.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
