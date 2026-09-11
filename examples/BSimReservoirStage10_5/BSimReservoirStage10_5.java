package BSimReservoirStage10_5;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.awt.Color;
import java.util.Vector;

import javax.vecmath.Vector3d;

import processing.core.PGraphics3D;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.particle.BSimBacterium;

/**
 * Stage 10.5: Isolated chemotaxis test.
 *
 * Minimal setup to test whether BSim's run-tumble chemotaxis actually
 * produces measurable drift toward (attractant) or away from (repellent)
 * a single continuous AC source.
 *
 * Stripped of: QS ODE, AHL field, growth, death, dose gating, sequences.
 * Just one AC producing continuously into one field, and bacteria
 * doing run-tumble with setGoal() or repellent override.
 *
 * Usage:
 *   java ... BSimReservoirStage10_5 attractant        (default)
 *   java ... BSimReservoirStage10_5 repellent
 *   java ... BSimReservoirStage10_5 attractant preview
 *   java ... BSimReservoirStage10_5 batch                 (5 seeds × 3 modes, prod=1e6)
 *   java ... BSimReservoirStage10_5 batch_lowfield        (prod scaled → field_at_ac≈38, Stage11 strength)
 *   java ... BSimReservoirStage10_5 none              (no chemotaxis -- null control)
 */
public class BSimReservoirStage10_5 {

    // Domain (same as all stages)
    static final double BOUND_X = 1000.0;
    static final double BOUND_Y = 500.0;
    static final double BOUND_Z = 10.0;

    // Grid
    static final int GRID_X = 200;
    static final int GRID_Y = 100;
    static final int GRID_Z = 1;

    // Signal field
    static final double DIFFUSIVITY = 100.0;   // um^2/s
    static final double DECAY_RATE  = 0.01;    // 1/s
    static final double PROD_RATE_FULL = 1e6;  // molecules/s — field_at_ac ~311 at t=600s
    /** Stage11-equivalent: peak_att_at_ac0≈38 vs ~311 at full rate (38/311 scaling). */
    static final double PROD_RATE_STAGE11 = PROD_RATE_FULL * 38.0 / 311.0;

    // AC position: near corner
    static final Vector3d AC_POS = new Vector3d(100, 100, 5);

    // Bacteria: seeded within gradient region (~200 um from AC)
    // Char length = sqrt(D/k) = 100 um, so gradient is detectable within ~200-300 um
    static final int NUM_BACTERIA = 500;
    static final double SEED_X = 250.0;   // center of seed region
    static final double SEED_Y = 250.0;
    static final double SEED_SPREAD = 50.0;  // +/- spread

    // Sim
    static final double DT = 0.05;
    static final double SIM_TIME = 600.0;  // 10 minutes

    // =========================================================
    //  Bacterium with optional repellent override
    // =========================================================
    static final Vector<TestBacterium> bacteria = new Vector<>();

    enum ChemoMode { ATTRACTANT, REPELLENT, NONE }

    static class TestBacterium extends BSimBacterium {
        final int id;
        final ChemoMode mode;
        final BSimChemicalField field;
        double[] repMemory;

        TestBacterium(BSim sim, Vector3d position, int id,
                      BSimChemicalField field, ChemoMode mode) {
            super(sim, position);
            this.id = id;
            this.mode = mode;
            this.field = field;

            if (mode == ChemoMode.ATTRACTANT) {
                setGoal(field);
            } else if (mode == ChemoMode.REPELLENT) {
                // Don't call setGoal -- we override pEndRun manually
                // But we DO need repellent memory
                int memLen = sim.timesteps(shortTermMemoryDuration + longTermMemoryDuration);
                repMemory = new double[memLen];
                double initConc = field.getConc(position);
                for (int i = 0; i < repMemory.length; i++) repMemory[i] = initConc;
            }
            // NONE: no setGoal, no repMemory -- pure random walk
        }

        public boolean movingUpRepellentGradient() {
            return repellentGradientDelta() > sensitivity;
        }

        /** Signed temporal contrast (same as Stage10 ReservoirBacterium). */
        private double repellentGradientDelta() {
            double shortTermCounter = 0, longTermCounter = 0;

            System.arraycopy(repMemory, 0, repMemory, 1, repMemory.length - 1);
            repMemory[0] = field.getConc(position);

            for (int i = 0; i < repMemory.length; i++) {
                if (i < shortTermMemoryLength) {
                    shortTermCounter += repMemory[i];
                } else {
                    longTermCounter += repMemory[i];
                }
            }
            double shortTermMean = shortTermCounter / shortTermMemoryLength;
            double longTermMean  = longTermCounter / longTermMemoryLength;

            return shortTermMean - longTermMean;
        }

        @Override
        public double pEndRun() {
            if (mode == ChemoMode.ATTRACTANT) {
                // Default BSimBacterium behavior via goal field
                return super.pEndRun();
            } else if (mode == ChemoMode.REPELLENT) {
                // Stage10 sign-flip plumbing: up-repellent shorten, down-repellent extend (flee)
                double repDelta = repellentGradientDelta();
                boolean upRepellent = repDelta > sensitivity;
                boolean downRepellent = repDelta < -sensitivity;
                if (upRepellent) return pEndRunElse;
                if (downRepellent) return pEndRunUp;
                return pEndRunElse;
            }
            // NONE: baseline tumble rate always
            return pEndRunElse;
        }

        @Override
        public void action() {
            super.action();
            // No growth, no QS, no death -- pure chemotaxis
        }

        @Override
        public void updatePosition() {
            super.updatePosition();
            // Reflective Y boundaries
            double boundY = sim.getBound().y;
            if (position.y < 0.0)         position.y = -position.y;
            else if (position.y > boundY) position.y = 2.0 * boundY - position.y;
            position.y = Math.max(0.0, Math.min(boundY, position.y));
            position.z = Math.max(0.0, Math.min(sim.getBound().z, position.z));
        }
    }

    static final int[] BATCH_SEEDS = {11, 22, 33, 44, 55};

    static final class DriftMetrics {
        final String mode;
        final int seed;
        final double meanDistT0, meanDistTend;
        final int nearT0, nearTend;

        DriftMetrics(String mode, int seed, double[] t0, double[] tend) {
            this.mode = mode;
            this.seed = seed;
            this.meanDistT0 = t0[0];
            this.meanDistTend = tend[0];
            this.nearT0 = (int) t0[1];
            this.nearTend = (int) tend[1];
        }
    }

    // =========================================================
    //  Main
    // =========================================================
    public static void main(String[] args) {

        boolean batch = false;
        boolean batchLowField = false;
        boolean preview = false;
        ChemoMode mode = ChemoMode.ATTRACTANT;
        int seed = 42;

        for (String arg : args) {
            if (arg.equalsIgnoreCase("batch")) batch = true;
            else if (arg.equalsIgnoreCase("batch_lowfield")) batchLowField = true;
            else if (arg.equalsIgnoreCase("repellent")) mode = ChemoMode.REPELLENT;
            else if (arg.equalsIgnoreCase("none")) mode = ChemoMode.NONE;
            else if (arg.equalsIgnoreCase("attractant")) mode = ChemoMode.ATTRACTANT;
            else if (arg.equalsIgnoreCase("preview")) preview = true;
            else if (arg.startsWith("seed=")) seed = Integer.parseInt(arg.substring(5));
        }

        double prodRate = batchLowField ? PROD_RATE_STAGE11 : PROD_RATE_FULL;
        String statsSuffix = batchLowField ? "_lowfield" : "";

        if (batch || batchLowField) {
            java.util.List<DriftMetrics> all = new java.util.ArrayList<>();
            for (ChemoMode m : ChemoMode.values()) {
                for (int s : BATCH_SEEDS) {
                    all.add(runOnce(m, s, false, prodRate, statsSuffix));
                }
            }
            writeAggregateStats(all, statsSuffix, prodRate);
        } else {
            runOnce(mode, seed, preview, prodRate, statsSuffix);
        }
    }

    /** Run one condition; returns drift metrics (null if preview). */
    static DriftMetrics runOnce(ChemoMode mode, int seed, boolean preview,
                                double prodRate, String statsSuffix) {
        bacteria.clear();

        final ChemoMode finalMode = mode;
        String modeStr = mode.name().toLowerCase();
        java.util.Random rng = new java.util.Random(seed);

        System.out.println("Stage 10.5: Isolated chemotaxis test");
        System.out.println("  Mode: " + modeStr + "  seed: " + seed);
        System.out.println("  AC at: (" + AC_POS.x + ", " + AC_POS.y + ")");
        System.out.println("  Bacteria seeded near: (" + SEED_X + ", " + SEED_Y + ")");
        System.out.println("  N=" + NUM_BACTERIA + ", sim=" + SIM_TIME + "s, prod_rate=" + prodRate);

        String timestamp = BSimUtils.timeStamp();
        String exportPath = "./results/" + timestamp + "_" + modeStr + "_s" + seed + "/";

        // Sim
        BSim sim = new BSim();
        sim.setDt(DT);
        sim.setSimulationTime(SIM_TIME);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);

        // Signal field
        final BSimChemicalField signalField = new BSimChemicalField(sim,
                new int[]{GRID_X, GRID_Y, GRID_Z}, DIFFUSIVITY, DECAY_RATE);

        // Seed bacteria
        for (int i = 0; i < NUM_BACTERIA; i++) {
            double bx = SEED_X + (rng.nextDouble() - 0.5) * 2 * SEED_SPREAD;
            double by = SEED_Y + (rng.nextDouble() - 0.5) * 2 * SEED_SPREAD;
            double bz = BOUND_Z / 2.0;
            TestBacterium b = new TestBacterium(sim,
                    new Vector3d(bx, by, bz), i, signalField, mode);
            bacteria.add(b);
        }

        // Ticker
        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                // AC produces continuously
                signalField.addQuantity(AC_POS, prodRate * sim.getDt());

                // Field update
                signalField.update();

                // Bacteria
                for (TestBacterium b : bacteria) {
                    b.action();
                    b.updatePosition();
                }
            }
        });

        // Drawer
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

                Color fieldColor = (finalMode == ChemoMode.REPELLENT)
                        ? new Color(255, 0, 255) : Color.CYAN;
                draw(signalField, fieldColor, (float)(255.0 / 540.0));

                // AC
                Color acColor = (finalMode == ChemoMode.REPELLENT) ? Color.RED : Color.GREEN;
                sphere(AC_POS, 15.0, acColor, 255);

                // Bacteria
                for (TestBacterium b : bacteria) {
                    sphere(b.getPosition(), 8.0, Color.YELLOW, 255);
                }
            }
        });

        // Exporters
        if (!preview) {
            BSimUtils.generateDirectoryPath(exportPath);

            // Distance logger (every 1s)
            BSimLogger distLogger = new BSimLogger(sim, exportPath + "distances.csv") {
                @Override
                public void before() {
                    super.before();
                    write("time_s,mean_dist,min_dist,max_dist,median_x,median_y,"
                            + "n_near_ac,n_far,field_at_ac,mode");
                }
                @Override
                public void during() {
                    double sumDist = 0;
                    double minDist = Double.MAX_VALUE;
                    double maxDist = 0;
                    int nearAC = 0;  // within 200 um of AC
                    int far = 0;     // > 500 um from AC
                    double[] xs = new double[bacteria.size()];
                    double[] ys = new double[bacteria.size()];

                    for (int i = 0; i < bacteria.size(); i++) {
                        Vector3d p = bacteria.get(i).getPosition();
                        double dx = p.x - AC_POS.x;
                        double dy = p.y - AC_POS.y;
                        double d = Math.sqrt(dx * dx + dy * dy);
                        sumDist += d;
                        if (d < minDist) minDist = d;
                        if (d > maxDist) maxDist = d;
                        if (d < 200) nearAC++;
                        if (d > 500) far++;
                        xs[i] = p.x;
                        ys[i] = p.y;
                    }
                    java.util.Arrays.sort(xs);
                    java.util.Arrays.sort(ys);
                    int n = bacteria.size();
                    double medX = xs[n / 2];
                    double medY = ys[n / 2];

                    write(sim.getFormattedTime()
                            + "," + String.format("%.1f", sumDist / n)
                            + "," + String.format("%.1f", minDist)
                            + "," + String.format("%.1f", maxDist)
                            + "," + String.format("%.1f", medX)
                            + "," + String.format("%.1f", medY)
                            + "," + nearAC
                            + "," + far
                            + "," + String.format("%.2f", signalField.getConc(AC_POS))
                            + "," + finalMode.name());
                }
            };
            distLogger.setDt(1.0);
            sim.addExporter(distLogger);

            // Per-bacterium snapshot (every 30s)
            BSimLogger bacLogger = new BSimLogger(sim, exportPath + "bacteria.csv") {
                @Override
                public void before() {
                    super.before();
                    write("time_s,id,x,y,dist_ac,field_conc");
                }
                @Override
                public void during() {
                    String t = sim.getFormattedTime();
                    for (TestBacterium b : bacteria) {
                        Vector3d p = b.getPosition();
                        double dx = p.x - AC_POS.x;
                        double dy = p.y - AC_POS.y;
                        double d = Math.sqrt(dx * dx + dy * dy);
                        write(t + "," + b.id
                                + "," + String.format("%.1f", p.x)
                                + "," + String.format("%.1f", p.y)
                                + "," + String.format("%.1f", d)
                                + "," + String.format("%.2f", signalField.getConc(p)));
                    }
                }
            };
            bacLogger.setDt(30.0);
            sim.addExporter(bacLogger);

            // Params
            BSimLogger paramLogger = new BSimLogger(sim, exportPath + "params.csv") {
                @Override
                public void before() {
                    super.before();
                    write("parameter,value");
                    write("mode," + finalMode.name());
                    write("rng_seed," + seed);
                    write("ac_x," + AC_POS.x);
                    write("ac_y," + AC_POS.y);
                    write("seed_x," + SEED_X);
                    write("seed_y," + SEED_Y);
                    write("seed_spread," + SEED_SPREAD);
                    write("num_bacteria," + NUM_BACTERIA);
                    write("diffusivity," + DIFFUSIVITY);
                    write("decay_rate," + DECAY_RATE);
                    write("prod_rate," + prodRate);
                    write("dt," + DT);
                    write("sim_time," + SIM_TIME);
                    write("sensitivity," + 1.0);
                    write("pEndRunUp," + (1.0/1.07));
                    write("pEndRunElse," + (1.0/0.86));
                    write("shortTermMemory_s," + 1.0);
                    write("longTermMemory_s," + 3.0);
                }
                @Override
                public void during() { }
            };
            sim.addExporter(paramLogger);

            System.out.println("Exporting to " + exportPath);
            sim.export();

            appendDriftSummary(exportPath + "distances.csv", modeStr, seed, statsSuffix);
            double[] t0 = readDriftRow(exportPath + "distances.csv", true);
            double[] tend = readDriftRow(exportPath + "distances.csv", false);
            return (t0 != null && tend != null) ? new DriftMetrics(modeStr, seed, t0, tend) : null;
        } else {
            System.out.println("Preview mode");
            sim.preview();
            return null;
        }
    }

    /** Append per-seed drift line to results/stage10_5_drift_summary.txt */
    static void appendDriftSummary(String distancesCsv, String modeStr, int seed, String statsSuffix) {
        java.io.File resultsDir = new java.io.File("results");
        if (!resultsDir.exists()) resultsDir.mkdirs();
        java.io.File summary = new java.io.File(resultsDir,
                "stage10_5_drift_summary" + statsSuffix + ".txt");

        double[] t0 = readDriftRow(distancesCsv, true);
        double[] tend = readDriftRow(distancesCsv, false);
        if (t0 == null || tend == null) {
            System.err.println("Could not parse " + distancesCsv + " for drift summary");
            return;
        }

        boolean append = summary.exists();
        try (FileWriter fw = new FileWriter(summary, append)) {
            if (!append) {
                fw.write("Stage10_5 chemotaxis drift (N=500, sim=600s, replicated seeds)\n");
                fw.write("Per-run: mode seed mean_dist_t0 mean_dist_tend near_ac_t0 near_ac_tend\n");
                fw.write("near_ac = within 200 um of AC\n\n");
            }
            fw.write(String.format("%s seed=%d: mean dist %.0f->%.0f, cells near AC %d->%d%n",
                    modeStr.toUpperCase(), seed, t0[0], tend[0], (int) t0[1], (int) tend[1]));
            System.out.printf("Drift summary: %s seed=%d mean dist %.0f->%.0f, near AC %d->%d%n",
                    modeStr.toUpperCase(), seed, t0[0], tend[0], (int) t0[1], (int) tend[1]);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /** Mean ± std across seeds; also delta vs NONE for tend metrics. */
    static void writeAggregateStats(java.util.List<DriftMetrics> all, String statsSuffix,
                                      double prodRate) {
        java.io.File out = new java.io.File("results/stage10_5_drift_stats" + statsSuffix + ".txt");
        java.io.File resultsDir = new java.io.File("results");
        if (!resultsDir.exists()) resultsDir.mkdirs();

        java.util.Map<String, java.util.List<DriftMetrics>> byMode = new java.util.LinkedHashMap<>();
        for (DriftMetrics d : all) {
            byMode.computeIfAbsent(d.mode, k -> new java.util.ArrayList<>()).add(d);
        }

        try (FileWriter fw = new FileWriter(out)) {
            fw.write(String.format("Stage10_5 replicated drift stats (5 seeds, N=500, sim=600s, prod=%.0f)%n%n",
                    prodRate));
            for (String mode : new String[]{"attractant", "repellent", "none"}) {
                java.util.List<DriftMetrics> runs = byMode.get(mode);
                if (runs == null || runs.isEmpty()) continue;
                double[] distTend = runs.stream().mapToDouble(d -> d.meanDistTend).toArray();
                double[] nearTend = runs.stream().mapToDouble(d -> d.nearTend).toArray();
                fw.write(String.format("%s (n=%d):%n", mode.toUpperCase(), runs.size()));
                fw.write(String.format("  mean_dist_tend: %.1f ± %.1f%n",
                        mean(distTend), std(distTend)));
                fw.write(String.format("  near_ac_tend:   %.1f ± %.1f%n",
                        mean(nearTend), std(nearTend)));
            }
            java.util.List<DriftMetrics> att = byMode.get("attractant");
            java.util.List<DriftMetrics> rep = byMode.get("repellent");
            java.util.List<DriftMetrics> none = byMode.get("none");
            if (att != null && none != null) {
                double dAtt = mean(att.stream().mapToDouble(d -> d.meanDistTend).toArray())
                        - mean(none.stream().mapToDouble(d -> d.meanDistTend).toArray());
                double dNear = mean(att.stream().mapToDouble(d -> d.nearTend).toArray())
                        - mean(none.stream().mapToDouble(d -> d.nearTend).toArray());
                fw.write(String.format("%nATTRACTANT vs NONE (mean_dist_tend delta): %.1f%n", dAtt));
                fw.write(String.format("ATTRACTANT vs NONE (near_ac_tend delta): %.1f%n", dNear));
            }
            if (rep != null && none != null) {
                double dRep = mean(rep.stream().mapToDouble(d -> d.meanDistTend).toArray())
                        - mean(none.stream().mapToDouble(d -> d.meanDistTend).toArray());
                double dNear = mean(rep.stream().mapToDouble(d -> d.nearTend).toArray())
                        - mean(none.stream().mapToDouble(d -> d.nearTend).toArray());
                fw.write(String.format("REPELLENT vs NONE (mean_dist_tend delta): %.1f%n", dRep));
                fw.write(String.format("REPELLENT vs NONE (near_ac_tend delta): %.1f%n", dNear));
            }
            System.out.println("Aggregate stats written to results/stage10_5_drift_stats"
                    + statsSuffix + ".txt");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static double mean(double[] v) {
        double s = 0;
        for (double x : v) s += x;
        return s / v.length;
    }

    static double std(double[] v) {
        double m = mean(v);
        double s2 = 0;
        for (double x : v) s2 += (x - m) * (x - m);
        return Math.sqrt(s2 / Math.max(1, v.length - 1));
    }

    /** Returns {mean_dist, n_near_ac} from first or last data row. */
    static double[] readDriftRow(String path, boolean first) {
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            br.readLine(); // header
            String line = null;
            String last = null;
            while ((line = br.readLine()) != null) {
                if (!line.trim().isEmpty()) {
                    if (first) {
                        last = line;
                        break;
                    }
                    last = line;
                }
            }
            if (last == null) return null;
            String[] cols = last.split(",");
            return new double[]{Double.parseDouble(cols[1]), Double.parseDouble(cols[6])};
        } catch (IOException e) {
            return null;
        }
    }
}
