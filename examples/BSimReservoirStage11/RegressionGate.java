package BSimReservoirStage11;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

/**
 * Stage 1 numerics-invariance regression gate.
 *
 * Arm A: dt=0.05 vs dt=0.01 (grid 50x25 fixed)
 * Arm B: grid 50x25 vs 200x100 (dt=0.05 fixed)
 * 3 seeds each; compare Total_Count, Mean_LuxI, Total_Deaths (mean ± sd overlap).
 */
public class RegressionGate {

    static final int[] SEEDS = {11, 22, 33};

    static final class RunMetrics {
        final String tag;
        final double meanCount;
        final double meanLuxI;
        final double totalDeaths;

        RunMetrics(String tag, double meanCount, double meanLuxI, double totalDeaths) {
            this.tag = tag;
            this.meanCount = meanCount;
            this.meanLuxI = meanLuxI;
            this.totalDeaths = totalDeaths;
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length > 0 && args[0].equals("summarize")) {
            summarizeFromDisk();
            return;
        }
        java.io.File regDir = new java.io.File("results/regression");
        regDir.mkdirs();

        List<RunMetrics> dt005 = new ArrayList<>();
        List<RunMetrics> dt001 = new ArrayList<>();
        List<RunMetrics> grid5025 = new ArrayList<>();
        List<RunMetrics> grid200100 = new ArrayList<>();

        for (int seed : SEEDS) {
            String tA = String.format("dt005_s%d", seed);
            writeConfig(regDir, tA, 0.05, 50, 25, seed);
            RunMetrics mA = runAndParse(tA);
            dt005.add(mA);
            grid5025.add(mA);

            String tB = String.format("dt001_s%d", seed);
            writeConfig(regDir, tB, 0.01, 50, 25, seed);
            dt001.add(runAndParse(tB));

            String tD = String.format("grid200100_s%d", seed);
            writeConfig(regDir, tD, 0.05, 200, 100, seed);
            grid200100.add(runAndParse(tD));
        }

        StringBuilder report = new StringBuilder();
        report.append("Stage 1 numerics-invariance regression gate\n");
        report.append("==========================================\n");
        report.append("Protocol: 5 analysis windows, warmup=900s, initial.pop=850\n");
        report.append("Metrics: mean Total_Count & Mean_LuxI over windows; final Total_Deaths\n\n");

        boolean armA = writeArmComparison(report, "Arm A (dt)", "0.05", dt005, "0.01", dt001);
        boolean armB = writeArmComparison(report, "Arm B (grid)", "50x25", grid5025, "200x100", grid200100);

        report.append(String.format("%nregression_gate=%s%n", (armA && armB) ? "PASS" : "FAIL"));

        // Population vs input correlation from full 40-window production run
        appendInputCorrelation(report);

        // Wall-time gate from production run
        appendWallTimeGate(report);

        java.nio.file.Files.write(
                new java.io.File("results/regression/regression_gate_summary.txt").toPath(),
                report.toString().getBytes(StandardCharsets.UTF_8));
        System.out.print(report);
    }

    static void writeConfig(java.io.File dir, String tag, double dt, int gx, int gy, int seed)
            throws IOException {
        String path = new java.io.File(dir, tag + ".properties").getPath();
        try (PrintWriter pw = new PrintWriter(path)) {
            pw.println("dt=" + dt);
            pw.println("grid.x=" + gx);
            pw.println("grid.y=" + gy);
            pw.println("warmup.s=900");
            pw.println("num.windows=5");
            pw.println("initial.pop=850");
            pw.println("headless=true");
            pw.println("regression.metrics.only=true");
            pw.println("rng.seed=" + seed);
            pw.println("output.metrics.tag=" + tag);
        }
    }

    static RunMetrics runAndParse(String tag) throws Exception {
        java.io.File cfg = new java.io.File("results/regression/" + tag + ".properties");
        long t0 = System.currentTimeMillis();
        BSimReservoirStage11.main(new String[]{cfg.getPath()});
        long wallMs = System.currentTimeMillis() - t0;
        System.out.printf("  [%s] wall=%.1fs%n", tag, wallMs / 1000.0);

        java.io.File mf = new java.io.File("results/regression/metrics_" + tag + ".txt");
        if (!mf.exists()) throw new IOException("Missing metrics: " + mf);
        Properties p = new Properties();
        try (FileInputStream fis = new FileInputStream(mf)) { p.load(fis); }
        return new RunMetrics(tag,
                Double.parseDouble(p.getProperty("mean_total_count")),
                Double.parseDouble(p.getProperty("mean_luxI")),
                Double.parseDouble(p.getProperty("total_deaths")));
    }

    static boolean writeArmComparison(StringBuilder sb, String label,
                                      String nameA, List<RunMetrics> a,
                                      String nameB, List<RunMetrics> b) {
        sb.append(label).append(": ").append(nameA).append(" vs ").append(nameB).append('\n');
        boolean pass = true;
        pass &= compareMetric(sb, "Total_Count", a, b, m -> m.meanCount);
        pass &= compareMetric(sb, "Mean_LuxI", a, b, m -> m.meanLuxI);
        pass &= compareMetric(sb, "Total_Deaths", a, b, m -> m.totalDeaths);
        sb.append(String.format("  arm_pass=%s%n%n", pass ? "PASS" : "FAIL"));
        return pass;
    }

    interface MetricGetter { double get(RunMetrics m); }

    static boolean compareMetric(StringBuilder sb, String name,
                                 List<RunMetrics> a, List<RunMetrics> b,
                                 MetricGetter g) {
        double[] va = a.stream().mapToDouble(g::get).toArray();
        double[] vb = b.stream().mapToDouble(g::get).toArray();
        double ma = mean(va), sa = std(va);
        double mb = mean(vb), sb_ = std(vb);
        boolean overlap = intervalsOverlap(ma, sa, mb, sb_);
        sb.append(String.format("  %s: %s=%.6g±%.6g  %s=%.6g±%.6g  overlap=%s%n",
                name, "A", ma, sa, "B", mb, sb_, overlap ? "YES" : "NO"));
        return overlap;
    }

    static boolean intervalsOverlap(double m1, double s1, double m2, double s2) {
        double scale = Math.max(Math.max(Math.abs(m1), Math.abs(m2)), 1e-12);
        double tol = Math.max(s1 + s2, 1e-4 * scale);
        return Math.abs(m1 - m2) <= tol;
    }

    /** Recompute summary from existing metrics_*.txt (no sim re-run). */
    static void summarizeFromDisk() throws Exception {
        java.io.File regDir = new java.io.File("results/regression");
        List<RunMetrics> dt005 = new ArrayList<>();
        List<RunMetrics> dt001 = new ArrayList<>();
        List<RunMetrics> grid200100 = new ArrayList<>();
        for (int seed : SEEDS) {
            dt005.add(loadMetrics("dt005_s" + seed));
            dt001.add(loadMetrics("dt001_s" + seed));
            grid200100.add(loadMetrics("grid200100_s" + seed));
        }
        List<RunMetrics> grid5025 = dt005;

        StringBuilder report = new StringBuilder();
        report.append("Stage 1 numerics-invariance regression gate\n");
        report.append("==========================================\n");
        report.append("Protocol: 5 analysis windows, warmup=900s, initial.pop=850\n");
        report.append("Metrics: mean Total_Count & Mean_LuxI over windows; final Total_Deaths\n");
        report.append("Overlap: mean±sd intervals, or |Δmean|≤1e-4·|mean| when sd≈0\n\n");

        boolean armA = writeArmComparison(report, "Arm A (dt)", "0.05", dt005, "0.01", dt001);
        boolean armB = writeArmComparison(report, "Arm B (grid)", "50x25", grid5025, "200x100", grid200100);
        report.append(String.format("%nregression_gate=%s%n", (armA && armB) ? "PASS" : "FAIL"));
        appendInputCorrelation(report);
        appendWallTimeGate(report);

        java.nio.file.Files.write(
                new java.io.File("results/regression/regression_gate_summary.txt").toPath(),
                report.toString().getBytes(StandardCharsets.UTF_8));
        System.out.print(report);
    }

    static RunMetrics loadMetrics(String tag) throws IOException {
        Properties p = new Properties();
        try (FileInputStream fis = new FileInputStream(
                new java.io.File("results/regression/metrics_" + tag + ".txt"))) {
            p.load(fis);
        }
        return new RunMetrics(tag,
                Double.parseDouble(p.getProperty("mean_total_count")),
                Double.parseDouble(p.getProperty("mean_luxI")),
                Double.parseDouble(p.getProperty("total_deaths")));
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

    static double pearson(double[] x, double[] y) {
        int n = Math.min(x.length, y.length);
        if (n < 2) return Double.NaN;
        double mx = mean(x), my = mean(y);
        double num = 0, dx = 0, dy = 0;
        for (int i = 0; i < n; i++) {
            num += (x[i] - mx) * (y[i] - my);
            dx += (x[i] - mx) * (x[i] - mx);
            dy += (y[i] - my) * (y[i] - my);
        }
        return (dx > 0 && dy > 0) ? num / Math.sqrt(dx * dy) : 0;
    }

    static void appendInputCorrelation(StringBuilder sb) throws IOException {
        sb.append("Population vs input (40-window production run)\n");
        sb.append("----------------------------------------------\n");
        int[] pops = readPopSeries(new java.io.File("results/gate_summary.txt"));
        if (pops == null || pops.length < 5) {
            sb.append("  (population series not found in gate_summary.txt)\n\n");
            return;
        }
        int n = pops.length;
        double[] u0 = readInputCol("input_AC0.txt", n);
        double[] u1 = readInputCol("input_AC1.txt", n);
        double[] u2 = readInputCol("input_AC2.txt", n);
        double[] np = new double[n];
        for (int i = 0; i < n; i++) np[i] = pops[i];
        sb.append(String.format("  r(N, input_AC0)=%.3f%n", pearson(np, u0)));
        sb.append(String.format("  r(N, input_AC1)=%.3f%n", pearson(np, u1)));
        sb.append(String.format("  r(N, input_AC2)=%.3f%n", pearson(np, u2)));
        sb.append('\n');
    }

    static int[] readPopSeries(java.io.File gateSummary) throws IOException {
        if (!gateSummary.exists()) return null;
        try (BufferedReader br = new BufferedReader(new FileReader(gateSummary))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("population_by_window_start:")) {
                    String[] parts = line.substring("population_by_window_start:".length()).trim().split("\\s+");
                    int[] pops = new int[parts.length];
                    for (int i = 0; i < parts.length; i++) pops[i] = Integer.parseInt(parts[i]);
                    return pops;
                }
            }
        }
        return null;
    }

    static double[] readInputCol(String path, int n) throws IOException {
        double[] v = new double[n];
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            for (int i = 0; i < n; i++) {
                String line = br.readLine();
                if (line == null) { v[i] = 0; continue; }
                v[i] = Double.parseDouble(line.trim());
            }
        }
        return v;
    }

    static void appendWallTimeGate(StringBuilder sb) throws IOException {
        java.io.File wt = new java.io.File("results/regression/production_wall_time.txt");
        sb.append("Production 40-window wall time\n");
        sb.append("------------------------------\n");
        if (wt.exists()) {
            sb.append(new String(java.nio.file.Files.readAllBytes(wt.toPath()), StandardCharsets.UTF_8));
        } else {
            sb.append("  (not recorded — re-run with record_wall_time=true)\n");
        }
        sb.append("  gate: wall < 30 min, voxels < 50 MB\n\n");
    }
}
