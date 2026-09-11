package BSimReservoirPlan;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Properties;

/** Three-seed Stage 1 field-grid regression. */
public final class RegressionGate {
    static final int[] SEEDS = {11, 22, 33};

    static final class Metrics {
        final double totalCount;
        final double meanLuxI;
        final double totalDeaths;

        Metrics(double totalCount, double meanLuxI, double totalDeaths) {
            this.totalCount = totalCount;
            this.meanLuxI = meanLuxI;
            this.totalDeaths = totalDeaths;
        }
    }

    interface Value {
        double get(Metrics metrics);
    }

    public static void main(String[] args) throws Exception {
        File regressionDir = new File("results/regression");
        if (!regressionDir.exists() && !regressionDir.mkdirs()) {
            throw new IOException("Cannot create " + regressionDir);
        }

        List<Metrics> oldGrid = new ArrayList<>();
        List<Metrics> newGrid = new ArrayList<>();
        for (int seed : SEEDS) {
            String oldTag = "old_200x100_seed" + seed;
            String newTag = "new_50x25_seed" + seed;
            writeConfig(new File(regressionDir, oldTag + ".properties"),
                    oldTag, 200, 100, seed);
            writeConfig(new File(regressionDir, newTag + ".properties"),
                    newTag, 50, 25, seed);
            oldGrid.add(run(oldTag));
            newGrid.add(run(newTag));
        }

        StringBuilder report = new StringBuilder();
        report.append("Plan Stage 1 three-seed regression gate\n");
        report.append("========================================\n");
        report.append("comparison=Stage10 field grid 200x100 vs Plan field grid 50x25\n");
        report.append("dt_s=0.050000 windows=40 total_sim_time_s=3600.000000 seeds=11,22,33\n");
        boolean countPass = compare(report, "Total_Count", oldGrid, newGrid,
                metrics -> metrics.totalCount, true);
        boolean luxOverlap = compare(report, "Mean_LuxI", oldGrid, newGrid,
                metrics -> metrics.meanLuxI, false);
        boolean deathsPass = compare(report, "Total_Deaths", oldGrid, newGrid,
                metrics -> metrics.totalDeaths, true);
        report.append("Mean_LuxI_blocking=0 caveat=may_be_degenerate_until_Plan_Stage3\n");
        report.append("Mean_LuxI_overlap_numeric=").append(luxOverlap ? 1 : 0).append('\n');
        report.append("n3_interpretation=no_detectable_difference_not_invariance\n");
        report.append("regression_gate=")
                .append(countPass && deathsPass ? "PASS" : "FAIL").append('\n');
        report.append("plan_stage3_started=0\n");

        File reportFile = new File(regressionDir, "regression_gate_summary.txt");
        Files.write(reportFile.toPath(), report.toString().getBytes(StandardCharsets.UTF_8));
        System.out.print(report);
    }

    static void writeConfig(File file, String tag, int gridX, int gridY, int seed)
            throws IOException {
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println("dt=0.05");
            writer.println("grid.x=" + gridX);
            writer.println("grid.y=" + gridY);
            writer.println("grid.z=1");
            writer.println("readout.grid.x=20");
            writer.println("readout.grid.y=10");
            writer.println("readout.grid.z=1");
            writer.println("readout.countgrid.x=4");
            writer.println("readout.countgrid.y=2");
            writer.println("readout.countgrid.z=1");
            writer.println("warmup.s=0");
            writer.println("window.duration.s=90");
            writer.println("pulse.duration.s=90");
            writer.println("sampling.duration.s=90");
            writer.println("sampling.interval.s=90");
            writer.println("num.windows=40");
            writer.println("initial.pop=100");
            writer.println("carrying.capacity=2000");
            writer.println("field.att.diff=100.0");
            writer.println("field.rep.diff=100.0");
            writer.println("field.ahl.diff=159.0");
            writer.println("field.att.decay=0.01");
            writer.println("field.rep.decay=0.01");
            writer.println("field.ahl.decay=0.000046");
            writer.println("headless=true");
            writer.println("regression.metrics.only=true");
            writer.println("rng.seed=" + seed);
            writer.println("output.metrics.tag=" + tag);
        }
    }

    static Metrics run(String tag) throws Exception {
        File config = new File("results/regression/" + tag + ".properties");
        long start = System.currentTimeMillis();
        BSimReservoirPlan.main(new String[]{config.getPath()});
        System.out.printf(Locale.US, "regression_run=%s wall_s=%.3f%n",
                tag, (System.currentTimeMillis() - start) / 1000.0);
        Properties p = new Properties();
        try (FileInputStream input = new FileInputStream(
                new File("results/regression/metrics_" + tag + ".properties"))) {
            p.load(input);
        }
        return new Metrics(
                Double.parseDouble(p.getProperty("mean_total_count")),
                Double.parseDouble(p.getProperty("mean_luxI")),
                Double.parseDouble(p.getProperty("total_deaths")));
    }

    static boolean compare(StringBuilder report, String name,
                           List<Metrics> oldValues, List<Metrics> newValues,
                           Value value, boolean blocking) {
        double[] oldData = oldValues.stream().mapToDouble(value::get).toArray();
        double[] newData = newValues.stream().mapToDouble(value::get).toArray();
        double oldMean = mean(oldData);
        double oldSd = sampleSd(oldData);
        double newMean = mean(newData);
        double newSd = sampleSd(newData);
        double oldLow = oldMean - oldSd;
        double oldHigh = oldMean + oldSd;
        double newLow = newMean - newSd;
        double newHigh = newMean + newSd;
        double tolerance = Math.max(oldSd + newSd,
                1e-4 * Math.max(Math.max(Math.abs(oldMean), Math.abs(newMean)), 1e-12));
        boolean overlap = Math.abs(oldMean - newMean) <= tolerance;
        report.append(String.format(Locale.US,
                "%s_old_mean=%.9g old_sd=%.9g old_interval=[%.9g,%.9g]%n",
                name, oldMean, oldSd, oldLow, oldHigh));
        report.append(String.format(Locale.US,
                "%s_new_mean=%.9g new_sd=%.9g new_interval=[%.9g,%.9g]%n",
                name, newMean, newSd, newLow, newHigh));
        report.append(String.format(Locale.US,
                "%s_abs_mean_difference=%.9g overlap_tolerance=%.9g overlap=%s blocking=%d%n",
                name, Math.abs(oldMean - newMean), tolerance,
                overlap ? "PASS" : "FAIL", blocking ? 1 : 0));
        return overlap;
    }

    static double mean(double[] values) {
        double sum = 0;
        for (double value : values) sum += value;
        return sum / values.length;
    }

    static double sampleSd(double[] values) {
        double mean = mean(values);
        double sum = 0;
        for (double value : values) sum += (value - mean) * (value - mean);
        return Math.sqrt(sum / Math.max(1, values.length - 1));
    }
}
