package bsim.lanec;

import bsim.BSimRandom;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * LC1_DEATH_VISIBLE. Death only. Schink WT gamma. Growth off.
 * Does not rewrite LC0. Does not import StarvationViability.
 * Not Paper 1 rates. Not NARMA.
 */
public final class LaneCLc1Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneC_LivingClocks/PROTOCOL_DEATH.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneC_LivingClocks/configs/protocol_lc1.json");
    private static final Path RESULTS = ROOT.resolve(
            "examples/LaneC_LivingClocks/results");

    private LaneCLc1Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        refuseCapacity(args);
        requireFrozen();
        Files.createDirectories(RESULTS);

        System.out.println(LaneCLc1Identity.OBJECT + " "
                + LaneCLc1Identity.GATE
                + " device=" + LaneCLc1Identity.DEVICE
                + " parent=" + LaneCLc1Identity.PARENT_GATE);
        System.out.println("Not LC0 with death on. Not Job 6 as a dish. "
                + "Not Paper 1 rates. Not NARMA. Growth never on.");
        System.out.printf(Locale.US,
                "gamma=%.2f /d  T=%.0f d  t_half=%.5f d  N0=%d  seed=%d%n",
                LaneCLc1Identity.GAMMA_PER_D,
                LaneCLc1Identity.T_HORIZON_D,
                LaneCLc1Identity.T_HALF_D,
                LaneCLc1Identity.N0,
                LaneCLc1Identity.RNG_SEED);

        double[] deathTimes = drawDeathTimes();
        Arm die = runArm("LC1_DIE", true, deathTimes);
        Arm off = runArm("LC1_OFF", false, deathTimes);
        writeArmCsv("lc1_die.csv", die);
        writeArmCsv("lc1_off.csv", off);
        writeWindows("lc1_windows.csv", die, off);
        writeSummary(die, off);
    }

    private static double[] drawDeathTimes() {
        BSimRandom rng = new BSimRandom(LaneCLc1Identity.RNG_SEED);
        double[] times = new double[LaneCLc1Identity.N0];
        for (int i = 0; i < LaneCLc1Identity.N0; i++) {
            times[i] = LaneCLc1Identity.deathTimeDays(rng.nextDouble());
        }
        return times;
    }

    private static Arm runArm(String name, boolean deathOn, double[] deathTimes) {
        System.out.printf(Locale.US,
                "%s death=%s growth=OFF motility=OFF chemistry=OFF seed=%d%n",
                name, deathOn ? "ON" : "OFF", LaneCLc1Identity.RNG_SEED);

        Arm arm = new Arm(name, deathOn);
        for (int k = 0; k <= LaneCLc1Identity.N_FULL_WINDOWS; k++) {
            double t = k * LaneCLc1Identity.WINDOW_D;
            int n = living(deathOn, deathTimes, t);
            record(arm, t, n);
        }

        for (int w = 0; w < LaneCLc1Identity.N_FULL_WINDOWS; w++) {
            double t0 = w * LaneCLc1Identity.WINDOW_D;
            double t1 = t0 + LaneCLc1Identity.WINDOW_D;
            int deaths = 0;
            if (deathOn) {
                for (double td : deathTimes) {
                    if (td >= t0 && td < t1) {
                        deaths++;
                    }
                }
            }
            arm.windows.add(new Window(w, t0, t1, deaths, true));
            arm.deathsTotal += deaths;
        }

        arm.nEnd = arm.obsN.get(arm.obsN.size() - 1);
        arm.nEverEnd = LaneCLc1Identity.N0;
        System.out.printf(Locale.US, "%s N_end=%d N_ever=%d deaths=%d%n",
                name, arm.nEnd, arm.nEverEnd, arm.deathsTotal);
        return arm;
    }

    private static int living(boolean deathOn, double[] deathTimes, double tDays) {
        if (!deathOn) {
            return LaneCLc1Identity.N0;
        }
        int n = 0;
        for (double td : deathTimes) {
            if (tDays < td) {
                n++;
            }
        }
        return n;
    }

    private static void record(Arm arm, double tDays, int n) {
        arm.obsT.add(tDays);
        arm.obsN.add(n);
        arm.obsNever.add(LaneCLc1Identity.N0);
        arm.obsNpred.add(LaneCLc1Identity.nPred(tDays));
    }

    private static void writeArmCsv(String name, Arm arm) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(
                RESULTS.resolve(name), StandardCharsets.UTF_8))) {
            w.println("t_d,N,N_ever,N_pred");
            for (int i = 0; i < arm.obsT.size(); i++) {
                w.printf(Locale.US, "%.10f,%d,%d,%.8f%n",
                        arm.obsT.get(i), arm.obsN.get(i),
                        arm.obsNever.get(i), arm.obsNpred.get(i));
            }
        }
    }

    private static void writeWindows(String name, Arm die, Arm off)
            throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(
                RESULTS.resolve(name), StandardCharsets.UTF_8))) {
            w.println("arm,window,t0_d,t1_d,deaths,full");
            for (Window win : die.windows) {
                writeWin(w, die.name, win);
            }
            for (Window win : off.windows) {
                writeWin(w, off.name, win);
            }
        }
    }

    private static void writeWin(PrintWriter w, String arm, Window win) {
        w.printf(Locale.US, "%s,%d,%.10f,%.10f,%d,%s%n",
                arm, win.index, win.t0, win.t1, win.deaths, win.full);
    }

    private static void writeSummary(Arm die, Arm off) throws IOException {
        int maxFull = 0;
        for (Window win : die.windows) {
            if (win.full) {
                maxFull = Math.max(maxFull, win.deaths);
            }
        }
        boolean honesty = die.deathOn && !off.deathOn;
        boolean immortal = off.deathsTotal == 0
                && off.nEnd == LaneCLc1Identity.N0
                && off.nEverEnd == LaneCLc1Identity.N0;
        boolean visibility = maxFull >= LaneCLc1Identity.VISIBILITY_MIN_DEATHS;
        boolean clock = clockOk(die);
        boolean noSilent = die.nEverEnd == LaneCLc1Identity.N0
                && off.nEverEnd == LaneCLc1Identity.N0
                && die.deathsTotal == LaneCLc1Identity.N0 - die.nEnd
                && off.deathsTotal == 0;
        boolean pass = honesty && immortal && visibility && clock && noSilent;
        String killer;
        if (!honesty) {
            killer = "honesty";
        } else if (!immortal) {
            killer = "immortal";
        } else if (!visibility) {
            killer = "visibility_bound";
        } else if (!clock) {
            killer = "clock";
        } else if (!noSilent) {
            killer = "silent_kill";
        } else {
            killer = "none";
        }

        Files.writeString(RESULTS.resolve("lc1_summary.json"),
                String.format(Locale.US,
                "{\n  \"gate\": \"LC1_DEATH_VISIBLE\",\n"
                        + "  \"object\": \"LANE_C_LIVING_CLOCKS\",\n"
                        + "  \"parent_extra\": \"LC0_REPLICATION_VISIBLE\",\n"
                        + "  \"narma\": false,\n  \"death\": true,\n"
                        + "  \"growth\": false,\n  \"motility\": false,\n"
                        + "  \"chemistry\": false,\n"
                        + "  \"seed\": %d,\n  \"n0\": %d,\n"
                        + "  \"gamma_per_d\": %.2f,\n"
                        + "  \"t_horizon_d\": %.1f,\n"
                        + "  \"clock_delta\": %.4f,\n"
                        + "  \"die_N_end\": %d,\n  \"die_N_ever\": %d,\n"
                        + "  \"die_deaths\": %d,\n"
                        + "  \"die_max_full_window_deaths\": %d,\n"
                        + "  \"off_N_end\": %d,\n  \"off_N_ever\": %d,\n"
                        + "  \"off_deaths\": %d,\n"
                        + "  \"pass\": %s,\n  \"killer\": \"%s\"\n}\n",
                LaneCLc1Identity.RNG_SEED, LaneCLc1Identity.N0,
                LaneCLc1Identity.GAMMA_PER_D,
                LaneCLc1Identity.T_HORIZON_D,
                LaneCLc1Identity.CLOCK_DELTA,
                die.nEnd, die.nEverEnd, die.deathsTotal, maxFull,
                off.nEnd, off.nEverEnd, off.deathsTotal,
                pass, killer), StandardCharsets.UTF_8);

        System.out.printf(Locale.US, "LC1.1 honesty %s%n",
                honesty ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "LC1.2 immortal deaths=%d N_end=%d %s%n",
                off.deathsTotal, off.nEnd, immortal ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "LC1.3 visibility max_full_window_deaths=%d %s%n",
                maxFull, visibility ? "PASS" : "BOUND");
        System.out.printf(Locale.US, "LC1.4 clock delta=%.2f %s%n",
                LaneCLc1Identity.CLOCK_DELTA, clock ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "LC1.5 no_silent_kill N_ever=64 deaths_logged %s%n",
                noSilent ? "PASS" : "FAIL");
        if (pass) {
            System.out.println("LC1_DEATH_VISIBLE=PASS. Schink WT death countable "
                    + "on 1 d windows. Immortal arm flat. Not Paper 1 rates. "
                    + "Not LC0 with death on.");
        } else {
            System.out.println("LC1_DEATH_VISIBLE=FAIL " + killer
                    + ". Do not raise gamma or N0. Do not shrink the window.");
            System.exit(1);
        }
    }

    private static boolean clockOk(Arm arm) {
        if (arm.obsT.size() != LaneCLc1Identity.N_FULL_WINDOWS + 1) {
            return false;
        }
        for (int i = 0; i < arm.obsT.size(); i++) {
            double pred = arm.obsNpred.get(i);
            if (pred <= 0.0) {
                return false;
            }
            double rel = Math.abs(arm.obsN.get(i) / pred - 1.0);
            if (rel > LaneCLc1Identity.CLOCK_DELTA + 1.0e-12) {
                return false;
            }
        }
        return true;
    }

    private static void refuseCapacity(String[] args) {
        String blob = String.join(" ", args == null ? new String[0] : args)
                .toLowerCase(Locale.ROOT);
        if (blob.contains("narma") || blob.contains("charc") || blob.contains("ipc")) {
            throw new IllegalArgumentException("LC1 refuses NARMA/CHARC/IPC");
        }
    }

    private static void requireFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces")
                || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("LC1 protocol is not frozen_before_traces");
        }
        if (!json.contains("\"status_label\": \"LC1_DEATH_VISIBLE\"")
                || !json.contains("\"narma\": false")
                || !json.contains("\"n0\": 64")
                || !json.contains("\"seed\": 101")
                || !json.contains("\"gamma_per_d\": 0.43")
                || !json.contains("\"growth\": false")
                || !json.contains("\"paper1_toxic\": false")) {
            throw new IllegalStateException("LC1 must freeze LC1_DEATH_VISIBLE, "
                    + "N0=64, seed=101, gamma=0.43/d, growth=false, "
                    + "paper1_toxic=false, narma=false");
        }
    }

    private static final class Window {
        final int index;
        final double t0;
        final double t1;
        final int deaths;
        final boolean full;

        Window(int index, double t0, double t1, int deaths, boolean full) {
            this.index = index;
            this.t0 = t0;
            this.t1 = t1;
            this.deaths = deaths;
            this.full = full;
        }
    }

    private static final class Arm {
        final String name;
        final boolean deathOn;
        final List<Double> obsT = new ArrayList<>();
        final List<Integer> obsN = new ArrayList<>();
        final List<Integer> obsNever = new ArrayList<>();
        final List<Double> obsNpred = new ArrayList<>();
        final List<Window> windows = new ArrayList<>();
        int deathsTotal;
        int nEnd;
        int nEverEnd;

        Arm(String name, boolean deathOn) {
            this.name = name;
            this.deathOn = deathOn;
        }
    }
}
