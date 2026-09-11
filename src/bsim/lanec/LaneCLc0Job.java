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
 * LC0_REPLICATION_VISIBLE. Replication only. No PDE, no Hertzian, no death.
 * Not NARMA. Not Lane A. Not Lane B occupation. Not P0 packing.
 */
public final class LaneCLc0Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneC_LivingClocks/PROTOCOL.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneC_LivingClocks/configs/protocol_lc0.json");
    private static final Path RESULTS = ROOT.resolve(
            "examples/LaneC_LivingClocks/results");

    private LaneCLc0Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        refuseCapacity(args);
        requireFrozen();
        Files.createDirectories(RESULTS);

        System.out.println(LaneCLc0Identity.OBJECT + " " + LaneCLc0Identity.GATE
                + " device=" + LaneCLc0Identity.DEVICE);
        System.out.println("Not Lane A. Not Lane B occupation. Not P0. Not NARMA.");
        System.out.printf(Locale.US, "T_div=%.10f s  N0=%d  seed=%d  birth_rate=ln2/T%n",
                LaneCLc0Identity.T_DIV_S, LaneCLc0Identity.N0, LaneCLc0Identity.RNG_SEED);

        Arm grow = runArm("LC0_GROW", true);
        Arm off = runArm("LC0_OFF", false);
        writeArmCsv("lc0_grow.csv", grow);
        writeArmCsv("lc0_off.csv", off);
        writeWindows("lc0_windows.csv", grow, off);
        writeSummary(grow, off);
    }

    private static Arm runArm(String name, boolean growOn) {
        System.out.printf(Locale.US,
                "%s growth=%s death=OFF motility=OFF chemistry=OFF seed=%d%n",
                name, growOn ? "ON" : "OFF", LaneCLc0Identity.RNG_SEED);

        BSimRandom rng = new BSimRandom(LaneCLc0Identity.RNG_SEED);
        List<Cell> cells = new ArrayList<>();
        for (int i = 0; i < LaneCLc0Identity.N0; i++) {
            double age = LaneCLc0Identity.stableAge(rng.nextDouble());
            cells.add(new Cell(
                    LaneCLc0Identity.founderX(i),
                    LaneCLc0Identity.founderY(i),
                    0.5 * LaneCLc0Identity.BOX_Z_UM,
                    age));
        }

        Arm arm = new Arm(name, growOn);
        arm.capFired = false;
        record(arm, 0.0, cells.size(), cells.size());

        int birthsInWindow = 0;
        int windowIdx = 0;
        double t = 0.0;
        int nEver = LaneCLc0Identity.N0;
        while (t < LaneCLc0Identity.T_HORIZON_S - 1.0e-12) {
            double dt = Math.min(LaneCLc0Identity.DT_S, LaneCLc0Identity.T_HORIZON_S - t);
            if (growOn) {
                int nBefore = cells.size();
                for (int i = 0; i < nBefore; i++) {
                    Cell c = cells.get(i);
                    c.age += dt;
                    c.length = LaneCLc0Identity.ELL0_UM
                            * Math.exp(LaneCLc0Identity.ALPHA_PER_S * c.age);
                    if (c.length + 1.0e-15 >= LaneCLc0Identity.ELL_DIV_UM) {
                        nEver++;
                        if (nEver > LaneCLc0Identity.N_EVER_COMPUTE_CAP) {
                            arm.capFired = true;
                            System.out.println("SCOPE_NOTE compute-cap N_ever>"
                                    + LaneCLc0Identity.N_EVER_COMPUTE_CAP
                                    + ". Not a death. Not PASS.");
                            arm.nEnd = cells.size();
                            arm.nEverEnd = nEver - 1;
                            return arm;
                        }
                        c.resetAfterFission();
                        cells.add(new Cell(c.x, c.y, c.z, 0.0));
                        birthsInWindow++;
                        arm.birthsTotal++;
                    }
                }
            }
            t += dt;

            double nextBoundary = (windowIdx + 1) * LaneCLc0Identity.WINDOW_S;
            if (windowIdx < LaneCLc0Identity.N_FULL_WINDOWS
                    && t + 1.0e-12 >= nextBoundary) {
                arm.windows.add(new Window(windowIdx,
                        windowIdx * LaneCLc0Identity.WINDOW_S,
                        nextBoundary, birthsInWindow, true));
                birthsInWindow = 0;
                windowIdx++;
            }

            if (isObservation(t)) {
                record(arm, t, cells.size(), nEver);
            }
        }

        if (windowIdx == LaneCLc0Identity.N_FULL_WINDOWS) {
            arm.windows.add(new Window(windowIdx,
                    LaneCLc0Identity.N_FULL_WINDOWS * LaneCLc0Identity.WINDOW_S,
                    LaneCLc0Identity.T_HORIZON_S, birthsInWindow, false));
        }

        if (arm.obsT.isEmpty()
                || Math.abs(arm.obsT.get(arm.obsT.size() - 1)
                - LaneCLc0Identity.T_HORIZON_S) > 1.0e-6) {
            record(arm, LaneCLc0Identity.T_HORIZON_S, cells.size(), nEver);
        }

        arm.nEnd = cells.size();
        arm.nEverEnd = nEver;
        System.out.printf(Locale.US, "%s N_end=%d N_ever=%d births=%d%n",
                name, arm.nEnd, arm.nEverEnd, arm.birthsTotal);
        return arm;
    }

    private static boolean isObservation(double t) {
        if (Math.abs(t - LaneCLc0Identity.T_HORIZON_S) < 1.0e-9) {
            return true;
        }
        double w = LaneCLc0Identity.WINDOW_S;
        if (t < 0 || t > LaneCLc0Identity.N_FULL_WINDOWS * w + 1.0e-9) {
            return false;
        }
        double k = t / w;
        return Math.abs(k - Math.round(k)) < 1.0e-9;
    }

    private static void record(Arm arm, double t, int n, int nEver) {
        arm.obsT.add(t);
        arm.obsN.add(n);
        arm.obsNever.add(nEver);
        arm.obsNpred.add(LaneCLc0Identity.nPred(t));
    }

    private static void writeArmCsv(String name, Arm arm) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(
                RESULTS.resolve(name), StandardCharsets.UTF_8))) {
            w.println("t,N,N_ever,N_pred");
            for (int i = 0; i < arm.obsT.size(); i++) {
                w.printf(Locale.US, "%.10f,%d,%d,%.8f%n",
                        arm.obsT.get(i), arm.obsN.get(i),
                        arm.obsNever.get(i), arm.obsNpred.get(i));
            }
        }
    }

    private static void writeWindows(String name, Arm grow, Arm off)
            throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(
                RESULTS.resolve(name), StandardCharsets.UTF_8))) {
            w.println("arm,window,t0,t1,births,full");
            for (Window win : grow.windows) {
                writeWin(w, grow.name, win);
            }
            for (Window win : off.windows) {
                writeWin(w, off.name, win);
            }
        }
    }

    private static void writeWin(PrintWriter w, String arm, Window win) {
        w.printf(Locale.US, "%s,%d,%.10f,%.10f,%d,%s%n",
                arm, win.index, win.t0, win.t1, win.births, win.full);
    }

    private static void writeSummary(Arm grow, Arm off) throws IOException {
        int maxFull = 0;
        for (Window win : grow.windows) {
            if (win.full) {
                maxFull = Math.max(maxFull, win.births);
            }
        }
        boolean honesty = grow.growOn && !off.growOn;
        boolean processOff = off.birthsTotal == 0
                && off.nEnd == LaneCLc0Identity.N0
                && off.nEverEnd == LaneCLc0Identity.N0;
        boolean visibility = maxFull >= LaneCLc0Identity.VISIBILITY_MIN_BIRTHS;
        boolean clock = clockOk(grow);
        boolean noKill = grow.nEnd == grow.nEverEnd && off.nEnd == off.nEverEnd;
        boolean cap = grow.capFired || off.capFired;
        boolean pass = honesty && processOff && visibility && clock && noKill && !cap;
        String killer;
        if (cap) {
            killer = "compute_cap";
        } else if (!honesty) {
            killer = "honesty";
        } else if (!processOff) {
            killer = "process_off";
        } else if (!visibility) {
            killer = "visibility_bound";
        } else if (!clock) {
            killer = "clock";
        } else if (!noKill) {
            killer = "silent_kill";
        } else {
            killer = "none";
        }

        Files.writeString(RESULTS.resolve("lc0_summary.json"), String.format(Locale.US,
                "{\n  \"gate\": \"LC0_REPLICATION_VISIBLE\",\n"
                        + "  \"object\": \"LANE_C_LIVING_CLOCKS\",\n"
                        + "  \"narma\": false,\n  \"death\": false,\n"
                        + "  \"motility\": false,\n  \"chemistry\": false,\n"
                        + "  \"seed\": %d,\n  \"n0\": %d,\n"
                        + "  \"t_div_s\": %.12f,\n  \"birth_rate_law\": \"ln2/T_div\",\n"
                        + "  \"clock_delta\": %.4f,\n"
                        + "  \"grow_N_end\": %d,\n  \"grow_N_ever\": %d,\n"
                        + "  \"grow_births\": %d,\n  \"grow_max_full_window_births\": %d,\n"
                        + "  \"off_N_end\": %d,\n  \"off_N_ever\": %d,\n"
                        + "  \"off_births\": %d,\n"
                        + "  \"cap_fired\": %s,\n  \"pass\": %s,\n  \"killer\": \"%s\"\n}\n",
                LaneCLc0Identity.RNG_SEED, LaneCLc0Identity.N0,
                LaneCLc0Identity.T_DIV_S, LaneCLc0Identity.CLOCK_DELTA,
                grow.nEnd, grow.nEverEnd, grow.birthsTotal, maxFull,
                off.nEnd, off.nEverEnd, off.birthsTotal,
                cap, pass, killer), StandardCharsets.UTF_8);

        System.out.printf(Locale.US, "LC0.1 honesty %s%n", honesty ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "LC0.2 process-off births=%d N_end=%d %s%n",
                off.birthsTotal, off.nEnd, processOff ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "LC0.3 visibility max_full_window_births=%d %s%n",
                maxFull, visibility ? "PASS" : (cap ? "SCOPE_NOTE" : "BOUND"));
        System.out.printf(Locale.US, "LC0.4 clock delta=%.2f %s%n",
                LaneCLc0Identity.CLOCK_DELTA, clock ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "LC0.5 no_silent_kill N_end=N_ever %s%n",
                noKill && !cap ? "PASS" : (cap ? "SCOPE_NOTE" : "FAIL"));
        if (pass) {
            System.out.println("LC0_REPLICATION_VISIBLE=PASS. Births countable on 300 s "
                    + "Warren clock. Growth-off flat. No death clamp. Not a chemostat.");
        } else if (cap) {
            System.out.println("LC0_REPLICATION_VISIBLE=SCOPE_NOTE compute_cap. Not PASS.");
            System.exit(2);
        } else {
            System.out.println("LC0_REPLICATION_VISIBLE=FAIL " + killer
                    + ". Do not raise lambda_S or N0.");
            System.exit(1);
        }
    }

    private static boolean clockOk(Arm arm) {
        if (arm.obsT.isEmpty()) {
            return false;
        }
        for (int i = 0; i < arm.obsT.size(); i++) {
            double pred = arm.obsNpred.get(i);
            if (pred <= 0.0) {
                return false;
            }
            double nRel = Math.abs(arm.obsN.get(i) / pred - 1.0);
            double eRel = Math.abs(arm.obsNever.get(i) / pred - 1.0);
            if (nRel > LaneCLc0Identity.CLOCK_DELTA + 1.0e-12
                    || eRel > LaneCLc0Identity.CLOCK_DELTA + 1.0e-12) {
                return false;
            }
        }
        return true;
    }

    private static void refuseCapacity(String[] args) {
        String blob = String.join(" ", args == null ? new String[0] : args)
                .toLowerCase(Locale.ROOT);
        if (blob.contains("narma") || blob.contains("charc") || blob.contains("ipc")) {
            throw new IllegalArgumentException("LC0 refuses NARMA/CHARC/IPC");
        }
    }

    private static void requireFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces")
                || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("LC0 protocol is not frozen_before_traces");
        }
        if (!json.contains("\"status_label\": \"LC0_REPLICATION_VISIBLE\"")
                || !json.contains("\"narma\": false")
                || !json.contains("\"n0\": 64")
                || !json.contains("\"seed\": 101")
                || !json.contains("\"birth_rate_law\": \"ln2/T_div\"")
                || !json.contains("\"death\": false")) {
            throw new IllegalStateException("LC0 must freeze LC0_REPLICATION_VISIBLE, "
                    + "N0=64, seed=101, ln2/T, death=false, narma=false");
        }
    }

    private static final class Cell {
        final double x;
        final double y;
        final double z;
        double age;
        double length;

        Cell(double x, double y, double z, double age) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.age = age;
            this.length = LaneCLc0Identity.ELL0_UM
                    * Math.exp(LaneCLc0Identity.ALPHA_PER_S * age);
        }

        void resetAfterFission() {
            age = 0.0;
            length = LaneCLc0Identity.ELL0_UM;
        }
    }

    private static final class Window {
        final int index;
        final double t0;
        final double t1;
        final int births;
        final boolean full;

        Window(int index, double t0, double t1, int births, boolean full) {
            this.index = index;
            this.t0 = t0;
            this.t1 = t1;
            this.births = births;
            this.full = full;
        }
    }

    private static final class Arm {
        final String name;
        final boolean growOn;
        final List<Double> obsT = new ArrayList<>();
        final List<Integer> obsN = new ArrayList<>();
        final List<Integer> obsNever = new ArrayList<>();
        final List<Double> obsNpred = new ArrayList<>();
        final List<Window> windows = new ArrayList<>();
        int birthsTotal;
        int nEnd;
        int nEverEnd;
        boolean capFired;

        Arm(String name, boolean growOn) {
            this.name = name;
            this.growOn = growOn;
        }
    }
}
