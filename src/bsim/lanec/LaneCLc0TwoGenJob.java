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
 * LC0_TWO_GENERATION. Replication only. Two generations on the LC0 clock.
 * Does not rewrite LC0. No PDE, no Hertzian, no death. Not NARMA.
 */
public final class LaneCLc0TwoGenJob {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneC_LivingClocks/PROTOCOL_TWO_GEN.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneC_LivingClocks/configs/protocol_lc0_two_gen.json");
    private static final Path RESULTS = ROOT.resolve(
            "examples/LaneC_LivingClocks/results");

    private LaneCLc0TwoGenJob() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        refuseCapacity(args);
        requireFrozen();
        Files.createDirectories(RESULTS);

        System.out.println(LaneCLc0TwoGenIdentity.OBJECT + " "
                + LaneCLc0TwoGenIdentity.GATE
                + " device=" + LaneCLc0TwoGenIdentity.DEVICE
                + " parent=" + LaneCLc0TwoGenIdentity.PARENT_GATE);
        System.out.println("Not a restage of LC0. Not Lane A. Not Lane B. "
                + "Not P0. Not NARMA. Not death.");
        System.out.printf(Locale.US,
                "T_div=%.10f s  T=%.10f s  N0=%d  seed=%d  birth_rate=ln2/T%n",
                LaneCLc0TwoGenIdentity.T_DIV_S,
                LaneCLc0TwoGenIdentity.T_HORIZON_S,
                LaneCLc0TwoGenIdentity.N0,
                LaneCLc0TwoGenIdentity.RNG_SEED);

        Arm grow = runArm("LC0_TWO_GROW", true);
        Arm off = runArm("LC0_TWO_OFF", false);
        writeArmCsv("lc0_two_gen_grow.csv", grow);
        writeArmCsv("lc0_two_gen_off.csv", off);
        writeWindows("lc0_two_gen_windows.csv", grow, off);
        writeSummary(grow, off);
    }

    private static Arm runArm(String name, boolean growOn) {
        System.out.printf(Locale.US,
                "%s growth=%s death=OFF motility=OFF chemistry=OFF seed=%d%n",
                name, growOn ? "ON" : "OFF", LaneCLc0TwoGenIdentity.RNG_SEED);

        BSimRandom rng = new BSimRandom(LaneCLc0TwoGenIdentity.RNG_SEED);
        List<Cell> cells = new ArrayList<>();
        for (int i = 0; i < LaneCLc0TwoGenIdentity.N0; i++) {
            double age = LaneCLc0TwoGenIdentity.stableAge(rng.nextDouble());
            cells.add(new Cell(
                    LaneCLc0TwoGenIdentity.founderX(i),
                    LaneCLc0TwoGenIdentity.founderY(i),
                    0.5 * LaneCLc0TwoGenIdentity.BOX_Z_UM,
                    age));
        }

        Arm arm = new Arm(name, growOn);
        arm.capFired = false;
        record(arm, 0.0, cells.size(), cells.size());

        int birthsInWindow = 0;
        int windowIdx = 0;
        double t = 0.0;
        int nEver = LaneCLc0TwoGenIdentity.N0;
        while (t < LaneCLc0TwoGenIdentity.T_HORIZON_S - 1.0e-12) {
            double dt = Math.min(LaneCLc0TwoGenIdentity.DT_S,
                    LaneCLc0TwoGenIdentity.T_HORIZON_S - t);
            if (growOn) {
                int nBefore = cells.size();
                for (int i = 0; i < nBefore; i++) {
                    Cell c = cells.get(i);
                    c.age += dt;
                    c.length = LaneCLc0TwoGenIdentity.ELL0_UM
                            * Math.exp(LaneCLc0TwoGenIdentity.ALPHA_PER_S * c.age);
                    if (c.length + 1.0e-15 >= LaneCLc0TwoGenIdentity.ELL_DIV_UM) {
                        nEver++;
                        if (nEver > LaneCLc0TwoGenIdentity.N_EVER_COMPUTE_CAP) {
                            arm.capFired = true;
                            System.out.println("SCOPE_NOTE compute-cap N_ever>"
                                    + LaneCLc0TwoGenIdentity.N_EVER_COMPUTE_CAP
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

            double nextBoundary = (windowIdx + 1) * LaneCLc0TwoGenIdentity.WINDOW_S;
            if (windowIdx < LaneCLc0TwoGenIdentity.N_FULL_WINDOWS
                    && t + 1.0e-12 >= nextBoundary) {
                arm.windows.add(new Window(windowIdx,
                        windowIdx * LaneCLc0TwoGenIdentity.WINDOW_S,
                        nextBoundary, birthsInWindow, true));
                birthsInWindow = 0;
                windowIdx++;
            }

            if (isObservation(t)) {
                record(arm, t, cells.size(), nEver);
            }
        }

        if (windowIdx == LaneCLc0TwoGenIdentity.N_FULL_WINDOWS) {
            arm.windows.add(new Window(windowIdx,
                    LaneCLc0TwoGenIdentity.N_FULL_WINDOWS
                            * LaneCLc0TwoGenIdentity.WINDOW_S,
                    LaneCLc0TwoGenIdentity.T_HORIZON_S, birthsInWindow, false));
        }

        if (arm.obsT.isEmpty()
                || Math.abs(arm.obsT.get(arm.obsT.size() - 1)
                - LaneCLc0TwoGenIdentity.T_HORIZON_S) > 1.0e-6) {
            record(arm, LaneCLc0TwoGenIdentity.T_HORIZON_S, cells.size(), nEver);
        }

        arm.nEnd = cells.size();
        arm.nEverEnd = nEver;
        System.out.printf(Locale.US, "%s N_end=%d N_ever=%d births=%d%n",
                name, arm.nEnd, arm.nEverEnd, arm.birthsTotal);
        return arm;
    }

    private static boolean isObservation(double t) {
        if (Math.abs(t - LaneCLc0TwoGenIdentity.T_HORIZON_S) < 1.0e-9) {
            return true;
        }
        double w = LaneCLc0TwoGenIdentity.WINDOW_S;
        if (t < 0 || t > LaneCLc0TwoGenIdentity.N_FULL_WINDOWS * w + 1.0e-9) {
            return false;
        }
        double k = t / w;
        return Math.abs(k - Math.round(k)) < 1.0e-9;
    }

    private static void record(Arm arm, double t, int n, int nEver) {
        arm.obsT.add(t);
        arm.obsN.add(n);
        arm.obsNever.add(nEver);
        arm.obsNpred.add(LaneCLc0TwoGenIdentity.nPred(t));
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
                && off.nEnd == LaneCLc0TwoGenIdentity.N0
                && off.nEverEnd == LaneCLc0TwoGenIdentity.N0;
        boolean visibility = maxFull >= LaneCLc0TwoGenIdentity.VISIBILITY_MIN_BIRTHS;
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

        Files.writeString(RESULTS.resolve("lc0_two_gen_summary.json"),
                String.format(Locale.US,
                "{\n  \"gate\": \"LC0_TWO_GENERATION\",\n"
                        + "  \"object\": \"LANE_C_LIVING_CLOCKS\",\n"
                        + "  \"parent_extra\": \"LC0_REPLICATION_VISIBLE\",\n"
                        + "  \"narma\": false,\n  \"death\": false,\n"
                        + "  \"motility\": false,\n  \"chemistry\": false,\n"
                        + "  \"seed\": %d,\n  \"n0\": %d,\n"
                        + "  \"t_div_s\": %.12f,\n  \"t_horizon_s\": %.12f,\n"
                        + "  \"birth_rate_law\": \"ln2/T_div\",\n"
                        + "  \"clock_delta\": %.4f,\n"
                        + "  \"grow_N_end\": %d,\n  \"grow_N_ever\": %d,\n"
                        + "  \"grow_births\": %d,\n  \"grow_max_full_window_births\": %d,\n"
                        + "  \"off_N_end\": %d,\n  \"off_N_ever\": %d,\n"
                        + "  \"off_births\": %d,\n"
                        + "  \"cap_fired\": %s,\n  \"pass\": %s,\n  \"killer\": \"%s\"\n}\n",
                LaneCLc0TwoGenIdentity.RNG_SEED, LaneCLc0TwoGenIdentity.N0,
                LaneCLc0TwoGenIdentity.T_DIV_S,
                LaneCLc0TwoGenIdentity.T_HORIZON_S,
                LaneCLc0TwoGenIdentity.CLOCK_DELTA,
                grow.nEnd, grow.nEverEnd, grow.birthsTotal, maxFull,
                off.nEnd, off.nEverEnd, off.birthsTotal,
                cap, pass, killer), StandardCharsets.UTF_8);

        System.out.printf(Locale.US, "LC0_TWO.1 honesty %s%n",
                honesty ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "LC0_TWO.2 process-off births=%d N_end=%d %s%n",
                off.birthsTotal, off.nEnd, processOff ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "LC0_TWO.3 visibility max_full_window_births=%d %s%n",
                maxFull, visibility ? "PASS" : (cap ? "SCOPE_NOTE" : "BOUND"));
        System.out.printf(Locale.US, "LC0_TWO.4 clock delta=%.2f %s%n",
                LaneCLc0TwoGenIdentity.CLOCK_DELTA, clock ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "LC0_TWO.5 no_silent_kill N_end=N_ever %s%n",
                noKill && !cap ? "PASS" : (cap ? "SCOPE_NOTE" : "FAIL"));
        if (pass) {
            System.out.println("LC0_TWO_GENERATION=PASS. Warren sizer tracks "
                    + "N0 2^{t/T} through two generations. Granddaughters divide. "
                    + "128 is not forced. Not a restage of LC0.");
        } else if (cap) {
            System.out.println("LC0_TWO_GENERATION=SCOPE_NOTE compute_cap. Not PASS.");
            System.exit(2);
        } else {
            System.out.println("LC0_TWO_GENERATION=FAIL " + killer
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
            if (nRel > LaneCLc0TwoGenIdentity.CLOCK_DELTA + 1.0e-12
                    || eRel > LaneCLc0TwoGenIdentity.CLOCK_DELTA + 1.0e-12) {
                return false;
            }
        }
        return true;
    }

    private static void refuseCapacity(String[] args) {
        String blob = String.join(" ", args == null ? new String[0] : args)
                .toLowerCase(Locale.ROOT);
        if (blob.contains("narma") || blob.contains("charc") || blob.contains("ipc")) {
            throw new IllegalArgumentException("LC0 two-gen refuses NARMA/CHARC/IPC");
        }
    }

    private static void requireFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces")
                || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException(
                    "LC0 two-gen protocol is not frozen_before_traces");
        }
        if (!json.contains("\"status_label\": \"LC0_TWO_GENERATION\"")
                || !json.contains("\"narma\": false")
                || !json.contains("\"n0\": 64")
                || !json.contains("\"seed\": 101")
                || !json.contains("\"birth_rate_law\": \"ln2/T_div\"")
                || !json.contains("\"death\": false")
                || !json.contains("\"n_full_windows\": 17")) {
            throw new IllegalStateException("LC0 two-gen must freeze "
                    + "LC0_TWO_GENERATION, N0=64, seed=101, ln2/T, death=false, "
                    + "narma=false, 17 full windows");
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
            this.length = LaneCLc0TwoGenIdentity.ELL0_UM
                    * Math.exp(LaneCLc0TwoGenIdentity.ALPHA_PER_S * age);
        }

        void resetAfterFission() {
            age = 0.0;
            length = LaneCLc0TwoGenIdentity.ELL0_UM;
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
