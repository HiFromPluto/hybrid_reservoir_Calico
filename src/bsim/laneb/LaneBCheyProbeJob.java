package bsim.laneb;

import bsim.BSim;
import bsim.BSimRandom;

import javax.vecmath.Vector3d;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;

/**
 * Single-cell CheY-P probe. B4_CHEY_PROBE. Not colony kappa. Not NARMA.
 * Parents B0–B3 FAILs kept.
 */
public final class LaneBCheyProbeJob {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/PROTOCOL_CHEY_PROBE.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/configs/protocol_chey.json");
    private static final Path RESULTS = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/results");

    private LaneBCheyProbeJob() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        refuseCapacity(args);
        requireFrozen();
        Files.createDirectories(RESULTS);

        System.out.println(LaneBCheyIdentity.OBJECT + " " + LaneBCheyIdentity.GATE
                + " device=" + LaneBCheyIdentity.DEVICE);
        System.out.println("Parents B0–B3 FAIL kept. Not colony kappa. Not NARMA.");

        Trace hold = runArm("B4_HOLD_0", 0.0, 0.0);
        Trace up = runArm("B4_STEP_UP", 0.0, LaneBCheyIdentity.DL);
        Trace down = runArm("B4_STEP_DOWN", LaneBCheyIdentity.DL, 0.0);
        writeCsv("b4_hold.csv", hold);
        writeCsv("b4_step_up.csv", up);
        writeCsv("b4_step_down.csv", down);

        double holdMean = meanY(hold, 3.0, LaneBCheyIdentity.T_TOTAL);
        double upMin = minY(up, LaneBCheyIdentity.WIN_LO, LaneBCheyIdentity.WIN_HI);
        double downMax = maxY(down, LaneBCheyIdentity.WIN_LO, LaneBCheyIdentity.WIN_HI);
        boolean holdPass = holdMean >= LaneBCheyIdentity.HOLD_LO
                && holdMean <= LaneBCheyIdentity.HOLD_HI;
        boolean upPass = upMin <= LaneBCheyIdentity.STEP_UP_MAX;
        boolean downPass = downMax >= LaneBCheyIdentity.STEP_DOWN_MIN;
        boolean pass = holdPass && upPass && downPass;
        String killer = pass ? "none" : (!holdPass ? "hold" : (!upPass ? "step_up" : "step_down"));

        System.out.printf(Locale.US, "B4.1 HOLD mean_y[3,15]=%.4f %s%n",
                holdMean, holdPass ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "B4.2 STEP_UP min_y[5.4,7]=%.4f %s%n",
                upMin, upPass ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "B4.3 STEP_DOWN max_y[5.4,7]=%.4f %s%n",
                downMax, downPass ? "PASS" : "FAIL");
        Files.writeString(RESULTS.resolve("b4_summary.json"), String.format(Locale.US,
                "{\n  \"gate\": \"B4_CHEY_PROBE\",\n  \"narma\": false,\n"
                        + "  \"pass\": %s,\n  \"killer\": \"%s\",\n"
                        + "  \"hold_mean_y\": %.8g,\n  \"step_up_min_y\": %.8g,\n"
                        + "  \"step_down_max_y\": %.8g\n}\n",
                pass, killer, holdMean, upMin, downMax), StandardCharsets.UTF_8);
        if (pass) {
            System.out.println("B4_CHEY_PROBE=PASS. ODE transduces dL=0.3. "
                    + "B3 colony FAIL is not a dead CheY. Do not start NARMA.");
        } else {
            System.out.println("B4_CHEY_PROBE=FAIL " + killer
                    + ". Port does not transduce dL=0.3. Do not raise J. Do not start NARMA.");
            System.exit(1);
        }
    }

    private static Trace runArm(String name, double l0, double l1) {
        BSim sim = new BSim();
        sim.setDt(LaneBIdentity.DT);
        sim.setBound(LaneBIdentity.LX, LaneBIdentity.LY, LaneBIdentity.LZ);
        sim.setSolid(true, true, true);
        sim.setTemperature(LaneBIdentity.T_K);
        sim.setVisc(LaneBIdentity.ETA);
        BSimRandom rng = new BSimRandom(LaneBCheyIdentity.RNG_SEED);
        sim.setRandom(rng);
        LaneBBacterium cell = new LaneBBacterium(
                sim, new Vector3d(500.0, 250.0, 5.0), rng, LaneBBacterium.Motility.MOTILE);
        int n = (int) Math.round(LaneBCheyIdentity.T_TOTAL / LaneBIdentity.DT);
        Trace tr = new Trace(n);
        for (int i = 0; i < n; i++) {
            double t = i * LaneBIdentity.DT;
            double L = t < LaneBCheyIdentity.T_ADAPT ? l0 : l1;
            cell.sampleLigand(Math.max(0.0, L), Math.max(0.0, -L));
            cell.action();
            tr.t[i] = t;
            tr.L[i] = L;
            tr.y[i] = cell.cheyP();
            tr.p[i] = cell.pEndRun();
        }
        System.out.printf(Locale.US, "%s steps=%d%n", name, n);
        return tr;
    }

    private static double meanY(Trace tr, double t0, double t1) {
        double s = 0.0;
        int n = 0;
        for (int i = 0; i < tr.t.length; i++) {
            if (tr.t[i] >= t0 && tr.t[i] <= t1) {
                s += tr.y[i];
                n++;
            }
        }
        return n == 0 ? Double.NaN : s / n;
    }

    private static double minY(Trace tr, double t0, double t1) {
        double m = Double.POSITIVE_INFINITY;
        for (int i = 0; i < tr.t.length; i++) {
            if (tr.t[i] >= t0 && tr.t[i] <= t1 && tr.y[i] < m) m = tr.y[i];
        }
        return m;
    }

    private static double maxY(Trace tr, double t0, double t1) {
        double m = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < tr.t.length; i++) {
            if (tr.t[i] >= t0 && tr.t[i] <= t1 && tr.y[i] > m) m = tr.y[i];
        }
        return m;
    }

    private static void writeCsv(String name, Trace tr) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(
                RESULTS.resolve(name), StandardCharsets.UTF_8))) {
            w.println("t,L,cheyP,p_end_run");
            for (int i = 0; i < tr.t.length; i++) {
                w.printf(Locale.US, "%.6g,%.6g,%.8g,%.8g%n", tr.t[i], tr.L[i], tr.y[i], tr.p[i]);
            }
        }
    }

    private static void refuseCapacity(String[] args) {
        String blob = String.join(" ", args == null ? new String[0] : args).toLowerCase(Locale.ROOT);
        if (blob.contains("narma") || blob.contains("charc") || blob.contains("ipc")) {
            throw new IllegalArgumentException("B4 refuses NARMA/CHARC/IPC");
        }
    }

    private static void requireFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("B4 protocol is not frozen_before_traces");
        }
        if (!json.contains("\"status_label\": \"B4_CHEY_PROBE\"") || !json.contains("\"narma\": false")) {
            throw new IllegalStateException("B4 must freeze B4_CHEY_PROBE and narma=false");
        }
        if (!json.contains("\"dL\": 0.3")) {
            throw new IllegalStateException("B4 must freeze dL=0.3");
        }
    }

    private static final class Trace {
        final double[] t;
        final double[] L;
        final double[] y;
        final double[] p;

        Trace(int n) {
            t = new double[n];
            L = new double[n];
            y = new double[n];
            p = new double[n];
        }
    }
}
