package bsim.lanea;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Extra LaneA_LORENZ: L3 k=10 AUTO_X vs silent, field, persist, AR10 reported.
 * Not L3 0.826. Does not reopen LANE_A_CARRIER or rewrite MG / NARMA.
 */
public final class LaneALorenzJob {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path OCC_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_OCCUPIED_MILLIMETRE_STANDING.md");
    private static final Path CARRIER_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_CARRIER_STANDING.md");
    private static final Path AUDIT_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_CARRIER_AUDIT_STANDING.md");
    private static final Path NARMA_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_NARMA10_STANDING.md");
    private static final Path MG_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_MG_STANDING.md");
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/PROTOCOL_LORENZ.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/configs/lorenz_protocol.json");
    private static final Path U_FILE = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/input_u_lorenz_k10.txt");
    private static final Path RESULTS = ROOT.resolve("examples/LaneA_OccupiedMillimetre/results");
    private static final Path SILENT_REUSE = RESULTS.resolve("java_LANE_A_NARMA10_SILENT.csv");
    private static final Path L3_SCOUT = ROOT.resolve(
            "examples/BSimReservoirPlanLorenzL3/results/LORENZ_L3_SCOUT.md");

    private LaneALorenzJob() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        refuseForbidden(args);
        requireParentsUnchanged();
        requireProtocolFrozen();
        requireL3UntouchedHint();
        LaneACli cli = LaneACli.parse(args);
        Path out = cli.outDir(RESULTS);
        if (!Files.exists(SILENT_REUSE)) {
            throw new IllegalStateException(
                    "NARMA silent maps missing; cannot reuse empty dish " + LaneAConstants.LORENZ_LABEL);
        }
        double[] u = loadFrozenU();
        Files.createDirectories(out);

        System.out.println(LaneAConstants.LORENZ_LABEL + " " + LaneAConstants.LORENZ_GATE);
        System.out.println("System L3 k=10 AUTO_X vs silent. Field, persist, AR10 reported. "
                + "TAKEN L3 u. Not L3 0.826. Not k=1. Not CHARC. Not carrier PASS. "
                + "NARMA/MG system PASS stand.");
        cli.announce(LaneAConstants.LORENZ_LABEL);

        LaneAOccupiedDish dish = new LaneAOccupiedDish();
        System.out.println("  LANE_A_LORENZ_DRIVEN " + LaneAConstants.LORENZ_LABEL);
        long t0 = System.nanoTime();
        LaneAOccupiedDish.Trajectory driven = dish.run(LaneAOccupiedDish.Arm.DRIVEN, true, u, cli.seed);
        writeMaps(out.resolve("java_LANE_A_LORENZ_DRIVEN.csv"), driven);
        double trainR = bandMean(driven, true, LaneAConstants.NARMA10_TRAIN_LO, LaneAConstants.NARMA10_TRAIN_HI);
        double expected = expectedMass(u);
        double cmdRel = rel(driven.commanded, driven.sourceAdded, expected);
        double n0Rel = Math.abs(driven.residual) / mStar(driven);
        System.out.printf(Locale.US,
                "    occupancy=%s test_mean_R=%.6g train_mean_R=%.6g mean_L=%.6g samples=%d "
                        + "cmd=%.6g src=%.6g remain=%.6g n0_rel=%.3e cmd_rel=%.3e (%.1fs) %s%n",
                driven.occupancyFlag, driven.meanROccupancy, trainR, driven.meanLOccupancy,
                driven.samples.size(), driven.commanded, driven.sourceAdded, driven.remaining,
                n0Rel, cmdRel, elapsed(t0), LaneAConstants.LORENZ_LABEL);
        System.out.println("  LANE_A_LORENZ_SILENT reused " + SILENT_REUSE.getFileName()
                + " (U=0 empty dish, u-independent) " + LaneAConstants.LORENZ_LABEL);

        boolean ledgerOk = cmdRel <= LaneAConstants.LEDGER_CMD_REL
                && n0Rel <= LaneAConstants.N0_MASS_REL
                && Math.abs(driven.commanded - expected) <= 1e-6 * Math.max(expected, 1.0);
        boolean alive = "ALIVE".equals(driven.occupancyFlag)
                && driven.meanROccupancy >= LaneAConstants.OCCUPANCY_ALIVE_MIN_R;

        if (!alive) {
            System.out.println("LANE_A_LORENZ=NOT_SCORED driven occupancy DEAD on this u. "
                    + "Do not raise J_max. Occupancy standing unchanged. Not L3 0.826.");
            writeStub("NOT_SCORED", driven, trainR, n0Rel, cmdRel, false);
            System.exit(2);
        }
        if (!ledgerOk) {
            throw new IllegalStateException(
                    "Lorenz ledger failed n0_rel=" + n0Rel + " " + LaneAConstants.LORENZ_LABEL);
        }
        writeStub("MAPS_READY", driven, trainR, n0Rel, cmdRel, true);
        System.out.println("LANE_A_LORENZ maps ready occupancy ALIVE. Ridge is the Python checker. "
                + "Not L3 0.826. Occupancy remains PASS. Carrier remains FAIL. "
                + "NARMA remains system PASS. MG remains system PASS.");
    }

    private static void refuseForbidden(String[] args) {
        for (String a : args) {
            String s = a.toLowerCase(Locale.ROOT);
            if (s.contains("charc") || s.contains("waveform") || s.contains("fig4b")
                    || s.contains("narma10b") || s.contains("cross_y") || s.contains("cross-y")
                    || s.contains("k=1") || s.contains("k=50") || s.contains("memory_capacity")) {
                throw new IllegalArgumentException(
                        "Lane A Lorenz refuses CHARC/waveform/Fig4b/Narma10b/CROSS_Y/k=1/k=50 "
                                + LaneAConstants.LORENZ_LABEL);
            }
        }
    }

    private static void requireParentsUnchanged() throws IOException {
        String occ = Files.readString(OCC_STANDING, StandardCharsets.UTF_8);
        String car = Files.readString(CARRIER_STANDING, StandardCharsets.UTF_8);
        String aud = Files.readString(AUDIT_STANDING, StandardCharsets.UTF_8);
        String nar = Files.readString(NARMA_STANDING, StandardCharsets.UTF_8);
        String mg = Files.readString(MG_STANDING, StandardCharsets.UTF_8);
        if (!occ.contains("**Status: PASS**") || !occ.contains("LANE_A_OCCUPIED_MILLIMETRE")) {
            throw new IllegalStateException("occupancy standing must remain PASS");
        }
        if (!car.contains("**Status: FAIL**") || !car.contains("LANE_A_CARRIER")) {
            throw new IllegalStateException("carrier standing must remain FAIL vs field");
        }
        if (!aud.contains("**Status: SCOPE_NOTE**") || !aud.contains("LANE_A_CARRIER_AUDIT")) {
            throw new IllegalStateException("audit standing must remain SCOPE_NOTE");
        }
        if (!nar.contains("**Status: PASS**") || !nar.contains("LANE_A_NARMA10")) {
            throw new IllegalStateException("NARMA standing must remain system PASS");
        }
        if (!mg.contains("**Status: PASS**") || !mg.contains("LANE_A_MG")) {
            throw new IllegalStateException("MG standing must remain system PASS");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("Lorenz PROTOCOL is not frozen_before_traces");
        }
        if (!json.contains("\"l3_nrmse_copy\": false") || !json.contains("\"parent_narma_rewrite\": false")) {
            throw new IllegalStateException("Lorenz must not copy L3 NRMSE or rewrite NARMA");
        }
        if (!json.contains("\"parent_mg_rewrite\": false") || !json.contains("\"cross_y\": false")) {
            throw new IllegalStateException("Lorenz must not rewrite MG or run CROSS_Y");
        }
        if (!json.contains(LaneAConstants.LORENZ_U_SHA256)) {
            throw new IllegalStateException("PROTOCOL u sha256 mismatch");
        }
        if (!json.contains("\"lorenz_k\": 10") || !json.contains("\"lorenz_k1\": false")) {
            throw new IllegalStateException("Lorenz must freeze k=10 AUTO_X, not k=1");
        }
        if (!json.contains("\"charc\": false") || !json.contains("\"memory_capacity\": false")) {
            throw new IllegalStateException("Lorenz must freeze not CHARC / MC");
        }
        if (!json.contains("\"warmup_18000_identity\": false") || !json.contains("\"J_max\": 128000000.0")) {
            throw new IllegalStateException("Lorenz must copy occupancy J_max and warmup 0");
        }
        if (!json.contains("\"occupied_mask_is_gate\": false")) {
            throw new IllegalStateException("occupied-mask must not be the system gate");
        }
        if (!json.contains("reuse_LANE_A_NARMA10_SILENT")) {
            throw new IllegalStateException("Lorenz must freeze silent reuse of NARMA empty dish");
        }
        if (!json.contains("\"u_one_per_line\": true")) {
            throw new IllegalStateException("Lorenz u file must be one value per line");
        }
        if (!json.contains("\"ar_m\": 10") || !json.contains("\"ar_m_search_on_test\": false")) {
            throw new IllegalStateException("AR m=10 must be frozen; no test search");
        }
    }

    private static void requireL3UntouchedHint() throws IOException {
        if (Files.exists(L3_SCOUT)) {
            String text = Files.readString(L3_SCOUT, StandardCharsets.UTF_8);
            if (!text.contains("0.8256") && !text.contains("0.826")) {
                throw new IllegalStateException("L3 scout missing living NRMSE markers; do not rewrite it");
            }
        }
    }

    private static double[] loadFrozenU() throws Exception {
        List<String> lines = Files.readAllLines(U_FILE, StandardCharsets.UTF_8);
        if (lines.size() < LaneAConstants.NARMA10_WINDOWS) {
            throw new IllegalStateException(
                    "u file looks like a one-line comma dump; write one u per line "
                            + LaneAConstants.LORENZ_LABEL);
        }
        List<Double> vals = new ArrayList<Double>();
        for (String line : lines) {
            String s = line.trim();
            if (s.isEmpty() || s.startsWith("#")) {
                continue;
            }
            if (s.contains(",")) {
                throw new IllegalStateException(
                        "u line contains a comma; one u per line " + LaneAConstants.LORENZ_LABEL);
            }
            vals.add(Double.parseDouble(s));
        }
        if (vals.size() != LaneAConstants.NARMA10_WINDOWS) {
            throw new IllegalStateException(
                    "u file has " + vals.size() + " values " + LaneAConstants.LORENZ_LABEL);
        }
        double[] u = new double[vals.size()];
        StringBuilder payload = new StringBuilder();
        for (int i = 0; i < vals.size(); i++) {
            u[i] = vals.get(i);
            if (u[i] < 0.0 || u[i] > 0.5) {
                throw new IllegalStateException("u escaped [0, 0.5] " + LaneAConstants.LORENZ_LABEL);
            }
            if (i > 0) {
                payload.append(',');
            }
            payload.append(String.format(Locale.US, "%.12f", u[i]));
        }
        String sha = sha256Hex(payload.toString().getBytes(StandardCharsets.US_ASCII));
        if (!LaneAConstants.LORENZ_U_SHA256.equals(sha)) {
            throw new IllegalStateException("u sha256 drifted: " + sha + " " + LaneAConstants.LORENZ_LABEL);
        }
        return u;
    }

    private static String sha256Hex(byte[] bytes) throws Exception {
        byte[] d = MessageDigest.getInstance("SHA-256").digest(bytes);
        StringBuilder sb = new StringBuilder(64);
        for (byte b : d) {
            sb.append(String.format("%02x", b));
        }
        return sb.toString();
    }

    private static double expectedMass(double[] u) {
        double s = 0.0;
        for (double v : u) {
            s += v;
        }
        return LaneAConstants.J_MAX * s * LaneAConstants.PULSE_S;
    }

    private static double bandMean(LaneAOccupiedDish.Trajectory traj, boolean rNotL, int lo, int hi) {
        double s = 0.0;
        int n = 0;
        for (LaneAOccupiedDish.Sample sample : traj.samples) {
            if (sample.window >= lo && sample.window <= hi) {
                s += rNotL ? sample.meanR : sample.meanL;
                n++;
            }
        }
        return n == 0 ? Double.NaN : s / n;
    }

    private static double rel(double a, double b, double floor) {
        return Math.abs(a - b) / Math.max(Math.max(Math.abs(a), Math.abs(b)), Math.max(floor, 1e-15));
    }

    private static double mStar(LaneAOccupiedDish.Trajectory traj) {
        return Math.max(Math.max(Math.abs(traj.commanded), Math.abs(traj.sourceAdded)),
                Math.max(Math.abs(traj.remaining), 1e-15));
    }

    private static double elapsed(long t0) {
        return (System.nanoTime() - t0) / 1e9;
    }

    private static void writeMaps(Path path, LaneAOccupiedDish.Trajectory traj) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            StringBuilder h = new StringBuilder(
                    "t,window,u,mean_R,mean_L,remaining,commanded,source_added,decay_loss,residual");
            for (int i = 0; i < LaneAConstants.READOUT_BINS; i++) {
                h.append(",AHL_").append(i);
            }
            for (int i = 0; i < LaneAConstants.READOUT_BINS; i++) {
                h.append(",R_").append(i);
            }
            for (int i = 0; i < LaneAConstants.READOUT_BINS; i++) {
                h.append(",L_").append(i);
            }
            h.append(",LANE_A_LORENZ");
            w.println(h);
            for (LaneAOccupiedDish.Sample s : traj.samples) {
                if (s.ahlReadout == null) {
                    throw new IllegalStateException("Lorenz maps missing " + LaneAConstants.LORENZ_LABEL);
                }
                w.printf(Locale.US, "%.6f,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e",
                        s.t, s.window, s.uCommand, s.meanR, s.meanL,
                        s.remaining, s.commanded, s.sourceAdded, s.decayLoss, s.residual);
                writeVec(w, s.ahlReadout);
                writeVec(w, s.rReadout);
                writeVec(w, s.lReadout);
                w.printf(Locale.US, ",%s%n", LaneAConstants.LORENZ_LABEL);
            }
        }
    }

    private static void writeVec(PrintWriter w, double[] v) {
        for (double x : v) {
            w.printf(Locale.US, ",%.16e", x);
        }
    }

    private static void writeStub(
            String flag,
            LaneAOccupiedDish.Trajectory driven,
            double trainR,
            double n0Rel,
            double cmdRel,
            boolean mapsReady) throws IOException {
        String json = "{\n"
                + "  \"gate\": \"" + LaneAConstants.LORENZ_GATE + "\",\n"
                + "  \"status_label\": \"" + LaneAConstants.LORENZ_LABEL + "\",\n"
                + "  \"maps\": \"" + flag + "\",\n"
                + "  \"occupancy_parent\": \"PASS\",\n"
                + "  \"carrier_parent\": \"FAIL\",\n"
                + "  \"audit_parent\": \"SCOPE_NOTE\",\n"
                + "  \"narma_parent\": \"PASS\",\n"
                + "  \"mg_parent\": \"PASS\",\n"
                + "  \"occupancy\": \"" + driven.occupancyFlag + "\",\n"
                + "  \"test_mean_R\": " + driven.meanROccupancy + ",\n"
                + "  \"train_mean_R\": " + trainR + ",\n"
                + "  \"silent_maps\": \"reuse_LANE_A_NARMA10_SILENT\",\n"
                + "  \"n0_rel\": " + n0Rel + ",\n"
                + "  \"cmd_rel\": " + cmdRel + ",\n"
                + "  \"l3_nrmse_copy\": false,\n"
                + "  \"maps_ready\": " + mapsReady + "\n"
                + "}\n";
        Files.writeString(RESULTS.resolve("lane_a_lorenz_maps.json"), json, StandardCharsets.UTF_8);
    }
}
