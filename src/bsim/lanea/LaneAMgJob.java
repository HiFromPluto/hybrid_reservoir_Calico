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
 * Extra LaneA_MG: Mackey-Glass one-step vs silent, field and persist reported.
 * Not BenchA 0.260. Does not reopen LANE_A_CARRIER or generalize NARMA note.
 */
public final class LaneAMgJob {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path OCC_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_OCCUPIED_MILLIMETRE_STANDING.md");
    private static final Path CARRIER_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_CARRIER_STANDING.md");
    private static final Path AUDIT_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_CARRIER_AUDIT_STANDING.md");
    private static final Path NARMA_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_NARMA10_STANDING.md");
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/PROTOCOL_MG.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/configs/mg_protocol.json");
    private static final Path U_FILE = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/input_u_mg200.txt");
    private static final Path RESULTS = ROOT.resolve("examples/LaneA_OccupiedMillimetre/results");
    private static final Path SILENT_REUSE = RESULTS.resolve("java_LANE_A_NARMA10_SILENT.csv");
    private static final Path BENCHA_EVIDENCE = ROOT.resolve(
            "examples/BSimReservoirPlanBenchA/results/GATE_EVIDENCE.md");

    private LaneAMgJob() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        refuseForbidden(args);
        requireParentsUnchanged();
        requireProtocolFrozen();
        requireBenchAUntouchedHint();
        LaneACli cli = LaneACli.parse(args);
        Path out = cli.outDir(RESULTS);
        if (!Files.exists(SILENT_REUSE)) {
            throw new IllegalStateException(
                    "NARMA silent maps missing; cannot reuse empty dish " + LaneAConstants.MG_LABEL);
        }
        double[] u = loadFrozenU();
        Files.createDirectories(out);

        System.out.println(LaneAConstants.MG_LABEL + " " + LaneAConstants.MG_GATE);
        System.out.println("System MG one-step vs silent. Field and persist reported. "
                + "TAKEN BenchA u map. Not BenchA 0.260. Not carrier PASS. NARMA system PASS stands.");
        cli.announce(LaneAConstants.MG_LABEL);

        LaneAOccupiedDish dish = new LaneAOccupiedDish();
        System.out.println("  LANE_A_MG_DRIVEN " + LaneAConstants.MG_LABEL);
        long t0 = System.nanoTime();
        LaneAOccupiedDish.Trajectory driven = dish.run(LaneAOccupiedDish.Arm.DRIVEN, true, u, cli.seed);
        writeMaps(out.resolve("java_LANE_A_MG_DRIVEN.csv"), driven);
        double trainR = bandMean(driven, true, LaneAConstants.NARMA10_TRAIN_LO, LaneAConstants.NARMA10_TRAIN_HI);
        double expected = expectedMass(u);
        double cmdRel = rel(driven.commanded, driven.sourceAdded, expected);
        double n0Rel = Math.abs(driven.residual) / mStar(driven);
        System.out.printf(Locale.US,
                "    occupancy=%s test_mean_R=%.6g train_mean_R=%.6g mean_L=%.6g samples=%d "
                        + "cmd=%.6g src=%.6g remain=%.6g n0_rel=%.3e cmd_rel=%.3e (%.1fs) %s%n",
                driven.occupancyFlag, driven.meanROccupancy, trainR, driven.meanLOccupancy,
                driven.samples.size(), driven.commanded, driven.sourceAdded, driven.remaining,
                n0Rel, cmdRel, elapsed(t0), LaneAConstants.MG_LABEL);
        System.out.println("  LANE_A_MG_SILENT reused " + SILENT_REUSE.getFileName()
                + " (U=0 empty dish, u-independent) " + LaneAConstants.MG_LABEL);

        boolean ledgerOk = cmdRel <= LaneAConstants.LEDGER_CMD_REL
                && n0Rel <= LaneAConstants.N0_MASS_REL
                && Math.abs(driven.commanded - expected) <= 1e-6 * Math.max(expected, 1.0);
        boolean alive = "ALIVE".equals(driven.occupancyFlag)
                && driven.meanROccupancy >= LaneAConstants.OCCUPANCY_ALIVE_MIN_R;

        if (!alive) {
            System.out.println("LANE_A_MG=NOT_SCORED driven occupancy DEAD on this u. "
                    + "Do not raise J_max. Occupancy standing unchanged. Not BenchA.");
            writeStub("NOT_SCORED", driven, trainR, n0Rel, cmdRel, false);
            System.exit(2);
        }
        if (!ledgerOk) {
            throw new IllegalStateException(
                    "MG ledger failed n0_rel=" + n0Rel + " " + LaneAConstants.MG_LABEL);
        }
        writeStub("MAPS_READY", driven, trainR, n0Rel, cmdRel, true);
        System.out.println("LANE_A_MG maps ready occupancy ALIVE. Ridge is the Python checker. "
                + "Not BenchA 0.260. Occupancy remains PASS. Carrier remains FAIL. "
                + "NARMA remains system PASS.");
    }

    private static void refuseForbidden(String[] args) {
        for (String a : args) {
            String s = a.toLowerCase(Locale.ROOT);
            if (s.contains("charc") || s.contains("lorenz") || s.contains("fig4b")
                    || s.contains("waveform") || s.contains("narma10b")) {
                throw new IllegalArgumentException(
                        "Lane A MG refuses CHARC/Lorenz/Fig4b/waveform/Narma10b copy "
                                + LaneAConstants.MG_LABEL);
            }
        }
    }

    private static void requireParentsUnchanged() throws IOException {
        String occ = Files.readString(OCC_STANDING, StandardCharsets.UTF_8);
        String car = Files.readString(CARRIER_STANDING, StandardCharsets.UTF_8);
        String aud = Files.readString(AUDIT_STANDING, StandardCharsets.UTF_8);
        String nar = Files.readString(NARMA_STANDING, StandardCharsets.UTF_8);
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
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("MG PROTOCOL is not frozen_before_traces");
        }
        if (!json.contains("\"bencha_nrmse_copy\": false") || !json.contains("\"parent_narma_rewrite\": false")) {
            throw new IllegalStateException("MG must not copy BenchA NRMSE or rewrite NARMA");
        }
        if (!json.contains(LaneAConstants.MG_U_SHA256)) {
            throw new IllegalStateException("PROTOCOL u sha256 mismatch");
        }
        if (!json.contains("\"warmup_18000_identity\": false") || !json.contains("\"J_max\": 128000000.0")) {
            throw new IllegalStateException("MG must copy occupancy J_max and warmup 0");
        }
        if (!json.contains("\"occupied_mask_is_gate\": false")) {
            throw new IllegalStateException("occupied-mask must not be the system gate");
        }
        if (!json.contains("reuse_LANE_A_NARMA10_SILENT")) {
            throw new IllegalStateException("MG must freeze silent reuse of NARMA empty dish");
        }
    }

    private static void requireBenchAUntouchedHint() throws IOException {
        if (Files.exists(BENCHA_EVIDENCE)) {
            String text = Files.readString(BENCHA_EVIDENCE, StandardCharsets.UTF_8);
            if (!text.contains("0.260") && !text.contains("Mackey")) {
                throw new IllegalStateException("BenchA GATE_EVIDENCE missing MG markers; do not rewrite it");
            }
        }
    }

    private static double[] loadFrozenU() throws Exception {
        List<Double> vals = new ArrayList<Double>();
        for (String line : Files.readAllLines(U_FILE, StandardCharsets.UTF_8)) {
            String s = line.trim();
            if (s.isEmpty() || s.startsWith("#")) {
                continue;
            }
            vals.add(Double.parseDouble(s));
        }
        if (vals.size() != LaneAConstants.NARMA10_WINDOWS) {
            throw new IllegalStateException("u file has " + vals.size() + " values " + LaneAConstants.MG_LABEL);
        }
        double[] u = new double[vals.size()];
        StringBuilder payload = new StringBuilder();
        for (int i = 0; i < vals.size(); i++) {
            u[i] = vals.get(i);
            if (u[i] < 0.0 || u[i] > 0.5) {
                throw new IllegalStateException("u escaped [0, 0.5] " + LaneAConstants.MG_LABEL);
            }
            if (i > 0) {
                payload.append(',');
            }
            payload.append(String.format(Locale.US, "%.12f", u[i]));
        }
        String sha = sha256Hex(payload.toString().getBytes(StandardCharsets.US_ASCII));
        if (!LaneAConstants.MG_U_SHA256.equals(sha)) {
            throw new IllegalStateException("u sha256 drifted: " + sha + " " + LaneAConstants.MG_LABEL);
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
            h.append(",LANE_A_MG");
            w.println(h);
            for (LaneAOccupiedDish.Sample s : traj.samples) {
                if (s.ahlReadout == null) {
                    throw new IllegalStateException("MG maps missing " + LaneAConstants.MG_LABEL);
                }
                w.printf(Locale.US, "%.6f,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e",
                        s.t, s.window, s.uCommand, s.meanR, s.meanL,
                        s.remaining, s.commanded, s.sourceAdded, s.decayLoss, s.residual);
                writeVec(w, s.ahlReadout);
                writeVec(w, s.rReadout);
                writeVec(w, s.lReadout);
                w.printf(Locale.US, ",%s%n", LaneAConstants.MG_LABEL);
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
                + "  \"gate\": \"" + LaneAConstants.MG_GATE + "\",\n"
                + "  \"status_label\": \"" + LaneAConstants.MG_LABEL + "\",\n"
                + "  \"maps\": \"" + flag + "\",\n"
                + "  \"occupancy_parent\": \"PASS\",\n"
                + "  \"carrier_parent\": \"FAIL\",\n"
                + "  \"audit_parent\": \"SCOPE_NOTE\",\n"
                + "  \"narma_parent\": \"PASS\",\n"
                + "  \"occupancy\": \"" + driven.occupancyFlag + "\",\n"
                + "  \"test_mean_R\": " + driven.meanROccupancy + ",\n"
                + "  \"train_mean_R\": " + trainR + ",\n"
                + "  \"silent_maps\": \"reuse_LANE_A_NARMA10_SILENT\",\n"
                + "  \"n0_rel\": " + n0Rel + ",\n"
                + "  \"cmd_rel\": " + cmdRel + ",\n"
                + "  \"bencha_nrmse_copy\": false,\n"
                + "  \"maps_ready\": " + mapsReady + "\n"
                + "}\n";
        Files.writeString(RESULTS.resolve("lane_a_mg_maps.json"), json, StandardCharsets.UTF_8);
    }
}
