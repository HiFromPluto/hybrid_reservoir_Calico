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
 * Extra LaneA_WAVEFORM: Waveform2c 3-class order vs silent, field and MOMENTS.
 * Not 0.835. 400 windows. Does not reuse NARMA silent. U-only must PASS first.
 */
public final class LaneAWaveformJob {

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
    private static final Path LORENZ_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_LORENZ_STANDING.md");
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/PROTOCOL_WAVEFORM.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/configs/waveform_protocol.json");
    private static final Path U_FILE = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/input_u_waveform2c_400.txt");
    private static final Path U_ONLY = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/results/lane_a_waveform_u_only.json");
    private static final Path RESULTS = ROOT.resolve("examples/LaneA_OccupiedMillimetre/results");
    private static final Path NARMA_SILENT = RESULTS.resolve("java_LANE_A_NARMA10_SILENT.csv");
    private static final Path W2C_SCOUT = ROOT.resolve(
            "examples/BSimReservoirPlanWaveform2c/results/WAVEFORM2C_SCOUT.md");

    private LaneAWaveformJob() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        refuseForbidden(args);
        requireParentsUnchanged();
        requireProtocolFrozen();
        requireUOnlyPass();
        requireW2cUntouchedHint();
        LaneACli cli = LaneACli.parse(args);
        Path out = cli.outDir(RESULTS);
        double[] u = loadFrozenU();
        Files.createDirectories(out);

        System.out.println(LaneAConstants.WAVEFORM_LABEL + " " + LaneAConstants.WAVEFORM_GATE);
        System.out.println("System Waveform2c order vs silent and MOMENTS. Field reported. "
                + "TAKEN Waveform2c u. 400 windows. Not 0.835. Not sine/square/triangle. "
                + "Not NARMA silent reuse. NARMA/MG/Lorenz stand.");
        cli.announce(LaneAConstants.WAVEFORM_LABEL);

        LaneAOccupiedDish dish = new LaneAOccupiedDish();
        System.out.println("  LANE_A_WAVEFORM_DRIVEN " + LaneAConstants.WAVEFORM_LABEL);
        long t0 = System.nanoTime();
        LaneAOccupiedDish.Trajectory driven = dish.run(LaneAOccupiedDish.Arm.DRIVEN, true, u, cli.seed);
        writeMaps(out.resolve("java_LANE_A_WAVEFORM_DRIVEN.csv"), driven);
        double trainR = bandMean(driven, true, 64, 271);
        double expected = expectedMass(u);
        double cmdRel = rel(driven.commanded, driven.sourceAdded, expected);
        double n0Rel = Math.abs(driven.residual) / mStar(driven);
        System.out.printf(Locale.US,
                "    occupancy=%s test_mean_R=%.6g train_mean_R=%.6g mean_L=%.6g samples=%d "
                        + "cmd=%.6g src=%.6g remain=%.6g n0_rel=%.3e cmd_rel=%.3e (%.1fs) %s%n",
                driven.occupancyFlag, driven.meanROccupancy, trainR, driven.meanLOccupancy,
                driven.samples.size(), driven.commanded, driven.sourceAdded, driven.remaining,
                n0Rel, cmdRel, elapsed(t0), LaneAConstants.WAVEFORM_LABEL);

        boolean ledgerOk = cmdRel <= LaneAConstants.LEDGER_CMD_REL
                && n0Rel <= LaneAConstants.N0_MASS_REL
                && Math.abs(driven.commanded - expected) <= 1e-6 * Math.max(expected, 1.0);
        boolean alive = "ALIVE".equals(driven.occupancyFlag)
                && driven.meanROccupancy >= LaneAConstants.OCCUPANCY_ALIVE_MIN_R;

        if (!alive) {
            System.out.println("LANE_A_WAVEFORM=NOT_SCORED driven occupancy DEAD on this u. "
                    + "Do not raise J_max. Occupancy standing unchanged. Not Waveform2c 0.835.");
            writeStub("NOT_SCORED", driven, trainR, n0Rel, cmdRel, false);
            System.exit(2);
        }
        if (!ledgerOk) {
            throw new IllegalStateException(
                    "Waveform ledger failed n0_rel=" + n0Rel + " " + LaneAConstants.WAVEFORM_LABEL);
        }

        if (cli.drivenOnly) {
            System.out.println("  LANE_A_WAVEFORM_SILENT skipped (--driven-only); reuse seed 101 silent "
                    + LaneAConstants.WAVEFORM_LABEL);
        } else {
            System.out.println("  LANE_A_WAVEFORM_SILENT 400-window U=0 empty dish "
                    + LaneAConstants.WAVEFORM_LABEL);
            long t1 = System.nanoTime();
            double[] silentLen = new double[LaneAConstants.WAVEFORM_WINDOWS];
            LaneAOccupiedDish.Trajectory silent = dish.run(
                    LaneAOccupiedDish.Arm.SILENT, true, silentLen, cli.seed);
            Path silentPath = out.resolve("java_LANE_A_WAVEFORM_SILENT.csv");
            if (silentPath.toAbsolutePath().equals(NARMA_SILENT.toAbsolutePath())) {
                throw new IllegalStateException("must not reuse NARMA 200-window silent");
            }
            writeMaps(silentPath, silent);
            System.out.printf(Locale.US, "    silent occupancy=%s test_mean_R=%.6g samples=%d (%.1fs) %s%n",
                    silent.occupancyFlag, silent.meanROccupancy, silent.samples.size(),
                    elapsed(t1), LaneAConstants.WAVEFORM_LABEL);
            if (silent.samples.size() != LaneAConstants.WAVEFORM_WINDOWS * LaneAConstants.SAMPLES_PER_WINDOW) {
                throw new IllegalStateException("silent maps are not 400 windows " + LaneAConstants.WAVEFORM_LABEL);
            }
        }

        writeStub("MAPS_READY", driven, trainR, n0Rel, cmdRel, true);
        System.out.println("LANE_A_WAVEFORM maps ready occupancy ALIVE. Ridge is the Python checker. "
                + "Not Waveform2c 0.835. Occupancy remains PASS. Carrier remains FAIL. "
                + "NARMA/MG/Lorenz remain system PASS.");
    }

    private static void refuseForbidden(String[] args) {
        for (String a : args) {
            String s = a.toLowerCase(Locale.ROOT);
            if (s.contains("charc") || s.contains("fig4b") || s.contains("narma10b")
                    || s.contains("sine") || s.contains("square") || s.contains("triangle")
                    || s.contains("ipc")) {
                throw new IllegalArgumentException(
                        "Lane A Waveform refuses CHARC/Fig4b/Narma10b/sine-square-triangle/IPC "
                                + LaneAConstants.WAVEFORM_LABEL);
            }
        }
    }

    private static void requireParentsUnchanged() throws IOException {
        String occ = Files.readString(OCC_STANDING, StandardCharsets.UTF_8);
        String car = Files.readString(CARRIER_STANDING, StandardCharsets.UTF_8);
        String aud = Files.readString(AUDIT_STANDING, StandardCharsets.UTF_8);
        String nar = Files.readString(NARMA_STANDING, StandardCharsets.UTF_8);
        String mg = Files.readString(MG_STANDING, StandardCharsets.UTF_8);
        String lor = Files.readString(LORENZ_STANDING, StandardCharsets.UTF_8);
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
        if (!lor.contains("**Status: PASS**") || !lor.contains("LANE_A_LORENZ")) {
            throw new IllegalStateException("Lorenz standing must remain system PASS");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("Waveform PROTOCOL is not frozen_before_traces");
        }
        if (!text.contains("PREHASH_GATES_PASS")) {
            throw new IllegalStateException("Waveform pre-hash gates not recorded");
        }
        if (!json.contains("\"waveform2c_auc_copy\": false") || !json.contains("\"sine_square_triangle\": false")) {
            throw new IllegalStateException("Waveform must not copy 0.835 or use sine/square/triangle");
        }
        if (!json.contains("\"reuse_LANE_A_NARMA10_SILENT\": false")) {
            throw new IllegalStateException("Waveform must not reuse NARMA 200-window silent");
        }
        if (!json.contains(LaneAConstants.WAVEFORM_U_SHA256) || !json.contains("\"num_windows\": 400")) {
            throw new IllegalStateException("PROTOCOL u sha256 / 400 windows mismatch");
        }
        if (!json.contains("\"u_only_before_java\": true") || !json.contains("\"u_one_per_line\": true")) {
            throw new IllegalStateException("Waveform must freeze u-only before Java and one u per line");
        }
        if (!json.contains("\"warmup_18000_identity\": false") || !json.contains("\"J_max\": 128000000.0")) {
            throw new IllegalStateException("Waveform must copy occupancy J_max and warmup 0");
        }
        if (!json.contains("\"occupied_mask_is_gate\": false")) {
            throw new IllegalStateException("occupied-mask must not be the system gate");
        }
    }

    private static void requireUOnlyPass() throws IOException {
        if (!Files.exists(U_ONLY)) {
            throw new IllegalStateException(
                    "u-only json missing; run check_lane_a_waveform.py --u-only first "
                            + LaneAConstants.WAVEFORM_LABEL);
        }
        String json = Files.readString(U_ONLY, StandardCharsets.UTF_8);
        if (!json.contains("\"u_only_pass\": true") || json.contains("\"void\": true")) {
            throw new IllegalStateException(
                    "u-only did not PASS; VOID / no Java " + LaneAConstants.WAVEFORM_LABEL);
        }
    }

    private static void requireW2cUntouchedHint() throws IOException {
        if (Files.exists(W2C_SCOUT)) {
            String text = Files.readString(W2C_SCOUT, StandardCharsets.UTF_8);
            if (!text.contains("0.835") && !text.contains("0.88")) {
                throw new IllegalStateException("Waveform2c scout missing AUC markers; do not rewrite it");
            }
        }
    }

    private static double[] loadFrozenU() throws Exception {
        List<String> lines = Files.readAllLines(U_FILE, StandardCharsets.UTF_8);
        if (lines.size() < LaneAConstants.WAVEFORM_WINDOWS) {
            throw new IllegalStateException(
                    "u file looks like a one-line comma dump; write one u per line "
                            + LaneAConstants.WAVEFORM_LABEL);
        }
        List<Double> vals = new ArrayList<Double>();
        for (String line : lines) {
            String s = line.trim();
            if (s.isEmpty() || s.startsWith("#")) {
                continue;
            }
            if (s.contains(",")) {
                throw new IllegalStateException(
                        "u line contains a comma; one u per line " + LaneAConstants.WAVEFORM_LABEL);
            }
            vals.add(Double.parseDouble(s));
        }
        if (vals.size() != LaneAConstants.WAVEFORM_WINDOWS) {
            throw new IllegalStateException(
                    "u file has " + vals.size() + " values " + LaneAConstants.WAVEFORM_LABEL);
        }
        double[] u = new double[vals.size()];
        StringBuilder payload = new StringBuilder();
        for (int i = 0; i < vals.size(); i++) {
            u[i] = vals.get(i);
            if (u[i] < 0.0 || u[i] > 0.5) {
                throw new IllegalStateException("u escaped [0, 0.5] " + LaneAConstants.WAVEFORM_LABEL);
            }
            if (i > 0) {
                payload.append(',');
            }
            payload.append(String.format(Locale.US, "%.12f", u[i]));
        }
        String sha = sha256Hex(payload.toString().getBytes(StandardCharsets.US_ASCII));
        if (!LaneAConstants.WAVEFORM_U_SHA256.equals(sha)) {
            throw new IllegalStateException("u sha256 drifted: " + sha + " " + LaneAConstants.WAVEFORM_LABEL);
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
            h.append(",LANE_A_WAVEFORM");
            w.println(h);
            for (LaneAOccupiedDish.Sample s : traj.samples) {
                if (s.ahlReadout == null) {
                    throw new IllegalStateException("Waveform maps missing " + LaneAConstants.WAVEFORM_LABEL);
                }
                w.printf(Locale.US, "%.6f,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e",
                        s.t, s.window, s.uCommand, s.meanR, s.meanL,
                        s.remaining, s.commanded, s.sourceAdded, s.decayLoss, s.residual);
                writeVec(w, s.ahlReadout);
                writeVec(w, s.rReadout);
                writeVec(w, s.lReadout);
                w.printf(Locale.US, ",%s%n", LaneAConstants.WAVEFORM_LABEL);
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
                + "  \"gate\": \"" + LaneAConstants.WAVEFORM_GATE + "\",\n"
                + "  \"status_label\": \"" + LaneAConstants.WAVEFORM_LABEL + "\",\n"
                + "  \"maps\": \"" + flag + "\",\n"
                + "  \"occupancy_parent\": \"PASS\",\n"
                + "  \"carrier_parent\": \"FAIL\",\n"
                + "  \"lorenz_parent\": \"PASS\",\n"
                + "  \"occupancy\": \"" + driven.occupancyFlag + "\",\n"
                + "  \"test_mean_R\": " + driven.meanROccupancy + ",\n"
                + "  \"train_mean_R\": " + trainR + ",\n"
                + "  \"silent_maps\": \"new_400_window_U0\",\n"
                + "  \"n0_rel\": " + n0Rel + ",\n"
                + "  \"cmd_rel\": " + cmdRel + ",\n"
                + "  \"waveform2c_auc_copy\": false,\n"
                + "  \"maps_ready\": " + mapsReady + "\n"
                + "}\n";
        Files.writeString(RESULTS.resolve("lane_a_waveform_maps.json"), json, StandardCharsets.UTF_8);
    }
}
