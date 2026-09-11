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
 * Extra LANE_A_FIELD_NULL_TASK: successive-bit XOR vs linear AHL.
 * 400 windows. New silent. Does not reuse NARMA or Waveform silent.
 * U-only must PASS first. No CHARC. No second K. No extra ACs. No two-way.
 */
public final class LaneAFieldNullJob {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path OCC_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_OCCUPIED_MILLIMETRE_STANDING.md");
    private static final Path CARRIER_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_CARRIER_STANDING.md");
    private static final Path NARMA_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_NARMA10_STANDING.md");
    private static final Path WAVEFORM_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_WAVEFORM_STANDING.md");
    private static final Path KR_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_KR_GR_STANDING.md");
    private static final Path METRIC_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_METRIC_LOCK_STANDING.md");
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/PROTOCOL_FIELD_NULL.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/configs/field_null_protocol.json");
    private static final Path U_FILE = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/input_u_field_null_400.txt");
    private static final Path U_ONLY = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/results/lane_a_field_null_u_only.json");
    private static final Path RESULTS = ROOT.resolve("examples/LaneA_OccupiedMillimetre/results");
    private static final Path NARMA_SILENT = RESULTS.resolve("java_LANE_A_NARMA10_SILENT.csv");
    private static final Path WAVEFORM_SILENT = RESULTS.resolve("java_LANE_A_WAVEFORM_SILENT.csv");

    private LaneAFieldNullJob() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        refuseForbidden(args);
        requireParentsUnchanged();
        requireProtocolFrozen();
        requireUOnlyPass();
        LaneACli cli = LaneACli.parse(args);
        Path out = cli.outDir(RESULTS);
        double[] u = loadFrozenU();
        Files.createDirectories(out);

        System.out.println(LaneAConstants.FIELD_NULL_LABEL + " " + LaneAConstants.FIELD_NULL_GATE);
        System.out.println("Successive-bit XOR field-null. 400 windows. New silent. "
                + "Not CHARC. Not extra ACs. Not two-way. Occupancy parent remains PASS.");
        cli.announce(LaneAConstants.FIELD_NULL_LABEL);

        LaneAOccupiedDish dish = new LaneAOccupiedDish();
        System.out.println("  LANE_A_FIELD_NULL_DRIVEN " + LaneAConstants.FIELD_NULL_LABEL);
        long t0 = System.nanoTime();
        LaneAOccupiedDish.Trajectory driven = dish.run(LaneAOccupiedDish.Arm.DRIVEN, true, u, cli.seed);
        writeMaps(out.resolve("java_LANE_A_FIELD_NULL_DRIVEN.csv"), driven);
        double trainR = bandMean(driven, true, 64, 271);
        double expected = expectedMass(u);
        double cmdRel = rel(driven.commanded, driven.sourceAdded, expected);
        double n0Rel = Math.abs(driven.residual) / mStar(driven);
        System.out.printf(Locale.US,
                "    occupancy=%s test_mean_R=%.6g train_mean_R=%.6g mean_L=%.6g samples=%d "
                        + "cmd=%.6g src=%.6g remain=%.6g n0_rel=%.3e cmd_rel=%.3e (%.1fs) %s%n",
                driven.occupancyFlag, driven.meanROccupancy, trainR, driven.meanLOccupancy,
                driven.samples.size(), driven.commanded, driven.sourceAdded, driven.remaining,
                n0Rel, cmdRel, elapsed(t0), LaneAConstants.FIELD_NULL_LABEL);

        boolean ledgerOk = cmdRel <= LaneAConstants.LEDGER_CMD_REL
                && n0Rel <= LaneAConstants.N0_MASS_REL
                && Math.abs(driven.commanded - expected) <= 1e-6 * Math.max(expected, 1.0);
        boolean alive = "ALIVE".equals(driven.occupancyFlag)
                && driven.meanROccupancy >= LaneAConstants.OCCUPANCY_ALIVE_MIN_R;

        if (!alive) {
            System.out.println("LANE_A_FIELD_NULL_TASK=NOT_SCORED driven occupancy DEAD on this u. "
                    + "Do not raise J_max. Occupancy standing unchanged.");
            writeStub("NOT_SCORED", driven, trainR, n0Rel, cmdRel, false);
            System.exit(2);
        }
        if (!ledgerOk) {
            throw new IllegalStateException(
                    "field-null ledger failed n0_rel=" + n0Rel + " " + LaneAConstants.FIELD_NULL_LABEL);
        }

        System.out.println("  LANE_A_FIELD_NULL_SILENT 400-window U=0 empty dish "
                + LaneAConstants.FIELD_NULL_LABEL);
        long t1 = System.nanoTime();
        double[] silentLen = new double[LaneAConstants.FIELD_NULL_WINDOWS];
        LaneAOccupiedDish.Trajectory silent = dish.run(
                LaneAOccupiedDish.Arm.SILENT, true, silentLen, cli.seed);
        Path silentPath = out.resolve("java_LANE_A_FIELD_NULL_SILENT.csv");
        if (silentPath.toAbsolutePath().equals(NARMA_SILENT.toAbsolutePath())
                || silentPath.toAbsolutePath().equals(WAVEFORM_SILENT.toAbsolutePath())) {
            throw new IllegalStateException("must not reuse NARMA or Waveform silent");
        }
        writeMaps(silentPath, silent);
        System.out.printf(Locale.US, "    silent occupancy=%s test_mean_R=%.6g samples=%d (%.1fs) %s%n",
                silent.occupancyFlag, silent.meanROccupancy, silent.samples.size(),
                elapsed(t1), LaneAConstants.FIELD_NULL_LABEL);
        if (silent.samples.size() != LaneAConstants.FIELD_NULL_WINDOWS * LaneAConstants.SAMPLES_PER_WINDOW) {
            throw new IllegalStateException("silent maps are not 400 windows " + LaneAConstants.FIELD_NULL_LABEL);
        }

        writeStub("MAPS_READY", driven, trainR, n0Rel, cmdRel, true);
        System.out.println("LANE_A_FIELD_NULL_TASK maps ready occupancy ALIVE. Ridge is the Python checker. "
                + "Occupancy remains PASS. Carrier remains FAIL. NARMA remains 0.895. "
                + "KR/GR remains SCORED not CHARC.");
    }

    private static void refuseForbidden(String[] args) {
        for (String a : args) {
            String s = a.toLowerCase(Locale.ROOT);
            if (s.contains("charc") || s.contains("fig4b") || s.contains("narma10b")
                    || s.contains("0.928") || s.contains("gamma") || s.contains("two-way")
                    || s.contains("twoway") || s.contains("jmax") || s.contains("j_max")) {
                throw new IllegalArgumentException(
                        "LANE_A_FIELD_NULL_TASK refuses CHARC / Fig4b / Narma10b / 0.928 / "
                                + "Lane D gamma / two-way / J_max " + LaneAConstants.FIELD_NULL_LABEL);
            }
        }
    }

    private static void requireParentsUnchanged() throws IOException {
        String occ = Files.readString(OCC_STANDING, StandardCharsets.UTF_8);
        String car = Files.readString(CARRIER_STANDING, StandardCharsets.UTF_8);
        String nar = Files.readString(NARMA_STANDING, StandardCharsets.UTF_8);
        String wav = Files.readString(WAVEFORM_STANDING, StandardCharsets.UTF_8);
        String kr = Files.readString(KR_STANDING, StandardCharsets.UTF_8);
        String met = Files.readString(METRIC_STANDING, StandardCharsets.UTF_8);
        if (!occ.contains("**Status: PASS**") || !occ.contains("LANE_A_OCCUPIED_MILLIMETRE")) {
            throw new IllegalStateException("occupancy standing must remain PASS");
        }
        if (!car.contains("**Status: FAIL**") || !car.contains("LANE_A_CARRIER")) {
            throw new IllegalStateException("carrier standing must remain FAIL vs field");
        }
        if (!nar.contains("**Status: PASS**") || !nar.contains("LANE_A_NARMA10")) {
            throw new IllegalStateException("NARMA standing must remain system PASS");
        }
        if (!wav.contains("**Status: PASS**") || !wav.contains("LANE_A_WAVEFORM")) {
            throw new IllegalStateException("waveform standing must remain system PASS");
        }
        if (!kr.contains("**Status: SCORED**") || !kr.contains("LANE_A_KR_GR")) {
            throw new IllegalStateException("KR/GR standing must remain SCORED");
        }
        if (!kr.contains("NOT CHARC")) {
            throw new IllegalStateException("KR/GR must remain not CHARC");
        }
        if (!met.contains("**Status: PASS**") || !met.contains("LANE_A_METRIC_LOCK")) {
            throw new IllegalStateException("metric-lock standing must remain PASS");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("FIELD_NULL PROTOCOL is not frozen_before_traces");
        }
        if (!json.contains("\"lane_a_charc_sweep_started\": false")) {
            throw new IllegalStateException("must not start LANE_A_CHARC_SWEEP");
        }
        if (!json.contains("\"extra_ACs\": false") || !json.contains("\"two_way\": false")) {
            throw new IllegalStateException("must not add ACs or two-way");
        }
        if (!json.contains("\"second_Hill_K\": false") || !json.contains("\"raise_J_max\": false")) {
            throw new IllegalStateException("must not add a second K or raise J_max");
        }
        if (!json.contains(LaneAConstants.FIELD_NULL_U_SHA256) || !json.contains("\"num_windows\": 400")) {
            throw new IllegalStateException("PROTOCOL u sha256 / 400 windows mismatch");
        }
        if (!json.contains("\"reuse_LANE_A_NARMA10_SILENT\": false")
                || !json.contains("\"reuse_LANE_A_WAVEFORM_SILENT\": false")) {
            throw new IllegalStateException("must not reuse NARMA or Waveform silent");
        }
        if (!json.contains("\"J_max\": 128000000.0") || !json.contains("\"warmup_s\": 0.0")) {
            throw new IllegalStateException("must copy occupancy J_max and warmup 0");
        }
        if (json.contains("\"retune_Hill\": true")) {
            throw new IllegalStateException("Hill must not be retuned");
        }
    }

    private static void requireUOnlyPass() throws IOException {
        if (!Files.exists(U_ONLY)) {
            throw new IllegalStateException(
                    "u-only json missing; run check_lane_a_field_null.py --u-only first "
                            + LaneAConstants.FIELD_NULL_LABEL);
        }
        String json = Files.readString(U_ONLY, StandardCharsets.UTF_8);
        if (!json.contains("\"verdict\": \"U_ONLY_PASS\"") || json.contains("\"verdict\": \"VOID\"")) {
            throw new IllegalStateException(
                    "u-only did not PASS; VOID / no Java " + LaneAConstants.FIELD_NULL_LABEL);
        }
    }

    private static double[] loadFrozenU() throws Exception {
        List<Double> vals = new ArrayList<Double>();
        for (String line : Files.readAllLines(U_FILE, StandardCharsets.UTF_8)) {
            String s = line.trim();
            if (s.isEmpty() || s.startsWith("#")) {
                continue;
            }
            if (s.contains(",")) {
                throw new IllegalStateException(
                        "u line contains a comma; one u per line " + LaneAConstants.FIELD_NULL_LABEL);
            }
            vals.add(Double.parseDouble(s));
        }
        if (vals.size() != LaneAConstants.FIELD_NULL_WINDOWS) {
            throw new IllegalStateException(
                    "u file has " + vals.size() + " values " + LaneAConstants.FIELD_NULL_LABEL);
        }
        double[] u = new double[vals.size()];
        StringBuilder payload = new StringBuilder();
        for (int i = 0; i < vals.size(); i++) {
            u[i] = vals.get(i);
            if (u[i] < 0.0 || u[i] > 0.5) {
                throw new IllegalStateException("u escaped [0, 0.5] " + LaneAConstants.FIELD_NULL_LABEL);
            }
            if (Math.abs(u[i]) > 1e-15 && Math.abs(u[i] - 0.5) > 1e-15) {
                throw new IllegalStateException("u is not a 0 / 0.5 bit command " + LaneAConstants.FIELD_NULL_LABEL);
            }
            if (i > 0) {
                payload.append(',');
            }
            payload.append(String.format(Locale.US, "%.12f", u[i]));
        }
        String sha = sha256Hex(payload.toString().getBytes(StandardCharsets.US_ASCII));
        if (!LaneAConstants.FIELD_NULL_U_SHA256.equals(sha)) {
            throw new IllegalStateException("u sha256 drifted: " + sha + " " + LaneAConstants.FIELD_NULL_LABEL);
        }
        if (LaneAConstants.NARMA10_U_SHA256.equals(sha)
                || LaneAConstants.NARMA10B_U_SHA256_FORBIDDEN.equals(sha)
                || LaneAConstants.WAVEFORM_U_SHA256.equals(sha)
                || LaneAConstants.KR_IID_U_SHA256.equals(sha)) {
            throw new IllegalStateException("reused a forbidden u hash " + LaneAConstants.FIELD_NULL_LABEL);
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
            h.append(",LANE_A_FIELD_NULL_TASK");
            w.println(h);
            for (LaneAOccupiedDish.Sample s : traj.samples) {
                if (s.ahlReadout == null) {
                    throw new IllegalStateException("field-null maps missing " + LaneAConstants.FIELD_NULL_LABEL);
                }
                w.printf(Locale.US, "%.6f,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e",
                        s.t, s.window, s.uCommand, s.meanR, s.meanL,
                        s.remaining, s.commanded, s.sourceAdded, s.decayLoss, s.residual);
                writeVec(w, s.ahlReadout);
                writeVec(w, s.rReadout);
                writeVec(w, s.lReadout);
                w.printf(Locale.US, ",%s%n", LaneAConstants.FIELD_NULL_LABEL);
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
                + "  \"gate\": \"" + LaneAConstants.FIELD_NULL_GATE + "\",\n"
                + "  \"status_label\": \"" + LaneAConstants.FIELD_NULL_LABEL + "\",\n"
                + "  \"maps\": \"" + flag + "\",\n"
                + "  \"occupancy_parent\": \"PASS\",\n"
                + "  \"carrier_parent\": \"FAIL\",\n"
                + "  \"narma_parent\": \"PASS\",\n"
                + "  \"kr_gr_parent\": \"SCORED\",\n"
                + "  \"occupancy\": \"" + driven.occupancyFlag + "\",\n"
                + "  \"test_mean_R\": " + driven.meanROccupancy + ",\n"
                + "  \"train_mean_R\": " + trainR + ",\n"
                + "  \"silent_maps\": \"new_400_window_U0\",\n"
                + "  \"n0_rel\": " + n0Rel + ",\n"
                + "  \"cmd_rel\": " + cmdRel + ",\n"
                + "  \"charc\": false,\n"
                + "  \"maps_ready\": " + mapsReady + "\n"
                + "}\n";
        Files.writeString(RESULTS.resolve("lane_a_field_null_maps.json"), json, StandardCharsets.UTF_8);
    }
}
