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
 * Extra LaneA_NARMA10: system NARMA-10 vs silent, field arm reported.
 * Not a rewrite of Narma10b 0.928. Does not reopen LANE_A_CARRIER.
 */
public final class LaneANarma10Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path OCC_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_OCCUPIED_MILLIMETRE_STANDING.md");
    private static final Path CARRIER_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_CARRIER_STANDING.md");
    private static final Path AUDIT_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_CARRIER_AUDIT_STANDING.md");
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/PROTOCOL_NARMA10.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/configs/narma10_protocol.json");
    private static final Path U_FILE = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/input_u_narma200.txt");
    private static final Path RESULTS = ROOT.resolve("examples/LaneA_OccupiedMillimetre/results");
    private static final Path HYBRID_OVERALL = ROOT.resolve(
            "examples/BSimReservoirPlanNarma10b/GATE_EVIDENCE.md");

    private LaneANarma10Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        refuseForbidden(args);
        requireParentsUnchanged();
        requireProtocolFrozen();
        requireHybridDishUntouchedHint();
        LaneACli cli = LaneACli.parse(args);
        Path out = cli.outDir(RESULTS);
        double[] u = loadFrozenU();
        Files.createDirectories(out);

        System.out.println(LaneAConstants.NARMA10_LABEL + " " + LaneAConstants.NARMA10_GATE);
        System.out.println("System NARMA-10 vs silent. Field arm reported. Not Narma10b. "
                + "Not carrier PASS. Occupancy parent remains PASS.");
        cli.announce(LaneAConstants.NARMA10_LABEL);

        LaneAOccupiedDish dish = new LaneAOccupiedDish();
        System.out.println("  LANE_A_NARMA10_DRIVEN " + LaneAConstants.NARMA10_LABEL);
        long t0 = System.nanoTime();
        LaneAOccupiedDish.Trajectory driven = dish.run(LaneAOccupiedDish.Arm.DRIVEN, true, u, cli.seed);
        writeMaps(out.resolve("java_LANE_A_NARMA10_DRIVEN.csv"), driven);
        double trainR = bandMean(driven, true, LaneAConstants.NARMA10_TRAIN_LO, LaneAConstants.NARMA10_TRAIN_HI);
        double expected = expectedMass(u);
        double cmdRel = rel(driven.commanded, driven.sourceAdded, expected);
        double n0Rel = Math.abs(driven.residual) / mStar(driven);
        System.out.printf(Locale.US,
                "    occupancy=%s test_mean_R=%.6g train_mean_R=%.6g mean_L=%.6g samples=%d "
                        + "cmd=%.6g src=%.6g remain=%.6g n0_rel=%.3e cmd_rel=%.3e (%.1fs) %s%n",
                driven.occupancyFlag, driven.meanROccupancy, trainR, driven.meanLOccupancy,
                driven.samples.size(), driven.commanded, driven.sourceAdded, driven.remaining,
                n0Rel, cmdRel, elapsed(t0), LaneAConstants.NARMA10_LABEL);

        LaneAOccupiedDish.Trajectory silent = null;
        if (cli.drivenOnly) {
            System.out.println("  LANE_A_NARMA10_SILENT skipped (--driven-only); reuse seed 101 silent "
                    + LaneAConstants.NARMA10_LABEL);
        } else {
            System.out.println("  LANE_A_NARMA10_SILENT " + LaneAConstants.NARMA10_LABEL);
            t0 = System.nanoTime();
            silent = dish.run(LaneAOccupiedDish.Arm.SILENT, true, u, cli.seed);
            writeMaps(out.resolve("java_LANE_A_NARMA10_SILENT.csv"), silent);
            System.out.printf(Locale.US,
                    "    occupancy=%s mean_R=%.6g samples=%d cmd=%.6g (%.1fs) %s%n",
                    silent.occupancyFlag, silent.meanROccupancy, silent.samples.size(),
                    silent.commanded, elapsed(t0), LaneAConstants.NARMA10_LABEL);
        }

        boolean ledgerOk = cmdRel <= LaneAConstants.LEDGER_CMD_REL
                && n0Rel <= LaneAConstants.N0_MASS_REL
                && Math.abs(driven.commanded - expected) <= 1e-6 * Math.max(expected, 1.0);
        boolean silentOk = cli.drivenOnly
                || (silent.meanROccupancy < LaneAConstants.OCCUPANCY_ALIVE_MIN_R
                && Math.abs(silent.commanded) <= 1e-15);
        boolean alive = "ALIVE".equals(driven.occupancyFlag)
                && driven.meanROccupancy >= LaneAConstants.OCCUPANCY_ALIVE_MIN_R;

        if (!alive) {
            System.out.println("LANE_A_NARMA10=NOT_SCORED driven occupancy DEAD on this u. "
                    + "Do not raise J_max. Occupancy standing unchanged. Not Narma10b.");
            writeStub("NOT_SCORED", driven, silent, trainR, n0Rel, cmdRel, false);
            System.exit(2);
        }
        if (!ledgerOk || !silentOk) {
            throw new IllegalStateException(
                    "NARMA-10 ledger/silent failed n0_rel=" + n0Rel
                            + " " + LaneAConstants.NARMA10_LABEL);
        }
        writeStub("MAPS_READY", driven, silent, trainR, n0Rel, cmdRel, true);
        System.out.println("LANE_A_NARMA10 maps ready occupancy ALIVE. Ridge is the Python checker. "
                + "Not Narma10b 0.928. Occupancy remains PASS. Carrier remains FAIL.");
    }

    private static void refuseForbidden(String[] args) {
        for (String a : args) {
            String s = a.toLowerCase(Locale.ROOT);
            if (s.contains("charc") || s.contains("lorenz") || s.contains("mackey")
                    || s.contains("fig4b") || s.contains("narma10b")) {
                throw new IllegalArgumentException(
                        "Lane A NARMA-10 refuses CHARC/Lorenz/MG/Fig4b/Narma10b copy "
                                + LaneAConstants.NARMA10_LABEL);
            }
        }
    }

    private static void requireParentsUnchanged() throws IOException {
        String occ = Files.readString(OCC_STANDING, StandardCharsets.UTF_8);
        String car = Files.readString(CARRIER_STANDING, StandardCharsets.UTF_8);
        String aud = Files.readString(AUDIT_STANDING, StandardCharsets.UTF_8);
        if (!occ.contains("**Status: PASS**") || !occ.contains("LANE_A_OCCUPIED_MILLIMETRE")) {
            throw new IllegalStateException("occupancy standing must remain PASS");
        }
        if (!car.contains("**Status: FAIL**") || !car.contains("LANE_A_CARRIER")) {
            throw new IllegalStateException("carrier standing must remain FAIL vs field");
        }
        if (!aud.contains("**Status: SCOPE_NOTE**") || !aud.contains("LANE_A_CARRIER_AUDIT")) {
            throw new IllegalStateException("audit standing must remain SCOPE_NOTE");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("NARMA-10 PROTOCOL is not frozen_before_traces");
        }
        if (!json.contains("\"narma10b_rewrite\": false") || !json.contains("\"parent_carrier_rewrite\": false")) {
            throw new IllegalStateException("NARMA-10 must not rewrite Narma10b or carrier");
        }
        if (!json.contains(LaneAConstants.NARMA10_U_SHA256)) {
            throw new IllegalStateException("PROTOCOL u sha256 mismatch");
        }
        if (!json.contains("\"warmup_18000_identity\": false") || !json.contains("\"J_max\": 128000000.0")) {
            throw new IllegalStateException("NARMA-10 must copy occupancy J_max and warmup 0");
        }
        if (!json.contains("\"occupied_mask_is_gate\": false")) {
            throw new IllegalStateException("occupied-mask must not be the system gate");
        }
        if (json.contains(LaneAConstants.NARMA10B_U_SHA256_FORBIDDEN)
                && json.contains("\"u_sha256\": \"" + LaneAConstants.NARMA10B_U_SHA256_FORBIDDEN)) {
            throw new IllegalStateException("do not copy Narma10b u hash as this identity");
        }
    }

    private static void requireHybridDishUntouchedHint() throws IOException {
        if (Files.exists(HYBRID_OVERALL)) {
            String text = Files.readString(HYBRID_OVERALL, StandardCharsets.UTF_8);
            if (!text.contains("0.928") && !text.contains("Narma10b")) {
                throw new IllegalStateException(
                        "Narma10b GATE_EVIDENCE missing frozen markers; do not rewrite it");
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
            throw new IllegalStateException(
                    "u file has " + vals.size() + " values " + LaneAConstants.NARMA10_LABEL);
        }
        double[] u = new double[vals.size()];
        StringBuilder payload = new StringBuilder();
        for (int i = 0; i < vals.size(); i++) {
            u[i] = vals.get(i);
            if (u[i] < 0.0 || u[i] > 0.5) {
                throw new IllegalStateException("u escaped [0, 0.5] " + LaneAConstants.NARMA10_LABEL);
            }
            if (i > 0) {
                payload.append(',');
            }
            payload.append(String.format(Locale.US, "%.12f", u[i]));
        }
        String sha = sha256Hex(payload.toString().getBytes(StandardCharsets.US_ASCII));
        if (!LaneAConstants.NARMA10_U_SHA256.equals(sha)) {
            throw new IllegalStateException("u sha256 drifted: " + sha + " " + LaneAConstants.NARMA10_LABEL);
        }
        if (LaneAConstants.NARMA10B_U_SHA256_FORBIDDEN.equals(sha)) {
            throw new IllegalStateException("copied Narma10b u hash " + LaneAConstants.NARMA10_LABEL);
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
            h.append(",LANE_A_NARMA10");
            w.println(h);
            for (LaneAOccupiedDish.Sample s : traj.samples) {
                if (s.ahlReadout == null) {
                    throw new IllegalStateException("NARMA maps missing " + LaneAConstants.NARMA10_LABEL);
                }
                w.printf(Locale.US, "%.6f,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e",
                        s.t, s.window, s.uCommand, s.meanR, s.meanL,
                        s.remaining, s.commanded, s.sourceAdded, s.decayLoss, s.residual);
                writeVec(w, s.ahlReadout);
                writeVec(w, s.rReadout);
                writeVec(w, s.lReadout);
                w.printf(Locale.US, ",%s%n", LaneAConstants.NARMA10_LABEL);
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
            LaneAOccupiedDish.Trajectory silent,
            double trainR,
            double n0Rel,
            double cmdRel,
            boolean mapsReady) throws IOException {
        String json = "{\n"
                + "  \"gate\": \"" + LaneAConstants.NARMA10_GATE + "\",\n"
                + "  \"status_label\": \"" + LaneAConstants.NARMA10_LABEL + "\",\n"
                + "  \"maps\": \"" + flag + "\",\n"
                + "  \"occupancy_parent\": \"PASS\",\n"
                + "  \"carrier_parent\": \"FAIL\",\n"
                + "  \"audit_parent\": \"SCOPE_NOTE\",\n"
                + "  \"occupancy\": \"" + driven.occupancyFlag + "\",\n"
                + "  \"test_mean_R\": " + driven.meanROccupancy + ",\n"
                + "  \"train_mean_R\": " + trainR + ",\n"
                + "  \"silent_mean_R\": " + (silent == null ? -1.0 : silent.meanROccupancy) + ",\n"
                + "  \"n0_rel\": " + n0Rel + ",\n"
                + "  \"cmd_rel\": " + cmdRel + ",\n"
                + "  \"narma10b_rewrite\": false,\n"
                + "  \"maps_ready\": " + mapsReady + "\n"
                + "}\n";
        Files.writeString(RESULTS.resolve("lane_a_narma10_maps.json"), json, StandardCharsets.UTF_8);
    }
}
