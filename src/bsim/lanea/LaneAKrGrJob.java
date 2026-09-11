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
 * Extra LANE_A_KR_GR: new i.i.d. and constant drives, T=840.
 * Not a rewrite of METRIC_LOCK / NARMA / IEEE draft. No parameter sweep.
 */
public final class LaneAKrGrJob {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path OCC_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_OCCUPIED_MILLIMETRE_STANDING.md");
    private static final Path CARRIER_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_CARRIER_STANDING.md");
    private static final Path NARMA_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_NARMA10_STANDING.md");
    private static final Path METRIC_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_METRIC_LOCK_STANDING.md");
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/PROTOCOL_KR_GR.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/configs/kr_gr_protocol.json");
    private static final Path KR_U_FILE = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/input_u_kr_iid840.txt");
    private static final Path GR_U_FILE = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/input_u_gr_const840.txt");
    private static final Path RESULTS = ROOT.resolve("examples/LaneA_OccupiedMillimetre/results");

    private LaneAKrGrJob() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        refuseForbidden(args);
        requireParentsUnchanged();
        requireProtocolFrozen();
        LaneACli cli = LaneACli.parse(args);
        Path out = cli.outDir(RESULTS);
        double[] uKr = loadFrozenU(KR_U_FILE, LaneAConstants.KR_IID_U_SHA256, "KR");
        double[] uGr = loadFrozenU(GR_U_FILE, LaneAConstants.GR_CONST_U_SHA256, "GR");
        for (double v : uGr) {
            if (Math.abs(v - 0.25) > 1e-15) {
                throw new IllegalStateException("GR u is not constant 0.25 " + LaneAConstants.KR_GR_LABEL);
            }
        }
        Files.createDirectories(out);

        System.out.println(LaneAConstants.KR_GR_LABEL + " " + LaneAConstants.KR_GR_GATE);
        System.out.println("New i.i.d. and constant drives T=840. Rank is the Python checker. "
                + "Not a NARMA-map rewrite. No parameter sweep. Occupancy parent remains PASS.");
        cli.announce(LaneAConstants.KR_GR_LABEL);

        LaneAOccupiedDish dish = new LaneAOccupiedDish();
        System.out.println("  LANE_A_KR_GR_IID_DRIVEN " + LaneAConstants.KR_GR_LABEL);
        long t0 = System.nanoTime();
        LaneAOccupiedDish.Trajectory iid = dish.run(LaneAOccupiedDish.Arm.DRIVEN, true, uKr, cli.seed);
        writeMaps(out.resolve("java_LANE_A_KR_GR_IID_DRIVEN.csv"), iid, "IID");
        double trainR = bandMean(iid, true, 40, 599);
        double expected = expectedMass(uKr);
        double cmdRel = rel(iid.commanded, iid.sourceAdded, expected);
        double n0Rel = Math.abs(iid.residual) / mStar(iid);
        System.out.printf(Locale.US,
                "    occupancy=%s test_mean_R=%.6g train_mean_R=%.6g mean_L=%.6g samples=%d "
                        + "cmd=%.6g src=%.6g remain=%.6g n0_rel=%.3e cmd_rel=%.3e (%.1fs) %s%n",
                iid.occupancyFlag, iid.meanROccupancy, trainR, iid.meanLOccupancy,
                iid.samples.size(), iid.commanded, iid.sourceAdded, iid.remaining,
                n0Rel, cmdRel, elapsed(t0), LaneAConstants.KR_GR_LABEL);

        boolean ledgerOk = cmdRel <= LaneAConstants.LEDGER_CMD_REL
                && n0Rel <= LaneAConstants.N0_MASS_REL
                && Math.abs(iid.commanded - expected) <= 1e-6 * Math.max(expected, 1.0);
        boolean alive = "ALIVE".equals(iid.occupancyFlag)
                && iid.meanROccupancy >= LaneAConstants.OCCUPANCY_ALIVE_MIN_R;

        if (!alive) {
            System.out.println("LANE_A_KR_GR=NOT_SCORED driven occupancy DEAD on the i.i.d. u. "
                    + "Do not raise J_max. Occupancy standing unchanged. No rank.");
            writeStub("NOT_SCORED", iid, null, null, trainR, n0Rel, cmdRel, false);
            System.exit(2);
        }
        if (!ledgerOk) {
            throw new IllegalStateException(
                    "KR i.i.d. ledger failed n0_rel=" + n0Rel + " " + LaneAConstants.KR_GR_LABEL);
        }

        System.out.println("  LANE_A_KR_GR_CONST_DRIVEN " + LaneAConstants.KR_GR_LABEL);
        t0 = System.nanoTime();
        LaneAOccupiedDish.Trajectory constant = dish.run(LaneAOccupiedDish.Arm.DRIVEN, true, uGr, cli.seed);
        writeMaps(out.resolve("java_LANE_A_KR_GR_CONST_DRIVEN.csv"), constant, "CONST");
        System.out.printf(Locale.US,
                "    occupancy=%s test_mean_R=%.6g samples=%d cmd=%.6g (%.1fs) %s%n",
                constant.occupancyFlag, constant.meanROccupancy, constant.samples.size(),
                constant.commanded, elapsed(t0), LaneAConstants.KR_GR_LABEL);

        System.out.println("  LANE_A_KR_GR_SILENT " + LaneAConstants.KR_GR_LABEL);
        t0 = System.nanoTime();
        LaneAOccupiedDish.Trajectory silent = dish.run(LaneAOccupiedDish.Arm.SILENT, true, uKr, cli.seed);
        writeMaps(out.resolve("java_LANE_A_KR_GR_SILENT.csv"), silent, "SILENT");
        System.out.printf(Locale.US,
                "    occupancy=%s mean_R=%.6g samples=%d cmd=%.6g (%.1fs) %s%n",
                silent.occupancyFlag, silent.meanROccupancy, silent.samples.size(),
                silent.commanded, elapsed(t0), LaneAConstants.KR_GR_LABEL);

        boolean silentOk = silent.meanROccupancy < LaneAConstants.OCCUPANCY_ALIVE_MIN_R
                && Math.abs(silent.commanded) <= 1e-15;
        if (!silentOk) {
            throw new IllegalStateException(
                    "silent failed mean_R=" + silent.meanROccupancy + " " + LaneAConstants.KR_GR_LABEL);
        }

        writeStub("MAPS_READY", iid, constant, silent, trainR, n0Rel, cmdRel, true);
        System.out.println("LANE_A_KR_GR maps ready occupancy ALIVE. Rank and Jaeger MC are the Python checker. "
                + "Not a NARMA-map rewrite. Occupancy remains PASS. Carrier remains FAIL.");
    }

    private static void refuseForbidden(String[] args) {
        for (String a : args) {
            String s = a.toLowerCase(Locale.ROOT);
            if (s.contains("charc") || s.contains("sweep") || s.contains("narma10b")
                    || s.contains("0.928") || s.contains("gamma") || s.contains("jmax")
                    || s.contains("j_max")) {
                throw new IllegalArgumentException(
                        "LANE_A_KR_GR refuses sweep / Narma10b / 0.928 / Lane D gamma / J_max "
                                + LaneAConstants.KR_GR_LABEL);
            }
        }
    }

    private static void requireParentsUnchanged() throws IOException {
        String occ = Files.readString(OCC_STANDING, StandardCharsets.UTF_8);
        String car = Files.readString(CARRIER_STANDING, StandardCharsets.UTF_8);
        String nar = Files.readString(NARMA_STANDING, StandardCharsets.UTF_8);
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
        if (!met.contains("**Status: PASS**") || !met.contains("LANE_A_METRIC_LOCK")) {
            throw new IllegalStateException("metric-lock standing must remain PASS");
        }
        if (!met.contains("forbidden")) {
            throw new IllegalStateException("metric lock must still forbid KR/GR on NARMA maps");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("KR_GR PROTOCOL is not frozen_before_traces");
        }
        if (!json.contains("\"lane_a_charc_sweep_started\": false")) {
            throw new IllegalStateException("must not start LANE_A_CHARC_SWEEP");
        }
        if (!json.contains("\"sweep_K_tau_J\": false") || !json.contains("\"raise_J_max\": false")) {
            throw new IllegalStateException("must not sweep K/tau/J or raise J_max");
        }
        if (!json.contains("\"reuse_narma_u_as_kr_stream\": false")) {
            throw new IllegalStateException("must not reuse NARMA u as the KR stream");
        }
        if (!json.contains(LaneAConstants.KR_IID_U_SHA256)
                || !json.contains(LaneAConstants.GR_CONST_U_SHA256)) {
            throw new IllegalStateException("PROTOCOL u sha256 mismatch");
        }
        if (!json.contains("\"J_max\": 128000000.0") || !json.contains("\"warmup_s\": 0.0")) {
            throw new IllegalStateException("must copy occupancy J_max and warmup 0");
        }
        if (!json.contains("\"T_windows\": 840") || !json.contains("\"S_post_washout\": 800")) {
            throw new IllegalStateException("T/S freeze drifted");
        }
        if (json.contains("\"retune_Hill\": true")) {
            throw new IllegalStateException("Hill must not be retuned");
        }
    }

    private static double[] loadFrozenU(Path file, String expectedSha, String tag) throws Exception {
        List<Double> vals = new ArrayList<Double>();
        for (String line : Files.readAllLines(file, StandardCharsets.UTF_8)) {
            String s = line.trim();
            if (s.isEmpty() || s.startsWith("#")) {
                continue;
            }
            if (s.contains(",")) {
                throw new IllegalStateException("u line contains a comma; one u per line "
                        + LaneAConstants.KR_GR_LABEL);
            }
            vals.add(Double.parseDouble(s));
        }
        if (vals.size() != LaneAConstants.KR_GR_WINDOWS) {
            throw new IllegalStateException(
                    tag + " u file has " + vals.size() + " values " + LaneAConstants.KR_GR_LABEL);
        }
        double[] u = new double[vals.size()];
        StringBuilder payload = new StringBuilder();
        for (int i = 0; i < vals.size(); i++) {
            u[i] = vals.get(i);
            if (u[i] < 0.0 || u[i] > 0.5) {
                throw new IllegalStateException(tag + " u escaped [0, 0.5] " + LaneAConstants.KR_GR_LABEL);
            }
            if (i > 0) {
                payload.append(',');
            }
            payload.append(String.format(Locale.US, "%.12f", u[i]));
        }
        String sha = sha256Hex(payload.toString().getBytes(StandardCharsets.US_ASCII));
        if (!expectedSha.equals(sha)) {
            throw new IllegalStateException(tag + " u sha256 drifted: " + sha + " " + LaneAConstants.KR_GR_LABEL);
        }
        if (LaneAConstants.NARMA10_U_SHA256.equals(sha)
                || LaneAConstants.NARMA10B_U_SHA256_FORBIDDEN.equals(sha)) {
            throw new IllegalStateException("reused a NARMA u hash " + LaneAConstants.KR_GR_LABEL);
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

    private static void writeMaps(Path path, LaneAOccupiedDish.Trajectory traj, String tag)
            throws IOException {
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
            h.append(",LANE_A_KR_GR");
            w.println(h);
            for (LaneAOccupiedDish.Sample s : traj.samples) {
                if (s.ahlReadout == null) {
                    throw new IllegalStateException(tag + " maps missing " + LaneAConstants.KR_GR_LABEL);
                }
                w.printf(Locale.US, "%.6f,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e",
                        s.t, s.window, s.uCommand, s.meanR, s.meanL,
                        s.remaining, s.commanded, s.sourceAdded, s.decayLoss, s.residual);
                writeVec(w, s.ahlReadout);
                writeVec(w, s.rReadout);
                writeVec(w, s.lReadout);
                w.printf(Locale.US, ",%s%n", LaneAConstants.KR_GR_LABEL);
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
            LaneAOccupiedDish.Trajectory iid,
            LaneAOccupiedDish.Trajectory constant,
            LaneAOccupiedDish.Trajectory silent,
            double trainR,
            double n0Rel,
            double cmdRel,
            boolean mapsReady) throws IOException {
        String json = "{\n"
                + "  \"gate\": \"" + LaneAConstants.KR_GR_GATE + "\",\n"
                + "  \"status_label\": \"" + LaneAConstants.KR_GR_LABEL + "\",\n"
                + "  \"maps\": \"" + flag + "\",\n"
                + "  \"occupancy_parent\": \"PASS\",\n"
                + "  \"carrier_parent\": \"FAIL\",\n"
                + "  \"narma_parent\": \"PASS\",\n"
                + "  \"metric_lock_parent\": \"PASS\",\n"
                + "  \"occupancy\": \"" + iid.occupancyFlag + "\",\n"
                + "  \"test_mean_R\": " + iid.meanROccupancy + ",\n"
                + "  \"train_mean_R\": " + trainR + ",\n"
                + "  \"const_mean_R\": " + (constant == null ? -1.0 : constant.meanROccupancy) + ",\n"
                + "  \"silent_mean_R\": " + (silent == null ? -1.0 : silent.meanROccupancy) + ",\n"
                + "  \"n0_rel\": " + n0Rel + ",\n"
                + "  \"cmd_rel\": " + cmdRel + ",\n"
                + "  \"narma_map_rewrite\": false,\n"
                + "  \"maps_ready\": " + mapsReady + "\n"
                + "}\n";
        Files.writeString(RESULTS.resolve("lane_a_kr_gr_maps.json"), json, StandardCharsets.UTF_8);
    }
}
