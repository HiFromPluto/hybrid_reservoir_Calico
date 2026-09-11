package bsim.lanea;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;

/**
 * Extra LaneA_CARRIER: reconstruct u from field vs R vs L.
 * Living-layer vs plume. Not NARMA. Not occupancy rewrite.
 */
public final class LaneACarrierJob {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path OCC_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_A_OCCUPIED_MILLIMETRE_STANDING.md");
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/PROTOCOL_CARRIER.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/configs/carrier_protocol.json");
    private static final Path RESULTS = ROOT.resolve("examples/LaneA_OccupiedMillimetre/results");

    private LaneACarrierJob() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        refuseNarma(args);
        requireOccupancyStillPass();
        requireProtocolFrozen();
        LaneACli cli = LaneACli.parse(args);
        Path out = cli.outDir(RESULTS);
        Files.createDirectories(out);

        System.out.println(LaneAConstants.CARRIER_LABEL + " " + LaneAConstants.CARRIER_GATE);
        System.out.println("Reconstruction of frozen u. Not forecasting. Not NARMA. Not occupancy rewrite.");
        System.out.println("Parent LANE_A_OCCUPIED_MILLIMETRE remains PASS. Not paper 1. Not C1c. Not Lane B.");
        cli.announce(LaneAConstants.CARRIER_LABEL);

        LaneAOccupiedDish dish = new LaneAOccupiedDish();
        System.out.println("  LANE_A_CARRIER_DRIVEN " + LaneAConstants.CARRIER_LABEL);
        LaneAOccupiedDish.Trajectory driven = dish.run(LaneAOccupiedDish.Arm.DRIVEN, true, null, cli.seed);
        writeMaps(out.resolve("java_LANE_A_CARRIER_DRIVEN.csv"), driven);
        System.out.printf(Locale.US,
                "    occupancy=%s mean_R=%.6g mean_L=%.6g samples=%d %s%n",
                driven.occupancyFlag, driven.meanROccupancy, driven.meanLOccupancy,
                driven.samples.size(), LaneAConstants.CARRIER_LABEL);

        if (cli.drivenOnly) {
            System.out.println("  LANE_A_CARRIER_SILENT skipped (--driven-only); reuse seed 101 silent "
                    + LaneAConstants.CARRIER_LABEL);
        } else {
            System.out.println("  LANE_A_CARRIER_SILENT " + LaneAConstants.CARRIER_LABEL);
            LaneAOccupiedDish.Trajectory silent = dish.run(LaneAOccupiedDish.Arm.SILENT, true, null, cli.seed);
            writeMaps(out.resolve("java_LANE_A_CARRIER_SILENT.csv"), silent);
            System.out.printf(Locale.US,
                    "    occupancy=%s mean_R=%.6g samples=%d %s%n",
                    silent.occupancyFlag, silent.meanROccupancy, silent.samples.size(),
                    LaneAConstants.CARRIER_LABEL);
        }

        if (!"ALIVE".equals(driven.occupancyFlag)
                || driven.meanROccupancy < LaneAConstants.OCCUPANCY_ALIVE_MIN_R) {
            System.out.println("LANE_A_CARRIER=NOT_SCORED driven occupancy DEAD. Do not raise J_max. "
                    + "Occupancy standing unchanged. Not NARMA.");
            writeStub(false, "NOT_SCORED");
            System.exit(2);
        }
        writeStub(true, "MAPS_READY");
        System.out.println("LANE_A_CARRIER maps ready occupancy ALIVE. Ridge is the Python checker. "
                + "Not NARMA. Occupancy remains PASS.");
    }

    private static void refuseNarma(String[] args) {
        for (String a : args) {
            String s = a.toLowerCase(Locale.ROOT);
            if (s.contains("narma") || s.contains("charc") || s.contains("lorenz")
                    || s.contains("mackey") || s.contains("fig4b")) {
                throw new IllegalArgumentException(
                        "Lane A carrier refuses NARMA/CHARC/Lorenz/MG/Fig4b " + LaneAConstants.CARRIER_LABEL);
            }
        }
    }

    private static void requireOccupancyStillPass() throws IOException {
        String text = Files.readString(OCC_STANDING, StandardCharsets.UTF_8);
        if (!text.contains("**Status: PASS**") || !text.contains("LANE_A_OCCUPIED_MILLIMETRE")) {
            throw new IllegalStateException("occupancy standing must remain PASS before carrier");
        }
        if (text.contains("living-layer PASS") && text.contains("Narma10b Overall")) {
            throw new IllegalStateException("do not convert occupancy standing into a task PASS");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("carrier PROTOCOL is not frozen_before_traces");
        }
        if (!json.contains("\"narma\": false") || !json.contains("\"forecasting\": false")) {
            throw new IllegalStateException("carrier PROTOCOL must freeze narma/forecasting false");
        }
        if (!json.contains("\"margin_delta\": 0.02") || !json.contains("\"warmup_18000_identity\": false")) {
            throw new IllegalStateException("carrier PROTOCOL must freeze delta=0.02 and no 18000 identity");
        }
        if (!json.contains("\"window_mean_u\": \"VOID\"")) {
            throw new IllegalStateException("carrier PROTOCOL must VOID window-mean u");
        }
        if (!json.contains("\"J_max\": 128000000.0")) {
            throw new IllegalStateException("carrier must copy occupancy J_max");
        }
    }

    private static void writeMaps(Path path, LaneAOccupiedDish.Trajectory traj) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            StringBuilder h = new StringBuilder("t,window,u");
            for (int i = 0; i < LaneAConstants.READOUT_BINS; i++) {
                h.append(",AHL_").append(i);
            }
            for (int i = 0; i < LaneAConstants.READOUT_BINS; i++) {
                h.append(",R_").append(i);
            }
            for (int i = 0; i < LaneAConstants.READOUT_BINS; i++) {
                h.append(",L_").append(i);
            }
            h.append(",LANE_A_CARRIER");
            w.println(h);
            for (LaneAOccupiedDish.Sample s : traj.samples) {
                if (s.ahlReadout == null) {
                    throw new IllegalStateException("carrier maps missing " + LaneAConstants.CARRIER_LABEL);
                }
                w.printf(Locale.US, "%.6f,%d,%.16e", s.t, s.window, s.uCommand);
                writeVec(w, s.ahlReadout);
                writeVec(w, s.rReadout);
                writeVec(w, s.lReadout);
                w.printf(Locale.US, ",%s%n", LaneAConstants.CARRIER_LABEL);
            }
        }
    }

    private static void writeVec(PrintWriter w, double[] v) {
        for (double x : v) {
            w.printf(Locale.US, ",%.16e", x);
        }
    }

    private static void writeStub(boolean mapsReady, String flag) throws IOException {
        String json = "{\n"
                + "  \"gate\": \"" + LaneAConstants.CARRIER_GATE + "\",\n"
                + "  \"status_label\": \"" + LaneAConstants.CARRIER_LABEL + "\",\n"
                + "  \"maps\": \"" + flag + "\",\n"
                + "  \"occupancy_parent\": \"PASS\",\n"
                + "  \"narma\": false,\n"
                + "  \"forecasting\": false,\n"
                + "  \"maps_ready\": " + mapsReady + "\n"
                + "}\n";
        Files.writeString(RESULTS.resolve("lane_a_carrier_maps.json"), json, StandardCharsets.UTF_8);
    }
}
