package bsim.lanea;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Gate LaneA_N0_OCCUPANCY. LANE_A_OCCUPIED_MILLIMETRE.
 * Not paper 1. Not HybridDish rewrite. Not NARMA. Not C1c. Not Fig. 4b. Not Lane B.
 */
public final class LaneAN0OccupancyJob {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path PROTOCOL = ROOT.resolve("examples/LaneA_OccupiedMillimetre/PROTOCOL.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneA_OccupiedMillimetre/configs/protocol.json");
    private static final Path RESULTS = ROOT.resolve("examples/LaneA_OccupiedMillimetre/results");
    private static final Path HYBRID_OVERALL = ROOT.resolve(
            "examples/BSimReservoirPlanNarma10b/GATE_EVIDENCE.md");

    private LaneAN0OccupancyJob() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        refuseNarma(args);
        requireProtocolFrozen();
        requireHybridDishUntouchedHint();
        Files.createDirectories(RESULTS);

        System.out.println(LaneAConstants.LABEL + " " + LaneAConstants.GATE);
        System.out.println("Not paper 1. Not C1c. Not Fig. 4b. Not Lane B. Not NARMA. Not calibrated AC.");
        System.out.printf(
                Locale.US,
                "%s dish=1000x500x10 CENTER=(%.1f,%.1f,%.1f) voxel=(%d,%d,%d) "
                        + "FLOW=0 J_max=%.6g K=%.1f n=%.0f tau_R=%.0f tau_L=%.0f motility=OFF%n",
                LaneAConstants.LABEL,
                LaneAConstants.AC_X, LaneAConstants.AC_Y, LaneAConstants.AC_Z,
                LaneAConstants.I_AC, LaneAConstants.J_AC, LaneAConstants.K_AC,
                LaneAConstants.J_MAX,
                LaneAConstants.HILL_K_UM, LaneAConstants.HILL_N,
                LaneAConstants.TAU_R_S, LaneAConstants.TAU_L_S);

        LaneAOccupiedDish dish = new LaneAOccupiedDish();
        List<ArmResult> arms = new ArrayList<ArmResult>();
        String want = args.length > 0 ? args[0].trim() : "";

        if (wantEmptyOr(want, "LANE_A_CLOSED_CONSERVATIVE")) {
            arms.add(runClosed(dish));
        }
        if (wantEmptyOr(want, "LANE_A_DRIVEN")) {
            arms.add(runLiving(dish, LaneAOccupiedDish.Arm.DRIVEN, "LANE_A_DRIVEN"));
        }
        if (wantEmptyOr(want, "LANE_A_SILENT")) {
            arms.add(runLiving(dish, LaneAOccupiedDish.Arm.SILENT, "LANE_A_SILENT"));
        }
        if (wantEmptyOr(want, "LANE_A_BROWNIAN")) {
            arms.add(runLiving(dish, LaneAOccupiedDish.Arm.BROWNIAN, "LANE_A_BROWNIAN"));
        }
        if (arms.isEmpty()) {
            throw new IllegalArgumentException("no Lane A arms selected: " + want);
        }

        ArmResult driven = find(arms, "LANE_A_DRIVEN");
        ArmResult silent = find(arms, "LANE_A_SILENT");
        if (driven != null && silent != null) {
            boolean silentOk = silent.meanR < LaneAConstants.OCCUPANCY_ALIVE_MIN_R
                    && driven.meanR > silent.meanR;
            if (!silentOk) {
                silent = silent.fail("silent occupied like driven or not below driven " + LaneAConstants.LABEL);
                replace(arms, silent);
            }
        }

        boolean pass = true;
        for (ArmResult a : arms) {
            pass = pass && a.pass;
        }

        writeSummary(arms, pass);
        System.out.println(pass
                ? "LANE_A_OCCUPIED_MILLIMETRE=PASS occupancy-first N0 millimetre dish. "
                + "HybridDish Overall unchanged. Paper 1 unchanged. Lane B not started. "
                + "Not Fig. 4b. Not C1c. Not NARMA."
                : "LANE_A_OCCUPIED_MILLIMETRE=FAIL. Do not raise J_max. Do not edit HybridDish "
                + "GATE_EVIDENCE. Do not start NARMA. Not Fig. 4b. Not C1c.");
        if (!pass) {
            System.exit(1);
        }
    }

    private static boolean wantEmptyOr(String want, String name) {
        return want.isEmpty() || want.equals(name);
    }

    private static void refuseNarma(String[] args) {
        for (String a : args) {
            String s = a.toLowerCase(Locale.ROOT);
            if (s.contains("narma") || s.contains("charc") || s.contains("fig4b")
                    || s.contains("c1c") || s.contains("time_adj")) {
                throw new IllegalArgumentException(
                        "Lane A occupancy refuses NARMA/CHARC/Fig4b/C1c/TIME_ADJ " + LaneAConstants.LABEL);
            }
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("Lane A PROTOCOL is not frozen_before_traces");
        }
        if (!text.contains("LANE_A_OCCUPIED_MILLIMETRE") || !json.contains("\"status_label\": \"LANE_A_OCCUPIED_MILLIMETRE\"")) {
            throw new IllegalStateException("Lane A PROTOCOL must freeze LANE_A_OCCUPIED_MILLIMETRE");
        }
        if (!json.contains("\"narma\": false") || !json.contains("\"J_max\": 128000000.0")) {
            throw new IllegalStateException("Lane A PROTOCOL must freeze narma=false and J_max=1.28e8");
        }
        if (!json.contains("\"occupancy_alive_min_R\": 0.05")) {
            throw new IllegalStateException("Lane A PROTOCOL must freeze occupancy threshold 0.05 before traces");
        }
        if (!json.contains("\"TIME_ADJ\": false") || !json.contains("\"FLOW\": 0.0")) {
            throw new IllegalStateException("Lane A PROTOCOL must freeze TIME_ADJ unused and FLOW=0");
        }
        if (json.contains("\"narma\": true") || text.toLowerCase(Locale.ROOT).contains("time_adj=60")) {
            throw new IllegalStateException("Lane A PROTOCOL forbids NARMA identity and TIME_ADJ=60");
        }
        if (Math.abs(LaneAConstants.EXPECTED_DRIVEN_MASS - 38400000000.0) > 1.0) {
            throw new IllegalStateException("driven commanded mass must stay 3.84e10 " + LaneAConstants.LABEL);
        }
    }

    private static void requireHybridDishUntouchedHint() throws IOException {
        if (Files.exists(HYBRID_OVERALL)) {
            String text = Files.readString(HYBRID_OVERALL, StandardCharsets.UTF_8);
            if (!text.contains("0.928") && !text.contains("PASS")) {
                throw new IllegalStateException(
                        "Narma10b GATE_EVIDENCE missing frozen Overall markers; do not rewrite it");
            }
        }
    }

    private static ArmResult runClosed(LaneAOccupiedDish dish) throws IOException {
        long t0 = System.nanoTime();
        System.out.println("  LANE_A_CLOSED_CONSERVATIVE " + LaneAConstants.LABEL + " decay=0 impulse");
        LaneAOccupiedDish.Trajectory traj = dish.run(LaneAOccupiedDish.Arm.CLOSED_CONSERVATIVE);
        writeCsv(RESULTS.resolve("java_LANE_A_CLOSED_CONSERVATIVE.csv"), traj);
        double cmdRel = rel(traj.commanded, traj.sourceAdded, LaneAConstants.CLOSED_DELTA_M);
        double massRel = Math.abs(traj.residual) / mStar(traj);
        boolean pass = cmdRel <= LaneAConstants.LEDGER_CMD_REL
                && massRel <= LaneAConstants.N0_MASS_REL
                && Math.abs(traj.commanded - LaneAConstants.CLOSED_DELTA_M)
                <= 1e-6 * LaneAConstants.CLOSED_DELTA_M;
        System.out.printf(
                Locale.US,
                "    cmd=%.6g src=%.6g remain=%.6g decay=%.6g residual=%.6g cmd_rel=%.3e n0_rel=%.3e %s (%.1fs)%n",
                traj.commanded, traj.sourceAdded, traj.remaining, traj.decayLoss, traj.residual,
                cmdRel, massRel, pass ? "PASS" : "FAIL", elapsed(t0));
        return new ArmResult(
                "LANE_A_CLOSED_CONSERVATIVE", traj, cmdRel, massRel, pass,
                pass ? "closed N0 band" : "closed ledger failed");
    }

    private static ArmResult runLiving(
            LaneAOccupiedDish dish, LaneAOccupiedDish.Arm arm, String name) throws IOException {
        long t0 = System.nanoTime();
        System.out.println("  " + name + " " + LaneAConstants.LABEL);
        LaneAOccupiedDish.Trajectory traj = dish.run(arm);
        writeCsv(RESULTS.resolve("java_" + name + ".csv"), traj);
        double cmdRel = rel(traj.commanded, traj.sourceAdded, LaneAConstants.EXPECTED_DRIVEN_MASS);
        double massRel = Math.abs(traj.residual) / mStar(traj);
        boolean pass;
        String note;
        if (arm == LaneAOccupiedDish.Arm.BROWNIAN) {
            pass = traj.meanROccupancy == 0.0
                    && "NO_RECEIVER".equals(traj.occupancyFlag)
                    && cmdRel <= LaneAConstants.LEDGER_CMD_REL
                    && massRel <= LaneAConstants.N0_MASS_REL;
            note = "mean_R=0 no Hill";
        } else if (arm == LaneAOccupiedDish.Arm.SILENT) {
            pass = "DEAD".equals(traj.occupancyFlag)
                    && traj.meanROccupancy < LaneAConstants.OCCUPANCY_ALIVE_MIN_R
                    && Math.abs(traj.commanded) <= 1e-15
                    && cmdRel <= LaneAConstants.LEDGER_CMD_REL
                    && massRel <= LaneAConstants.N0_MASS_REL;
            note = "mean_R=" + traj.meanROccupancy + " " + traj.occupancyFlag;
        } else {
            pass = "ALIVE".equals(traj.occupancyFlag)
                    && traj.meanROccupancy >= LaneAConstants.OCCUPANCY_ALIVE_MIN_R
                    && cmdRel <= LaneAConstants.LEDGER_CMD_REL
                    && massRel <= LaneAConstants.N0_MASS_REL
                    && Math.abs(traj.commanded - LaneAConstants.EXPECTED_DRIVEN_MASS)
                    <= 1e-6 * LaneAConstants.EXPECTED_DRIVEN_MASS;
            note = "mean_R=" + traj.meanROccupancy + " " + traj.occupancyFlag;
        }
        System.out.printf(
                Locale.US,
                "    occupancy=%s mean_R=%.6g mean_L=%.6g n_occ=%d cmd=%.6g src=%.6g "
                        + "remain=%.6g decay=%.6g residual=%.6g cmd_rel=%.3e n0_rel=%.3e %s (%.1fs)%n",
                traj.occupancyFlag, traj.meanROccupancy, traj.meanLOccupancy, traj.occupancySamples,
                traj.commanded, traj.sourceAdded, traj.remaining, traj.decayLoss, traj.residual,
                cmdRel, massRel, pass ? "PASS" : "FAIL", elapsed(t0));
        return new ArmResult(name, traj, cmdRel, massRel, pass, note);
    }

    private static void writeCsv(Path path, LaneAOccupiedDish.Trajectory traj) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.println("t,window,n_bodies,mean_R,mean_L,frac_R_gt_0_5,ahl_uM_mean,ahl_uM_center,"
                    + "remaining,commanded,source_added,decay_loss,residual,LANE_A_OCCUPIED_MILLIMETRE");
            for (LaneAOccupiedDish.Sample s : traj.samples) {
                w.printf(Locale.US,
                        "%.6f,%d,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%s%n",
                        s.t, s.window, s.nBodies, s.meanR, s.meanL, s.fracRgt05,
                        s.ahlUmMean, s.ahlUmCenter, s.remaining, s.commanded,
                        s.sourceAdded, s.decayLoss, s.residual, LaneAConstants.LABEL);
            }
        }
    }

    private static void writeSummary(List<ArmResult> arms, boolean pass) throws IOException {
        StringBuilder json = new StringBuilder();
        json.append("{\n");
        json.append("  \"gate\": \"").append(LaneAConstants.GATE).append("\",\n");
        json.append("  \"status_label\": \"").append(LaneAConstants.LABEL).append("\",\n");
        json.append("  \"LANE_A_OCCUPIED_MILLIMETRE\": \"").append(pass ? "PASS" : "FAIL").append("\",\n");
        json.append("  \"not_paper_1\": true,\n");
        json.append("  \"not_c1c\": true,\n");
        json.append("  \"not_fig4b\": true,\n");
        json.append("  \"not_lane_b\": true,\n");
        json.append("  \"narma\": false,\n");
        json.append("  \"TIME_ADJ\": false,\n");
        json.append("  \"J_max\": ").append(LaneAConstants.J_MAX).append(",\n");
        json.append("  \"occupancy_alive_min_R\": ").append(LaneAConstants.OCCUPANCY_ALIVE_MIN_R).append(",\n");
        json.append("  \"hybriddish_overall_rewrite\": false,\n");
        json.append("  \"paper1_rewrite\": false,\n");
        json.append("  \"lane_b_started\": false,\n");
        json.append("  \"arms\": [\n");
        for (int i = 0; i < arms.size(); i++) {
            ArmResult a = arms.get(i);
            json.append("    {\"name\":\"").append(a.name).append("\",");
            json.append("\"occupancy\":\"").append(a.traj.occupancyFlag).append("\",");
            json.append("\"mean_R\":").append(num(a.meanR)).append(",");
            json.append("\"mean_L\":").append(num(a.traj.meanLOccupancy)).append(",");
            json.append("\"cmd_mass\":").append(num(a.traj.commanded)).append(",");
            json.append("\"src_mass\":").append(num(a.traj.sourceAdded)).append(",");
            json.append("\"remaining\":").append(num(a.traj.remaining)).append(",");
            json.append("\"decay\":").append(num(a.traj.decayLoss)).append(",");
            json.append("\"residual\":").append(num(a.traj.residual)).append(",");
            json.append("\"cmd_rel\":").append(num(a.cmdRel)).append(",");
            json.append("\"n0_rel\":").append(num(a.massRel)).append(",");
            json.append("\"pass\":").append(a.pass).append(",");
            json.append("\"LANE_A_OCCUPIED_MILLIMETRE\": true}");
            json.append(i + 1 < arms.size() ? ",\n" : "\n");
        }
        json.append("  ]\n}\n");
        Files.writeString(RESULTS.resolve("lane_a_summary.json"), json.toString(), StandardCharsets.UTF_8);
    }

    private static ArmResult find(List<ArmResult> arms, String name) {
        for (ArmResult a : arms) {
            if (name.equals(a.name)) {
                return a;
            }
        }
        return null;
    }

    private static void replace(List<ArmResult> arms, ArmResult next) {
        for (int i = 0; i < arms.size(); i++) {
            if (arms.get(i).name.equals(next.name)) {
                arms.set(i, next);
                return;
            }
        }
    }

    private static double rel(double a, double b, double floor) {
        return Math.abs(a - b) / Math.max(Math.max(Math.abs(a), Math.abs(b)), Math.max(floor, 1e-15));
    }

    private static double mStar(LaneAOccupiedDish.Trajectory traj) {
        return Math.max(Math.max(Math.abs(traj.commanded), Math.abs(traj.sourceAdded)),
                Math.max(Math.abs(traj.remaining), 1e-15));
    }

    private static String num(double x) {
        return Double.isNaN(x) ? "null" : Double.toString(x);
    }

    private static double elapsed(long t0) {
        return (System.nanoTime() - t0) / 1e9;
    }

    private static final class ArmResult {
        final String name;
        final LaneAOccupiedDish.Trajectory traj;
        final double cmdRel;
        final double massRel;
        final double meanR;
        final boolean pass;
        final String note;

        ArmResult(String name, LaneAOccupiedDish.Trajectory traj, double cmdRel, double massRel,
                  boolean pass, String note) {
            this.name = name;
            this.traj = traj;
            this.cmdRel = cmdRel;
            this.massRel = massRel;
            this.meanR = traj.meanROccupancy;
            this.pass = pass;
            this.note = note;
        }

        ArmResult fail(String note) {
            return new ArmResult(name, traj, cmdRel, massRel, false, note);
        }
    }
}
