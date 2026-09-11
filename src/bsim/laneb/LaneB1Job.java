package bsim.laneb;

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
 * Gate B1: paint occupancy and hold after field fade.
 * B1_PAINT_HOLD. Parent B0_MOTILITY_MEMORY FAIL is kept.
 * Not Fig. 4b. Not Lane A. Not C1c. No NARMA.
 */
public final class LaneB1Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/PROTOCOL_PAINT_HOLD.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/configs/protocol_b1.json");
    private static final Path B0_PROTOCOL = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/PROTOCOL.md");
    private static final Path RESULTS = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/results");

    private LaneB1Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        refuseCapacityWords(args);
        requireProtocolFrozen();
        Files.createDirectories(RESULTS);

        System.out.println(LaneB1Identity.OBJECT + " " + LaneB1Identity.GATE
                + " device=" + LaneB1Identity.DEVICE + " " + LaneB1Identity.NOT_FIG4B);
        System.out.println("Parent FAIL kept: " + LaneB1Identity.PARENT_FAIL
                + ". Not a B0 rewrite. Not Lane A. Not NARMA/CHARC.");
        System.out.printf(Locale.US,
                "t_off=%.3g s t_occ=%.3g s T_hold=%.3g s T_mix=%.3g s stripe=%.3g um J_max=%.6g%n",
                LaneB1Identity.T_OFF, LaneB1Identity.T_OCC, LaneB1Identity.T_HOLD,
                LaneB1Identity.T_MIX, LaneB1Identity.STRIPE_WAVELENGTH_UM, LaneBIdentity.J_MAX);

        BSimRandom rng = new BSimRandom(LaneB1Identity.RNG_SEED);
        LaneBDish dish = new LaneBDish(rng);
        String want = args.length > 0 ? args[0].trim() : "";
        List<LaneBDish.Result> runs = new ArrayList<LaneBDish.Result>();

        LaneBDish.Result mix = null;
        LaneBDish.Result paintMotile = null;
        LaneBDish.Result paintOff = null;
        LaneBDish.Result paintField = null;
        LaneBDish.Result paintSilent = null;

        if (wantEmptyOr(want, "B1_MIX_STRIPE")) {
            mix = dish.run(LaneBDish.Arm.MIX_STRIPE, LaneB1Identity.T_MIX);
            runs.add(mix);
            System.out.println(LaneBDish.line(mix));
        }

        double tauMix = mix == null ? Double.NaN : LaneBDish.tauMix(mix);
        boolean mixPass = gateMix(mix, tauMix);
        System.out.printf(Locale.US, "B1.1 stripe mix tau_mix=%.6g s (T_mix=%.3g) %s%n",
                tauMix, LaneB1Identity.T_MIX, mixPass ? "PASS" : "FAIL");

        boolean runPaint = wantEmptyOr(want, "B1_PAINT_MOTILE")
                || wantEmptyOr(want, "B1_PAINT_OFF")
                || wantEmptyOr(want, "B1_PAINT_FIELD")
                || wantEmptyOr(want, "B1_PAINT_SILENT");
        if (want.isEmpty() && !mixPass) {
            System.out.println("B1_PAINT_HOLD=FAIL mix. STOP. No paint interpretation. Do not start NARMA.");
            writeMixCsv(mix);
            writeSummary(false, mixPass, false, false, false, mix, null, null, null, null,
                    tauMix, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                    Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                    Double.NaN, false, "mix", runs);
            System.exit(1);
            return;
        }

        double tPaint = LaneB1Identity.paintDuration();
        if (runPaint && (want.isEmpty() || want.startsWith("B1_PAINT"))) {
            if (wantEmptyOr(want, "B1_PAINT_MOTILE")) {
                paintMotile = dish.run(LaneBDish.Arm.PAINT_MOTILE, tPaint, LaneB1Identity.T_OFF);
                runs.add(paintMotile);
                System.out.println(LaneBDish.line(paintMotile));
            }
            if (wantEmptyOr(want, "B1_PAINT_OFF")) {
                paintOff = dish.run(LaneBDish.Arm.PAINT_OFF, tPaint, LaneB1Identity.T_OFF);
                runs.add(paintOff);
                System.out.println(LaneBDish.line(paintOff));
            }
            if (wantEmptyOr(want, "B1_PAINT_FIELD")) {
                paintField = dish.run(LaneBDish.Arm.PAINT_FIELD, tPaint, LaneB1Identity.T_OFF);
                runs.add(paintField);
                System.out.println(LaneBDish.line(paintField));
            }
            if (wantEmptyOr(want, "B1_PAINT_SILENT")) {
                paintSilent = dish.run(LaneBDish.Arm.PAINT_SILENT, tPaint, LaneB1Identity.T_OFF);
                runs.add(paintSilent);
                System.out.println(LaneBDish.line(paintSilent));
            }
        }

        if ("LANE_B_MM_DISH_FLOW8".equals(want) || "B1_FLOW8".equals(want)) {
            throw new IllegalArgumentException(
                    "FLOW=8 is a named extra, not B1 identity. Not run. B1_PAINT_HOLD");
        }
        if (runs.isEmpty()) {
            throw new IllegalArgumentException("no B1 arms selected: " + want);
        }

        double kMotOcc = atTime(paintMotile, true, LaneB1Identity.T_OCC);
        double kOffOcc = atTime(paintOff, true, LaneB1Identity.T_OCC);
        double kFieldOcc = atTime(paintField, false, LaneB1Identity.T_OCC);
        double kSilentOcc = atTime(paintSilent, true, LaneB1Identity.T_OCC);
        boolean occPass = gateOccupation(kMotOcc, kOffOcc, kSilentOcc);

        double kMotHold = atTime(paintMotile, true, LaneB1Identity.T_HOLD);
        double kOffHold = atTime(paintOff, true, LaneB1Identity.T_HOLD);
        double kFieldHold = atTime(paintField, false, LaneB1Identity.T_HOLD);
        double kSilentHold = atTime(paintSilent, true, LaneB1Identity.T_HOLD);
        boolean fieldFadeOk = Double.isFinite(kFieldHold)
                && Math.abs(kFieldHold) < LaneB1Identity.KAPPA_FIELD_HOLD_MAX;
        boolean holdPass = occPass && gateHold(kMotHold, kOffHold, kFieldHold);
        boolean ledgerPass = gateLedger(paintMotile, paintOff, paintField, paintSilent);

        double pulseL = paintField != null ? paintField.pulseMeanAbsL
                : (paintMotile == null ? Double.NaN : paintMotile.pulseMeanAbsL);

        String killer = "none";
        boolean pass = mixPass && occPass && holdPass && ledgerPass;
        if (!mixPass) killer = "mix";
        else if (!occPass) killer = "occupation";
        else if (!holdPass) killer = "hold";
        else if (!ledgerPass) killer = "ledger";

        System.out.printf(Locale.US,
                "B1.2 occupation t=%.3g kappa_motile=%.6g kappa_off=%.6g kappa_silent=%.6g kappa_field=%.6g %s%n",
                LaneB1Identity.T_OCC, kMotOcc, kOffOcc, kSilentOcc, kFieldOcc,
                occPass ? "PASS" : "FAIL");
        if (occPass) {
            System.out.printf(Locale.US,
                    "B1.3 hold T_hold=%.3g kappa_motile=%.6g kappa_off=%.6g kappa_field=%.6g kappa_silent=%.6g %s%n",
                    LaneB1Identity.T_HOLD, kMotHold, kOffHold, kFieldHold, kSilentHold,
                    holdPass ? "PASS" : "FAIL");
            if (!fieldFadeOk) {
                System.out.printf(Locale.US,
                        "SCOPE_NOTE scout vs Java: |kappa_field|(T_hold)=%.6g not < 0.15. Do not retune D.%n",
                        kFieldHold);
            }
        } else {
            System.out.println("B1.3 hold NOT_SCORED (occupation FAIL). Do not dress hold as a maybe.");
        }
        System.out.printf(Locale.US, "B1.4 N0 ledger %s%n", ledgerPass ? "PASS" : "FAIL");
        System.out.printf(Locale.US,
                "REPORT mean_|L|_pulse=%.6g (not a gate; do not raise J_max)%n", pulseL);

        writeMixCsv(mix);
        writePaintCsv("b1_paint_motile.csv", paintMotile);
        writePaintCsv("b1_paint_off.csv", paintOff);
        writePaintCsv("b1_paint_field.csv", paintField);
        writePaintCsv("b1_paint_silent.csv", paintSilent);
        writeSummary(pass, mixPass, occPass, holdPass, ledgerPass, mix,
                paintMotile, paintOff, paintField, paintSilent,
                tauMix, kMotOcc, kOffOcc, kFieldOcc, kSilentOcc,
                kMotHold, kOffHold, kFieldHold, kSilentHold, pulseL,
                fieldFadeOk, killer, runs);

        if (pass) {
            System.out.println("B1_PAINT_HOLD=PASS on LANE_B_MM_DISH_FLOW0. "
                    + "B0 FAIL is not rewritten. Not NARMA. Not Fig. 4b.");
        } else {
            System.out.println("B1_PAINT_HOLD=FAIL " + killer
                    + ". Do not start NARMA/CHARC. Do not raise J_max. B0 FAIL stays.");
            System.exit(1);
        }
    }

    private static boolean wantEmptyOr(String want, String name) {
        return want.isEmpty() || want.equals(name);
    }

    private static void refuseCapacityWords(String[] argv) {
        String blob = String.join(" ", argv == null ? new String[0] : argv).toLowerCase(Locale.ROOT);
        if (blob.contains("charc") || blob.contains("ipc") || blob.contains("narma")) {
            throw new IllegalArgumentException("B1 cannot see a capacity/CHARC/NARMA target.");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("B1 protocol is not frozen_before_traces");
        }
        if (!json.contains("\"status_label\": \"B1_PAINT_HOLD\"")
                || !json.contains("\"device\": \"LANE_B_MM_DISH_FLOW0\"")
                || !json.contains("\"parent_fail\": \"B0_MOTILITY_MEMORY\"")) {
            throw new IllegalStateException("B1 protocol must freeze B1_PAINT_HOLD parent B0 FAIL");
        }
        if (!json.contains("\"narma\": false") || json.contains("\"narma\": true")) {
            throw new IllegalStateException("B1 protocol must freeze NARMA off");
        }
        if (!json.contains("\"t_off_s\": 60.0") || !json.contains("\"J_max\": 20000.0")) {
            throw new IllegalStateException("B1 protocol must freeze t_off=60 and J_max=20000");
        }
        if (!json.contains("\"T_hold_s\": 275.3")) {
            throw new IllegalStateException("B1 protocol must freeze scout T_hold=275.3 before motile traces");
        }
        if (!json.contains("\"T_hold_status\": \"SCOUT_SET_BEFORE_MOTILE_TRACES\"")) {
            throw new IllegalStateException("B1 T_hold must be scout-set before motile traces");
        }
        if (Files.exists(B0_PROTOCOL)) {
            String b0 = Files.readString(B0_PROTOCOL, StandardCharsets.UTF_8);
            if (!b0.contains("B0_MOTILITY_MEMORY") || !b0.contains("**10 s**")) {
                throw new IllegalStateException("B0 PROTOCOL.md must remain the 10 s parent; do not edit it");
            }
        }
    }

    private static double atTime(LaneBDish.Result r, boolean cell, double tq) {
        if (r == null) return Double.NaN;
        return LaneBDish.interpolate(r.time, cell ? r.kappaCell : r.kappaField, tq);
    }

    private static boolean gateMix(LaneBDish.Result mix, double tauMix) {
        return mix != null && Double.isFinite(tauMix) && tauMix > 0.0
                && tauMix < LaneB1Identity.T_MIX;
    }

    private static boolean gateOccupation(double kMot, double kOff, double kSilent) {
        if (!Double.isFinite(kMot) || !Double.isFinite(kOff) || !Double.isFinite(kSilent)) {
            return false;
        }
        return Math.abs(kMot) >= LaneB1Identity.KAPPA_MOTILE_MIN
                && Math.abs(kMot) - Math.abs(kOff) >= LaneB1Identity.KAPPA_MINUS_OFF
                && Math.abs(kSilent) < LaneB1Identity.KAPPA_SILENT_MAX;
    }

    private static boolean gateHold(double kMot, double kOff, double kField) {
        if (!Double.isFinite(kMot) || !Double.isFinite(kOff) || !Double.isFinite(kField)) {
            return false;
        }
        return Math.abs(kMot) - Math.abs(kField) >= LaneB1Identity.KAPPA_MINUS_FIELD
                && Math.abs(kMot) - Math.abs(kOff) >= LaneB1Identity.KAPPA_MINUS_OFF;
    }

    private static boolean gateLedger(LaneBDish.Result... chemistry) {
        for (LaneBDish.Result r : chemistry) {
            if (r == null) return false;
            if (!r.ledgerPass) return false;
        }
        return true;
    }

    private static void writeMixCsv(LaneBDish.Result r) throws IOException {
        if (r == null) return;
        Path p = RESULTS.resolve("b1_mix.csv");
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(p, StandardCharsets.UTF_8))) {
            w.println("t,auto_corr,kappa_cell");
            for (int i = 0; i < r.time.length; i++) {
                if (i > 0 && r.time[i] < r.time[i - 1]) break;
                w.printf(Locale.US, "%.6g,%.8g,%.8g%n", r.time[i], r.autoCorr[i], r.kappaCell[i]);
            }
        }
    }

    private static void writePaintCsv(String name, LaneBDish.Result r) throws IOException {
        if (r == null) return;
        Path p = RESULTS.resolve(name);
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(p, StandardCharsets.UTF_8))) {
            w.println("t,kappa_cell,kappa_field,mean_abs_L");
            for (int i = 0; i < r.time.length; i++) {
                if (i > 0 && r.time[i] < r.time[i - 1]) break;
                double ell = r.meanAbsL == null ? Double.NaN : r.meanAbsL[i];
                w.printf(Locale.US, "%.6g,%.8g,%.8g,%.8g%n",
                        r.time[i], r.kappaCell[i], r.kappaField[i], ell);
            }
        }
    }

    private static void writeSummary(
            boolean pass, boolean mixPass, boolean occPass, boolean holdPass, boolean ledgerPass,
            LaneBDish.Result mix, LaneBDish.Result paintMotile, LaneBDish.Result paintOff,
            LaneBDish.Result paintField, LaneBDish.Result paintSilent,
            double tauMix, double kMotOcc, double kOffOcc, double kFieldOcc, double kSilentOcc,
            double kMotHold, double kOffHold, double kFieldHold, double kSilentHold,
            double pulseL, boolean fieldFadeOk, String killer,
            List<LaneBDish.Result> runs) throws IOException {
        Path p = RESULTS.resolve("b1_summary.json");
        StringBuilder sb = new StringBuilder();
        sb.append("{\n");
        sb.append("  \"object\": \"").append(LaneB1Identity.OBJECT).append("\",\n");
        sb.append("  \"gate\": \"").append(LaneB1Identity.GATE).append("\",\n");
        sb.append("  \"device\": \"").append(LaneB1Identity.DEVICE).append("\",\n");
        sb.append("  \"parent_fail\": \"").append(LaneB1Identity.PARENT_FAIL).append("\",\n");
        sb.append("  \"NOT_FIG4B\": true,\n");
        sb.append("  \"not_lane_a\": true,\n");
        sb.append("  \"narma\": false,\n");
        sb.append("  \"b0_overall_rewrite\": false,\n");
        sb.append("  \"pass\": ").append(pass).append(",\n");
        sb.append("  \"B1.1_mix\": ").append(mixPass).append(",\n");
        sb.append("  \"B1.2_occupation\": ").append(occPass).append(",\n");
        sb.append("  \"B1.3_hold\": ").append(holdPass).append(",\n");
        sb.append("  \"B1.4_ledger\": ").append(ledgerPass).append(",\n");
        sb.append("  \"killer\": \"").append(killer).append("\",\n");
        sb.append(String.format(Locale.US, "  \"tau_mix_s\": %.8g,%n", tauMix));
        sb.append(String.format(Locale.US, "  \"t_off_s\": %.8g,%n", LaneB1Identity.T_OFF));
        sb.append(String.format(Locale.US, "  \"T_hold_s\": %.8g,%n", LaneB1Identity.T_HOLD));
        sb.append(String.format(Locale.US, "  \"kappa_motile_occ\": %.8g,%n", kMotOcc));
        sb.append(String.format(Locale.US, "  \"kappa_off_occ\": %.8g,%n", kOffOcc));
        sb.append(String.format(Locale.US, "  \"kappa_field_occ\": %.8g,%n", kFieldOcc));
        sb.append(String.format(Locale.US, "  \"kappa_silent_occ\": %.8g,%n", kSilentOcc));
        sb.append(String.format(Locale.US, "  \"kappa_motile_hold\": %.8g,%n", kMotHold));
        sb.append(String.format(Locale.US, "  \"kappa_off_hold\": %.8g,%n", kOffHold));
        sb.append(String.format(Locale.US, "  \"kappa_field_hold\": %.8g,%n", kFieldHold));
        sb.append(String.format(Locale.US, "  \"kappa_silent_hold\": %.8g,%n", kSilentHold));
        sb.append(String.format(Locale.US, "  \"mean_abs_L_pulse\": %.8g,%n", pulseL));
        sb.append("  \"java_field_fade_ok\": ").append(fieldFadeOk).append(",\n");
        sb.append("  \"J_max_raised\": false,\n");
        sb.append("  \"capacity_may_start\": ").append(pass).append(",\n");
        sb.append("  \"arms\": [");
        for (int i = 0; i < runs.size(); i++) {
            if (i > 0) sb.append(", ");
            sb.append("\"").append(runs.get(i).name).append("\"");
        }
        sb.append("]\n}\n");
        Files.writeString(p, sb.toString(), StandardCharsets.UTF_8);
    }
}
