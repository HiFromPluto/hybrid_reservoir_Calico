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
 * Gate B2: stripe occupation at the 200 µm mix length.
 * B2_STRIPE_OCCUPY. Parent B0 and B1 FAILs are kept.
 * Not Fig. 4b. Not Lane A. Not C1c. No NARMA. Not a B1 J_max hunt.
 */
public final class LaneB2Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/PROTOCOL_STRIPE_OCCUPY.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/configs/protocol_b2.json");
    private static final Path B0_PROTOCOL = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/PROTOCOL.md");
    private static final Path B1_PROTOCOL = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/PROTOCOL_PAINT_HOLD.md");
    private static final Path RESULTS = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/results");

    private LaneB2Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        refuseCapacityWords(args);
        requireProtocolFrozen();
        Files.createDirectories(RESULTS);

        System.out.println(LaneB2Identity.OBJECT + " " + LaneB2Identity.GATE
                + " device=" + LaneB2Identity.DEVICE + " " + LaneB2Identity.NOT_FIG4B);
        System.out.println("Parent FAILs kept: " + LaneB2Identity.PARENT_FAIL_B0
                + ", " + LaneB2Identity.PARENT_FAIL_B1
                + ". Not a B0/B1 rewrite. Not half-dish kappa. Not NARMA/CHARC.");
        System.out.printf(Locale.US,
                "L_slab=%.6g t_occ=%.3g s T_hold=%.3g s lambda=%.3g um ell=%.3g um%n",
                LaneB2Identity.L_SLAB, LaneB2Identity.T_OCC, LaneB2Identity.T_HOLD,
                LaneB2Identity.LAMBDA_UM, LaneB2Identity.ELL_RUN_UM);

        BSimRandom rng = new BSimRandom(LaneB2Identity.RNG_SEED);
        LaneBDish dish = new LaneBDish(rng);
        String want = args.length > 0 ? args[0].trim() : "";
        List<LaneBDish.Result> runs = new ArrayList<LaneBDish.Result>();

        LaneBDish.Result stripeMotile = null;
        LaneBDish.Result stripeOff = null;
        LaneBDish.Result stripeField = null;
        LaneBDish.Result stripeSilent = null;

        if ("LANE_B_MM_DISH_FLOW8".equals(want) || "B2_FLOW8".equals(want)) {
            throw new IllegalArgumentException(
                    "FLOW=8 is a named extra, not B2 identity. Not run. B2_STRIPE_OCCUPY");
        }

        double tPaint = LaneB2Identity.paintDuration();
        if (wantEmptyOr(want, "B2_STRIPE_MOTILE")) {
            stripeMotile = dish.run(LaneBDish.Arm.STRIPE_MOTILE, tPaint, 0.0);
            runs.add(stripeMotile);
            System.out.println(LaneBDish.line(stripeMotile));
        }
        if (wantEmptyOr(want, "B2_STRIPE_OFF")) {
            stripeOff = dish.run(LaneBDish.Arm.STRIPE_OFF, tPaint, 0.0);
            runs.add(stripeOff);
            System.out.println(LaneBDish.line(stripeOff));
        }
        if (wantEmptyOr(want, "B2_STRIPE_FIELD")) {
            stripeField = dish.run(LaneBDish.Arm.STRIPE_FIELD, tPaint, 0.0);
            runs.add(stripeField);
            System.out.println(LaneBDish.line(stripeField));
        }
        if (wantEmptyOr(want, "B2_STRIPE_SILENT")) {
            stripeSilent = dish.run(LaneBDish.Arm.STRIPE_SILENT, tPaint, 0.0);
            runs.add(stripeSilent);
            System.out.println(LaneBDish.line(stripeSilent));
        }
        if (runs.isEmpty()) {
            throw new IllegalArgumentException("no B2 arms selected: " + want);
        }

        double kMotOcc = atTime(stripeMotile, true, LaneB2Identity.T_OCC);
        double kOffOcc = atTime(stripeOff, true, LaneB2Identity.T_OCC);
        double kFieldOcc = atTime(stripeField, false, LaneB2Identity.T_OCC);
        double kSilentOcc = atTime(stripeSilent, true, LaneB2Identity.T_OCC);
        boolean occPass = gateOccupation(kMotOcc, kOffOcc, kSilentOcc);

        double kMotHold = atTime(stripeMotile, true, LaneB2Identity.T_HOLD);
        double kOffHold = atTime(stripeOff, true, LaneB2Identity.T_HOLD);
        double kFieldHold = atTime(stripeField, false, LaneB2Identity.T_HOLD);
        double kSilentHold = atTime(stripeSilent, true, LaneB2Identity.T_HOLD);
        boolean fieldFadeOk = Double.isFinite(kFieldHold)
                && Math.abs(kFieldHold) < LaneB2Identity.KAPPA_FIELD_HOLD_MAX;
        boolean holdPass = occPass && gateHold(kMotHold, kOffHold, kFieldHold);
        boolean ledgerPass = gateLedger(stripeMotile, stripeOff, stripeField, stripeSilent);

        double dLOcc = deltaLAt(stripeField, LaneB2Identity.T_OCC);
        if (!Double.isFinite(dLOcc)) dLOcc = deltaLAt(stripeMotile, LaneB2Identity.T_OCC);

        String killer = "none";
        boolean pass = occPass && holdPass && ledgerPass;
        if (!occPass) killer = "occupation";
        else if (!holdPass) killer = "hold";
        else if (!ledgerPass) killer = "ledger";

        System.out.printf(Locale.US,
                "B2.1 occupation t=%.3g kappa_motile=%.6g kappa_off=%.6g kappa_silent=%.6g kappa_field=%.6g %s%n",
                LaneB2Identity.T_OCC, kMotOcc, kOffOcc, kSilentOcc, kFieldOcc,
                occPass ? "PASS" : "FAIL");
        if (occPass) {
            System.out.printf(Locale.US,
                    "B2.2 hold T_hold=%.3g kappa_motile=%.6g kappa_off=%.6g kappa_field=%.6g kappa_silent=%.6g %s%n",
                    LaneB2Identity.T_HOLD, kMotHold, kOffHold, kFieldHold, kSilentHold,
                    holdPass ? "PASS" : "FAIL");
            if (!fieldFadeOk) {
                System.out.printf(Locale.US,
                        "SCOPE_NOTE scout vs Java: |kappa_lambda_field|(T_hold)=%.6g not < 0.15. Do not retune D.%n",
                        kFieldHold);
            }
        } else {
            System.out.println("B2.2 hold NOT_SCORED (occupation FAIL). Do not dress hold as a maybe.");
        }
        System.out.printf(Locale.US, "B2.3 N0 ledger %s%n", ledgerPass ? "PASS" : "FAIL");
        System.out.printf(Locale.US,
                "REPORT DeltaL(ell) at t_occ=%.6g (Java field; not a gate)%n", dLOcc);

        writePaintCsv("b2_stripe_motile.csv", stripeMotile);
        writePaintCsv("b2_stripe_off.csv", stripeOff);
        writePaintCsv("b2_stripe_field.csv", stripeField);
        writePaintCsv("b2_stripe_silent.csv", stripeSilent);
        writeSummary(pass, occPass, holdPass, ledgerPass,
                stripeMotile, stripeOff, stripeField, stripeSilent,
                kMotOcc, kOffOcc, kFieldOcc, kSilentOcc,
                kMotHold, kOffHold, kFieldHold, kSilentHold, dLOcc,
                fieldFadeOk, killer, runs);

        if (pass) {
            System.out.println("B2_STRIPE_OCCUPY=PASS on LANE_B_MM_DISH_FLOW0. "
                    + "B0 FAIL and B1 FAIL occupation are not rewritten. Not NARMA. Not Fig. 4b.");
        } else {
            System.out.println("B2_STRIPE_OCCUPY=FAIL " + killer
                    + ". Do not start NARMA/CHARC. Do not raise L_slab. Do not drop D. "
                    + "B0 FAIL and B1 FAIL stay.");
            System.exit(1);
        }
    }

    private static boolean wantEmptyOr(String want, String name) {
        return want.isEmpty() || want.equals(name);
    }

    private static void refuseCapacityWords(String[] argv) {
        String blob = String.join(" ", argv == null ? new String[0] : argv).toLowerCase(Locale.ROOT);
        if (blob.contains("charc") || blob.contains("ipc") || blob.contains("narma")) {
            throw new IllegalArgumentException("B2 cannot see a capacity/CHARC/NARMA target.");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("B2 protocol is not frozen_before_traces");
        }
        if (!json.contains("\"status_label\": \"B2_STRIPE_OCCUPY\"")
                || !json.contains("\"device\": \"LANE_B_MM_DISH_FLOW0\"")
                || !json.contains("\"B0_MOTILITY_MEMORY\"")
                || !json.contains("\"B1_PAINT_HOLD\"")) {
            throw new IllegalStateException("B2 protocol must freeze B2_STRIPE_OCCUPY parent B0/B1 FAILs");
        }
        if (!json.contains("\"narma\": false") || json.contains("\"narma\": true")) {
            throw new IllegalStateException("B2 protocol must freeze NARMA off");
        }
        if (!json.contains("\"T_hold_s\": 3.1") || !json.contains("\"t_occ_s\": 1.1")) {
            throw new IllegalStateException("B2 protocol must freeze scout t_occ=1.1 and T_hold=3.1");
        }
        if (!json.contains("\"T_hold_status\": \"SCOUT_SET_BEFORE_MOTILE_TRACES\"")) {
            throw new IllegalStateException("B2 T_hold must be scout-set before motile traces");
        }
        if (!json.contains("\"D_att_um2_s\": 800.0") || json.contains("\"D_att_um2_s\": 50")) {
            throw new IllegalStateException("B2 must keep D=800");
        }
        if (Files.exists(B0_PROTOCOL)) {
            String b0 = Files.readString(B0_PROTOCOL, StandardCharsets.UTF_8);
            if (!b0.contains("B0_MOTILITY_MEMORY") || !b0.contains("**10 s**")) {
                throw new IllegalStateException("B0 PROTOCOL.md must remain the 10 s parent; do not edit it");
            }
        }
        if (Files.exists(B1_PROTOCOL)) {
            String b1 = Files.readString(B1_PROTOCOL, StandardCharsets.UTF_8);
            if (!b1.contains("B1_PAINT_HOLD") || !b1.contains("**60** s")
                    || !b1.contains("**275.3 s**")) {
                throw new IllegalStateException("B1 PROTOCOL_PAINT_HOLD.md must remain the 60 s / 275.3 s parent");
            }
        }
    }

    private static double atTime(LaneBDish.Result r, boolean cell, double tq) {
        if (r == null) return Double.NaN;
        return LaneBDish.interpolate(r.time, cell ? r.kappaCell : r.kappaField, tq);
    }

    private static double deltaLAt(LaneBDish.Result r, double tq) {
        if (r == null || r.deltaLEll == null) return Double.NaN;
        return LaneBDish.interpolate(r.time, r.deltaLEll, tq);
    }

    private static boolean gateOccupation(double kMot, double kOff, double kSilent) {
        if (!Double.isFinite(kMot) || !Double.isFinite(kOff) || !Double.isFinite(kSilent)) {
            return false;
        }
        return Math.abs(kMot) >= LaneB2Identity.KAPPA_MOTILE_MIN
                && Math.abs(kMot) - Math.abs(kOff) >= LaneB2Identity.KAPPA_MINUS_OFF
                && Math.abs(kSilent) < LaneB2Identity.KAPPA_SILENT_MAX;
    }

    private static boolean gateHold(double kMot, double kOff, double kField) {
        if (!Double.isFinite(kMot) || !Double.isFinite(kOff) || !Double.isFinite(kField)) {
            return false;
        }
        return Math.abs(kMot) - Math.abs(kField) >= LaneB2Identity.KAPPA_MINUS_FIELD
                && Math.abs(kMot) - Math.abs(kOff) >= LaneB2Identity.KAPPA_MINUS_OFF;
    }

    private static boolean gateLedger(LaneBDish.Result... chemistry) {
        for (LaneBDish.Result r : chemistry) {
            if (r == null) return false;
            if (!r.ledgerPass) return false;
        }
        return true;
    }

    private static void writePaintCsv(String name, LaneBDish.Result r) throws IOException {
        if (r == null) return;
        Path p = RESULTS.resolve(name);
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(p, StandardCharsets.UTF_8))) {
            w.println("t,kappa_cell,kappa_field,mean_abs_L,median_abs_dL_ell");
            for (int i = 0; i < r.time.length; i++) {
                if (i > 0 && r.time[i] < r.time[i - 1]) break;
                double ell = r.meanAbsL == null ? Double.NaN : r.meanAbsL[i];
                double dL = r.deltaLEll == null ? Double.NaN : r.deltaLEll[i];
                w.printf(Locale.US, "%.6g,%.8g,%.8g,%.8g,%.8g%n",
                        r.time[i], r.kappaCell[i], r.kappaField[i], ell, dL);
            }
        }
    }

    private static void writeSummary(
            boolean pass, boolean occPass, boolean holdPass, boolean ledgerPass,
            LaneBDish.Result stripeMotile, LaneBDish.Result stripeOff,
            LaneBDish.Result stripeField, LaneBDish.Result stripeSilent,
            double kMotOcc, double kOffOcc, double kFieldOcc, double kSilentOcc,
            double kMotHold, double kOffHold, double kFieldHold, double kSilentHold,
            double dLOcc, boolean fieldFadeOk, String killer,
            List<LaneBDish.Result> runs) throws IOException {
        Path p = RESULTS.resolve("b2_summary.json");
        StringBuilder sb = new StringBuilder();
        sb.append("{\n");
        sb.append("  \"object\": \"").append(LaneB2Identity.OBJECT).append("\",\n");
        sb.append("  \"gate\": \"").append(LaneB2Identity.GATE).append("\",\n");
        sb.append("  \"device\": \"").append(LaneB2Identity.DEVICE).append("\",\n");
        sb.append("  \"parent_fail_b0\": \"").append(LaneB2Identity.PARENT_FAIL_B0).append("\",\n");
        sb.append("  \"parent_fail_b1\": \"").append(LaneB2Identity.PARENT_FAIL_B1).append("\",\n");
        sb.append("  \"NOT_FIG4B\": true,\n");
        sb.append("  \"not_lane_a\": true,\n");
        sb.append("  \"narma\": false,\n");
        sb.append("  \"b0_overall_rewrite\": false,\n");
        sb.append("  \"b1_overall_rewrite\": false,\n");
        sb.append("  \"pass\": ").append(pass).append(",\n");
        sb.append("  \"B2.1_occupation\": ").append(occPass).append(",\n");
        sb.append("  \"B2.2_hold\": ").append(holdPass).append(",\n");
        sb.append("  \"B2.3_ledger\": ").append(ledgerPass).append(",\n");
        sb.append("  \"killer\": \"").append(killer).append("\",\n");
        sb.append(String.format(Locale.US, "  \"L_slab\": %.8g,%n", LaneB2Identity.L_SLAB));
        sb.append(String.format(Locale.US, "  \"t_occ_s\": %.8g,%n", LaneB2Identity.T_OCC));
        sb.append(String.format(Locale.US, "  \"T_hold_s\": %.8g,%n", LaneB2Identity.T_HOLD));
        sb.append(String.format(Locale.US, "  \"kappa_motile_occ\": %.8g,%n", kMotOcc));
        sb.append(String.format(Locale.US, "  \"kappa_off_occ\": %.8g,%n", kOffOcc));
        sb.append(String.format(Locale.US, "  \"kappa_field_occ\": %.8g,%n", kFieldOcc));
        sb.append(String.format(Locale.US, "  \"kappa_silent_occ\": %.8g,%n", kSilentOcc));
        sb.append(String.format(Locale.US, "  \"kappa_motile_hold\": %.8g,%n", kMotHold));
        sb.append(String.format(Locale.US, "  \"kappa_off_hold\": %.8g,%n", kOffHold));
        sb.append(String.format(Locale.US, "  \"kappa_field_hold\": %.8g,%n", kFieldHold));
        sb.append(String.format(Locale.US, "  \"kappa_silent_hold\": %.8g,%n", kSilentHold));
        sb.append(String.format(Locale.US, "  \"delta_L_ell_t_occ\": %.8g,%n", dLOcc));
        sb.append("  \"java_field_fade_ok\": ").append(fieldFadeOk).append(",\n");
        sb.append("  \"L_slab_raised\": false,\n");
        sb.append("  \"D_dropped\": false,\n");
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
