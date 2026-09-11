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
 * Gate B3: sustained 200 µm stripe for a 25 s swim.
 * B3_SUSTAINED_STRIPE. Parent B0, B1, B2 FAILs are kept.
 * Not Fig. 4b. Not Lane A. Not C1c. No NARMA. Not a B2 L_slab raise.
 */
public final class LaneB3Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/PROTOCOL_SUSTAINED_STRIPE.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/configs/protocol_b3.json");
    private static final Path B0_PROTOCOL = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/PROTOCOL.md");
    private static final Path B1_PROTOCOL = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/PROTOCOL_PAINT_HOLD.md");
    private static final Path B2_PROTOCOL = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/PROTOCOL_STRIPE_OCCUPY.md");
    private static final Path RESULTS = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/results");

    private LaneB3Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        refuseCapacityWords(args);
        requireProtocolFrozen();
        Files.createDirectories(RESULTS);

        System.out.println(LaneB3Identity.OBJECT + " " + LaneB3Identity.GATE
                + " device=" + LaneB3Identity.DEVICE + " " + LaneB3Identity.NOT_FIG4B);
        System.out.println("Parent FAILs kept: " + LaneB3Identity.PARENT_FAIL_B0
                + ", " + LaneB3Identity.PARENT_FAIL_B1
                + ", " + LaneB3Identity.PARENT_FAIL_B2
                + ". Not a B0/B1/B2 rewrite. Not half-dish kappa. Not NARMA/CHARC.");
        System.out.printf(Locale.US,
                "J_slab=%.6g t_occ=%.3g s T_hold=%.3g s T_live=%.3g s lambda=%.3g um%n",
                LaneB3Identity.J_SLAB, LaneB3Identity.T_OCC, LaneB3Identity.T_HOLD,
                LaneB3Identity.T_LIVE, LaneB3Identity.LAMBDA_UM);

        BSimRandom rng = new BSimRandom(LaneB3Identity.RNG_SEED);
        LaneBDish dish = new LaneBDish(rng);
        String want = args.length > 0 ? args[0].trim() : "";
        List<LaneBDish.Result> runs = new ArrayList<LaneBDish.Result>();

        LaneBDish.Result stripeMotile = null;
        LaneBDish.Result stripeOff = null;
        LaneBDish.Result stripeField = null;
        LaneBDish.Result stripeSilent = null;

        if ("LANE_B_MM_DISH_FLOW8".equals(want) || "B3_FLOW8".equals(want)) {
            throw new IllegalArgumentException(
                    "FLOW=8 is a named extra, not B3 identity. Not run. B3_SUSTAINED_STRIPE");
        }

        double tPaint = LaneB3Identity.paintDuration();
        if (wantEmptyOr(want, "B3_STRIPE_MOTILE")) {
            stripeMotile = dish.run(LaneBDish.Arm.SUSTAIN_MOTILE, tPaint, LaneB3Identity.T_LIVE);
            runs.add(stripeMotile);
            System.out.println(LaneBDish.line(stripeMotile));
        }
        if (wantEmptyOr(want, "B3_STRIPE_OFF")) {
            stripeOff = dish.run(LaneBDish.Arm.SUSTAIN_OFF, tPaint, LaneB3Identity.T_LIVE);
            runs.add(stripeOff);
            System.out.println(LaneBDish.line(stripeOff));
        }
        if (wantEmptyOr(want, "B3_STRIPE_FIELD")) {
            stripeField = dish.run(LaneBDish.Arm.SUSTAIN_FIELD, tPaint, LaneB3Identity.T_LIVE);
            runs.add(stripeField);
            System.out.println(LaneBDish.line(stripeField));
        }
        if (wantEmptyOr(want, "B3_STRIPE_SILENT")) {
            stripeSilent = dish.run(LaneBDish.Arm.SUSTAIN_SILENT, tPaint, LaneB3Identity.T_LIVE);
            runs.add(stripeSilent);
            System.out.println(LaneBDish.line(stripeSilent));
        }
        if (runs.isEmpty()) {
            throw new IllegalArgumentException("no B3 arms selected: " + want);
        }

        double kMotOcc = atTime(stripeMotile, true, LaneB3Identity.T_OCC);
        double kOffOcc = atTime(stripeOff, true, LaneB3Identity.T_OCC);
        double kFieldOcc = atTime(stripeField, false, LaneB3Identity.T_OCC);
        double kSilentOcc = atTime(stripeSilent, true, LaneB3Identity.T_OCC);
        boolean occPass = gateOccupation(kMotOcc, kOffOcc, kSilentOcc);

        double kMotHold = atTime(stripeMotile, true, LaneB3Identity.T_HOLD);
        double kOffHold = atTime(stripeOff, true, LaneB3Identity.T_HOLD);
        double kFieldHold = atTime(stripeField, false, LaneB3Identity.T_HOLD);
        double kSilentHold = atTime(stripeSilent, true, LaneB3Identity.T_HOLD);
        boolean fieldFadeOk = Double.isFinite(kFieldHold)
                && Math.abs(kFieldHold) < LaneB3Identity.KAPPA_FIELD_HOLD_MAX;
        boolean holdPass = occPass && gateHold(kMotHold, kOffHold, kFieldHold);
        boolean ledgerPass = gateLedger(stripeMotile, stripeOff, stripeField, stripeSilent);

        double dLOcc = deltaLAt(stripeField, LaneB3Identity.T_OCC);
        if (!Double.isFinite(dLOcc)) dLOcc = deltaLAt(stripeMotile, LaneB3Identity.T_OCC);

        String killer = "none";
        boolean pass = occPass && holdPass && ledgerPass;
        if (!occPass) killer = "occupation";
        else if (!holdPass) killer = "hold";
        else if (!ledgerPass) killer = "ledger";

        System.out.printf(Locale.US,
                "B3.1 occupation t=%.3g kappa_motile=%.6g kappa_off=%.6g kappa_silent=%.6g kappa_field=%.6g %s%n",
                LaneB3Identity.T_OCC, kMotOcc, kOffOcc, kSilentOcc, kFieldOcc,
                occPass ? "PASS" : "FAIL");
        if (occPass) {
            System.out.printf(Locale.US,
                    "B3.2 hold T_hold=%.3g kappa_motile=%.6g kappa_off=%.6g kappa_field=%.6g kappa_silent=%.6g %s%n",
                    LaneB3Identity.T_HOLD, kMotHold, kOffHold, kFieldHold, kSilentHold,
                    holdPass ? "PASS" : "FAIL");
            if (!fieldFadeOk) {
                System.out.printf(Locale.US,
                        "SCOPE_NOTE scout vs Java: |kappa_lambda_field|(T_hold)=%.6g not < 0.15. Do not retune D.%n",
                        kFieldHold);
            }
        } else {
            System.out.println("B3.2 hold NOT_SCORED (occupation FAIL). Do not dress hold as a maybe.");
        }
        System.out.printf(Locale.US, "B3.3 N0 ledger %s%n", ledgerPass ? "PASS" : "FAIL");
        System.out.printf(Locale.US,
                "REPORT DeltaL(ell) at t=25=%.6g (Java field; not a gate)%n", dLOcc);

        writePaintCsv("b3_stripe_motile.csv", stripeMotile);
        writePaintCsv("b3_stripe_off.csv", stripeOff);
        writePaintCsv("b3_stripe_field.csv", stripeField);
        writePaintCsv("b3_stripe_silent.csv", stripeSilent);
        writeSummary(pass, occPass, holdPass, ledgerPass,
                stripeMotile, stripeOff, stripeField, stripeSilent,
                kMotOcc, kOffOcc, kFieldOcc, kSilentOcc,
                kMotHold, kOffHold, kFieldHold, kSilentHold, dLOcc,
                fieldFadeOk, killer, runs);

        if (pass) {
            System.out.println("B3_SUSTAINED_STRIPE=PASS on LANE_B_MM_DISH_FLOW0. "
                    + "B0/B1/B2 FAILs are not rewritten. Not NARMA. Not Fig. 4b.");
        } else {
            System.out.println("B3_SUSTAINED_STRIPE=FAIL " + killer
                    + ". Do not start NARMA/CHARC. Do not raise J_slab. Do not drop D. "
                    + "B0, B1, B2 FAILs stay.");
            System.exit(1);
        }
    }

    private static boolean wantEmptyOr(String want, String name) {
        return want.isEmpty() || want.equals(name);
    }

    private static void refuseCapacityWords(String[] argv) {
        String blob = String.join(" ", argv == null ? new String[0] : argv).toLowerCase(Locale.ROOT);
        if (blob.contains("charc") || blob.contains("ipc") || blob.contains("narma")) {
            throw new IllegalArgumentException("B3 cannot see a capacity/CHARC/NARMA target.");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("B3 protocol is not frozen_before_traces");
        }
        if (!json.contains("\"status_label\": \"B3_SUSTAINED_STRIPE\"")
                || !json.contains("\"device\": \"LANE_B_MM_DISH_FLOW0\"")
                || !json.contains("\"B0_MOTILITY_MEMORY\"")
                || !json.contains("\"B1_PAINT_HOLD\"")
                || !json.contains("\"B2_STRIPE_OCCUPY\"")) {
            throw new IllegalStateException("B3 protocol must freeze B3_SUSTAINED_STRIPE parent B0/B1/B2 FAILs");
        }
        if (!json.contains("\"narma\": false") || json.contains("\"narma\": true")) {
            throw new IllegalStateException("B3 protocol must freeze NARMA off");
        }
        if (!json.contains("\"T_hold_s\": 25.1") || !json.contains("\"t_occ_s\": 25.0")
                || !json.contains("\"T_live_s\": 25.0")) {
            throw new IllegalStateException("B3 protocol must freeze scout t_occ=25 and T_hold=25.1");
        }
        if (!json.contains("\"T_hold_status\": \"SCOUT_SET_BEFORE_MOTILE_TRACES\"")) {
            throw new IllegalStateException("B3 T_hold must be scout-set before motile traces");
        }
        if (!json.contains("\"D_att_um2_s\": 800.0")) {
            throw new IllegalStateException("B3 must keep D=800");
        }
        if (Files.exists(B0_PROTOCOL)) {
            String b0 = Files.readString(B0_PROTOCOL, StandardCharsets.UTF_8);
            if (!b0.contains("B0_MOTILITY_MEMORY") || !b0.contains("**10 s**")) {
                throw new IllegalStateException("B0 PROTOCOL.md must remain the 10 s parent; do not edit it");
            }
        }
        if (Files.exists(B1_PROTOCOL)) {
            String b1 = Files.readString(B1_PROTOCOL, StandardCharsets.UTF_8);
            if (!b1.contains("B1_PAINT_HOLD") || !b1.contains("**275.3 s**")) {
                throw new IllegalStateException("B1 PROTOCOL_PAINT_HOLD.md must remain the 275.3 s parent");
            }
        }
        if (Files.exists(B2_PROTOCOL)) {
            String b2 = Files.readString(B2_PROTOCOL, StandardCharsets.UTF_8);
            if (!b2.contains("B2_STRIPE_OCCUPY") || !b2.contains("**1.1** s")
                    || !b2.contains("**3.1** s")) {
                throw new IllegalStateException("B2 PROTOCOL_STRIPE_OCCUPY.md must remain the 1.1 / 3.1 s parent");
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
        return Math.abs(kMot) >= LaneB3Identity.KAPPA_MOTILE_MIN
                && Math.abs(kMot) - Math.abs(kOff) >= LaneB3Identity.KAPPA_MINUS_OFF
                && Math.abs(kSilent) < LaneB3Identity.KAPPA_SILENT_MAX;
    }

    private static boolean gateHold(double kMot, double kOff, double kField) {
        if (!Double.isFinite(kMot) || !Double.isFinite(kOff) || !Double.isFinite(kField)) {
            return false;
        }
        return Math.abs(kMot) - Math.abs(kField) >= LaneB3Identity.KAPPA_MINUS_FIELD
                && Math.abs(kMot) - Math.abs(kOff) >= LaneB3Identity.KAPPA_MINUS_OFF;
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
        Path p = RESULTS.resolve("b3_summary.json");
        StringBuilder sb = new StringBuilder();
        sb.append("{\n");
        sb.append("  \"object\": \"").append(LaneB3Identity.OBJECT).append("\",\n");
        sb.append("  \"gate\": \"").append(LaneB3Identity.GATE).append("\",\n");
        sb.append("  \"device\": \"").append(LaneB3Identity.DEVICE).append("\",\n");
        sb.append("  \"parent_fail_b0\": \"").append(LaneB3Identity.PARENT_FAIL_B0).append("\",\n");
        sb.append("  \"parent_fail_b1\": \"").append(LaneB3Identity.PARENT_FAIL_B1).append("\",\n");
        sb.append("  \"parent_fail_b2\": \"").append(LaneB3Identity.PARENT_FAIL_B2).append("\",\n");
        sb.append("  \"NOT_FIG4B\": true,\n");
        sb.append("  \"not_lane_a\": true,\n");
        sb.append("  \"narma\": false,\n");
        sb.append("  \"b0_overall_rewrite\": false,\n");
        sb.append("  \"b1_overall_rewrite\": false,\n");
        sb.append("  \"b2_overall_rewrite\": false,\n");
        sb.append("  \"pass\": ").append(pass).append(",\n");
        sb.append("  \"B3.1_occupation\": ").append(occPass).append(",\n");
        sb.append("  \"B3.2_hold\": ").append(holdPass).append(",\n");
        sb.append("  \"B3.3_ledger\": ").append(ledgerPass).append(",\n");
        sb.append("  \"killer\": \"").append(killer).append("\",\n");
        sb.append(String.format(Locale.US, "  \"J_slab\": %.8g,%n", LaneB3Identity.J_SLAB));
        sb.append(String.format(Locale.US, "  \"t_occ_s\": %.8g,%n", LaneB3Identity.T_OCC));
        sb.append(String.format(Locale.US, "  \"T_hold_s\": %.8g,%n", LaneB3Identity.T_HOLD));
        sb.append(String.format(Locale.US, "  \"kappa_motile_occ\": %.8g,%n", kMotOcc));
        sb.append(String.format(Locale.US, "  \"kappa_off_occ\": %.8g,%n", kOffOcc));
        sb.append(String.format(Locale.US, "  \"kappa_field_occ\": %.8g,%n", kFieldOcc));
        sb.append(String.format(Locale.US, "  \"kappa_silent_occ\": %.8g,%n", kSilentOcc));
        sb.append(String.format(Locale.US, "  \"kappa_motile_hold\": %.8g,%n", kMotHold));
        sb.append(String.format(Locale.US, "  \"kappa_off_hold\": %.8g,%n", kOffHold));
        sb.append(String.format(Locale.US, "  \"kappa_field_hold\": %.8g,%n", kFieldHold));
        sb.append(String.format(Locale.US, "  \"kappa_silent_hold\": %.8g,%n", kSilentHold));
        sb.append(String.format(Locale.US, "  \"delta_L_ell_t25\": %.8g,%n", dLOcc));
        sb.append("  \"java_field_fade_ok\": ").append(fieldFadeOk).append(",\n");
        sb.append("  \"J_slab_raised\": false,\n");
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
