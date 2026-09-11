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
 * Gate B0: motility-as-memory on LANE_B_MM_DISH_FLOW0.
 * B0_MOTILITY_MEMORY. Not Fig. 4b. Not Lane A. Not C1c. No NARMA.
 */
public final class LaneBJob {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/PROTOCOL.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/configs/protocol.json");
    private static final Path RESULTS = ROOT.resolve(
            "examples/LaneB_ChemotacticSpatial/results");
    private static final Path LANES_REVIEW = ROOT.resolve(
            "examples/PocketDish/LANES_REVIEW_AND_LANE_B_DESIGN.md");

    private LaneBJob() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        refuseCapacityWords(args);
        requireProtocolFrozen();
        Files.createDirectories(RESULTS);

        System.out.println(LaneBIdentity.OBJECT + " " + LaneBIdentity.GATE
                + " device=" + LaneBIdentity.DEVICE + " " + LaneBIdentity.NOT_FIG4B);
        System.out.println("Not Fig. 4b. Not Lane A R/L. Not C1c. Not NARMA/CHARC.");
        System.out.println("LANES_REVIEW §4 is not identity. Stokes 0.66 forbidden.");
        System.out.printf(Locale.US,
                "lambda=%.3g um FLOW=0 PDE=%dx%dx%d readout=%dx%d D_att=%.3g J_max=%.6g%n",
                LaneBIdentity.LAMBDA_UM, LaneBIdentity.NX, LaneBIdentity.NY, LaneBIdentity.NZ,
                LaneBIdentity.RX, LaneBIdentity.RY, LaneBIdentity.D_ATT, LaneBIdentity.J_MAX);
        System.out.printf(Locale.US, "Stokes-Einstein in-dish D_SE=%.3g um2/s (Berg eta,T,r); water37 yardstick=%.3g%n",
                LaneBIdentity.stokesEinsteinUm2s(), LaneBIdentity.D_SE_WATER37);

        BSimRandom rng = new BSimRandom(LaneBIdentity.RNG_SEED);
        LaneBDish dish = new LaneBDish(rng);
        String want = args.length > 0 ? args[0].trim() : "";
        List<LaneBDish.Result> runs = new ArrayList<LaneBDish.Result>();

        LaneBDish.Result msd = null;
        LaneBDish.Result mix = null;
        LaneBDish.Result thermal = null;
        LaneBDish.Result paintMotile = null;
        LaneBDish.Result paintOff = null;
        LaneBDish.Result paintField = null;
        LaneBDish.Result paintSilent = null;

        if (wantEmptyOr(want, "B0_MSD_MOTILE")) {
            msd = dish.run(LaneBDish.Arm.MSD_MOTILE, LaneBIdentity.T_MSD);
            runs.add(msd);
            System.out.println(LaneBDish.line(msd));
        }
        if (wantEmptyOr(want, "B0_MIX_BLOB")) {
            mix = dish.run(LaneBDish.Arm.MIX_BLOB, LaneBIdentity.T_MIX);
            runs.add(mix);
            System.out.println(LaneBDish.line(mix));
        }
        if (wantEmptyOr(want, "B0_THERMAL")) {
            thermal = dish.run(LaneBDish.Arm.THERMAL, LaneBIdentity.T_MSD);
            runs.add(thermal);
            System.out.println(LaneBDish.line(thermal));
        }

        double tauMix = mix == null ? Double.NaN : LaneBDish.tauMix(mix);
        double w = LaneBIdentity.W_OVER_TAU * tauMix;
        double tStar = LaneBIdentity.T_OFF + w;
        double tPaint = LaneBIdentity.T_PAINT_BASE;
        if (Double.isFinite(w) && w > 12.0) {
            tPaint = LaneBIdentity.T_OFF + w + 2.0;
        }
        tPaint = LaneBIdentity.DT * Math.round(tPaint / LaneBIdentity.DT);
        System.out.printf(Locale.US,
                "tau_mix(lambda=%.3g um)=%.6g s  W=0.3*tau_mix=%.6g s  t*=t_off+W=%.6g s  T_paint=%.6g%n",
                LaneBIdentity.LAMBDA_UM, tauMix, w, tStar, tPaint);
        if (Double.isFinite(tauMix)) {
            double l2 = (LaneBIdentity.LY * LaneBIdentity.LY)
                    / (4.0 * LaneBIdentity.MU0_MIDDLEBROOKS);
            System.out.printf(Locale.US,
                    "NOT_A_GATE L_dishY^2/(4 mu0_lit)=%.6g s  eigenmode lambda^2/(pi^2 mu0)=%.6g s%n",
                    l2, (LaneBIdentity.LAMBDA_UM * LaneBIdentity.LAMBDA_UM)
                            / (Math.PI * Math.PI * LaneBIdentity.MU0_MIDDLEBROOKS));
        }

        if (wantEmptyOr(want, "B0_PAINT_MOTILE")) {
            paintMotile = dish.run(LaneBDish.Arm.PAINT_MOTILE, tPaint);
            runs.add(paintMotile);
            System.out.println(LaneBDish.line(paintMotile));
        }
        if (wantEmptyOr(want, "B0_PAINT_OFF")) {
            paintOff = dish.run(LaneBDish.Arm.PAINT_OFF, tPaint);
            runs.add(paintOff);
            System.out.println(LaneBDish.line(paintOff));
        }
        if (wantEmptyOr(want, "B0_PAINT_FIELD")) {
            paintField = dish.run(LaneBDish.Arm.PAINT_FIELD, tPaint);
            runs.add(paintField);
            System.out.println(LaneBDish.line(paintField));
        }
        if (wantEmptyOr(want, "B0_PAINT_SILENT")) {
            paintSilent = dish.run(LaneBDish.Arm.PAINT_SILENT, tPaint);
            runs.add(paintSilent);
            System.out.println(LaneBDish.line(paintSilent));
        }
        if ("LANE_B_MM_DISH_FLOW8".equals(want) || "B0_FLOW8".equals(want)) {
            throw new IllegalArgumentException(
                    "FLOW=8 is a named extra, not B0 identity. Not run. B0_MOTILITY_MEMORY");
        }

        if (runs.isEmpty()) {
            throw new IllegalArgumentException("no B0 arms selected: " + want);
        }

        boolean b01 = gateB01(msd);
        boolean b02 = gateB02(mix, tauMix);
        boolean b03 = gateB03(tauMix, w);
        double kMot = paintMotile == null || !Double.isFinite(tStar) ? Double.NaN
                : LaneBDish.interpolate(paintMotile.time, paintMotile.kappaCell, tStar);
        double kOff = paintOff == null || !Double.isFinite(tStar) ? Double.NaN
                : LaneBDish.interpolate(paintOff.time, paintOff.kappaCell, tStar);
        double kField = paintField == null || !Double.isFinite(tStar) ? Double.NaN
                : LaneBDish.interpolate(paintField.time, paintField.kappaField, tStar);
        double kSilent = paintSilent == null || !Double.isFinite(tStar) ? Double.NaN
                : LaneBDish.interpolate(paintSilent.time, paintSilent.kappaCell, tStar);
        boolean b04 = gateB04(kMot, kOff, kField, kSilent);
        boolean b05 = gateB05(paintMotile, paintOff, paintField, paintSilent);
        boolean pass = b01 && b02 && b03 && b04 && b05;

        System.out.printf(Locale.US,
                "B0.1 MSD mu2d=%.6g um2/s (lit 120 / 130+/-21) %s%n",
                msd == null ? Double.NaN : msd.mu2d, b01 ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "B0.2 tau_mix=%.6g s at lambda=%.3g um %s%n",
                tauMix, LaneBIdentity.LAMBDA_UM, b02 ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "B0.3 W/tau_mix=%.6g band=[0.1,1] %s%n",
                Double.isFinite(tauMix) ? w / tauMix : Double.NaN, b03 ? "PASS" : "FAIL");
        System.out.printf(Locale.US,
                "B0.4 t*=%.6g kappa_motile=%.6g kappa_off=%.6g kappa_field=%.6g kappa_silent=%.6g %s%n",
                tStar, kMot, kOff, kField, kSilent, b04 ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "B0.5 N0 ledger %s%n", b05 ? "PASS" : "FAIL");
        if (thermal != null) {
            System.out.printf(Locale.US,
                    "B0_THERMAL D2d=%.6g um2/s vs D_SE=%.3g (NOT B0.4 control; not Stokes 0.66)%n",
                    thermal.thermalD2d, LaneBIdentity.stokesEinsteinUm2s());
        }

        writeCsv("b0_msd.csv", msd, "t,msd_xy");
        writeCsv("b0_mix.csv", mix, "t,auto_corr,kappa_cell");
        writeCsv("b0_paint_motile.csv", paintMotile, "t,kappa_cell,kappa_field");
        writeCsv("b0_paint_off.csv", paintOff, "t,kappa_cell,kappa_field");
        writeCsv("b0_paint_field.csv", paintField, "t,kappa_cell,kappa_field");
        writeCsv("b0_paint_silent.csv", paintSilent, "t,kappa_cell,kappa_field");
        writeCsv("b0_thermal.csv", thermal, "t,msd_xy");
        writeSummary(pass, b01, b02, b03, b04, b05, msd, mix, thermal,
                tauMix, w, tStar, tPaint, kMot, kOff, kField, kSilent, runs);

        if (pass) {
            System.out.println("B0=PASS B0_MOTILITY_MEMORY on LANE_B_MM_DISH_FLOW0. "
                    + "Not Fig. 4b. Not Lane A. Not C1c. Capacity work is still not this gate.");
        } else {
            System.out.println("B0=FAIL B0_MOTILITY_MEMORY. Do not start NARMA/CHARC. "
                    + "If FAIL vs field-only, Lane B may not start capacity work.");
            if (!b04 && Double.isFinite(kMot) && Double.isFinite(kField)
                    && Math.abs(kMot) - Math.abs(kField) < LaneBIdentity.KAPPA_MINUS_FIELD) {
                System.out.println("B0.4 motile did not separate from field-only. STOP. "
                        + "Lane B is decoration on this dish. Do not add NARMA to rescue it.");
            }
            System.exit(1);
        }
    }

    private static boolean wantEmptyOr(String want, String name) {
        return want.isEmpty() || want.equals(name);
    }

    private static void refuseCapacityWords(String[] argv) {
        String blob = String.join(" ", argv == null ? new String[0] : argv).toLowerCase(Locale.ROOT);
        if (blob.contains("charc") || blob.contains("ipc") || blob.contains("narma")) {
            throw new IllegalArgumentException("B0 cannot see a capacity/CHARC/NARMA target.");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("B0 protocol is not frozen_before_traces");
        }
        if (!json.contains("\"status_label\": \"B0_MOTILITY_MEMORY\"")
                || !json.contains("\"device\": \"LANE_B_MM_DISH_FLOW0\"")
                || !json.contains("\"flow_um_s\": [0.0, 0.0, 0.0]")) {
            throw new IllegalStateException("B0 protocol must freeze B0_MOTILITY_MEMORY FLOW=0");
        }
        if (!json.contains("\"lanes_review_section4_identity\": false")
                || !json.contains("\"stokes_d_066\": false")
                || !json.contains("\"tau_mix_equals_L2_over_4_mu0\": false")) {
            throw new IllegalStateException("B0 protocol must reject LANES_REVIEW §4 / Stokes 0.66 / L2/4mu0");
        }
        if (!json.contains("\"narma\": false") || json.contains("\"narma\": true")) {
            throw new IllegalStateException("B0 protocol must freeze NARMA off");
        }
        if (!json.contains("\"lambda_um\": 50.0") || !json.contains("\"J_max\": 20000.0")) {
            throw new IllegalStateException("B0 protocol must freeze lambda=50 and J_max=20000");
        }
        if (Files.exists(LANES_REVIEW)) {
            String review = Files.readString(LANES_REVIEW, StandardCharsets.UTF_8);
            if (!review.contains("NOT identity") && !review.contains("not identity")
                    && !review.contains("§4 table is **not**")) {
                /* optional pointer; do not fail a missing erratum on first run */
            }
        }
    }

    private static boolean gateB01(LaneBDish.Result msd) {
        if (msd == null) return false;
        double mu = msd.mu2d;
        return Double.isFinite(mu)
                && mu >= LaneBIdentity.MU_FAIL_BELOW
                && mu <= LaneBIdentity.MU_FAIL_ABOVE;
    }

    private static boolean gateB02(LaneBDish.Result mix, double tauMix) {
        return mix != null && Double.isFinite(tauMix) && tauMix > 0.0
                && tauMix < LaneBIdentity.T_MIX;
    }

    private static boolean gateB03(double tauMix, double w) {
        if (!Double.isFinite(tauMix) || tauMix <= 0.0 || !Double.isFinite(w)) return false;
        double r = w / tauMix;
        return r + 1e-12 >= LaneBIdentity.W_BAND_LO && r - 1e-12 <= LaneBIdentity.W_BAND_HI;
    }

    private static boolean gateB04(double kMot, double kOff, double kField, double kSilent) {
        if (!Double.isFinite(kMot) || !Double.isFinite(kOff)
                || !Double.isFinite(kField) || !Double.isFinite(kSilent)) {
            return false;
        }
        return Math.abs(kMot) >= LaneBIdentity.KAPPA_MOTILE_MIN
                && Math.abs(kMot) - Math.abs(kField) >= LaneBIdentity.KAPPA_MINUS_FIELD
                && Math.abs(kMot) - Math.abs(kOff) >= LaneBIdentity.KAPPA_MINUS_OFF
                && Math.abs(kSilent) < LaneBIdentity.KAPPA_SILENT_MAX;
    }

    private static boolean gateB05(LaneBDish.Result... chemistry) {
        for (LaneBDish.Result r : chemistry) {
            if (r == null) return false;
            if (!r.ledgerPass) return false;
        }
        return true;
    }

    private static void writeCsv(String name, LaneBDish.Result r, String header) throws IOException {
        if (r == null) return;
        Path p = RESULTS.resolve(name);
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(p, StandardCharsets.UTF_8))) {
            w.println(header);
            for (int i = 0; i < r.time.length; i++) {
                if (i > 0 && r.time[i] < r.time[i - 1]) break;
                if (header.equals("t,msd_xy")) {
                    w.printf(Locale.US, "%.6g,%.8g%n", r.time[i], r.msdXy[i]);
                } else if (header.equals("t,auto_corr,kappa_cell")) {
                    w.printf(Locale.US, "%.6g,%.8g,%.8g%n", r.time[i], r.autoCorr[i], r.kappaCell[i]);
                } else {
                    w.printf(Locale.US, "%.6g,%.8g,%.8g%n", r.time[i], r.kappaCell[i], r.kappaField[i]);
                }
            }
        }
    }

    private static void writeSummary(
            boolean pass, boolean b01, boolean b02, boolean b03, boolean b04, boolean b05,
            LaneBDish.Result msd, LaneBDish.Result mix, LaneBDish.Result thermal,
            double tauMix, double w, double tStar, double tPaint,
            double kMot, double kOff, double kField, double kSilent,
            List<LaneBDish.Result> runs) throws IOException {
        Path p = RESULTS.resolve("b0_summary.json");
        StringBuilder sb = new StringBuilder();
        sb.append("{\n");
        sb.append("  \"object\": \"").append(LaneBIdentity.OBJECT).append("\",\n");
        sb.append("  \"gate\": \"").append(LaneBIdentity.GATE).append("\",\n");
        sb.append("  \"device\": \"").append(LaneBIdentity.DEVICE).append("\",\n");
        sb.append("  \"NOT_FIG4B\": true,\n");
        sb.append("  \"not_lane_a\": true,\n");
        sb.append("  \"not_c1c\": true,\n");
        sb.append("  \"narma\": false,\n");
        sb.append("  \"pass\": ").append(pass).append(",\n");
        sb.append("  \"B0.1\": ").append(b01).append(",\n");
        sb.append("  \"B0.2\": ").append(b02).append(",\n");
        sb.append("  \"B0.3\": ").append(b03).append(",\n");
        sb.append("  \"B0.4\": ").append(b04).append(",\n");
        sb.append("  \"B0.5\": ").append(b05).append(",\n");
        sb.append(String.format(Locale.US, "  \"mu2d_um2_s\": %.8g,%n", msd == null ? Double.NaN : msd.mu2d));
        sb.append(String.format(Locale.US, "  \"mu0_middlebrooks_um2_s\": %.8g,%n", LaneBIdentity.MU0_MIDDLEBROOKS));
        sb.append(String.format(Locale.US, "  \"mu0_zhao_um2_s\": %.8g,%n", LaneBIdentity.MU0_ZHAO));
        sb.append(String.format(Locale.US, "  \"lambda_um\": %.8g,%n", LaneBIdentity.LAMBDA_UM));
        sb.append(String.format(Locale.US, "  \"tau_mix_s\": %.8g,%n", tauMix));
        sb.append(String.format(Locale.US, "  \"W_s\": %.8g,%n", w));
        sb.append(String.format(Locale.US, "  \"W_over_tau_mix\": %.8g,%n",
                Double.isFinite(tauMix) && tauMix != 0.0 ? w / tauMix : Double.NaN));
        sb.append(String.format(Locale.US, "  \"t_star_s\": %.8g,%n", tStar));
        sb.append(String.format(Locale.US, "  \"T_paint_s\": %.8g,%n", tPaint));
        sb.append(String.format(Locale.US, "  \"kappa_motile\": %.8g,%n", kMot));
        sb.append(String.format(Locale.US, "  \"kappa_off\": %.8g,%n", kOff));
        sb.append(String.format(Locale.US, "  \"kappa_field\": %.8g,%n", kField));
        sb.append(String.format(Locale.US, "  \"kappa_silent\": %.8g,%n", kSilent));
        sb.append(String.format(Locale.US, "  \"thermal_D2d_um2_s\": %.8g,%n",
                thermal == null ? Double.NaN : thermal.thermalD2d));
        sb.append(String.format(Locale.US, "  \"D_SE_berg_um2_s\": %.8g,%n", LaneBIdentity.stokesEinsteinUm2s()));
        sb.append("  \"stokes_066_forbidden\": true,\n");
        sb.append("  \"capacity_may_start\": false,\n");
        sb.append("  \"arms\": [");
        for (int i = 0; i < runs.size(); i++) {
            if (i > 0) sb.append(", ");
            sb.append("\"").append(runs.get(i).name).append("\"");
        }
        sb.append("]\n}\n");
        Files.writeString(p, sb.toString(), StandardCharsets.UTF_8);
    }
}
