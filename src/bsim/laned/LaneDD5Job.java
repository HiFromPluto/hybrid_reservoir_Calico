package bsim.laned;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;

/**
 * D5_GATES. Regression is the checker. This job: memory sweep,
 * feast-famine point, refuse-concurrency. Not a new D0 feast.
 */
public final class LaneDD5Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneD_LifeCycle/PROTOCOL_D5.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneD_LifeCycle/configs/protocol_d5.json");
    private static final Path D4_STANDING = ROOT.resolve(
            "examples/PocketDish/LANE_D_D4_THREE_PHASE_STANDING.md");
    private static final Path RESULTS = ROOT.resolve(
            "examples/LaneD_LifeCycle/results");

    private LaneDD5Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        refuseCapacity(args);
        requireFrozen();
        Files.createDirectories(RESULTS);

        System.out.println(LaneDD5Identity.OBJECT + " " + LaneDD5Identity.GATE
                + " device=" + LaneDD5Identity.DEVICE);
        System.out.println("phi_X=phi_Rb  not C_up/C_down/C_flat");
        System.out.println("Not a C_s sweep. Not invented EV1 rows.");
        System.out.println("Not LC1 restage. Not NARMA.");

        boolean refuseReplete = probeRefuse(
                LaneDD0Identity.C_S_MM, true, true);
        boolean allowCut = probeRefuse(LaneDD0Identity.C_CUT_MM, true, false);
        boolean allowBelow = probeRefuse(0.5 * LaneDD0Identity.C_CUT_MM, true, false);

        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(
                RESULTS.resolve("d5_memory.csv"), StandardCharsets.UTF_8))) {
            w.println("mu_per_h,phi_Rb,gamma_law,N_end,gamma_meas");
            for (double mu : LaneDD5Identity.SWEEP_MU) {
                double phi = LaneDD5Identity.phiRbFromMu(mu);
                double gamma = LaneDD4Identity.gammaPerDay(mu);
                double n = famine(LaneDD5Identity.N0_FAMINE, gamma);
                double gMeas = -Math.log(n / LaneDD5Identity.N0_FAMINE)
                        / LaneDD4Identity.T_FAMINE_D;
                w.printf(Locale.US, "%.10f,%.16e,%.16e,%.16e,%.16e%n",
                        mu, phi, gamma, n, gMeas);
                System.out.printf(Locale.US,
                        "D5_MEMORY mu=%.4f phi_Rb=%.6f gamma=%.6f "
                                + "gamma_meas=%.6f%n",
                        mu, phi, gamma, gMeas);
            }
        }

        double muStar = LaneDD5Identity.optimumMu(
                LaneDD5Identity.T_PLUS_H, LaneDD5Identity.T_MINUS_D);
        double muNum = LaneDD5Identity.numericalOptimum(
                LaneDD5Identity.T_PLUS_H, LaneDD5Identity.T_MINUS_D);
        double gStar = LaneDD4Identity.gammaPerDay(muStar);
        double muOld = LaneDD5Identity.optimumMu(
                LaneDD5Identity.T_PLUS_H, LaneDD5Identity.OLD_T_MINUS_D);
        double gOld = LaneDD4Identity.gammaPerDay(muOld);

        System.out.printf(Locale.US,
                "D5_OPTIMUM 3h/6d mu*=%.6f gamma=%.6f numeric=%.6f%n",
                muStar, gStar, muNum);
        System.out.printf(Locale.US,
                "D5_OPTIMUM 3h/3d mu*=%.6f gamma=%.6f (must miss)%n",
                muOld, gOld);
        System.out.printf(Locale.US,
                "D5_REFUSE replete_threw=%s cut_allowed=%s%n",
                refuseReplete, allowCut && allowBelow);

        Files.writeString(RESULTS.resolve("d5_summary.json"), String.format(
                Locale.US,
                "{\n  \"gate\": \"D5_GATES\",\n"
                        + "  \"object\": \"LANE_D_LIFE_CYCLE\",\n"
                        + "  \"device\": \"LANE_D_CLOSED_BATH\",\n"
                        + "  \"phi_x\": \"phi_Rb\",\n"
                        + "  \"narma\": false,\n"
                        + "  \"invent_ev1_rows\": false,\n"
                        + "  \"refuse_replete_threw\": %s,\n"
                        + "  \"refuse_cut_allowed\": %s,\n"
                        + "  \"mu_star\": %.16e,\n"
                        + "  \"mu_star_numeric\": %.16e,\n"
                        + "  \"gamma_star\": %.16e,\n"
                        + "  \"mu_old\": %.16e,\n"
                        + "  \"gamma_old\": %.16e\n}\n",
                refuseReplete, allowCut && allowBelow,
                muStar, muNum, gStar, muOld, gOld),
                StandardCharsets.UTF_8);
    }

    private static boolean probeRefuse(double cMm, boolean deathOn,
            boolean expectThrow) {
        LaneDClosedBath bath = new LaneDClosedBath(cMm, false);
        try {
            bath.refuseDeath(deathOn);
            if (expectThrow) {
                System.out.println("D5_REFUSE FAIL: replete death-on did not throw");
                return false;
            }
            return true;
        } catch (IllegalStateException ex) {
            boolean ok = expectThrow && ex.getMessage() != null
                    && ex.getMessage().contains("refuses death-on while C > C_cut");
            if (!ok) {
                System.out.println("D5_REFUSE unexpected: " + ex.getMessage());
            }
            return ok;
        }
    }

    private static double famine(double n0, double gamma) {
        double n = n0;
        double t = 0.0;
        while (t < LaneDD4Identity.T_FAMINE_D - 1.0e-15) {
            double dt = Math.min(LaneDD4Identity.FAMINE_DT_D,
                    LaneDD4Identity.T_FAMINE_D - t);
            double k1 = -gamma * n;
            double k2 = -gamma * (n + 0.5 * dt * k1);
            double k3 = -gamma * (n + 0.5 * dt * k2);
            double k4 = -gamma * (n + dt * k3);
            n += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            t += dt;
        }
        return n;
    }

    private static void refuseCapacity(String[] args) {
        String blob = String.join(" ", args == null ? new String[0] : args)
                .toLowerCase(Locale.ROOT);
        if (blob.contains("narma") || blob.contains("charc") || blob.contains("ipc")) {
            throw new IllegalArgumentException("D5 refuses NARMA/CHARC/IPC");
        }
    }

    private static void requireFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        String d4 = Files.readString(D4_STANDING, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces")
                || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("D5 protocol is not frozen");
        }
        if (!json.contains("\"status_label\": \"D5_GATES\"")
                || !json.contains("\"phi_x\": \"phi_Rb\"")
                || !json.contains("\"invent_ev1_rows\": false")
                || !json.contains("\"lc1_restage\": false")) {
            throw new IllegalStateException("D5 freeze incomplete");
        }
        if (!d4.contains("**Status: PASS**")
                || !d4.contains("D4_THREE_PHASE")) {
            throw new IllegalStateException("D5 requires D4 PASS");
        }
    }
}
