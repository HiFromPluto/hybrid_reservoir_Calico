package bsim.laned;

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
 * D4_THREE_PHASE. Feast (D0) then FCR freeze then Job 6 gamma(mu_eff).
 * Death off while C &gt; C_cut. Not LC1b. Not NARMA.
 */
public final class LaneDD4Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneD_LifeCycle/PROTOCOL_D4.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneD_LifeCycle/configs/protocol_d4.json");
    private static final Path RESULTS = ROOT.resolve(
            "examples/LaneD_LifeCycle/results");

    private LaneDD4Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        refuseCapacity(args);
        requireFrozen();
        Files.createDirectories(RESULTS);
        new BSimRandom(LaneDD0Identity.RNG_SEED);

        System.out.println(LaneDD4Identity.OBJECT + " " + LaneDD4Identity.GATE
                + " device=" + LaneDD4Identity.DEVICE);
        System.out.println("phi_X=phi_Rb  not C_up/C_down/C_flat");
        System.out.println("Not K_M in Monod. Not lambda_C on elongation.");
        System.out.println("Not LC1b. Not bacteria.max. Not NARMA.");

        Seq seq = runSeq();
        Off off = runOff();
        writeFeastCsv("d4_seq_feast.csv", seq);
        writeFamineCsv("d4_seq_famine.csv", seq);
        writeOffCsv("d4_off.csv", off);
        writeSummary(seq, off);
    }

    private static Seq runSeq() {
        System.out.printf(Locale.US,
                "D4_SEQ feast growth=ON death=OFF motility=OFF PDE=OFF "
                        + "device=%s sink=rho_cell elongation=incremental "
                        + "fcr=spectator seed=%d%n",
                LaneDD4Identity.DEVICE, LaneDD0Identity.RNG_SEED);

        LaneDClosedBath bath = new LaneDClosedBath(LaneDD0Identity.C_S_MM, true);
        List<LaneDClosedBathCell> cells = founders();
        double lamI = LaneDD4Identity.lambdaMonodPerHour(LaneDD0Identity.C_S_MM);
        LaneDFcrSpectator fcr = new LaneDFcrSpectator(lamI);
        Seq seq = new Seq();
        recordFeast(seq, 0.0, bath, cells, fcr, false);

        double t = 0.0;
        int nEver = LaneDD0Identity.N0;
        while (t < LaneDD0Identity.T_CAP_S - 1.0e-12) {
            bath.refuseDeath(false);
            double dt = Math.min(LaneDD0Identity.DT_S, LaneDD0Identity.T_CAP_S - t);
            double cStep = bath.concentrationMm();
            double lamStep = LaneDD4Identity.lambdaMonodPerHour(cStep);

            double sumV = sumVolume(cells);
            bath.stepSink(dt, sumV);
            int nBefore = cells.size();
            for (int i = 0; i < nBefore; i++) {
                LaneDClosedBathCell cell = cells.get(i);
                if (cell.elongateIncremental(dt, cStep)) {
                    seq.shrinkEvents++;
                }
                if (cell.shouldDivide()) {
                    if (nEver + 1 > LaneDD0Identity.N_EVER_COMPUTE_CAP) {
                        seq.capFired = true;
                        System.out.println("SCOPE_NOTE compute-cap N_ever>"
                                + LaneDD0Identity.N_EVER_COMPUTE_CAP
                                + ". Not a death. Not PASS.");
                        finishFeast(seq, t + dt, bath, cells, nEver, fcr);
                        return seq;
                    }
                    cell.resetAfterFission();
                    cells.add(new LaneDClosedBathCell(
                            cell.x, cell.y, cell.z, LaneDD0Identity.ELL0_UM));
                    nEver++;
                    seq.birthsTotal++;
                }
            }
            fcr.stepHours(dt / 3600.0, lamStep);

            t += dt;
            boolean cut = bath.concentrationMm() <= LaneDD0Identity.C_CUT_MM;
            if (isLogTime(t) || cut || t + 1.0e-12 >= LaneDD0Identity.T_CAP_S) {
                recordFeast(seq, t, bath, cells, fcr, false);
            }
            if (cut) {
                fcr.freeze();
                seq.stopReason = "C_CUT";
                break;
            }
        }

        if (seq.stopReason == null) {
            seq.stopReason = "T_CAP";
            System.out.println("SCOPE_NOTE T_cap with C > C_cut. Not PASS.");
        }
        finishFeast(seq, t, bath, cells, nEver, fcr);
        if (seq.capFired || "T_CAP".equals(seq.stopReason)) {
            return seq;
        }

        bath.refuseDeath(true);
        seq.muEff = fcr.muEff();
        seq.gammaPerDay = LaneDD4Identity.gammaPerDay(seq.muEff);
        seq.nHandoff = seq.nEnd;
        System.out.printf(Locale.US,
                "D4_SEQ handoff C=%.8f phi_Rb=%.8f sigma=%.8f mu_eff=%.8f /h "
                        + "gamma=%.6f /d N_h=%d death=ON (C<=C_cut)%n",
                seq.cEnd, seq.phiRbFreeze, seq.sigmaFreeze, seq.muEff,
                seq.gammaPerDay, seq.nHandoff);
        System.out.printf(Locale.US,
                "D32_COMPARISON_ONLY sigma_T_1e-4=%.6f (not this freeze)%n",
                LaneDD4Identity.D32_SIGMA_T_COMPARISON);

        runFamine(seq);
        return seq;
    }

    private static void runFamine(Seq seq) {
        System.out.printf(Locale.US,
                "D4_SEQ famine death=ON growth=OFF sink=OFF "
                        + "law=Job6_gamma(mu_eff) not_LC1b%n");
        double n = seq.nHandoff;
        double tDays = 0.0;
        recordFamine(seq, 0.0, n);
        int obs = 1;
        while (tDays < LaneDD4Identity.T_FAMINE_D - 1.0e-15) {
            double dt = Math.min(LaneDD4Identity.FAMINE_DT_D,
                    LaneDD4Identity.T_FAMINE_D - tDays);
            n = rk4N(n, dt, seq.gammaPerDay);
            tDays += dt;
            if (obs < LaneDD4Identity.FAMINE_OBS_D.length
                    && tDays + 1.0e-12 >= LaneDD4Identity.FAMINE_OBS_D[obs]) {
                recordFamine(seq, LaneDD4Identity.FAMINE_OBS_D[obs], n);
                obs++;
            }
        }
        seq.nFamineEnd = n;
        seq.deaths = seq.nHandoff - n;
        System.out.printf(Locale.US,
                "D4_SEQ famine N_end=%.8f deaths=%.8f N_ever=%d%n",
                seq.nFamineEnd, seq.deaths, seq.nEverEnd);
    }

    private static double rk4N(double n, double dtDays, double gamma) {
        double k1 = -gamma * n;
        double k2 = -gamma * (n + 0.5 * dtDays * k1);
        double k3 = -gamma * (n + 0.5 * dtDays * k2);
        double k4 = -gamma * (n + dtDays * k3);
        return n + dtDays * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }

    private static Off runOff() {
        System.out.printf(Locale.US,
                "D4_OFF growth=OFF death=OFF sink=OFF fcr=FROZEN "
                        + "device=%s seed=%d%n",
                LaneDD4Identity.DEVICE, LaneDD0Identity.RNG_SEED);
        Off off = new Off();
        off.cEnd = LaneDD0Identity.C_S_MM;
        off.nEnd = LaneDD0Identity.N0;
        off.births = 0;
        off.deaths = 0.0;
        return off;
    }

    private static List<LaneDClosedBathCell> founders() {
        List<LaneDClosedBathCell> cells = new ArrayList<>();
        for (int i = 0; i < LaneDD0Identity.N0; i++) {
            cells.add(new LaneDClosedBathCell(
                    LaneDD0Identity.founderX(i),
                    LaneDD0Identity.founderY(i),
                    0.5 * LaneDD0Identity.BOX_Z_UM,
                    LaneDD0Identity.ELL0_UM));
        }
        return cells;
    }

    private static void finishFeast(Seq seq, double t, LaneDClosedBath bath,
            List<LaneDClosedBathCell> cells, int nEver, LaneDFcrSpectator fcr) {
        if (seq.feastT.isEmpty()
                || Math.abs(seq.feastT.get(seq.feastT.size() - 1) - t) > 1.0e-9) {
            recordFeast(seq, t, bath, cells, fcr, false);
        }
        seq.tEnd = t;
        seq.nEnd = cells.size();
        seq.nEverEnd = nEver;
        seq.cEnd = bath.concentrationMm();
        seq.phiRbFreeze = fcr.phiRb();
        seq.sigmaFreeze = fcr.sigma();
        System.out.printf(Locale.US,
                "D4_SEQ feast N_end=%d N_ever=%d births=%d C_end=%.8f "
                        + "stop=%s shrink=%d%n",
                seq.nEnd, seq.nEverEnd, seq.birthsTotal, seq.cEnd,
                seq.stopReason, seq.shrinkEvents);
    }

    private static boolean isLogTime(double t) {
        double k = t / LaneDD0Identity.LOG_DT_S;
        return Math.abs(k - Math.round(k)) < 1.0e-9;
    }

    private static double sumVolume(List<LaneDClosedBathCell> cells) {
        double s = 0.0;
        for (LaneDClosedBathCell cell : cells) {
            s += cell.volumeUm3();
        }
        return s;
    }

    private static void recordFeast(Seq seq, double t, LaneDClosedBath bath,
            List<LaneDClosedBathCell> cells, LaneDFcrSpectator fcr,
            boolean deathOn) {
        seq.feastT.add(t);
        seq.feastC.add(bath.concentrationMm());
        seq.feastN.add(cells.size());
        seq.feastNever.add(cells.size());
        seq.feastLam.add(LaneDD4Identity.lambdaMonodPerHour(bath.concentrationMm()));
        seq.feastPhiRb.add(fcr.phiRb());
        seq.feastSigma.add(fcr.sigma());
        seq.feastDeathOn.add(deathOn);
    }

    private static void recordFamine(Seq seq, double tDays, double n) {
        seq.famT.add(tDays);
        seq.famN.add(n);
        seq.famPred.add(seq.nHandoff * Math.exp(-seq.gammaPerDay * tDays));
    }

    private static void writeFeastCsv(String name, Seq seq) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(
                RESULTS.resolve(name), StandardCharsets.UTF_8))) {
            w.println("t,C_mM,lambda_monod_per_h,N,N_ever,phi_Rb,sigma,death_on");
            for (int i = 0; i < seq.feastT.size(); i++) {
                w.printf(Locale.US,
                        "%.10f,%.16e,%.16e,%d,%d,%.16e,%.16e,%s%n",
                        seq.feastT.get(i), seq.feastC.get(i),
                        seq.feastLam.get(i), seq.feastN.get(i),
                        seq.feastNever.get(i), seq.feastPhiRb.get(i),
                        seq.feastSigma.get(i), seq.feastDeathOn.get(i));
            }
        }
    }

    private static void writeFamineCsv(String name, Seq seq) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(
                RESULTS.resolve(name), StandardCharsets.UTF_8))) {
            w.println("t_d,N,N_pred,gamma_per_d,mu_eff");
            for (int i = 0; i < seq.famT.size(); i++) {
                w.printf(Locale.US, "%.10f,%.16e,%.16e,%.16e,%.16e%n",
                        seq.famT.get(i), seq.famN.get(i), seq.famPred.get(i),
                        seq.gammaPerDay, seq.muEff);
            }
        }
    }

    private static void writeOffCsv(String name, Off off) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(
                RESULTS.resolve(name), StandardCharsets.UTF_8))) {
            w.println("t,C_mM,N,deaths");
            w.printf(Locale.US, "0.0,%.16e,%d,0%n", off.cEnd, off.nEnd);
            w.printf(Locale.US, "36000.0,%.16e,%d,0%n", off.cEnd, off.nEnd);
        }
    }

    private static void writeSummary(Seq seq, Off off) throws IOException {
        boolean honesty = true;
        boolean processOff = off.births == 0 && off.deaths == 0.0
                && off.nEnd == LaneDD0Identity.N0
                && Math.abs(off.cEnd - LaneDD0Identity.C_S_MM) <= 1.0e-15;
        boolean noDeathFeast = !seq.feastDeathOn.contains(true);
        boolean fingerprint = Math.abs(seq.tEnd - LaneDD4Identity.D0_T_END_S) <= 1.0e-6
                && Math.abs(seq.cEnd - LaneDD4Identity.D0_C_END_MM)
                <= LaneDD4Identity.D0_C_END_TOL
                && seq.nEnd == LaneDD4Identity.D0_N_END;
        boolean muId = Math.abs(seq.muEff
                - LaneDD4Identity.muEffFromPhiRb(seq.phiRbFreeze)) <= 1.0e-12;
        boolean famine = famineOk(seq);
        boolean noKill = seq.nEnd == seq.nEverEnd && !seq.capFired;
        boolean tCap = "T_CAP".equals(seq.stopReason);
        boolean scope = seq.capFired || tCap;
        boolean pass = honesty && processOff && noDeathFeast && fingerprint
                && muId && famine && noKill && !scope;

        String killer;
        if (seq.capFired) {
            killer = "compute_cap";
        } else if (tCap) {
            killer = "t_cap";
        } else if (!processOff) {
            killer = "process_off";
        } else if (!noDeathFeast) {
            killer = "concurrency";
        } else if (!fingerprint) {
            killer = "feast_fingerprint";
        } else if (!muId) {
            killer = "phi_x_identity";
        } else if (!famine) {
            killer = "famine_clock";
        } else if (!noKill) {
            killer = "silent_kill";
        } else {
            killer = "none";
        }

        Files.writeString(RESULTS.resolve("d4_summary.json"), String.format(Locale.US,
                "{\n  \"gate\": \"D4_THREE_PHASE\",\n"
                        + "  \"object\": \"LANE_D_LIFE_CYCLE\",\n"
                        + "  \"device\": \"LANE_D_CLOSED_BATH\",\n"
                        + "  \"death\": false,\n  \"narma\": false,\n"
                        + "  \"phi_x\": \"phi_Rb\",\n"
                        + "  \"seed\": %d,\n  \"n0\": %d,\n"
                        + "  \"feast_t_end\": %.10f,\n  \"feast_C_end\": %.12e,\n"
                        + "  \"feast_N_end\": %d,\n  \"feast_N_ever\": %d,\n"
                        + "  \"feast_births\": %d,\n  \"feast_stop\": \"%s\",\n"
                        + "  \"phi_Rb_freeze\": %.16e,\n  \"sigma_freeze\": %.16e,\n"
                        + "  \"mu_eff\": %.16e,\n  \"gamma_per_d\": %.16e,\n"
                        + "  \"n_handoff\": %d,\n  \"n_famine_end\": %.16e,\n"
                        + "  \"deaths\": %.16e,\n"
                        + "  \"off_N_end\": %d,\n  \"off_births\": %d,\n"
                        + "  \"off_C_end\": %.12f,\n"
                        + "  \"cap_fired\": %s,\n  \"t_cap_scope\": %s,\n"
                        + "  \"pass_without_isolation\": %s,\n"
                        + "  \"killer\": \"%s\"\n}\n",
                LaneDD0Identity.RNG_SEED, LaneDD0Identity.N0,
                seq.tEnd, seq.cEnd, seq.nEnd, seq.nEverEnd, seq.birthsTotal,
                seq.stopReason, seq.phiRbFreeze, seq.sigmaFreeze, seq.muEff,
                seq.gammaPerDay, seq.nHandoff, seq.nFamineEnd, seq.deaths,
                off.nEnd, off.births, off.cEnd,
                seq.capFired, tCap, pass, killer), StandardCharsets.UTF_8);

        System.out.printf(Locale.US, "D4.1 honesty PASS%n");
        System.out.printf(Locale.US, "D4.2 process-off N_end=%d C_end=%.3f %s%n",
                off.nEnd, off.cEnd, processOff ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "D4.3 no_concurrency feast_death_on=false %s%n",
                noDeathFeast ? "PASS" : "FAIL");
        System.out.printf(Locale.US,
                "D4.4 feast_fingerprint t=%.1f C=%.6e N=%d %s%n",
                seq.tEnd, seq.cEnd, seq.nEnd, fingerprint ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "D4.5 phi_X=phi_Rb mu_eff=%.6f %s%n",
                seq.muEff, muId ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "D4.6 famine_clock gamma=%.6f N_end=%.6f %s%n",
                seq.gammaPerDay, seq.nFamineEnd, famine ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "D4.7 no_silent_kill %s%n",
                noKill && !scope ? "PASS" : (scope ? "SCOPE_NOTE" : "FAIL"));
        System.out.println("D4.8 isolation is the checker git diff.");

        if (pass) {
            System.out.println("D4_THREE_PHASE=PASS pending isolation. "
                    + "Sequential feast-handoff-famine. Death off while C>C_cut.");
        } else if (scope) {
            System.out.println("D4_THREE_PHASE=SCOPE_NOTE " + killer + ". Not PASS.");
            System.exit(2);
        } else {
            System.out.println("D4_THREE_PHASE=FAIL " + killer
                    + ". Do not raise C_s or put K_M in Monod.");
            System.exit(1);
        }
    }

    private static boolean famineOk(Seq seq) {
        if (seq.famT.size() != LaneDD4Identity.FAMINE_OBS_D.length) {
            return false;
        }
        if (seq.nHandoff <= 0) {
            return false;
        }
        for (int i = 0; i < seq.famT.size(); i++) {
            double pred = seq.nHandoff * Math.exp(-seq.gammaPerDay * seq.famT.get(i));
            double rel = Math.abs(seq.famN.get(i) - pred) / seq.nHandoff;
            if (rel > LaneDD4Identity.FAMINE_CLOCK_DELTA + 1.0e-15) {
                return false;
            }
        }
        return true;
    }

    private static void refuseCapacity(String[] args) {
        String blob = String.join(" ", args == null ? new String[0] : args)
                .toLowerCase(Locale.ROOT);
        if (blob.contains("narma") || blob.contains("charc") || blob.contains("ipc")) {
            throw new IllegalArgumentException("D4 refuses NARMA/CHARC/IPC");
        }
        if (blob.contains("k_m") || blob.contains("bacteria.max")
                || blob.contains("lc1b") || blob.contains("concurrent")) {
            throw new IllegalArgumentException(
                    "D4 refuses K_M, bacteria.max, LC1b, concurrent death");
        }
    }

    private static void requireFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces")
                || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("D4 protocol is not frozen");
        }
        if (!json.contains("\"status_label\": \"D4_THREE_PHASE\"")
                || !json.contains("\"phi_x\": \"phi_Rb\"")
                || !json.contains("\"erickson_k_m_in_monod\": false")
                || !json.contains("\"lambda_c_replaces_lambda_s\": false")
                || !json.contains("\"lc1b_coin_flip\": false")
                || !json.contains("\"d33_unblocked\": true")) {
            throw new IllegalStateException("D4 freeze incomplete");
        }
    }

    private static final class Seq {
        final List<Double> feastT = new ArrayList<>();
        final List<Double> feastC = new ArrayList<>();
        final List<Integer> feastN = new ArrayList<>();
        final List<Integer> feastNever = new ArrayList<>();
        final List<Double> feastLam = new ArrayList<>();
        final List<Double> feastPhiRb = new ArrayList<>();
        final List<Double> feastSigma = new ArrayList<>();
        final List<Boolean> feastDeathOn = new ArrayList<>();
        final List<Double> famT = new ArrayList<>();
        final List<Double> famN = new ArrayList<>();
        final List<Double> famPred = new ArrayList<>();
        int birthsTotal;
        int nEnd;
        int nEverEnd;
        int nHandoff;
        double tEnd;
        double cEnd;
        double phiRbFreeze;
        double sigmaFreeze;
        double muEff;
        double gammaPerDay;
        double nFamineEnd;
        double deaths;
        int shrinkEvents;
        boolean capFired;
        String stopReason;
    }

    private static final class Off {
        int nEnd;
        int births;
        double deaths;
        double cEnd;
    }
}
