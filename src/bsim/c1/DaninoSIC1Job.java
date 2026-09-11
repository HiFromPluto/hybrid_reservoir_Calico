package bsim.c1;

import bsim.BSimRandom;
import bsim.circuit.DaninoSIPeriodCheck;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

/**
 * Gate C1: packed spatial Object B. NOT_FIG4B. Not Fig. 4b. Not Object A.
 * Does not start A0/W0 on failure. Does not retune Object B / D1_spatial / d.
 */
public final class DaninoSIC1Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path CLAIM_FREEZE = ROOT.resolve("examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md");
    private static final Path PROTOCOL = ROOT.resolve("examples/BSimReservoirPlanDaninoPocketC1/PROTOCOL.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/BSimReservoirPlanDaninoPocketC1/configs/protocol.json");
    private static final Path RESULTS = ROOT.resolve("examples/BSimReservoirPlanDaninoPocketC1/results");
    private static final Path P0_CENTRES = ROOT.resolve(
            "examples/BSimReservoirPlanDaninoPocketP0/results/p0_seed101/p0_end_centres.csv");

    private static final double RK_DT = 0.001;
    private static final double RK_DT_LORES = 0.002;
    private static final double SAMPLE_DT = 0.5;
    private static final double T_END = 1000.0;
    private static final double PY_T_040 = 63.583333333333336;
    private static final double PERIOD_REL_BAND = 0.20;
    private static final double LEDGER_REL = 1e-3;
    private static final double SYNC_FRAC_OSC = 0.90;
    private static final double SYNC_WINDOW = 5.0;
    private static final double SYNC_FRAC_WINDOW = 0.90;
    private static final double HE_RANGE_MIN = 1e-6;

    private DaninoSIC1Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        DaninoSIPeriodCheck.refuseNarma(args);
        requireClaimFreeze();
        requireProtocolFrozen();
        Files.createDirectories(RESULTS);

        BSimRandom rng = new BSimRandom(0L);
        DaninoSIC1PackedSpatial device = new DaninoSIC1PackedSpatial(rng);
        if (device.unusedRng() != rng) {
            throw new IllegalStateException("RNG injection contract broken NOT_FIG4B");
        }

        System.out.println("C1 packed spatial Object B. DaninoSI_OccupiedDF. NOT_FIG4B.");
        System.out.println("Object A Danino2010_Fig4b_bulk_twin is not ported. Not Fig. 4b.");
        System.out.println("D0 remains FAIL. D0b remains FAIL_NO_IDENTITY.");
        System.out.println("D1 remains PASS on circuit only. P0 remains PASS on mechanics only.");
        System.out.println("C0 remains PASS on well-mixed Object B only.");
        System.out.println("Mechanics OFF. SI-scaled clock. TIME_ADJ unused. NOT_FIG4B.");
        System.out.printf(
                Locale.US,
                "held-pack NOT_FIG4B N=%d v_cell=%.6g Ve_local=%.6g d=%.16g D1_spatial=%.6g %s mask=%s%n",
                device.nCells(),
                device.vCell(),
                device.vELocal(),
                device.d(),
                device.d1Spatial(),
                DaninoSIC1PackedSpatial.D1_LABEL,
                device.maskSha256());

        DaninoSIPeriodCheck.Protocol peaks = DaninoSIPeriodCheck.Protocol.d0b();
        List<ArmResult> arms = new ArrayList<>();
        String want = args.length > 0 ? args[0].trim() : "";

        if (wantEmptyOr(want, "C1_HELD_PACK_D05_MU040")) {
            arms.add(runIdentity(device, peaks, "C1_HELD_PACK_D05_MU040",
                    "AHL_KICK_005", 0.40, RK_DT, true, true));
        }
        if (wantEmptyOr(want, "C1_HELD_PACK_WRONG_KICK")) {
            arms.add(runControl(device, peaks, "C1_HELD_PACK_WRONG_KICK",
                    "SI_BASAL_PERTURB", 0.40, "NO_PERIOD"));
        }
        if (wantEmptyOr(want, "C1_HELD_PACK_MU150")) {
            arms.add(runControl(device, peaks, "C1_HELD_PACK_MU150",
                    "AHL_KICK_005", 1.50, "NO_PERIOD"));
        }
        if (wantEmptyOr(want, "C1_HELD_PACK_MU0_LEDGER")) {
            arms.add(runLedger(device, peaks));
        }
        if (wantEmptyOr(want, "C1_HELD_PACK_MU040_LORES")) {
            arms.add(runIdentity(device, peaks, "C1_HELD_PACK_MU040_LORES",
                    "AHL_KICK_005", 0.40, RK_DT_LORES, true, false));
        }
        if (wantEmptyOr(want, "C1_P0_SNAPSHOT_SMOKE")) {
            arms.add(runP0SnapshotSmoke());
        }

        if (arms.isEmpty()) {
            throw new IllegalArgumentException("no C1 arms selected: " + want);
        }

        boolean pass = true;
        for (ArmResult a : arms) {
            pass = pass && a.pass;
        }
        if (want.isEmpty()) {
            pass = pass && loresOscUnchanged(arms);
        }

        writeSummary(device, arms, pass, want.isEmpty());
        System.out.println(pass
                ? "C1=PASS NOT_FIG4B. A0 may start (ideal AC source, still NOT_FIG4B, still not Fig. 4b). W0 later."
                : "C1=FAIL NOT_FIG4B. Do not retune Object B / D1_spatial / d. Do not start A0/W0. Not Fig. 4b.");
        if (!pass) {
            System.exit(1);
        }
    }

    private static boolean wantEmptyOr(String want, String name) {
        return want.isEmpty() || want.equals(name);
    }

    private static void requireClaimFreeze() throws IOException {
        if (!Files.exists(CLAIM_FREEZE)) {
            throw new IllegalStateException("CLAIM_FREEZE_DANINO_SI.md must exist before C1.");
        }
        String freeze = Files.readString(CLAIM_FREEZE, StandardCharsets.UTF_8);
        if (!freeze.contains("DaninoSI_OccupiedDF") || !freeze.contains("NOT_FIG4B")) {
            throw new IllegalStateException("claim freeze must name Object B as NOT_FIG4B");
        }
        if (!freeze.contains("Danino2010_Fig4b_bulk_twin")) {
            throw new IllegalStateException("claim freeze must name Object A");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        if (!Files.exists(PROTOCOL) || !Files.exists(PROTOCOL_JSON)) {
            throw new IllegalStateException("C1 PROTOCOL must be frozen before traces");
        }
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("C1 protocol is not frozen_before_traces");
        }
        if (!text.contains("NOT_FIG4B") || !json.contains("NOT_FIG4B")) {
            throw new IllegalStateException("C1 protocol must say NOT_FIG4B");
        }
        if (!text.contains("TIME_ADJ") || json.contains("\"TIME_ADJ\": true")) {
            throw new IllegalStateException("C1 protocol must freeze TIME_ADJ unused");
        }
        if (!text.contains("C1_HELD_PACK_D05_MU040") || !text.contains("SYNC_PEAK_WINDOW")) {
            throw new IllegalStateException("C1 protocol must freeze identity arm and sync metric");
        }
    }

    private static ArmResult runIdentity(
            DaninoSIC1PackedSpatial device,
            DaninoSIPeriodCheck.Protocol peaks,
            String name,
            String ensemble,
            double mu,
            double rkDt,
            boolean requireOsc,
            boolean requireSync) throws IOException {
        long t0 = System.nanoTime();
        System.out.printf("  %s NOT_FIG4B starting rk_dt=%.4g mu=%.2f%n", name, rkDt, mu);
        DaninoSIC1PackedSpatial.Trajectory traj = device.integrate(mu, ensemble, T_END, SAMPLE_DT, rkDt);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        if (traj.negativeState) {
            per = new DaninoSIPeriodCheck.Result(
                    "NEGATIVE_STATE", per.nPeaks, per.nPersistPeaks, Double.NaN,
                    per.peakTimes, per.relAmplitude, per.ampPersist, per.meanI, per.p5I, per.p95I);
        }
        SyncStats sync = syncStats(traj, peaks);
        double relT = Double.NaN;
        boolean periodOk = !requireOsc;
        boolean flagOk = requireOsc ? "OSC".equals(per.flag) : true;
        if (requireOsc && "OSC".equals(per.flag) && !Double.isNaN(per.period)) {
            relT = Math.abs(per.period - PY_T_040) / PY_T_040;
            periodOk = relT <= PERIOD_REL_BAND;
        }
        boolean heOk = traj.heFluidRangeMax >= HE_RANGE_MIN;
        boolean syncOk = !requireSync || (sync.fracOsc >= SYNC_FRAC_OSC
                && sync.fracInWindow >= SYNC_FRAC_WINDOW);
        boolean armPass = flagOk && periodOk && heOk && syncOk && !traj.negativeState;
        writeMeanCsv(RESULTS.resolve("java_" + name + ".csv"), traj);
        if (requireSync) {
            writeSyncCsv(RESULTS.resolve("sync_" + name + ".csv"), sync);
        }
        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  %s NOT_FIG4B mu=%.2f flag=%s T=%s relT=%s sync_osc=%.3f sync_win=%.3f "
                        + "He_fluid_range=%.3g He_patch_range=%.3g residual_rel=%.3g wall_s=%.1f gate=%s%n",
                name, mu, per.flag, fmt(per.period), fmt(relT),
                sync.fracOsc, sync.fracInWindow,
                traj.heFluidRangeMax, traj.hePatchRangeMax,
                traj.relativeResidual(), sec,
                armPass ? "PASS" : "FAIL");
        return new ArmResult(
                name, ensemble, mu, rkDt, per.flag, per.period, relT,
                sync.fracOsc, sync.fracInWindow,
                traj.heFluidRangeMax, traj.hePatchRangeMax,
                traj.relativeResidual(), traj.ledger.membraneCancel(),
                traj.extra.getDecayLoss(), traj.ledger.gammaHLoss(), traj.ledger.synthSource(),
                traj.extra.getOutletLoss(), traj.extra.getBoundaryLoss(),
                armPass, requireOsc, requireSync,
                "packed spatial Object B NOT_FIG4B");
    }

    private static ArmResult runControl(
            DaninoSIC1PackedSpatial device,
            DaninoSIPeriodCheck.Protocol peaks,
            String name,
            String ensemble,
            double mu,
            String requiredFlag) throws IOException {
        long t0 = System.nanoTime();
        System.out.printf("  %s NOT_FIG4B starting mu=%.2f ensemble=%s%n", name, mu, ensemble);
        DaninoSIC1PackedSpatial.Trajectory traj = device.integrate(mu, ensemble, T_END, SAMPLE_DT, RK_DT);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        boolean armPass = requiredFlag.equals(per.flag);
        writeMeanCsv(RESULTS.resolve("java_" + name + ".csv"), traj);
        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  %s NOT_FIG4B mu=%.2f flag=%s required=%s residual_rel=%.3g wall_s=%.1f gate=%s%n",
                name, mu, per.flag, requiredFlag, traj.relativeResidual(), sec,
                armPass ? "PASS" : "FAIL");
        return new ArmResult(
                name, ensemble, mu, RK_DT, per.flag, per.period, Double.NaN,
                Double.NaN, Double.NaN,
                traj.heFluidRangeMax, traj.hePatchRangeMax,
                traj.relativeResidual(), traj.ledger.membraneCancel(),
                traj.extra.getDecayLoss(), traj.ledger.gammaHLoss(), traj.ledger.synthSource(),
                traj.extra.getOutletLoss(), traj.extra.getBoundaryLoss(),
                armPass, false, false,
                "control NOT_FIG4B");
    }

    private static ArmResult runLedger(
            DaninoSIC1PackedSpatial device,
            DaninoSIPeriodCheck.Protocol peaks) throws IOException {
        long t0 = System.nanoTime();
        System.out.println("  C1_HELD_PACK_MU0_LEDGER NOT_FIG4B starting mu=0");
        DaninoSIC1PackedSpatial.Trajectory traj = device.integrate(
                0.0, "AHL_KICK_005", T_END, SAMPLE_DT, RK_DT);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        boolean ledgerOk = traj.relativeResidual() <= LEDGER_REL
                && Math.abs(traj.extra.getDecayLoss()) <= 1e-18
                && Math.abs(traj.ledger.membraneCancel()) <= 1e-12;
        writeMeanCsv(RESULTS.resolve("java_C1_HELD_PACK_MU0_LEDGER.csv"), traj);
        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  C1_HELD_PACK_MU0_LEDGER NOT_FIG4B mu=0 flag=%s residual_rel=%.3g mu_loss=%.3g "
                        + "gammaH_loss=%.6g synth=%.6g cancel=%.3g outlet=%.3g boundary=%.3g wall_s=%.1f gate=%s%n",
                per.flag,
                traj.relativeResidual(),
                traj.extra.getDecayLoss(),
                traj.ledger.gammaHLoss(),
                traj.ledger.synthSource(),
                traj.ledger.membraneCancel(),
                traj.extra.getOutletLoss(),
                traj.extra.getBoundaryLoss(),
                sec,
                ledgerOk ? "PASS" : "FAIL");
        return new ArmResult(
                "C1_HELD_PACK_MU0_LEDGER", "AHL_KICK_005", 0.0, RK_DT, per.flag, per.period, Double.NaN,
                Double.NaN, Double.NaN,
                traj.heFluidRangeMax, traj.hePatchRangeMax,
                traj.relativeResidual(), traj.ledger.membraneCancel(),
                traj.extra.getDecayLoss(), traj.ledger.gammaHLoss(), traj.ledger.synthSource(),
                traj.extra.getOutletLoss(), traj.extra.getBoundaryLoss(),
                ledgerOk, false, false,
                "mu=0 ledger split gammaH vs remaining NOT_FIG4B");
    }

    private static ArmResult runP0SnapshotSmoke() {
        boolean present = Files.exists(P0_CENTRES);
        System.out.printf(
                "  C1_P0_SNAPSHOT_SMOKE NOT_FIG4B not identity skipped=%s (P0 end centres %s)%n",
                !present,
                present ? "present" : "absent");
        if (!present) {
            return new ArmResult(
                    "C1_P0_SNAPSHOT_SMOKE", "SKIPPED", Double.NaN, Double.NaN, "SKIPPED",
                    Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                    Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                    true, false, false,
                    "optional smoke; P0 did not write end centres; not identity NOT_FIG4B");
        }
        return new ArmResult(
                "C1_P0_SNAPSHOT_SMOKE", "P0_END", Double.NaN, Double.NaN, "NOT_RUN",
                Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                true, false, false,
                "optional smoke present but not required for PASS NOT_FIG4B");
    }

    private static boolean loresOscUnchanged(List<ArmResult> arms) {
        ArmResult id = find(arms, "C1_HELD_PACK_D05_MU040");
        ArmResult lo = find(arms, "C1_HELD_PACK_MU040_LORES");
        if (id == null || lo == null) {
            return false;
        }
        boolean ok = "OSC".equals(id.flag) && "OSC".equals(lo.flag);
        System.out.printf(
                "  LORES OSC-call unchanged NOT_FIG4B identity=%s lores=%s gate=%s%n",
                id.flag, lo.flag, ok ? "PASS" : "FAIL");
        return ok;
    }

    private static ArmResult find(List<ArmResult> arms, String name) {
        for (ArmResult a : arms) {
            if (name.equals(a.name)) {
                return a;
            }
        }
        return null;
    }

    private static SyncStats syncStats(
            DaninoSIC1PackedSpatial.Trajectory traj,
            DaninoSIPeriodCheck.Protocol peaks) {
        int n = traj.I[0].length;
        int nOsc = 0;
        double[] lastPeak = new double[n];
        String[] flags = new String[n];
        double[] periods = new double[n];
        Arrays.fill(lastPeak, Double.NaN);
        for (int c = 0; c < n; c++) {
            double[] ic = new double[traj.t.length];
            for (int s = 0; s < traj.t.length; s++) {
                ic[s] = traj.I[s][c];
            }
            DaninoSIPeriodCheck.Result r = DaninoSIPeriodCheck.extract(traj.t, ic, peaks);
            flags[c] = r.flag;
            periods[c] = r.period;
            if ("OSC".equals(r.flag) && !r.peakTimes.isEmpty()) {
                nOsc++;
                lastPeak[c] = r.peakTimes.get(r.peakTimes.size() - 1);
            }
        }
        double fracOsc = nOsc / (double) n;
        double[] oscPeaks = new double[nOsc];
        int k = 0;
        for (int c = 0; c < n; c++) {
            if (!Double.isNaN(lastPeak[c])) {
                oscPeaks[k++] = lastPeak[c];
            }
        }
        Arrays.sort(oscPeaks);
        double median = nOsc == 0 ? Double.NaN : oscPeaks[nOsc / 2];
        int nWin = 0;
        for (double lp : oscPeaks) {
            if (Math.abs(lp - median) <= SYNC_WINDOW) {
                nWin++;
            }
        }
        double fracWin = nOsc == 0 ? 0.0 : nWin / (double) nOsc;
        return new SyncStats(fracOsc, fracWin, nOsc, median, flags, periods, lastPeak);
    }

    private static void writeMeanCsv(Path path, DaninoSIC1PackedSpatial.Trajectory traj) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.println("t,A_mean,I_mean,Hi_mean,He_patch_mean,He_patch_min,He_patch_max,"
                    + "He_fluid_min,He_fluid_max,NOT_FIG4B");
            for (int i = 0; i < traj.t.length; i++) {
                w.printf(Locale.US, "%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,NOT_FIG4B%n",
                        traj.t[i], traj.Y[i][0], traj.Y[i][1], traj.Y[i][2],
                        traj.hePatchMean[i], traj.hePatchMin[i], traj.hePatchMax[i],
                        traj.heFluidMin[i], traj.heFluidMax[i]);
            }
        }
    }

    private static void writeSyncCsv(Path path, SyncStats sync) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.println("cell,flag,period,last_peak_t,NOT_FIG4B");
            for (int c = 0; c < sync.flags.length; c++) {
                w.printf(Locale.US, "%d,%s,%s,%s,NOT_FIG4B%n",
                        c, sync.flags[c], fmt(sync.periods[c]), fmt(sync.lastPeak[c]));
            }
        }
    }

    private static void writeSummary(
            DaninoSIC1PackedSpatial device,
            List<ArmResult> arms,
            boolean pass,
            boolean fullSuite) throws IOException {
        StringBuilder json = new StringBuilder();
        json.append("{\n");
        json.append("  \"object\": \"DaninoSI_OccupiedDF\",\n");
        json.append("  \"object_status\": \"OCCUPIED_NOT_FIG4B\",\n");
        json.append("  \"literature_twin\": \"Danino2010_Fig4b_bulk_twin\",\n");
        json.append("  \"literature_twin_status\": \"FAIL_NO_IDENTITY\",\n");
        json.append("  \"D0\": \"FAIL\",\n");
        json.append("  \"D0b\": \"FAIL_NO_IDENTITY\",\n");
        json.append("  \"D1\": \"PASS\",\n");
        json.append("  \"P0\": \"PASS\",\n");
        json.append("  \"C0\": \"PASS\",\n");
        json.append("  \"C1\": \"").append(pass ? "PASS" : "FAIL").append("\",\n");
        json.append("  \"NOT_FIG4B\": true,\n");
        json.append("  \"a0_may_start\": ").append(pass && fullSuite).append(",\n");
        json.append("  \"w0_may_start\": false,\n");
        json.append("  \"n_cells\": ").append(device.nCells()).append(",\n");
        json.append("  \"v_cell_um3\": ").append(device.vCell()).append(",\n");
        json.append("  \"Ve_local_um3\": ").append(device.vELocal()).append(",\n");
        json.append("  \"d\": ").append(device.d()).append(",\n");
        json.append("  \"D1_spatial\": ").append(device.d1Spatial()).append(",\n");
        json.append("  \"D1_spatial_label\": \"").append(DaninoSIC1PackedSpatial.D1_LABEL).append("\",\n");
        json.append("  \"mask_sha256\": \"").append(device.maskSha256()).append("\",\n");
        json.append("  \"mechanics\": \"OFF\",\n");
        json.append("  \"TIME_ADJ\": false,\n");
        json.append("  \"clock\": \"SI_SCALED\",\n");
        json.append("  \"p0_box_used_as_Ve\": false,\n");
        json.append("  \"sync_metric\": \"SYNC_PEAK_WINDOW\",\n");
        json.append("  \"arms\": [\n");
        for (int i = 0; i < arms.size(); i++) {
            ArmResult a = arms.get(i);
            json.append("    {");
            json.append("\"name\":\"").append(a.name).append("\",");
            json.append("\"ensemble\":\"").append(a.ensemble).append("\",");
            json.append("\"mu\":").append(num(a.mu)).append(",");
            json.append("\"rk_dt\":").append(num(a.rkDt)).append(",");
            json.append("\"flag\":\"").append(a.flag).append("\",");
            json.append("\"period\":").append(num(a.period)).append(",");
            json.append("\"period_rel_err\":").append(num(a.relPeriod)).append(",");
            json.append("\"sync_frac_osc\":").append(num(a.fracOsc)).append(",");
            json.append("\"sync_frac_window\":").append(num(a.fracWindow)).append(",");
            json.append("\"he_fluid_range\":").append(num(a.heFluidRange)).append(",");
            json.append("\"he_patch_range\":").append(num(a.hePatchRange)).append(",");
            json.append("\"ledger_rel\":").append(num(a.ledgerRel)).append(",");
            json.append("\"membrane_cancel\":").append(num(a.memCancel)).append(",");
            json.append("\"mu_loss\":").append(num(a.muLoss)).append(",");
            json.append("\"gammaH_loss\":").append(num(a.gammaH)).append(",");
            json.append("\"synth\":").append(num(a.synth)).append(",");
            json.append("\"outlet\":").append(num(a.outlet)).append(",");
            json.append("\"boundary\":").append(num(a.boundary)).append(",");
            json.append("\"pass\":").append(a.pass).append(",");
            json.append("\"NOT_FIG4B\": true");
            json.append("}");
            json.append(i + 1 < arms.size() ? ",\n" : "\n");
        }
        json.append("  ]\n}\n");
        Files.writeString(RESULTS.resolve("c1_summary.json"), json.toString(), StandardCharsets.UTF_8);

        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(
                RESULTS.resolve("c1_arms.csv"), StandardCharsets.UTF_8))) {
            w.println("name,ensemble,mu,rk_dt,flag,period,period_rel_err,sync_frac_osc,sync_frac_window,"
                    + "he_fluid_range,he_patch_range,ledger_rel,pass,NOT_FIG4B,note");
            for (ArmResult a : arms) {
                w.printf(Locale.US, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,NOT_FIG4B,%s%n",
                        a.name, a.ensemble, fmt(a.mu), fmt(a.rkDt), a.flag, fmt(a.period), fmt(a.relPeriod),
                        fmt(a.fracOsc), fmt(a.fracWindow), fmt(a.heFluidRange), fmt(a.hePatchRange),
                        fmt(a.ledgerRel), a.pass, a.note.replace(',', ';'));
            }
        }
    }

    private static String fmt(double x) {
        if (Double.isNaN(x)) {
            return "nan";
        }
        return String.format(Locale.US, "%.6g", x);
    }

    private static String num(double x) {
        if (Double.isNaN(x)) {
            return "null";
        }
        return Double.toString(x);
    }

    private static final class SyncStats {
        final double fracOsc;
        final double fracInWindow;
        final int nOsc;
        final double medianLastPeak;
        final String[] flags;
        final double[] periods;
        final double[] lastPeak;

        SyncStats(
                double fracOsc,
                double fracInWindow,
                int nOsc,
                double medianLastPeak,
                String[] flags,
                double[] periods,
                double[] lastPeak) {
            this.fracOsc = fracOsc;
            this.fracInWindow = fracInWindow;
            this.nOsc = nOsc;
            this.medianLastPeak = medianLastPeak;
            this.flags = flags;
            this.periods = periods;
            this.lastPeak = lastPeak;
        }
    }

    private static final class ArmResult {
        final String name;
        final String ensemble;
        final double mu;
        final double rkDt;
        final String flag;
        final double period;
        final double relPeriod;
        final double fracOsc;
        final double fracWindow;
        final double heFluidRange;
        final double hePatchRange;
        final double ledgerRel;
        final double memCancel;
        final double muLoss;
        final double gammaH;
        final double synth;
        final double outlet;
        final double boundary;
        final boolean pass;
        final boolean requireOsc;
        final boolean requireSync;
        final String note;

        ArmResult(
                String name,
                String ensemble,
                double mu,
                double rkDt,
                String flag,
                double period,
                double relPeriod,
                double fracOsc,
                double fracWindow,
                double heFluidRange,
                double hePatchRange,
                double ledgerRel,
                double memCancel,
                double muLoss,
                double gammaH,
                double synth,
                double outlet,
                double boundary,
                boolean pass,
                boolean requireOsc,
                boolean requireSync,
                String note) {
            this.name = name;
            this.ensemble = ensemble;
            this.mu = mu;
            this.rkDt = rkDt;
            this.flag = flag;
            this.period = period;
            this.relPeriod = relPeriod;
            this.fracOsc = fracOsc;
            this.fracWindow = fracWindow;
            this.heFluidRange = heFluidRange;
            this.hePatchRange = hePatchRange;
            this.ledgerRel = ledgerRel;
            this.memCancel = memCancel;
            this.muLoss = muLoss;
            this.gammaH = gammaH;
            this.synth = synth;
            this.outlet = outlet;
            this.boundary = boundary;
            this.pass = pass;
            this.requireOsc = requireOsc;
            this.requireSync = requireSync;
            this.note = note;
        }
    }
}
