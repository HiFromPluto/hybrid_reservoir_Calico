package bsim.c0;

import bsim.BSimRandom;
import bsim.circuit.DaninoSIOccupiedDF;
import bsim.circuit.DaninoSIPeriodCheck;
import bsim.d1.NpzReader;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Gate C0: well-mixed Object B coupling. NOT_FIG4B. Not Fig. 4b.
 * Does not port Object A. Does not start C1 on failure.
 */
public final class DaninoSIC0Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path CLAIM_FREEZE = ROOT.resolve("examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md");
    private static final Path PROTOCOL = ROOT.resolve("examples/BSimReservoirPlanDaninoPocketC0/PROTOCOL.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/BSimReservoirPlanDaninoPocketC0/configs/protocol.json");
    private static final Path RESULTS = ROOT.resolve("examples/BSimReservoirPlanDaninoPocketC0/results");
    private static final Path D1_FIXTURES = ROOT.resolve("examples/BSimReservoirPlanPocketOscSI_D1/fixtures");
    private static final Path D0B_FIXTURES = ROOT.resolve(
            "examples/BSimReservoirPlanPocketOscSI_D0b/results/fixtures");

    private static final double RK_DT = 0.001;
    private static final double SAMPLE_DT = 0.5;
    private static final double T_END = 1000.0;
    private static final double T_END_CLOSED = 50.0;
    private static final double TRAJ_RTOL = 0.02;
    private static final double TRAJ_ATOL = 1e-4;
    private static final double FIRST_RTOL = 1e-3;
    private static final double FIRST_ATOL = 1e-6;
    private static final double PERIOD_REL = 0.02;
    private static final double LEDGER_REL = 1e-3;
    private static final double IDENTICAL_SPREAD = 1e-10;
    private static final double HETERO_HE_MOVE = 1e-6;
    private static final double PY_T_040 = 63.583333333333336;
    private static final double PY_T_032 = 56.30769230769231;

    private DaninoSIC0Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        DaninoSIPeriodCheck.refuseNarma(args);
        requireClaimFreeze();
        requireProtocolFrozen();
        Files.createDirectories(RESULTS);

        BSimRandom rng = new BSimRandom(0L);
        DaninoSIC0WellMixed device = new DaninoSIC0WellMixed(rng);
        if (device.unusedRng() != rng) {
            throw new IllegalStateException("RNG injection contract broken NOT_FIG4B");
        }

        System.out.println("C0 well-mixed Object B coupling. DaninoSI_OccupiedDF. NOT_FIG4B.");
        System.out.println("Object A Danino2010_Fig4b_bulk_twin is not ported.");
        System.out.println("D0 remains FAIL. D0b remains FAIL_NO_IDENTITY. D1 remains PASS on circuit only.");
        System.out.println("P0 remains PASS on mechanics only. Mechanics OFF. Not Fig. 4b.");
        System.out.printf(
                Locale.US,
                "volumes NOT_FIG4B N=%d v_cell=%.6g V_e=%.6g d=%.16g d/(1-d)=%.16g%n",
                device.nCells(),
                device.vCell(),
                device.vE(),
                device.d(),
                device.dOverOneMinusD());

        algebraCheck(device);

        DaninoSIPeriodCheck.Protocol peaks = DaninoSIPeriodCheck.Protocol.d0b();
        List<ArmResult> arms = new ArrayList<>();

        arms.add(runIdentity(
                device, peaks, "C0_IDENTICAL_D05_MU040", "AHL_KICK_005", 0.40,
                "OSC", PY_T_040, "primary_mu_0p40", true));
        arms.add(runIdentity(
                device, peaks, "C0_IDENTICAL_D05_MU032", "AHL_KICK_005", 0.32,
                "OSC", PY_T_032, "primary_osc", false));
        arms.add(runIdentity(
                device, peaks, "C0_IDENTICAL_D05_MU150", "AHL_KICK_005", 1.50,
                "NO_PERIOD", Double.NaN, null, false));
        arms.add(runIdentity(
                device, peaks, "C0_WRONG_KICK_BASAL", "SI_BASAL_PERTURB", 0.40,
                "NO_PERIOD", Double.NaN, "si_basal_perturb_mu_0p40", false));
        arms.add(runLedger(rng, peaks));
        arms.add(runClosedMembrane(rng));
        arms.add(runHetero(rng, peaks));

        boolean pass = true;
        for (ArmResult a : arms) {
            pass = pass && a.pass;
        }

        writeSummary(device, arms, pass);
        System.out.println(pass
                ? "C0=PASS NOT_FIG4B. C1 may start (packed spatial Object B, still NOT_FIG4B)."
                : "C0=FAIL NOT_FIG4B. Do not retune Object B. Do not start C1.");
        if (!pass) {
            System.exit(1);
        }
    }

    private static void requireClaimFreeze() throws IOException {
        if (!Files.exists(CLAIM_FREEZE)) {
            throw new IllegalStateException("CLAIM_FREEZE_DANINO_SI.md must exist before C0.");
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
            throw new IllegalStateException("C0 PROTOCOL must be frozen before traces");
        }
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("C0 protocol is not frozen_before_traces");
        }
        if (!text.contains("NOT_FIG4B") || !json.contains("NOT_FIG4B")) {
            throw new IllegalStateException("C0 protocol must say NOT_FIG4B");
        }
    }

    private static void algebraCheck(DaninoSIC0WellMixed device) {
        double[] y = {0.1, 0.2, 0.3, 0.4};
        device.assertReducesToObjectB(y, 0.25, 0.40);
        double[] kick = {0.0, 0.0, 0.05, 0.05};
        device.assertReducesToObjectB(kick, 0.05, 0.40);
        System.out.println("  algebra_reduction NOT_FIG4B d=0.5 identical He RHS matches Object B gate=PASS");
    }

    private static ArmResult runIdentity(
            DaninoSIC0WellMixed device,
            DaninoSIPeriodCheck.Protocol peaks,
            String name,
            String ensemble,
            double mu,
            String requiredFlag,
            double pythonPeriod,
            String fixtureName,
            boolean compareTraj) throws IOException {
        double[] y0 = ensemble.equals("SI_BASAL_PERTURB")
                ? new double[] {0.0, 1.0, 0.0, 0.0}
                : new double[] {0.0, 0.0, 0.05, 0.05};
        double hist = y0[2];
        double[][] cells = identicalCells(device.nCells(), y0);
        double[] hiHist = fill(device.nCells(), hist);

        DaninoSIC0WellMixed.Trajectory traj = device.integrate(
                mu, cells, hiHist, y0[3], T_END, SAMPLE_DT, RK_DT);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        if (traj.negativeState) {
            per = new DaninoSIPeriodCheck.Result(
                    "NEGATIVE_STATE", per.nPeaks, per.nPersistPeaks, Double.NaN,
                    per.peakTimes, per.relAmplitude, per.ampPersist, per.meanI, per.p5I, per.p95I);
        }

        TrajStats full = TrajStats.skipped();
        TrajStats first = TrajStats.skipped();
        if (compareTraj || fixtureName != null) {
            Path npz = resolveFixture(fixtureName);
            if (npz != null) {
                Map<String, NpzReader.Array> arrays = NpzReader.loadNpz(npz);
                double[] tPy = arrays.get("t").vector();
                double[][] yPy = arrays.get("Y").matrix2();
                full = compare(tPy, yPy, traj, TRAJ_ATOL, TRAJ_RTOL, 0.0, T_END);
                first = compare(tPy, yPy, traj, FIRST_ATOL, FIRST_RTOL, 0.0, DaninoSIOccupiedDF.TAU);
            }
        }

        boolean flagOk = requiredFlag.equals(per.flag);
        boolean periodOk = true;
        double relT = Double.NaN;
        if ("OSC".equals(requiredFlag)) {
            if (Double.isNaN(per.period) || Double.isNaN(pythonPeriod) || pythonPeriod == 0.0) {
                periodOk = false;
            } else {
                relT = Math.abs(per.period - pythonPeriod) / pythonPeriod;
                periodOk = relT < PERIOD_REL;
            }
        }
        boolean spreadOk = traj.maxHiSpread <= IDENTICAL_SPREAD;
        boolean trajOk = !compareTraj || (full.allInside && first.allInside);
        boolean armPass = trajOk && flagOk && periodOk && spreadOk;

        writeMeanCsv(RESULTS.resolve("java_" + name + ".csv"), traj);
        System.out.printf(
                "  %s NOT_FIG4B mu=%.2f flag=%s T=%s relT=%s traj_max=%s first_max=%s spread=%.3g residual_rel=%.3g gate=%s%n",
                name,
                mu,
                per.flag,
                fmt(per.period),
                fmt(relT),
                compareTraj ? String.format(Locale.US, "%.3g", full.maxNorm) : "n/a",
                compareTraj ? String.format(Locale.US, "%.3g", first.maxNorm) : "n/a",
                traj.maxHiSpread,
                traj.relativeResidual(),
                armPass ? "PASS" : "FAIL");
        if (compareTraj && !full.allInside) {
            System.out.printf("    traj violations=%d / %d worst=%s%n", full.violations, full.n, full.worst);
        }
        if (!spreadOk) {
            System.out.printf("    identical Hi spread %.3g exceeds %.3g NOT_FIG4B%n",
                    traj.maxHiSpread, IDENTICAL_SPREAD);
        }
        return new ArmResult(
                name, ensemble, mu, per.flag, per.period, relT,
                trajOk, first.allInside || !compareTraj, full.maxNorm, first.maxNorm,
                traj.relativeResidual(), traj.maxHiSpread, Double.NaN, armPass,
                "identical cells d=0.5 Object B NOT_FIG4B");
    }

    private static ArmResult runLedger(BSimRandom rng, DaninoSIPeriodCheck.Protocol peaks)
            throws IOException {
        DaninoSIC0WellMixed device = new DaninoSIC0WellMixed(rng);
        double[] y0 = {0.0, 0.0, 0.05, 0.05};
        DaninoSIC0WellMixed.Trajectory traj = device.integrate(
                0.0,
                identicalCells(device.nCells(), y0),
                fill(device.nCells(), 0.05),
                0.05,
                T_END,
                SAMPLE_DT,
                RK_DT);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        boolean ledgerOk = traj.relativeResidual() <= LEDGER_REL
                && Math.abs(traj.ledger.muLoss()) <= 1e-18
                && Math.abs(traj.ledger.membraneCancel()) == 0.0;
        writeMeanCsv(RESULTS.resolve("java_C0_MU0_LEDGER.csv"), traj);
        System.out.printf(
                "  C0_MU0_LEDGER NOT_FIG4B mu=0 flag=%s residual_rel=%.3g mu_loss=%.3g gammaH_loss=%.6g synth=%.6g cancel=%.3g gate=%s%n",
                per.flag,
                traj.relativeResidual(),
                traj.ledger.muLoss(),
                traj.ledger.gammaHLoss(),
                traj.ledger.synthSource(),
                traj.ledger.membraneCancel(),
                ledgerOk ? "PASS" : "FAIL");
        return new ArmResult(
                "C0_MU0_LEDGER", "AHL_KICK_005", 0.0, per.flag, per.period, Double.NaN,
                true, true, 0.0, 0.0,
                traj.relativeResidual(), traj.maxHiSpread, Double.NaN, ledgerOk,
                "mu=0 ledger split gammaH vs remaining NOT_FIG4B");
    }

    private static ArmResult runClosedMembrane(BSimRandom rng) throws IOException {
        DaninoSIC0WellMixed device = DaninoSIC0WellMixed.closedMembrane(rng);
        double[] y0 = {0.0, 0.0, 0.05};
        double[][] cells = new double[device.nCells()][3];
        double[] hist = new double[device.nCells()];
        for (int c = 0; c < device.nCells(); c++) {
            cells[c][0] = 0.0;
            cells[c][1] = 0.0;
            cells[c][2] = 0.05;
            hist[c] = 0.05;
        }
        DaninoSIC0WellMixed.Trajectory traj = device.integrate(
                0.0, cells, hist, 0.0, T_END_CLOSED, SAMPLE_DT, RK_DT);
        boolean ok = traj.relativeResidual() <= LEDGER_REL
                && Math.abs(traj.ledger.synthSource()) <= 1e-18
                && Math.abs(traj.ledger.gammaHLoss()) <= 1e-18
                && Math.abs(traj.ledger.muLoss()) <= 1e-18
                && Math.abs(traj.ledger.membraneCancel()) == 0.0;
        writeMeanCsv(RESULTS.resolve("java_C0_MU0_CLOSED_MEMBRANE.csv"), traj);
        System.out.printf(
                "  C0_MU0_CLOSED_MEMBRANE NOT_FIG4B residual_rel=%.3g remaining=%.16g cancel=%.3g gate=%s%n",
                traj.relativeResidual(),
                traj.remainingMass,
                traj.ledger.membraneCancel(),
                ok ? "PASS" : "FAIL");
        return new ArmResult(
                "C0_MU0_CLOSED_MEMBRANE", "MEMBRANE_ONLY", 0.0, "n/a", Double.NaN, Double.NaN,
                true, true, 0.0, 0.0,
                traj.relativeResidual(), traj.maxHiSpread, Double.NaN, ok,
                "b=gammaH=mu=0 conservative Hi+He NOT_FIG4B");
    }

    private static ArmResult runHetero(BSimRandom rng, DaninoSIPeriodCheck.Protocol peaks)
            throws IOException {
        DaninoSIC0WellMixed device = new DaninoSIC0WellMixed(rng);
        int n = device.nCells();
        double[][] cells = new double[n][3];
        double[] hist = new double[n];
        cells[0][0] = 0.0;
        cells[0][1] = 0.0;
        cells[0][2] = 0.05;
        hist[0] = 0.05;
        for (int c = 1; c < n; c++) {
            cells[c][0] = 0.0;
            cells[c][1] = 0.0;
            cells[c][2] = 0.0;
            hist[c] = 0.0;
        }
        DaninoSIC0WellMixed.Trajectory traj = device.integrate(
                0.40, cells, hist, 0.05, T_END, SAMPLE_DT, RK_DT);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        double he0 = traj.Y[0][3];
        double maxMove = 0.0;
        for (int k = 0; k < traj.t.length; k++) {
            maxMove = Math.max(maxMove, Math.abs(traj.Y[k][3] - he0));
        }
        boolean moved = maxMove > HETERO_HE_MOVE;
        writeHeteroCsv(RESULTS.resolve("java_C0_HETERO_ONE_KICK.csv"), traj);
        System.out.printf(
                "  C0_HETERO_ONE_KICK NOT_FIG4B not identity not Fig4e flag=%s max_|dHe|=%.6g He_moved=%s gate=%s%n",
                per.flag,
                maxMove,
                moved,
                moved ? "PASS" : "FAIL");
        return new ArmResult(
                "C0_HETERO_ONE_KICK", "ONE_KICK_OTHERS_BASAL", 0.40, per.flag, per.period, Double.NaN,
                true, true, 0.0, 0.0,
                traj.relativeResidual(), traj.maxHiSpread, maxMove, moved,
                "heterogeneous smoke shared He must move NOT_FIG4B not Fig. 4e");
    }

    private static Path resolveFixture(String name) {
        if (name == null) {
            return null;
        }
        Path npz = D1_FIXTURES.resolve(name + ".npz");
        if (Files.exists(npz)) {
            return npz;
        }
        npz = D0B_FIXTURES.resolve(name + ".npz");
        return Files.exists(npz) ? npz : null;
    }

    private static double[][] identicalCells(int n, double[] y0) {
        double[][] cells = new double[n][3];
        for (int c = 0; c < n; c++) {
            cells[c][0] = y0[0];
            cells[c][1] = y0[1];
            cells[c][2] = y0[2];
        }
        return cells;
    }

    private static double[] fill(int n, double v) {
        double[] x = new double[n];
        java.util.Arrays.fill(x, v);
        return x;
    }

    private static TrajStats compare(
            double[] tPy,
            double[][] yPy,
            DaninoSIC0WellMixed.Trajectory traj,
            double atol,
            double rtol,
            double tLo,
            double tHi) {
        if (tPy.length != traj.t.length) {
            throw new IllegalStateException(
                    "sample grid length mismatch python=" + tPy.length + " c0=" + traj.t.length);
        }
        int n = 0;
        int viol = 0;
        double maxNorm = 0.0;
        String worst = "";
        String[] names = {"A", "I", "Hi", "He"};
        for (int i = 0; i < tPy.length; i++) {
            double t = tPy[i];
            if (t < tLo - 1e-12 || t > tHi + 1e-12) {
                continue;
            }
            if (Math.abs(t - traj.t[i]) > 1e-9) {
                throw new IllegalStateException("sample time mismatch at " + i);
            }
            for (int j = 0; j < 4; j++) {
                n++;
                double py = yPy[i][j];
                double jv = traj.Y[i][j];
                double allowed = atol + rtol * Math.abs(py);
                double err = Math.abs(jv - py);
                double norm = allowed > 0.0 ? err / allowed : err;
                if (norm > maxNorm) {
                    maxNorm = norm;
                    worst = String.format(Locale.US, "t=%.1f %s py=%.6g c0=%.6g err=%.3g allow=%.3g",
                            t, names[j], py, jv, err, allowed);
                }
                if (err > allowed) {
                    viol++;
                }
            }
        }
        return new TrajStats(n, viol, maxNorm, worst, viol == 0);
    }

    private static void writeMeanCsv(Path path, DaninoSIC0WellMixed.Trajectory traj) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.println("t,A,I,Hi,He,NOT_FIG4B");
            for (int i = 0; i < traj.t.length; i++) {
                w.printf(Locale.US, "%.6f,%.16e,%.16e,%.16e,%.16e,NOT_FIG4B%n",
                        traj.t[i], traj.Y[i][0], traj.Y[i][1], traj.Y[i][2], traj.Y[i][3]);
            }
        }
    }

    private static void writeHeteroCsv(Path path, DaninoSIC0WellMixed.Trajectory traj) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.print("t,A_mean,I_mean,Hi_mean,He");
            for (int c = 0; c < traj.Hi[0].length; c++) {
                w.print(",Hi_" + c);
            }
            w.println(",NOT_FIG4B");
            for (int i = 0; i < traj.t.length; i++) {
                w.printf(Locale.US, "%.6f,%.16e,%.16e,%.16e,%.16e",
                        traj.t[i], traj.Y[i][0], traj.Y[i][1], traj.Y[i][2], traj.Y[i][3]);
                for (int c = 0; c < traj.Hi[i].length; c++) {
                    w.printf(Locale.US, ",%.16e", traj.Hi[i][c]);
                }
                w.println(",NOT_FIG4B");
            }
        }
    }

    private static void writeSummary(DaninoSIC0WellMixed device, List<ArmResult> arms, boolean pass)
            throws IOException {
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
        json.append("  \"C0\": \"").append(pass ? "PASS" : "FAIL").append("\",\n");
        json.append("  \"NOT_FIG4B\": true,\n");
        json.append("  \"c1_may_start\": ").append(pass).append(",\n");
        json.append("  \"n_cells\": ").append(device.nCells()).append(",\n");
        json.append("  \"v_cell_um3\": ").append(device.vCell()).append(",\n");
        json.append("  \"V_e_um3\": ").append(device.vE()).append(",\n");
        json.append("  \"d\": ").append(device.d()).append(",\n");
        json.append("  \"mechanics\": \"OFF\",\n");
        json.append("  \"TIME_ADJ\": false,\n");
        json.append("  \"arms\": [\n");
        for (int i = 0; i < arms.size(); i++) {
            ArmResult a = arms.get(i);
            json.append("    {");
            json.append("\"name\":\"").append(a.name).append("\",");
            json.append("\"ensemble\":\"").append(a.ensemble).append("\",");
            json.append("\"mu\":").append(a.mu).append(",");
            json.append("\"flag\":\"").append(a.flag).append("\",");
            json.append("\"period\":").append(Double.isNaN(a.period) ? "null" : a.period).append(",");
            json.append("\"period_rel_err\":").append(Double.isNaN(a.relPeriod) ? "null" : a.relPeriod).append(",");
            json.append("\"traj_ok\":").append(a.trajOk).append(",");
            json.append("\"first_ok\":").append(a.firstOk).append(",");
            json.append("\"traj_max_norm\":").append(Double.isNaN(a.trajMaxNorm) ? "null" : a.trajMaxNorm).append(",");
            json.append("\"ledger_rel\":").append(a.ledgerRel).append(",");
            json.append("\"hi_spread\":").append(a.hiSpread).append(",");
            json.append("\"he_move\":").append(Double.isNaN(a.heMove) ? "null" : a.heMove).append(",");
            json.append("\"pass\":").append(a.pass).append(",");
            json.append("\"NOT_FIG4B\": true");
            json.append("}");
            json.append(i + 1 < arms.size() ? ",\n" : "\n");
        }
        json.append("  ]\n}\n");
        Files.writeString(RESULTS.resolve("c0_summary.json"), json.toString(), StandardCharsets.UTF_8);

        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(RESULTS.resolve("c0_arms.csv"), StandardCharsets.UTF_8))) {
            w.println("name,ensemble,mu,flag,period,period_rel_err,traj_ok,ledger_rel,hi_spread,he_move,pass,NOT_FIG4B,note");
            for (ArmResult a : arms) {
                w.printf(Locale.US, "%s,%s,%.4f,%s,%s,%s,%s,%.6g,%.6g,%s,%s,NOT_FIG4B,%s%n",
                        a.name, a.ensemble, a.mu, a.flag, fmt(a.period), fmt(a.relPeriod),
                        a.trajOk, a.ledgerRel, a.hiSpread,
                        Double.isNaN(a.heMove) ? "nan" : fmt(a.heMove),
                        a.pass, a.note.replace(',', ';'));
            }
        }
    }

    private static String fmt(double x) {
        if (Double.isNaN(x)) {
            return "nan";
        }
        return String.format(Locale.US, "%.6g", x);
    }

    private static final class TrajStats {
        final int n;
        final int violations;
        final double maxNorm;
        final String worst;
        final boolean allInside;

        TrajStats(int n, int violations, double maxNorm, String worst, boolean allInside) {
            this.n = n;
            this.violations = violations;
            this.maxNorm = maxNorm;
            this.worst = worst;
            this.allInside = allInside;
        }

        static TrajStats skipped() {
            return new TrajStats(0, 0, Double.NaN, "", true);
        }
    }

    private static final class ArmResult {
        final String name;
        final String ensemble;
        final double mu;
        final String flag;
        final double period;
        final double relPeriod;
        final boolean trajOk;
        final boolean firstOk;
        final double trajMaxNorm;
        final double firstMaxNorm;
        final double ledgerRel;
        final double hiSpread;
        final double heMove;
        final boolean pass;
        final String note;

        ArmResult(
                String name,
                String ensemble,
                double mu,
                String flag,
                double period,
                double relPeriod,
                boolean trajOk,
                boolean firstOk,
                double trajMaxNorm,
                double firstMaxNorm,
                double ledgerRel,
                double hiSpread,
                double heMove,
                boolean pass,
                String note) {
            this.name = name;
            this.ensemble = ensemble;
            this.mu = mu;
            this.flag = flag;
            this.period = period;
            this.relPeriod = relPeriod;
            this.trajOk = trajOk;
            this.firstOk = firstOk;
            this.trajMaxNorm = trajMaxNorm;
            this.firstMaxNorm = firstMaxNorm;
            this.ledgerRel = ledgerRel;
            this.hiSpread = hiSpread;
            this.heMove = heMove;
            this.pass = pass;
            this.note = note;
        }
    }
}
