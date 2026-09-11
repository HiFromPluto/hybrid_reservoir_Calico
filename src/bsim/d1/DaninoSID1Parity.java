package bsim.d1;

import bsim.BSimRandom;
import bsim.circuit.DaninoSIOccupiedDF;
import bsim.circuit.DaninoSIPeriodCheck;

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
 * Gate D1: Java vs frozen D0b Python fixtures of Object B.
 * NOT_FIG4B. Does not port Object A (Fig. 4b twin).
 */
public final class DaninoSID1Parity {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path CLAIM_FREEZE = ROOT.resolve("examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md");
    private static final Path FIXTURES = ROOT.resolve("examples/BSimReservoirPlanPocketOscSI_D1/fixtures");
    private static final Path RESULTS = ROOT.resolve("examples/BSimReservoirPlanPocketOscSI_D1/results");
    private static final Path D0B_FIXTURES = ROOT.resolve(
            "examples/BSimReservoirPlanPocketOscSI_D0b/results/fixtures");

    private static final double RK_DT = 0.001;
    private static final double SAMPLE_DT = 0.5;
    private static final double T_END = 1000.0;
    private static final double TRAJ_RTOL = 0.02;
    private static final double TRAJ_ATOL = 1e-4;
    private static final double FIRST_RTOL = 1e-3;
    private static final double FIRST_ATOL = 1e-6;
    private static final double PERIOD_REL = 0.02;
    private static final double PY_T_040 = 63.583333333333336;
    private static final double PY_T_032 = 56.30769230769231;

    private DaninoSID1Parity() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        DaninoSIPeriodCheck.refuseNarma(args);
        if (!Files.exists(CLAIM_FREEZE)) {
            throw new IllegalStateException("CLAIM_FREEZE_DANINO_SI.md must exist before D1 Java.");
        }
        String freeze = Files.readString(CLAIM_FREEZE, StandardCharsets.UTF_8);
        if (!freeze.contains("DaninoSI_OccupiedDF") || !freeze.contains("NOT_FIG4B")) {
            throw new IllegalStateException("claim freeze must name Object B as NOT_FIG4B");
        }
        if (!freeze.contains("Danino2010_Fig4b_bulk_twin")) {
            throw new IllegalStateException("claim freeze must name Object A");
        }
        Files.createDirectories(RESULTS);

        BSimRandom rng = new BSimRandom(0L);
        DaninoSIOccupiedDF model = new DaninoSIOccupiedDF(rng);
        if (model.unusedRng() != rng) {
            throw new IllegalStateException("RNG injection contract broken");
        }

        System.out.println("D1 Java circuit parity. Object B DaninoSI_OccupiedDF. NOT_FIG4B.");
        System.out.println("Object A Danino2010_Fig4b_bulk_twin is not ported.");
        System.out.println("D0 remains FAIL. D0b remains FAIL_NO_IDENTITY.");

        DaninoSIPeriodCheck.Protocol peaks = DaninoSIPeriodCheck.Protocol.d0b();
        List<ArmResult> arms = new ArrayList<>();

        arms.add(runFixture(model, peaks, "primary_mu_0p40", "AHL_KICK_005", true, "OSC", PY_T_040));
        arms.add(runFixture(model, peaks, "primary_osc", "AHL_KICK_005", true, "OSC", PY_T_032));
        arms.add(runFixture(model, peaks, "si_basal_perturb_mu_0p40", "SI_BASAL_PERTURB", true, "NO_PERIOD", Double.NaN));
        arms.add(runFixture(model, peaks, "failed_prior_mu_1p5", "FAILED_PRIOR_IVP", true, "NO_PERIOD", Double.NaN));

        double[] yPrimary = {0.0, 0.0, 0.05, 0.05};
        DaninoSIOccupiedDF.Trajectory t15 = model.integrate(0.5, 1.5, yPrimary, 0.05, T_END, SAMPLE_DT, RK_DT);
        DaninoSIPeriodCheck.Result p15 = DaninoSIPeriodCheck.extract(t15.t, t15.luxI(), peaks);
        if (t15.negativeState) {
            p15 = DaninoSIPeriodCheck.Result.noPeriod(p15.meanI);
        }
        boolean flag15 = "NO_PERIOD".equals(p15.flag);
        arms.add(new ArmResult(
                "primary_mu_1p5",
                "AHL_KICK_005",
                1.5,
                p15.flag,
                p15.period,
                Double.NaN,
                true,
                true,
                0.0,
                0.0,
                flag15,
                "NO_PERIOD expected for AHL_KICK_005 at mu=1.5 (D0b extra OFF/collapse)"));
        writeCsv(RESULTS.resolve("java_primary_mu_1p5.csv"), t15);
        System.out.printf(
                "  primary_mu_1p5 NOT_FIG4B flag=%s T=%s meanI=%.4g gate=%s%n",
                p15.flag,
                fmt(p15.period),
                p15.meanI,
                flag15 ? "PASS" : "FAIL");

        boolean pass = true;
        for (ArmResult a : arms) {
            pass = pass && a.pass;
        }

        writeSummary(arms, pass);
        System.out.println(pass ? "D1=PASS NOT_FIG4B" : "D1=FAIL NOT_FIG4B");
        if (!pass) {
            System.exit(1);
        }
    }

    private static ArmResult runFixture(
            DaninoSIOccupiedDF model,
            DaninoSIPeriodCheck.Protocol peaks,
            String name,
            String ensemble,
            boolean compareTraj,
            String requiredFlag,
            double pythonPeriod) throws IOException {
        Path npz = FIXTURES.resolve(name + ".npz");
        if (!Files.exists(npz)) {
            npz = D0B_FIXTURES.resolve(name + ".npz");
        }
        Map<String, NpzReader.Array> arrays = NpzReader.loadNpz(npz);
        double mu = arrays.get("mu").scalar();
        double d = arrays.get("d").scalar();
        double hiHist = arrays.get("H_i_history").scalar();
        double[] y0 = arrays.get("y0").vector();
        double[] tPy = arrays.get("t").vector();
        double[][] yPy = arrays.get("Y").matrix2();

        DaninoSIOccupiedDF.Trajectory traj = model.integrate(d, mu, y0, hiHist, T_END, SAMPLE_DT, RK_DT);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxI(), peaks);
        if (traj.negativeState) {
            per = new DaninoSIPeriodCheck.Result(
                    "NEGATIVE_STATE", per.nPeaks, per.nPersistPeaks, Double.NaN,
                    per.peakTimes, per.relAmplitude, per.ampPersist, per.meanI, per.p5I, per.p95I);
        }

        TrajStats full = compare(tPy, yPy, traj, TRAJ_ATOL, TRAJ_RTOL, 0.0, T_END);
        TrajStats first = compare(tPy, yPy, traj, FIRST_ATOL, FIRST_RTOL, 0.0, DaninoSIOccupiedDF.TAU);
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
        boolean armPass = (!compareTraj || (full.allInside && first.allInside)) && flagOk && periodOk;
        writeCsv(RESULTS.resolve("java_" + name + ".csv"), traj);
        System.out.printf(
                "  %s NOT_FIG4B mu=%.2f flag=%s T=%s relT=%s traj_max=%.3g first_max=%.3g gate=%s%n",
                name,
                mu,
                per.flag,
                fmt(per.period),
                fmt(relT),
                full.maxNorm,
                first.maxNorm,
                armPass ? "PASS" : "FAIL");
        if (!full.allInside) {
            System.out.printf("    traj violations=%d / %d worst=%s%n", full.violations, full.n, full.worst);
        }
        if (!first.allInside) {
            System.out.printf("    first-interval violations=%d / %d worst=%s%n", first.violations, first.n, first.worst);
        }
        return new ArmResult(
                name, ensemble, mu, per.flag, per.period, relT,
                full.allInside, first.allInside, full.maxNorm, first.maxNorm, armPass,
                "fixture " + name + " Object B NOT_FIG4B");
    }

    private static TrajStats compare(
            double[] tPy,
            double[][] yPy,
            DaninoSIOccupiedDF.Trajectory traj,
            double atol,
            double rtol,
            double tLo,
            double tHi) {
        if (tPy.length != traj.t.length) {
            throw new IllegalStateException("sample grid length mismatch python=" + tPy.length + " java=" + traj.t.length);
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
                    worst = String.format(Locale.US, "t=%.1f %s py=%.6g java=%.6g err=%.3g allow=%.3g",
                            t, names[j], py, jv, err, allowed);
                }
                if (err > allowed) {
                    viol++;
                }
            }
        }
        return new TrajStats(n, viol, maxNorm, worst, viol == 0);
    }

    private static void writeCsv(Path path, DaninoSIOccupiedDF.Trajectory traj) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.println("t,A,I,Hi,He");
            for (int i = 0; i < traj.t.length; i++) {
                w.printf(Locale.US, "%.6f,%.16e,%.16e,%.16e,%.16e%n",
                        traj.t[i], traj.Y[i][0], traj.Y[i][1], traj.Y[i][2], traj.Y[i][3]);
            }
        }
    }

    private static void writeSummary(List<ArmResult> arms, boolean pass) throws IOException {
        StringBuilder json = new StringBuilder();
        json.append("{\n");
        json.append("  \"object\": \"DaninoSI_OccupiedDF\",\n");
        json.append("  \"object_status\": \"OCCUPIED_NOT_FIG4B\",\n");
        json.append("  \"literature_twin\": \"Danino2010_Fig4b_bulk_twin\",\n");
        json.append("  \"literature_twin_status\": \"FAIL_NO_IDENTITY\",\n");
        json.append("  \"D0\": \"FAIL\",\n");
        json.append("  \"D0b\": \"FAIL_NO_IDENTITY\",\n");
        json.append("  \"D1\": \"").append(pass ? "PASS" : "FAIL").append("\",\n");
        json.append("  \"NOT_FIG4B\": true,\n");
        json.append("  \"d1_may_start_c0\": ").append(pass).append(",\n");
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
            json.append("\"traj_max_norm\":").append(a.trajMaxNorm).append(",");
            json.append("\"first_max_norm\":").append(a.firstMaxNorm).append(",");
            json.append("\"pass\":").append(a.pass);
            json.append("}");
            json.append(i + 1 < arms.size() ? ",\n" : "\n");
        }
        json.append("  ]\n}\n");
        Files.writeString(RESULTS.resolve("d1_summary.json"), json.toString(), StandardCharsets.UTF_8);

        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(RESULTS.resolve("d1_parity.csv"), StandardCharsets.UTF_8))) {
            w.println("name,ensemble,mu,flag,period,period_rel_err,traj_ok,first_ok,traj_max_norm,first_max_norm,pass,note");
            for (ArmResult a : arms) {
                w.printf(Locale.US, "%s,%s,%.4f,%s,%s,%s,%s,%s,%.6g,%.6g,%s,%s%n",
                        a.name, a.ensemble, a.mu, a.flag, fmt(a.period), fmt(a.relPeriod),
                        a.trajOk, a.firstOk, a.trajMaxNorm, a.firstMaxNorm, a.pass, a.note.replace(',', ';'));
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
            this.pass = pass;
            this.note = note;
        }
    }
}
