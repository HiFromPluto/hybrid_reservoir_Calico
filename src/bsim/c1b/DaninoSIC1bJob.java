package bsim.c1b;

import bsim.BSimRandom;
import bsim.c0.DaninoSIC0WellMixed;
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
 * Gate C1b: island C0-equivalent occupancy. NOT_FIG4B.
 * C1 open-chip standing remains FAIL. Does not start A0/W0 on identity
 * PASS if the open chip is still OFF.
 */
public final class DaninoSIC1bJob {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path CLAIM_FREEZE = ROOT.resolve("examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md");
    private static final Path C1_STANDING = ROOT.resolve("examples/PocketDish/C1_PACKED_SPATIAL_STANDING.md");
    private static final Path PROTOCOL = ROOT.resolve("examples/BSimReservoirPlanDaninoPocketC1b/PROTOCOL.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/BSimReservoirPlanDaninoPocketC1b/configs/protocol.json");
    private static final Path RESULTS = ROOT.resolve("examples/BSimReservoirPlanDaninoPocketC1b/results");
    private static final Path C1_SUMMARY = ROOT.resolve(
            "examples/BSimReservoirPlanDaninoPocketC1/results/c1_summary.json");
    private static final Path D1_FIXTURES = ROOT.resolve("examples/BSimReservoirPlanPocketOscSI_D1/fixtures");

    private static final double RK_DT = 0.001;
    private static final double SAMPLE_DT = 0.5;
    private static final double T_END = 1000.0;
    private static final double TRAJ_RTOL = 0.02;
    private static final double TRAJ_ATOL = 1e-4;
    private static final double FIRST_RTOL = 1e-3;
    private static final double FIRST_ATOL = 1e-6;
    private static final double PERIOD_REL = 0.02;
    private static final double LEDGER_REL = 1e-3;
    private static final double IDENTICAL_SPREAD = 1e-10;
    private static final double PY_T_040 = 63.583333333333336;

    private DaninoSIC1bJob() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        DaninoSIPeriodCheck.refuseNarma(args);
        requireClaimFreeze();
        requireC1StillFail();
        requireProtocolFrozen();
        Files.createDirectories(RESULTS);

        BSimRandom rng = new BSimRandom(0L);
        DaninoSIC1bIsland island = new DaninoSIC1bIsland(rng);
        DaninoSIC0WellMixed device = island.wellMixed();
        if (island.unusedRng() != rng || device.unusedRng() != rng) {
            throw new IllegalStateException("RNG injection contract broken NOT_FIG4B");
        }

        System.out.println("C1b island C0-equivalent Object B. DaninoSI_OccupiedDF. NOT_FIG4B.");
        System.out.println("C1 open dilute remains FAIL. Not Fig. 4b. Object A not ported.");
        System.out.printf(
                Locale.US,
                "island identity NOT_FIG4B N=%d v_cell=%.6g Ve=%.6g voxel=%.6g d=%.16g D1=%s n_fluid=%d%n",
                island.nCells(), island.vCell(), island.vE(), island.boxVolume(),
                island.d(), island.d1Spatial(), island.nFluidVoxels());

        double[] y = {0.1, 0.2, 0.3, 0.4};
        device.assertReducesToObjectB(y, 0.25, 0.40);
        device.assertReducesToObjectB(new double[] {0.0, 0.0, 0.05, 0.05}, 0.05, 0.40);
        System.out.println("  algebra_reduction NOT_FIG4B d=0.5 island He RHS matches Object B gate=PASS");

        DaninoSIPeriodCheck.Protocol peaks = DaninoSIPeriodCheck.Protocol.d0b();
        List<ArmResult> arms = new ArrayList<>();

        arms.add(runIsland(island, peaks, "C1b_ISLAND_D05_MU040", "AHL_KICK_005", 0.40,
                "OSC", PY_T_040, "primary_mu_0p40", true));
        boolean identityPass = arms.get(0).pass;

        arms.add(runIsland(island, peaks, "C1b_ISLAND_WRONG_KICK", "SI_BASAL_PERTURB", 0.40,
                "NO_PERIOD", Double.NaN, null, false));
        arms.add(runIsland(island, peaks, "C1b_ISLAND_MU150", "AHL_KICK_005", 1.50,
                "NO_PERIOD", Double.NaN, null, false));
        arms.add(runLedger(island, peaks));
        arms.add(reportOpenDilute());

        if (identityPass) {
            arms.add(runIslandD1800(rng, peaks));
        } else {
            System.out.println("  C1b_ISLAND_D1_800 NOT_FIG4B skipped: identity FAIL. Do not open D1_spatial.");
            arms.add(new ArmResult(
                    "C1b_ISLAND_D1_800", "SKIPPED", 0.40, "SKIPPED", Double.NaN, Double.NaN,
                    true, true, Double.NaN, 0.0, 0.0, true,
                    "not opened after identity FAIL NOT_FIG4B"));
        }

        boolean pass = true;
        for (ArmResult a : arms) {
            if (!"C1_OPEN_DILUTE".equals(a.name) && !"C1b_ISLAND_D1_800".equals(a.name)) {
                pass = pass && a.pass;
            }
        }
        pass = pass && identityPass;

        writeSummary(island, arms, pass);
        System.out.println(pass
                ? "C1b=PASS NOT_FIG4B on island identity. C1 open chip remains FAIL. A0 is not automatic. W0 later."
                : "C1b=FAIL NOT_FIG4B. Do not retune Object B. Do not start A0/W0. C1 FAIL unchanged.");
        if (!pass) {
            System.exit(1);
        }
    }

    private static void requireClaimFreeze() throws IOException {
        String freeze = Files.readString(CLAIM_FREEZE, StandardCharsets.UTF_8);
        if (!freeze.contains("DaninoSI_OccupiedDF") || !freeze.contains("NOT_FIG4B")) {
            throw new IllegalStateException("claim freeze must name Object B as NOT_FIG4B");
        }
        if (!freeze.contains("Danino2010_Fig4b_bulk_twin")) {
            throw new IllegalStateException("claim freeze must name Object A");
        }
    }

    private static void requireC1StillFail() throws IOException {
        String standing = Files.readString(C1_STANDING, StandardCharsets.UTF_8);
        if (!standing.contains("**Status: FAIL**")) {
            throw new IllegalStateException("C1 standing must remain FAIL; do not rewrite as PASS");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("C1b protocol is not frozen_before_traces");
        }
        if (!text.contains("C1b_ISLAND") || !text.contains("NOT_FIG4B")) {
            throw new IllegalStateException("C1b protocol must freeze island identity NOT_FIG4B");
        }
        if (!text.contains("C0 volumes") && !text.contains("C0's 8-cell") && !json.contains("\"n_cells\": 8")) {
            throw new IllegalStateException("C1b identity must freeze C0 volumes");
        }
    }

    private static ArmResult runIsland(
            DaninoSIC1bIsland island,
            DaninoSIPeriodCheck.Protocol peaks,
            String name,
            String ensemble,
            double mu,
            String requiredFlag,
            double pythonPeriod,
            String fixtureName,
            boolean compareTraj) throws IOException {
        DaninoSIC1bIsland.Trajectory traj = island.integrate(
                mu, ensemble, T_END, SAMPLE_DT, RK_DT);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        TrajStats full = TrajStats.skipped();
        TrajStats first = TrajStats.skipped();
        if (compareTraj && fixtureName != null) {
            Path npz = D1_FIXTURES.resolve(fixtureName + ".npz");
            Map<String, NpzReader.Array> arrays = NpzReader.loadNpz(npz);
            full = compare(arrays.get("t").vector(), arrays.get("Y").matrix2(), traj,
                    TRAJ_ATOL, TRAJ_RTOL, 0.0, T_END);
            first = compare(arrays.get("t").vector(), arrays.get("Y").matrix2(), traj,
                    FIRST_ATOL, FIRST_RTOL, 0.0, DaninoSIOccupiedDF.TAU);
        }
        boolean flagOk = requiredFlag.equals(per.flag);
        boolean periodOk = true;
        double relT = Double.NaN;
        if ("OSC".equals(requiredFlag)) {
            relT = Math.abs(per.period - pythonPeriod) / pythonPeriod;
            periodOk = "OSC".equals(per.flag) && relT < PERIOD_REL;
        }
        boolean spreadOk = traj.maxHiSpread <= IDENTICAL_SPREAD;
        boolean trajOk = !compareTraj || (full.allInside && first.allInside);
        boolean armPass = trajOk && flagOk && periodOk && spreadOk;
        writeMeanCsv(RESULTS.resolve("java_" + name + ".csv"), traj);
        System.out.printf(
                Locale.US,
                "  %s NOT_FIG4B mu=%.2f D1=0 flag=%s T=%s relT=%s traj_max=%s spread=%.3g residual_rel=%.3g gate=%s%n",
                name, mu, per.flag, fmt(per.period), fmt(relT),
                compareTraj ? String.format(Locale.US, "%.3g", full.maxNorm) : "n/a",
                traj.maxHiSpread, traj.relativeResidual(),
                armPass ? "PASS" : "FAIL");
        return new ArmResult(
                name, ensemble, mu, per.flag, per.period, relT,
                trajOk, first.allInside || !compareTraj, full.maxNorm,
                traj.relativeResidual(), traj.maxHiSpread, armPass,
                "island C0-equivalent D1=0 NOT_FIG4B");
    }

    private static ArmResult runLedger(
            DaninoSIC1bIsland island,
            DaninoSIPeriodCheck.Protocol peaks) throws IOException {
        DaninoSIC1bIsland.Trajectory traj = island.integrate(
                0.0, "AHL_KICK_005", T_END, SAMPLE_DT, RK_DT);
        boolean ok = traj.relativeResidual() <= LEDGER_REL
                && Math.abs(traj.ledger.muLoss()) <= 1e-18
                && Math.abs(traj.ledger.membraneCancel()) == 0.0;
        writeMeanCsv(RESULTS.resolve("java_C1b_ISLAND_MU0_LEDGER.csv"), traj);
        System.out.printf(
                Locale.US,
                "  C1b_ISLAND_MU0_LEDGER NOT_FIG4B residual_rel=%.3g mu_loss=%.3g gammaH=%.6g synth=%.6g gate=%s%n",
                traj.relativeResidual(), traj.ledger.muLoss(), traj.ledger.gammaHLoss(),
                traj.ledger.synthSource(), ok ? "PASS" : "FAIL");
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        return new ArmResult(
                "C1b_ISLAND_MU0_LEDGER", "AHL_KICK_005", 0.0, per.flag, per.period, Double.NaN,
                true, true, Double.NaN, traj.relativeResidual(), traj.maxHiSpread, ok,
                "mu=0 island ledger NOT_FIG4B");
    }

    private static ArmResult runIslandD1800(BSimRandom rng, DaninoSIPeriodCheck.Protocol peaks)
            throws IOException {
        System.out.println("  C1b_ISLAND_D1_800 NOT_FIG4B extra starting (not identity)");
        DaninoSIC1bIslandSpatial extra = new DaninoSIC1bIslandSpatial(rng);
        DaninoSIC1bIslandSpatial.Trajectory traj = extra.integrate(
                0.40, "AHL_KICK_005", T_END, SAMPLE_DT, RK_DT);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        writeSpatialCsv(RESULTS.resolve("java_C1b_ISLAND_D1_800.csv"), traj);
        double relT = Double.NaN;
        if ("OSC".equals(per.flag)) {
            relT = Math.abs(per.period - PY_T_040) / PY_T_040;
        }
        System.out.printf(
                Locale.US,
                "  C1b_ISLAND_D1_800 NOT_FIG4B not identity N=%d d=%.2f D1=800 flag=%s T=%s relT=%s "
                        + "He_range=%.3g residual_rel=%.3g reported-only%n",
                extra.nCells(), extra.d(), per.flag, fmt(per.period), fmt(relT),
                traj.heFluidRangeMax, traj.relativeResidual());
        return new ArmResult(
                "C1b_ISLAND_D1_800", "AHL_KICK_005", 0.40, per.flag, per.period, relT,
                true, true, Double.NaN, traj.relativeResidual(), Double.NaN, true,
                "optional island D1=800 report-only NOT_FIG4B");
    }

    private static ArmResult reportOpenDilute() throws IOException {
        String flag = "MISSING";
        boolean ok = false;
        if (Files.exists(C1_SUMMARY)) {
            String json = Files.readString(C1_SUMMARY, StandardCharsets.UTF_8);
            if (json.contains("\"C1\": \"FAIL\"") && json.contains("NO_PERIOD")) {
                flag = "NO_PERIOD";
                ok = true;
            }
        }
        System.out.printf(
                "  C1_OPEN_DILUTE NOT_FIG4B report-only flag=%s (C1 standing FAIL, not re-run) gate=%s%n",
                flag, ok ? "PASS" : "FAIL");
        return new ArmResult(
                "C1_OPEN_DILUTE", "AHL_KICK_005", 0.40, flag, Double.NaN, Double.NaN,
                true, true, Double.NaN, Double.NaN, Double.NaN, ok,
                "predeclared NO_PERIOD report-only; do not retune NOT_FIG4B");
    }

    private static TrajStats compare(
            double[] tPy, double[][] yPy, DaninoSIC1bIsland.Trajectory traj,
            double atol, double rtol, double tLo, double tHi) {
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
            for (int j = 0; j < 4; j++) {
                n++;
                double py = yPy[i][j];
                double jv = traj.Y[i][j];
                double allowed = atol + rtol * Math.abs(py);
                double err = Math.abs(jv - py);
                double norm = allowed > 0.0 ? err / allowed : err;
                if (norm > maxNorm) {
                    maxNorm = norm;
                    worst = String.format(Locale.US, "t=%.1f %s py=%.6g c1b=%.6g", t, names[j], py, jv);
                }
                if (err > allowed) {
                    viol++;
                }
            }
        }
        return new TrajStats(n, viol, maxNorm, worst, viol == 0);
    }

    private static void writeMeanCsv(Path path, DaninoSIC1bIsland.Trajectory traj) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.println("t,A,I,Hi,He,NOT_FIG4B");
            for (int i = 0; i < traj.t.length; i++) {
                w.printf(Locale.US, "%.6f,%.16e,%.16e,%.16e,%.16e,NOT_FIG4B%n",
                        traj.t[i], traj.Y[i][0], traj.Y[i][1], traj.Y[i][2], traj.Y[i][3]);
            }
        }
    }

    private static void writeSpatialCsv(Path path, DaninoSIC1bIslandSpatial.Trajectory traj)
            throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.println("t,A_mean,I_mean,Hi_mean,He_patch_mean,NOT_FIG4B");
            for (int i = 0; i < traj.t.length; i++) {
                w.printf(Locale.US, "%.6f,%.16e,%.16e,%.16e,%.16e,NOT_FIG4B%n",
                        traj.t[i], traj.Y[i][0], traj.Y[i][1], traj.Y[i][2], traj.Y[i][3]);
            }
        }
    }

    private static void writeSummary(DaninoSIC1bIsland island, List<ArmResult> arms, boolean pass)
            throws IOException {
        StringBuilder json = new StringBuilder();
        json.append("{\n");
        json.append("  \"object\": \"DaninoSI_OccupiedDF\",\n");
        json.append("  \"C1\": \"FAIL\",\n");
        json.append("  \"C1b\": \"").append(pass ? "PASS" : "FAIL").append("\",\n");
        json.append("  \"NOT_FIG4B\": true,\n");
        json.append("  \"open_chip_still_off\": true,\n");
        json.append("  \"a0_automatic\": false,\n");
        json.append("  \"w0_may_start\": false,\n");
        json.append("  \"n_cells\": ").append(island.nCells()).append(",\n");
        json.append("  \"V_e_um3\": ").append(island.vE()).append(",\n");
        json.append("  \"d\": ").append(island.d()).append(",\n");
        json.append("  \"D1_spatial_identity\": 0.0,\n");
        json.append("  \"arms\": [\n");
        for (int i = 0; i < arms.size(); i++) {
            ArmResult a = arms.get(i);
            json.append("    {\"name\":\"").append(a.name).append("\",");
            json.append("\"flag\":\"").append(a.flag).append("\",");
            json.append("\"period\":").append(Double.isNaN(a.period) ? "null" : a.period).append(",");
            json.append("\"period_rel_err\":").append(Double.isNaN(a.relPeriod) ? "null" : a.relPeriod).append(",");
            json.append("\"ledger_rel\":").append(Double.isNaN(a.ledgerRel) ? "null" : a.ledgerRel).append(",");
            json.append("\"pass\":").append(a.pass).append(",\"NOT_FIG4B\": true}");
            json.append(i + 1 < arms.size() ? ",\n" : "\n");
        }
        json.append("  ]\n}\n");
        Files.writeString(RESULTS.resolve("c1b_summary.json"), json.toString(), StandardCharsets.UTF_8);
    }

    private static String fmt(double x) {
        return Double.isNaN(x) ? "nan" : String.format(Locale.US, "%.6g", x);
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
        final double ledgerRel;
        final double hiSpread;
        final boolean pass;
        final String note;

        ArmResult(String name, String ensemble, double mu, String flag, double period, double relPeriod,
                  boolean trajOk, boolean firstOk, double trajMaxNorm, double ledgerRel, double hiSpread,
                  boolean pass, String note) {
            this.name = name;
            this.ensemble = ensemble;
            this.mu = mu;
            this.flag = flag;
            this.period = period;
            this.relPeriod = relPeriod;
            this.trajOk = trajOk;
            this.firstOk = firstOk;
            this.trajMaxNorm = trajMaxNorm;
            this.ledgerRel = ledgerRel;
            this.hiSpread = hiSpread;
            this.pass = pass;
            this.note = note;
        }
    }
}
