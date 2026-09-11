package bsim.c1c;

import bsim.BSimRandom;
import bsim.circuit.DaninoSIPeriodCheck;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Gate C1c: filled-pocket spatial occupancy of Object B. NOT_FIG4B.
 * C1 open-chip standing remains FAIL. C1b island PASS is unchanged.
 * Does not start A0/W0 on FAIL. Does not retune Object B onto C1_OPEN_DILUTE.
 */
public final class DaninoSIC1cJob {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path CLAIM_FREEZE = ROOT.resolve("examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md");
    private static final Path C1_STANDING = ROOT.resolve("examples/PocketDish/C1_PACKED_SPATIAL_STANDING.md");
    private static final Path C1B_STANDING = ROOT.resolve("examples/PocketDish/C1B_ISLAND_STANDING.md");
    private static final Path PROTOCOL = ROOT.resolve("examples/BSimReservoirPlanDaninoPocketC1c/PROTOCOL.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/BSimReservoirPlanDaninoPocketC1c/configs/protocol.json");
    private static final Path RESULTS = ROOT.resolve("examples/BSimReservoirPlanDaninoPocketC1c/results");
    private static final Path C1_SUMMARY = ROOT.resolve(
            "examples/BSimReservoirPlanDaninoPocketC1/results/c1_summary.json");

    private static final double RK_DT = 0.001;
    private static final double SAMPLE_DT = 0.5;
    private static final double T_END = 1000.0;
    private static final double PY_T_040 = 63.583333333333336;
    private static final double PERIOD_REL_BAND = 0.20;
    private static final double LEDGER_REL = 1e-3;
    private static final double HE_RANGE_MIN = 1e-6;

    private DaninoSIC1cJob() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        DaninoSIPeriodCheck.refuseNarma(args);
        requireClaimFreeze();
        requireC1StillFail();
        requireC1bIslandPass();
        requireProtocolFrozen();
        Files.createDirectories(RESULTS);

        BSimRandom rng = new BSimRandom(0L);
        DaninoSIC1cFilledPocket device = new DaninoSIC1cFilledPocket(rng);
        if (device.unusedRng() != rng) {
            throw new IllegalStateException("RNG injection contract broken NOT_FIG4B");
        }

        System.out.println("C1c filled-pocket Object B. DaninoSI_OccupiedDF. NOT_FIG4B.");
        System.out.println("C1 open dilute remains FAIL. C1b island PASS unchanged.");
        System.out.println("Not Fig. 4b. Object A not ported. Not A0. Not W0.");
        System.out.printf(
                Locale.US,
                "filled pocket NOT_FIG4B N=%d v_cell=%.6g Ve_pocket=%.6g n_vox=%d d=%.16g "
                        + "D1=%s D1_spatial=%.6g open_edge=%s bus=CHEMICAL_SOLIDS%n",
                device.nCells(), device.vCell(), device.vEPocket(), device.nVoxels(),
                device.d(), DaninoSIC1cFilledPocket.D1_LABEL, device.d1Spatial(),
                device.openEdgeBc());

        DaninoSIPeriodCheck.Protocol peaks = DaninoSIPeriodCheck.Protocol.d0b();
        List<ArmResult> arms = new ArrayList<>();
        String want = args.length > 0 ? args[0].trim() : "";

        if (wantEmptyOr(want, "C1c_FILLED_POCKET_D05_MU040")) {
            arms.add(runIdentity(device, peaks));
        }
        if (wantEmptyOr(want, "C1c_WRONG_KICK")) {
            arms.add(runControl(device, peaks, "C1c_WRONG_KICK", "SI_BASAL_PERTURB", 0.40, "NO_PERIOD"));
        }
        if (wantEmptyOr(want, "C1c_MU150")) {
            arms.add(runControl(device, peaks, "C1c_MU150", "AHL_KICK_005", 1.50, "NO_PERIOD"));
        }
        if (wantEmptyOr(want, "C1c_MU0_LEDGER")) {
            arms.add(runLedger(device, peaks));
        }
        if (wantEmptyOr(want, "C1_OPEN_DILUTE") || want.isEmpty()) {
            arms.add(reportOpenDilute());
        }
        if ("C1c_BUS_CONNECTED".equals(want)) {
            arms.add(skipBusConnected("invoked explicitly; not implemented as a retune target NOT_FIG4B"));
        } else if (want.isEmpty()) {
            arms.add(skipBusConnected("optional extra; not required for PASS; predeclared NO_PERIOD NOT_FIG4B"));
        }

        if (arms.isEmpty()) {
            throw new IllegalArgumentException("no C1c arms selected: " + want);
        }

        boolean pass = true;
        for (ArmResult a : arms) {
            if ("C1_OPEN_DILUTE".equals(a.name) || "C1c_BUS_CONNECTED".equals(a.name)) {
                continue;
            }
            pass = pass && a.pass;
        }

        writeSummary(device, arms, pass);
        System.out.println(pass
                ? "C1c=PASS NOT_FIG4B on filled-pocket identity. C1 open chip remains FAIL. "
                + "A0 may start only as an AC on C1c_FILLED_POCKET, still NOT_FIG4B, not on C1_OPEN_DILUTE. W0 later."
                : "C1c=FAIL NOT_FIG4B. Do not retune Object B / D1_spatial / mu / d. "
                + "Do not start A0/W0. Do not occupy C1_OPEN_DILUTE. C1 FAIL unchanged.");
        if (!pass) {
            System.exit(1);
        }
    }

    private static boolean wantEmptyOr(String want, String name) {
        return want.isEmpty() || want.equals(name);
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

    private static void requireC1bIslandPass() throws IOException {
        String standing = Files.readString(C1B_STANDING, StandardCharsets.UTF_8);
        if (!standing.contains("**Status: PASS**")) {
            throw new IllegalStateException("C1b island standing must remain PASS");
        }
        if (standing.contains("C1 packed-spatial standing remains FAIL") == false
                && standing.contains("**C1 packed-spatial standing remains FAIL.**") == false) {
            throw new IllegalStateException("C1b standing must still record C1 FAIL");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("C1c protocol is not frozen_before_traces");
        }
        if (!text.contains("C1c_FILLED_POCKET") || !text.contains("NOT_FIG4B")) {
            throw new IllegalStateException("C1c protocol must freeze filled-pocket identity NOT_FIG4B");
        }
        if (!text.contains("CHEMICAL") && !json.contains("CHEMICAL_SOLIDS")) {
            throw new IllegalStateException("C1c protocol must freeze bus as chemical solids");
        }
        if (json.contains("\"TIME_ADJ\": true") || !text.contains("TIME_ADJ")) {
            throw new IllegalStateException("C1c protocol must freeze TIME_ADJ unused");
        }
        if (!text.contains("D1_spatial") || !json.contains("\"D1_spatial\": 800.0")) {
            throw new IllegalStateException("C1c identity D1_spatial must be frozen at 800");
        }
    }

    private static ArmResult runIdentity(
            DaninoSIC1cFilledPocket device,
            DaninoSIPeriodCheck.Protocol peaks) throws IOException {
        long t0 = System.nanoTime();
        System.out.println("  C1c_FILLED_POCKET_D05_MU040 NOT_FIG4B starting mu=0.40 D1=800");
        DaninoSIC1cFilledPocket.Trajectory traj = device.integrate(
                0.40, "AHL_KICK_005", T_END, SAMPLE_DT, RK_DT);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        if (traj.negativeState) {
            per = new DaninoSIPeriodCheck.Result(
                    "NEGATIVE_STATE", per.nPeaks, per.nPersistPeaks, Double.NaN,
                    per.peakTimes, per.relAmplitude, per.ampPersist, per.meanI, per.p5I, per.p95I);
        }
        double relT = Double.NaN;
        boolean periodOk = false;
        if ("OSC".equals(per.flag) && !Double.isNaN(per.period)) {
            relT = Math.abs(per.period - PY_T_040) / PY_T_040;
            periodOk = relT <= PERIOD_REL_BAND;
        }
        boolean heOk = traj.hePocketRangeMax >= HE_RANGE_MIN;
        boolean armPass = periodOk && heOk && !traj.negativeState;
        writeMeanCsv(RESULTS.resolve("java_C1c_FILLED_POCKET_D05_MU040.csv"), traj);
        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  C1c_FILLED_POCKET_D05_MU040 NOT_FIG4B mu=0.40 flag=%s T=%s relT=%s "
                        + "He_pocket_range=%.3g residual_rel=%.3g wall_s=%.1f gate=%s%n",
                per.flag, fmt(per.period), fmt(relT),
                traj.hePocketRangeMax, traj.relativeResidual(), sec,
                armPass ? "PASS" : "FAIL");
        return new ArmResult(
                "C1c_FILLED_POCKET_D05_MU040", "AHL_KICK_005", 0.40, per.flag, per.period, relT,
                traj.hePocketRangeMax, traj.relativeResidual(),
                traj.ledger.membraneCancel(), traj.extra.getDecayLoss(),
                traj.ledger.gammaHLoss(), traj.ledger.synthSource(),
                traj.extra.getOutletLoss(), traj.extra.getBoundaryLoss(),
                armPass, "filled-pocket identity D1=800 NOT_FIG4B");
    }

    private static ArmResult runControl(
            DaninoSIC1cFilledPocket device,
            DaninoSIPeriodCheck.Protocol peaks,
            String name,
            String ensemble,
            double mu,
            String requiredFlag) throws IOException {
        long t0 = System.nanoTime();
        System.out.printf("  %s NOT_FIG4B starting mu=%.2f ensemble=%s%n", name, mu, ensemble);
        DaninoSIC1cFilledPocket.Trajectory traj = device.integrate(
                mu, ensemble, T_END, SAMPLE_DT, RK_DT);
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
                name, ensemble, mu, per.flag, per.period, Double.NaN,
                traj.hePocketRangeMax, traj.relativeResidual(),
                traj.ledger.membraneCancel(), traj.extra.getDecayLoss(),
                traj.ledger.gammaHLoss(), traj.ledger.synthSource(),
                traj.extra.getOutletLoss(), traj.extra.getBoundaryLoss(),
                armPass, "control NOT_FIG4B");
    }

    private static ArmResult runLedger(
            DaninoSIC1cFilledPocket device,
            DaninoSIPeriodCheck.Protocol peaks) throws IOException {
        long t0 = System.nanoTime();
        System.out.println("  C1c_MU0_LEDGER NOT_FIG4B starting mu=0");
        DaninoSIC1cFilledPocket.Trajectory traj = device.integrate(
                0.0, "AHL_KICK_005", T_END, SAMPLE_DT, RK_DT);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        boolean ok = traj.relativeResidual() <= LEDGER_REL
                && Math.abs(traj.extra.getDecayLoss()) <= 1e-18
                && Math.abs(traj.ledger.membraneCancel()) <= 1e-12;
        writeMeanCsv(RESULTS.resolve("java_C1c_MU0_LEDGER.csv"), traj);
        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  C1c_MU0_LEDGER NOT_FIG4B residual_rel=%.3g mu_loss=%.3g gammaH=%.6g synth=%.6g "
                        + "cancel=%.3g outlet=%.3g boundary=%.3g wall_s=%.1f gate=%s%n",
                traj.relativeResidual(), traj.extra.getDecayLoss(), traj.ledger.gammaHLoss(),
                traj.ledger.synthSource(), traj.ledger.membraneCancel(),
                traj.extra.getOutletLoss(), traj.extra.getBoundaryLoss(), sec,
                ok ? "PASS" : "FAIL");
        return new ArmResult(
                "C1c_MU0_LEDGER", "AHL_KICK_005", 0.0, per.flag, per.period, Double.NaN,
                traj.hePocketRangeMax, traj.relativeResidual(),
                traj.ledger.membraneCancel(), traj.extra.getDecayLoss(),
                traj.ledger.gammaHLoss(), traj.ledger.synthSource(),
                traj.extra.getOutletLoss(), traj.extra.getBoundaryLoss(),
                ok, "mu=0 filled-pocket ledger NOT_FIG4B");
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
                Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                Double.NaN, Double.NaN, ok,
                "predeclared NO_PERIOD report-only; do not retune NOT_FIG4B");
    }

    private static ArmResult skipBusConnected(String note) {
        System.out.printf("  C1c_BUS_CONNECTED NOT_FIG4B skipped (%s)%n", note);
        return new ArmResult(
                "C1c_BUS_CONNECTED", "SKIPPED", 0.40, "SKIPPED", Double.NaN, Double.NaN,
                Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                Double.NaN, Double.NaN, true,
                note);
    }

    private static void writeMeanCsv(Path path, DaninoSIC1cFilledPocket.Trajectory traj)
            throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.println("t,A_mean,I_mean,Hi_mean,He_pocket_mean,He_pocket_min,He_pocket_max,NOT_FIG4B");
            for (int i = 0; i < traj.t.length; i++) {
                w.printf(Locale.US, "%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,NOT_FIG4B%n",
                        traj.t[i], traj.Y[i][0], traj.Y[i][1], traj.Y[i][2],
                        traj.Y[i][3], traj.hePocketMin[i], traj.hePocketMax[i]);
            }
        }
    }

    private static void writeSummary(
            DaninoSIC1cFilledPocket device,
            List<ArmResult> arms,
            boolean pass) throws IOException {
        StringBuilder json = new StringBuilder();
        json.append("{\n");
        json.append("  \"object\": \"DaninoSI_OccupiedDF\",\n");
        json.append("  \"object_status\": \"OCCUPIED_NOT_FIG4B\",\n");
        json.append("  \"C1\": \"FAIL\",\n");
        json.append("  \"C1b\": \"PASS_ISLAND_ONLY\",\n");
        json.append("  \"C1c\": \"").append(pass ? "PASS" : "FAIL").append("\",\n");
        json.append("  \"NOT_FIG4B\": true,\n");
        json.append("  \"open_chip_still_off\": true,\n");
        json.append("  \"a0_on_open_dilute\": false,\n");
        json.append("  \"a0_may_start_on_c1c_filled_pocket\": ").append(pass).append(",\n");
        json.append("  \"w0_may_start\": false,\n");
        json.append("  \"n_cells\": ").append(device.nCells()).append(",\n");
        json.append("  \"n_voxels\": ").append(device.nVoxels()).append(",\n");
        json.append("  \"Ve_pocket_um3\": ").append(device.vEPocket()).append(",\n");
        json.append("  \"d\": ").append(device.d()).append(",\n");
        json.append("  \"D1_spatial\": ").append(device.d1Spatial()).append(",\n");
        json.append("  \"open_edge_bc\": \"").append(device.openEdgeBc()).append("\",\n");
        json.append("  \"bus_on_identity\": \"CHEMICAL_SOLIDS\",\n");
        json.append("  \"TIME_ADJ\": false,\n");
        json.append("  \"mechanics\": \"OFF\",\n");
        json.append("  \"arms\": [\n");
        for (int i = 0; i < arms.size(); i++) {
            ArmResult a = arms.get(i);
            json.append("    {\"name\":\"").append(a.name).append("\",");
            json.append("\"flag\":\"").append(a.flag).append("\",");
            json.append("\"period\":").append(num(a.period)).append(",");
            json.append("\"period_rel_err\":").append(num(a.relPeriod)).append(",");
            json.append("\"he_pocket_range\":").append(num(a.heRange)).append(",");
            json.append("\"ledger_rel\":").append(num(a.ledgerRel)).append(",");
            json.append("\"mu_loss\":").append(num(a.muLoss)).append(",");
            json.append("\"gammaH_loss\":").append(num(a.gammaH)).append(",");
            json.append("\"synth\":").append(num(a.synth)).append(",");
            json.append("\"pass\":").append(a.pass).append(",\"NOT_FIG4B\": true}");
            json.append(i + 1 < arms.size() ? ",\n" : "\n");
        }
        json.append("  ]\n}\n");
        Files.writeString(RESULTS.resolve("c1c_summary.json"), json.toString(), StandardCharsets.UTF_8);
    }

    private static String fmt(double x) {
        return Double.isNaN(x) ? "nan" : String.format(Locale.US, "%.6g", x);
    }

    private static String num(double x) {
        return Double.isNaN(x) ? "null" : Double.toString(x);
    }

    private static final class ArmResult {
        final String name;
        final String ensemble;
        final double mu;
        final String flag;
        final double period;
        final double relPeriod;
        final double heRange;
        final double ledgerRel;
        final double memCancel;
        final double muLoss;
        final double gammaH;
        final double synth;
        final double outlet;
        final double boundary;
        final boolean pass;
        final String note;

        ArmResult(String name, String ensemble, double mu, String flag, double period, double relPeriod,
                  double heRange, double ledgerRel, double memCancel, double muLoss, double gammaH,
                  double synth, double outlet, double boundary, boolean pass, String note) {
            this.name = name;
            this.ensemble = ensemble;
            this.mu = mu;
            this.flag = flag;
            this.period = period;
            this.relPeriod = relPeriod;
            this.heRange = heRange;
            this.ledgerRel = ledgerRel;
            this.memCancel = memCancel;
            this.muLoss = muLoss;
            this.gammaH = gammaH;
            this.synth = synth;
            this.outlet = outlet;
            this.boundary = boundary;
            this.pass = pass;
            this.note = note;
        }
    }
}
