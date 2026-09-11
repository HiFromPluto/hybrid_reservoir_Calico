package bsim.a1;

import bsim.BSimRandom;
import bsim.a0.A0IdealSource;
import bsim.c1c.DaninoSIC1cFilledPocket;
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
 * Gate A1: bounded transducer on C1c_FILLED_POCKET.
 * A1_BOUNDED_TRANSDUCER. HYPOTHETICAL_DESIGN_ENVELOPE. NOT_FIG4B.
 * Not calibrated. Not Lentini TX-TL. Not A2/W0. Not C1_OPEN_DILUTE.
 */
public final class DaninoSIA1Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path CLAIM_FREEZE = ROOT.resolve("examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md");
    private static final Path C1_STANDING = ROOT.resolve("examples/PocketDish/C1_PACKED_SPATIAL_STANDING.md");
    private static final Path C1C_STANDING = ROOT.resolve("examples/PocketDish/C1C_FILLED_POCKET_STANDING.md");
    private static final Path A0_STANDING = ROOT.resolve("examples/PocketDish/A0_IDEAL_SOURCE_STANDING.md");
    private static final Path PROTOCOL = ROOT.resolve("examples/BSimReservoirPlanDaninoPocketA1/PROTOCOL.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/BSimReservoirPlanDaninoPocketA1/configs/protocol.json");
    private static final Path RESULTS = ROOT.resolve("examples/BSimReservoirPlanDaninoPocketA1/results");
    private static final Path A0_PULSE_CSV = ROOT.resolve(
            "examples/BSimReservoirPlanDaninoPocketA0/results/java_A0_J_PULSE_BACTERIA.csv");

    private static final double PY_T_040 = 63.583333333333336;
    private static final double C1C_T_040 = 61.25;
    private static final double PERIOD_REL_BAND = 0.20;
    private static final double LEDGER_REL = 1e-3;
    private static final double A0_CONTRAST = 0.10;

    private DaninoSIA1Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        DaninoSIPeriodCheck.refuseNarma(args);
        requireFreeze();
        Files.createDirectories(RESULTS);
        A1BoundedTransducer.freezeCheck();

        BSimRandom rng = new BSimRandom(0L);
        DaninoSIA1Field fieldDev = new DaninoSIA1Field(rng);
        DaninoSIC1cFilledPocket device = new DaninoSIC1cFilledPocket(rng);
        if (device.unusedRng() != rng || fieldDev.unusedRng() != rng) {
            throw new IllegalStateException("RNG injection contract broken A1_BOUNDED_TRANSDUCER NOT_FIG4B");
        }

        System.out.println("A1_BOUNDED_TRANSDUCER HYPOTHETICAL_DESIGN_ENVELOPE on C1c_FILLED_POCKET. NOT_FIG4B.");
        System.out.println("Not calibrated. Not Lentini TX-TL. Not Fig. 4b. Not A2/W0. Not C1_OPEN_DILUTE.");
        System.out.printf(
                Locale.US,
                "A1_BOUNDED_TRANSDUCER NOT_FIG4B tau_AC=%.3g n=%.3g K=%.3g J_max=%.6g J_leak=%.6g "
                        + "M0=%.6g M_half=%.6g payload=extracellular_AHL%n",
                A1BoundedTransducer.TAU_AC, A1BoundedTransducer.N_HILL, A1BoundedTransducer.K,
                A1BoundedTransducer.J_MAX, A1BoundedTransducer.J_LEAK,
                A1BoundedTransducer.M0, A1BoundedTransducer.M_HALF);

        DaninoSIPeriodCheck.Protocol peaks = DaninoSIPeriodCheck.Protocol.d0b();
        List<ArmResult> arms = new ArrayList<>();
        String want = args.length > 0 ? args[0].trim() : "";
        DaninoSIC1cFilledPocket.Trajectory j0 = null;

        if (wantEmptyOr(want, "A1_FIELD_STEP_U")) {
            arms.add(runField(fieldDev, A1BoundedTransducer.Command.U_STEP, "A1_FIELD_STEP_U",
                    A1BoundedTransducer.T_END_FIELD, true));
        }
        if (wantEmptyOr(want, "A1_FIELD_PULSE_U")) {
            arms.add(runField(fieldDev, A1BoundedTransducer.Command.U_PULSE, "A1_FIELD_PULSE_U",
                    A1BoundedTransducer.T_END_FIELD, true));
        }
        if (wantEmptyOr(want, "A1_U0_LEAK")) {
            arms.add(runLeak(fieldDev));
        }
        if (wantEmptyOr(want, "A1_PAYLOAD_EXHAUST")) {
            arms.add(runExhaust(fieldDev));
        }
        if (wantEmptyOr(want, "A1_J0_BACTERIA")) {
            ArmResult r = runJ0Bacteria(device, peaks);
            arms.add(r);
            j0 = r.traj;
        }
        if (wantEmptyOr(want, "A1_PULSE_BACTERIA")) {
            if (j0 == null) {
                System.out.println("  A1_PULSE_BACTERIA running A1_J0_BACTERIA first A1_BOUNDED_TRANSDUCER NOT_FIG4B");
                j0 = runJ0Bacteria(device, peaks).traj;
            }
            arms.add(runPulseBacteria(device, peaks, j0));
        }
        if (wantEmptyOr(want, "C1_OPEN_DILUTE") || want.isEmpty()) {
            arms.add(forbidOpenDilute());
        }
        if ("A1_MATCHED_MASS".equals(want)) {
            throw new IllegalArgumentException(
                    "A1_MATCHED_MASS is optional extra, not identity, skipped A1_BOUNDED_TRANSDUCER NOT_FIG4B");
        }

        if (arms.isEmpty()) {
            throw new IllegalArgumentException("no A1 arms selected: " + want);
        }

        boolean pass = true;
        for (ArmResult a : arms) {
            if ("C1_OPEN_DILUTE".equals(a.name)) {
                continue;
            }
            pass = pass && a.pass;
        }
        writeSummary(arms, pass);
        System.out.println(pass
                ? "A1=PASS A1_BOUNDED_TRANSDUCER HYPOTHETICAL_DESIGN_ENVELOPE NOT_FIG4B on C1c_FILLED_POCKET. "
                + "Not calibrated. A1C may not start (no lab u->J). A2 may not start. W0 later. "
                + "Not Fig. 4b. Not C1_OPEN_DILUTE."
                : "A1=FAIL A1_BOUNDED_TRANSDUCER HYPOTHETICAL_DESIGN_ENVELOPE NOT_FIG4B. "
                + "Do not retune Object B / tau_AC / K / n / J_max. Do not start A1C/A2/W0. "
                + "Do not place the AC on C1_OPEN_DILUTE.");
        if (!pass) {
            System.exit(1);
        }
    }

    private static boolean wantEmptyOr(String want, String name) {
        return want.isEmpty() || want.equals(name);
    }

    private static void requireFreeze() throws IOException {
        String freeze = Files.readString(CLAIM_FREEZE, StandardCharsets.UTF_8);
        if (!freeze.contains("NOT_FIG4B") || !freeze.contains("DaninoSI_OccupiedDF")) {
            throw new IllegalStateException("claim freeze must name Object B NOT_FIG4B");
        }
        if (!Files.readString(C1_STANDING, StandardCharsets.UTF_8).contains("**Status: FAIL**")) {
            throw new IllegalStateException("C1 standing must remain FAIL");
        }
        if (!Files.readString(C1C_STANDING, StandardCharsets.UTF_8).contains("**Status: PASS**")) {
            throw new IllegalStateException("C1c standing must remain PASS before A1");
        }
        if (!Files.readString(A0_STANDING, StandardCharsets.UTF_8).contains("**Status: PASS**")) {
            throw new IllegalStateException("A0 standing must remain PASS before A1");
        }
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("A1 protocol is not frozen_before_traces");
        }
        if (!text.contains("A1_BOUNDED_TRANSDUCER") || !text.contains("HYPOTHETICAL_DESIGN_ENVELOPE")) {
            throw new IllegalStateException("A1 protocol must freeze envelope label");
        }
        if (!json.contains("\"calibrated\": false") || json.contains("\"narma\": true")) {
            throw new IllegalStateException("A1 must not be calibrated or NARMA");
        }
        if (!json.contains("\"J_max\": 825.0") || !json.contains("\"tau_AC\": 2.0")) {
            throw new IllegalStateException("A1 J_max and tau_AC freeze missing");
        }
    }

    private static ArmResult runField(
            DaninoSIA1Field fieldDev,
            A1BoundedTransducer.Command command,
            String name,
            double tEnd,
            boolean contrast) throws IOException {
        long t0 = System.nanoTime();
        System.out.printf("  %s A1_BOUNDED_TRANSDUCER HYPOTHETICAL_DESIGN_ENVELOPE NOT_FIG4B command=%s%n",
                name, command.name());
        DaninoSIA1Field.Trajectory traj = fieldDev.run(
                command, tEnd, A1BoundedTransducer.SAMPLE_DT_FIELD, A1BoundedTransducer.RK_DT);
        writeFieldCsv(RESULTS.resolve("java_" + name + ".csv"), traj);
        boolean ledgerOk = traj.massReplayRel() <= LEDGER_REL
                && traj.payloadRel() <= LEDGER_REL
                && traj.transportRel() <= LEDGER_REL;
        boolean contrastOk = !contrast || (traj.a0ContrastRel() >= A0_CONTRAST && delayed(traj));
        boolean hotter = !contrast || acHotter(traj, command == A1BoundedTransducer.Command.U_STEP
                ? A1BoundedTransducer.T_OFF_STEP : A1BoundedTransducer.T_OFF_PULSE);
        boolean tail = command != A1BoundedTransducer.Command.U_PULSE || delayTail(traj);
        boolean armPass = ledgerOk && contrastOk && hotter && tail;
        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  %s A1_BOUNDED_TRANSDUCER NOT_FIG4B released=%.6g src=%.6g payload_drop=%.6g "
                        + "mass_rel=%.3g payload_rel=%.3g A0_contrast=%.3g J_end=%.6g wall_s=%.1f gate=%s%n",
                name, traj.ac.released(), traj.ledger.getSourceAdded(), traj.ac.payloadDrop(),
                traj.massReplayRel(), traj.payloadRel(), traj.a0ContrastRel(),
                traj.jS[traj.jS.length - 1], sec, armPass ? "PASS" : "FAIL");
        return new ArmResult(
                name, command.name(), "OFF", "FIELD", Double.NaN, Double.NaN,
                traj.massReplayRel(), traj.payloadRel(), traj.a0ContrastRel(),
                traj.ac.released(), traj.ledger.getSourceAdded(), traj.ac.payload(),
                Double.NaN, Double.NaN, armPass,
                "field A1_BOUNDED_TRANSDUCER HYPOTHETICAL_DESIGN_ENVELOPE NOT_FIG4B", null);
    }

    private static boolean delayed(DaninoSIA1Field.Trajectory traj) {
        for (int i = 1; i < traj.t.length; i++) {
            if (traj.t[i] > 0.0) {
                return traj.jS[i] < traj.jA0[i] - 1e-12;
            }
        }
        return false;
    }

    private static boolean delayTail(DaninoSIA1Field.Trajectory traj) {
        for (int i = 0; i < traj.t.length; i++) {
            if (traj.t[i] > A1BoundedTransducer.T_OFF_PULSE + 1e-12
                    && traj.jS[i] > 10.0 * A1BoundedTransducer.J_LEAK) {
                return true;
            }
        }
        return false;
    }

    private static boolean acHotter(DaninoSIA1Field.Trajectory traj, double tOff) {
        for (int i = 1; i < traj.t.length; i++) {
            if (traj.t[i] > tOff + 1e-12) {
                break;
            }
            if (!(traj.heAc[i] > traj.heFar[i])) {
                return false;
            }
        }
        return true;
    }

    private static ArmResult runLeak(DaninoSIA1Field fieldDev) throws IOException {
        long t0 = System.nanoTime();
        String name = "A1_U0_LEAK";
        System.out.printf("  %s A1_BOUNDED_TRANSDUCER HYPOTHETICAL_DESIGN_ENVELOPE NOT_FIG4B u=0%n", name);
        DaninoSIA1Field.Trajectory traj = fieldDev.run(
                A1BoundedTransducer.Command.U_ZERO,
                A1BoundedTransducer.T_END_FIELD,
                A1BoundedTransducer.SAMPLE_DT_FIELD,
                A1BoundedTransducer.RK_DT);
        writeFieldCsv(RESULTS.resolve("java_" + name + ".csv"), traj);
        double meanJ = 0.0;
        int n = 0;
        for (int i = 1; i < traj.jS.length; i++) {
            meanJ += traj.jS[i];
            n++;
        }
        meanJ /= n;
        boolean jOk = Math.abs(meanJ - A1BoundedTransducer.J_LEAK) <= 0.01;
        double expectDrop = A1BoundedTransducer.J_LEAK * A1BoundedTransducer.T_END_FIELD;
        boolean dropOk = Math.abs(traj.ac.payloadDrop() - expectDrop) / Math.max(expectDrop, 1e-15) <= LEDGER_REL;
        boolean ledgerOk = traj.massReplayRel() <= LEDGER_REL && traj.payloadRel() <= LEDGER_REL;
        boolean armPass = jOk && dropOk && ledgerOk;
        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  %s A1_BOUNDED_TRANSDUCER NOT_FIG4B meanJ=%.6g J_leak=%.6g drop=%.6g expect=%.6g "
                        + "mass_rel=%.3g wall_s=%.1f gate=%s%n",
                name, meanJ, A1BoundedTransducer.J_LEAK, traj.ac.payloadDrop(), expectDrop,
                traj.massReplayRel(), sec, armPass ? "PASS" : "FAIL");
        return new ArmResult(
                name, "U_ZERO", "OFF", "LEAK", Double.NaN, Double.NaN,
                traj.massReplayRel(), traj.payloadRel(), Double.NaN,
                traj.ac.released(), traj.ledger.getSourceAdded(), traj.ac.payload(),
                Double.NaN, Double.NaN, armPass,
                "u=0 basal leak A1_BOUNDED_TRANSDUCER NOT_FIG4B", null);
    }

    private static ArmResult runExhaust(DaninoSIA1Field fieldDev) throws IOException {
        long t0 = System.nanoTime();
        String name = "A1_PAYLOAD_EXHAUST";
        System.out.printf("  %s A1_BOUNDED_TRANSDUCER HYPOTHETICAL_DESIGN_ENVELOPE NOT_FIG4B long u=1%n", name);
        DaninoSIA1Field.Trajectory traj = fieldDev.run(
                A1BoundedTransducer.Command.U_EXHAUST,
                A1BoundedTransducer.T_END_EXHAUST,
                A1BoundedTransducer.SAMPLE_DT_FIELD,
                A1BoundedTransducer.RK_DT);
        writeFieldCsv(RESULTS.resolve("java_" + name + ".csv"), traj);
        double mFrac = traj.ac.payload() / A1BoundedTransducer.M0;
        double jMaxSeen = 0.0;
        for (double v : traj.jS) {
            jMaxSeen = Math.max(jMaxSeen, v);
        }
        double jEnd = traj.jS[traj.jS.length - 1];
        boolean empty = mFrac <= 0.10;
        boolean fell = jEnd <= 0.15 * jMaxSeen;
        boolean ledgerOk = traj.massReplayRel() <= LEDGER_REL && traj.payloadRel() <= LEDGER_REL;
        boolean noRefill = traj.ac.payload() <= traj.payload[0] + 1e-12;
        boolean armPass = empty && fell && ledgerOk && noRefill;
        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  %s A1_BOUNDED_TRANSDUCER NOT_FIG4B M_frac=%.4g J_end=%.6g J_max_seen=%.6g "
                        + "released=%.6g mass_rel=%.3g wall_s=%.1f gate=%s%n",
                name, mFrac, jEnd, jMaxSeen, traj.ac.released(), traj.massReplayRel(), sec,
                armPass ? "PASS" : "FAIL");
        return new ArmResult(
                name, "U_EXHAUST", "OFF", "EXHAUST", Double.NaN, Double.NaN,
                traj.massReplayRel(), traj.payloadRel(), Double.NaN,
                traj.ac.released(), traj.ledger.getSourceAdded(), traj.ac.payload(),
                Double.NaN, Double.NaN, armPass,
                "payload exhaust A1_BOUNDED_TRANSDUCER NOT_FIG4B", null);
    }

    private static ArmResult runJ0Bacteria(
            DaninoSIC1cFilledPocket device,
            DaninoSIPeriodCheck.Protocol peaks) throws IOException {
        long t0 = System.nanoTime();
        String name = "A1_J0_BACTERIA";
        System.out.println("  " + name + " A1_BOUNDED_TRANSDUCER NOT_FIG4B bacteria=ON AC=off");
        DaninoSIC1cFilledPocket.Probe probe = new DaninoSIC1cFilledPocket.Probe(
                A0IdealSource.I_AC, A0IdealSource.J_AC, A0IdealSource.K_AC,
                A0IdealSource.I_FAR, A0IdealSource.J_FAR, A0IdealSource.K_FAR);
        DaninoSIC1cFilledPocket.Trajectory traj = device.integrate(
                0.40, "AHL_KICK_005",
                A1BoundedTransducer.T_END_BACTERIA, A1BoundedTransducer.SAMPLE_DT_BACTERIA,
                A1BoundedTransducer.RK_DT,
                null, probe);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        double relT = Double.NaN;
        double relTC1c = Double.NaN;
        boolean periodOk = false;
        if ("OSC".equals(per.flag) && !Double.isNaN(per.period)) {
            relT = Math.abs(per.period - PY_T_040) / PY_T_040;
            relTC1c = Math.abs(per.period - C1C_T_040) / C1C_T_040;
            periodOk = relT <= PERIOD_REL_BAND;
        }
        double ahlRel = ahlResidualRel(traj);
        boolean armPass = periodOk && ahlRel <= LEDGER_REL && !traj.negativeState
                && Math.abs(traj.extra.getSourceAdded()) <= 1e-18;
        writeBacteriaCsv(RESULTS.resolve("java_" + name + ".csv"), traj, 0.0, A1BoundedTransducer.M0);
        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  %s A1_BOUNDED_TRANSDUCER NOT_FIG4B flag=%s T=%s relT_B=%s relT_C1c=%s "
                        + "ahl_rel=%.3g src=%.6g wall_s=%.1f gate=%s%n",
                name, per.flag, fmt(per.period), fmt(relT), fmt(relTC1c),
                ahlRel, traj.extra.getSourceAdded(), sec, armPass ? "PASS" : "FAIL");
        return new ArmResult(
                name, "AC_OFF", "ON", per.flag, per.period, relT,
                0.0, ahlRel, Double.NaN,
                0.0, traj.extra.getSourceAdded(), A1BoundedTransducer.M0,
                mean(traj.luxIMean()), mean(traj.iAc), armPass,
                "AC off C1c class A1_BOUNDED_TRANSDUCER NOT_FIG4B", traj);
    }

    private static ArmResult runPulseBacteria(
            DaninoSIC1cFilledPocket device,
            DaninoSIPeriodCheck.Protocol peaks,
            DaninoSIC1cFilledPocket.Trajectory j0) throws IOException {
        long t0 = System.nanoTime();
        String name = "A1_PULSE_BACTERIA";
        System.out.println("  " + name + " A1_BOUNDED_TRANSDUCER NOT_FIG4B bacteria=ON U_PULSE");
        A1BoundedTransducer ac = new A1BoundedTransducer(A1BoundedTransducer.Command.U_PULSE, false);
        DaninoSIC1cFilledPocket.Probe probe = new DaninoSIC1cFilledPocket.Probe(
                A0IdealSource.I_AC, A0IdealSource.J_AC, A0IdealSource.K_AC,
                A0IdealSource.I_FAR, A0IdealSource.J_FAR, A0IdealSource.K_FAR);
        DaninoSIC1cFilledPocket.Trajectory traj = device.integrate(
                0.40, "AHL_KICK_005",
                A1BoundedTransducer.T_END_BACTERIA, A1BoundedTransducer.SAMPLE_DT_BACTERIA,
                A1BoundedTransducer.RK_DT,
                ac::deposit, probe);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        int iOn = firstIndexAtLeast(traj.t, A1BoundedTransducer.T_OFF_PULSE);
        boolean heRise = iOn >= 0 && traj.heAc[iOn] > j0.heAc[iOn];
        double massRel = rel(ac.released(), traj.extra.getSourceAdded());
        double payloadRel = rel(ac.payloadDrop(), ac.released());
        double ahlRel = ahlResidualRel(traj);
        boolean ledgerOk = massRel <= LEDGER_REL && payloadRel <= LEDGER_REL && ahlRel <= LEDGER_REL;
        boolean armPass = heRise && ledgerOk && !traj.negativeState;
        writeBacteriaCsv(RESULTS.resolve("java_" + name + ".csv"), traj, ac.released(), ac.payload());
        double a0I = Double.NaN;
        double a0He = Double.NaN;
        if (Files.exists(A0_PULSE_CSV)) {
            List<String> lines = Files.readAllLines(A0_PULSE_CSV, StandardCharsets.UTF_8);
            if (lines.size() > 3) {
                String[] row = lines.get(3).split(",");
                a0He = Double.parseDouble(row[7]);
            }
            double s = 0.0;
            int n = 0;
            for (int i = 1; i < lines.size(); i++) {
                s += Double.parseDouble(lines.get(i).split(",")[2]);
                n++;
            }
            a0I = s / n;
        }
        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  %s A1_BOUNDED_TRANSDUCER NOT_FIG4B flag=%s T=%s He_ac_rise=%s "
                        + "I_mean=%.6g I_ac=%.6g vs_A0_I_mean=%s vs_A0_He_ac(t=1)=%s "
                        + "released=%.6g mass_rel=%.3g wall_s=%.1f gate=%s%n",
                name, per.flag, fmt(per.period), heRise ? "YES" : "NO",
                mean(traj.luxIMean()), mean(traj.iAc), fmt(a0I), fmt(a0He),
                ac.released(), massRel, sec, armPass ? "PASS" : "FAIL");
        return new ArmResult(
                name, "U_PULSE", "ON", per.flag, per.period, Double.NaN,
                massRel, ahlRel, Double.NaN,
                ac.released(), traj.extra.getSourceAdded(), ac.payload(),
                mean(traj.luxIMean()), mean(traj.iAc), armPass,
                "report I/He vs A0; not a task score A1_BOUNDED_TRANSDUCER NOT_FIG4B", traj);
    }

    private static ArmResult forbidOpenDilute() {
        System.out.println(
                "  C1_OPEN_DILUTE A1_BOUNDED_TRANSDUCER NOT_FIG4B forbidden; AC not placed");
        return new ArmResult(
                "C1_OPEN_DILUTE", "FORBIDDEN", "NA", "FORBIDDEN", Double.NaN, Double.NaN,
                Double.NaN, Double.NaN, Double.NaN,
                Double.NaN, Double.NaN, Double.NaN,
                Double.NaN, Double.NaN, true,
                "forbidden A1 device A1_BOUNDED_TRANSDUCER NOT_FIG4B", null);
    }

    private static double ahlResidualRel(DaninoSIC1cFilledPocket.Trajectory traj) {
        double r = traj.ledger.residual(traj.remainingMass, traj.extra) + traj.extra.getSourceAdded();
        double mStar = Math.max(
                Math.max(Math.abs(traj.remainingMass), Math.abs(traj.ledger.synthSource())),
                Math.max(Math.abs(traj.extra.getSourceAdded()), 1e-15));
        return Math.abs(r) / mStar;
    }

    private static double rel(double a, double b) {
        double mStar = Math.max(Math.max(Math.abs(a), Math.abs(b)),
                Math.max(A0IdealSource.DELTA_M_KICK, 1e-15));
        return Math.abs(a - b) / mStar;
    }

    private static int firstIndexAtLeast(double[] t, double value) {
        for (int i = 0; i < t.length; i++) {
            if (t[i] + 1e-12 >= value) {
                return i;
            }
        }
        return -1;
    }

    private static double mean(double[] x) {
        if (x == null || x.length == 0) {
            return Double.NaN;
        }
        double s = 0.0;
        for (double v : x) {
            s += v;
        }
        return s / x.length;
    }

    private static void writeFieldCsv(Path path, DaninoSIA1Field.Trajectory traj) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.println("t,J_S,J_A0,x1,M,He_ac,He_far,released,src,NOT_FIG4B,A1_BOUNDED_TRANSDUCER,HYPOTHETICAL_DESIGN_ENVELOPE");
            for (int i = 0; i < traj.t.length; i++) {
                w.printf(Locale.US,
                        "%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,NOT_FIG4B,A1_BOUNDED_TRANSDUCER,HYPOTHETICAL_DESIGN_ENVELOPE%n",
                        traj.t[i], traj.jS[i], traj.jA0[i], traj.x1[i], traj.payload[i],
                        traj.heAc[i], traj.heFar[i], traj.released[i], traj.src[i]);
            }
        }
    }

    private static void writeBacteriaCsv(
            Path path,
            DaninoSIC1cFilledPocket.Trajectory traj,
            double released,
            double payload) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.println("t,A_mean,I_mean,Hi_mean,He_pocket_mean,He_pocket_min,He_pocket_max,"
                    + "He_ac,He_far,I_ac,released,M,NOT_FIG4B,A1_BOUNDED_TRANSDUCER,HYPOTHETICAL_DESIGN_ENVELOPE");
            for (int i = 0; i < traj.t.length; i++) {
                w.printf(Locale.US,
                        "%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,"
                                + "NOT_FIG4B,A1_BOUNDED_TRANSDUCER,HYPOTHETICAL_DESIGN_ENVELOPE%n",
                        traj.t[i], traj.Y[i][0], traj.Y[i][1], traj.Y[i][2],
                        traj.Y[i][3], traj.hePocketMin[i], traj.hePocketMax[i],
                        traj.heAc[i], traj.heFar[i], traj.iAc[i], released, payload);
            }
        }
    }

    private static void writeSummary(List<ArmResult> arms, boolean pass) throws IOException {
        StringBuilder json = new StringBuilder();
        json.append("{\n");
        json.append("  \"gate\": \"A1\",\n");
        json.append("  \"status_label\": \"A1_BOUNDED_TRANSDUCER\",\n");
        json.append("  \"provenance\": \"HYPOTHETICAL_DESIGN_ENVELOPE\",\n");
        json.append("  \"A1\": \"").append(pass ? "PASS" : "FAIL").append("\",\n");
        json.append("  \"NOT_FIG4B\": true,\n");
        json.append("  \"calibrated\": false,\n");
        json.append("  \"device\": \"C1c_FILLED_POCKET\",\n");
        json.append("  \"C1\": \"FAIL\",\n");
        json.append("  \"A0\": \"PASS\",\n");
        json.append("  \"a1_on_open_dilute\": false,\n");
        json.append("  \"a1c_may_start\": false,\n");
        json.append("  \"a2_may_start\": false,\n");
        json.append("  \"w0_may_start\": false,\n");
        json.append("  \"J_max\": ").append(A1BoundedTransducer.J_MAX).append(",\n");
        json.append("  \"tau_AC\": ").append(A1BoundedTransducer.TAU_AC).append(",\n");
        json.append("  \"J_leak\": ").append(A1BoundedTransducer.J_LEAK).append(",\n");
        json.append("  \"M0\": ").append(A1BoundedTransducer.M0).append(",\n");
        json.append("  \"arms\": [\n");
        for (int i = 0; i < arms.size(); i++) {
            ArmResult a = arms.get(i);
            json.append("    {\"name\":\"").append(a.name).append("\",");
            json.append("\"flag\":\"").append(a.flag).append("\",");
            json.append("\"period\":").append(num(a.period)).append(",");
            json.append("\"period_rel_err\":").append(num(a.relPeriod)).append(",");
            json.append("\"mass_rel\":").append(num(a.massRel)).append(",");
            json.append("\"payload_rel\":").append(num(a.payloadRel)).append(",");
            json.append("\"a0_contrast\":").append(num(a.contrast)).append(",");
            json.append("\"released\":").append(num(a.released)).append(",");
            json.append("\"src_mass\":").append(num(a.src)).append(",");
            json.append("\"M\":").append(num(a.payload)).append(",");
            json.append("\"I_mean\":").append(num(a.iMean)).append(",");
            json.append("\"I_ac\":").append(num(a.iAc)).append(",");
            json.append("\"pass\":").append(a.pass).append(",");
            json.append("\"NOT_FIG4B\": true,\"A1_BOUNDED_TRANSDUCER\": true,");
            json.append("\"HYPOTHETICAL_DESIGN_ENVELOPE\": true}");
            json.append(i + 1 < arms.size() ? ",\n" : "\n");
        }
        json.append("  ]\n}\n");
        Files.writeString(RESULTS.resolve("a1_summary.json"), json.toString(), StandardCharsets.UTF_8);
    }

    private static String fmt(double x) {
        return Double.isNaN(x) ? "nan" : String.format(Locale.US, "%.6g", x);
    }

    private static String num(double x) {
        return Double.isNaN(x) ? "null" : Double.toString(x);
    }

    private static final class ArmResult {
        final String name;
        final String command;
        final String bacteria;
        final String flag;
        final double period;
        final double relPeriod;
        final double massRel;
        final double payloadRel;
        final double contrast;
        final double released;
        final double src;
        final double payload;
        final double iMean;
        final double iAc;
        final boolean pass;
        final String note;
        final DaninoSIC1cFilledPocket.Trajectory traj;

        ArmResult(String name, String command, String bacteria, String flag, double period, double relPeriod,
                  double massRel, double payloadRel, double contrast,
                  double released, double src, double payload,
                  double iMean, double iAc, boolean pass, String note,
                  DaninoSIC1cFilledPocket.Trajectory traj) {
            this.name = name;
            this.command = command;
            this.bacteria = bacteria;
            this.flag = flag;
            this.period = period;
            this.relPeriod = relPeriod;
            this.massRel = massRel;
            this.payloadRel = payloadRel;
            this.contrast = contrast;
            this.released = released;
            this.src = src;
            this.payload = payload;
            this.iMean = iMean;
            this.iAc = iAc;
            this.pass = pass;
            this.note = note;
            this.traj = traj;
        }
    }
}
