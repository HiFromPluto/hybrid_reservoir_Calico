package bsim.a0;

import bsim.BSimRandom;
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
 * Gate A0: ideal AC source on C1c_FILLED_POCKET. A0_IDEAL_SOURCE. NOT_FIG4B.
 * Not Fig. 4b. Not C1_OPEN_DILUTE. Not A1/A2/W0.
 */
public final class DaninoSIA0Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path CLAIM_FREEZE = ROOT.resolve("examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md");
    private static final Path C1_STANDING = ROOT.resolve("examples/PocketDish/C1_PACKED_SPATIAL_STANDING.md");
    private static final Path C1B_STANDING = ROOT.resolve("examples/PocketDish/C1B_ISLAND_STANDING.md");
    private static final Path C1C_STANDING = ROOT.resolve("examples/PocketDish/C1C_FILLED_POCKET_STANDING.md");
    private static final Path PROTOCOL = ROOT.resolve("examples/BSimReservoirPlanDaninoPocketA0/PROTOCOL.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/BSimReservoirPlanDaninoPocketA0/configs/protocol.json");
    private static final Path RESULTS = ROOT.resolve("examples/BSimReservoirPlanDaninoPocketA0/results");

    private static final double PY_T_040 = 63.583333333333336;
    private static final double C1C_T_040 = 61.25;
    private static final double PERIOD_REL_BAND = 0.20;
    private static final double LEDGER_REL = 1e-3;
    private static final double FIELD_HE_REL = 1e-9;
    private static final double CENTROID_ABS = 1e-9;

    private DaninoSIA0Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        DaninoSIPeriodCheck.refuseNarma(args);
        requireClaimFreeze();
        requireC1StillFail();
        requireC1bIslandPass();
        requireC1cFilledPass();
        requireProtocolFrozen();
        Files.createDirectories(RESULTS);

        BSimRandom rng = new BSimRandom(0L);
        DaninoSIA0Field fieldDev = new DaninoSIA0Field(rng);
        DaninoSIC1cFilledPocket device = new DaninoSIC1cFilledPocket(rng);
        if (device.unusedRng() != rng || fieldDev.unusedRng() != rng) {
            throw new IllegalStateException("RNG injection contract broken A0_IDEAL_SOURCE NOT_FIG4B");
        }

        System.out.println("A0_IDEAL_SOURCE on C1c_FILLED_POCKET. DaninoSI_OccupiedDF. NOT_FIG4B.");
        System.out.println("Not Fig. 4b. Not C1_OPEN_DILUTE. Not A1/A2/W0. Not Lentini TX-TL.");
        System.out.printf(
                Locale.US,
                "A0_IDEAL_SOURCE NOT_FIG4B device=C1c_FILLED_POCKET AC_xyz=(%.3f,%.3f,%.3f) "
                        + "AC_ijk=(%d,%d,%d) J_max=%.6g payload=extracellular_AHL%n",
                A0IdealSource.X_UM, A0IdealSource.Y_UM, A0IdealSource.Z_UM,
                A0IdealSource.I_AC, A0IdealSource.J_AC, A0IdealSource.K_AC,
                A0IdealSource.J_MAX);

        DaninoSIPeriodCheck.Protocol peaks = DaninoSIPeriodCheck.Protocol.d0b();
        List<ArmResult> arms = new ArrayList<>();
        String want = args.length > 0 ? args[0].trim() : "";
        DaninoSIC1cFilledPocket.Trajectory j0 = null;

        if (wantEmptyOr(want, "A0_FIELD_IMPULSE")) {
            arms.add(runFieldArm(
                    fieldDev, A0IdealSource.Command.U_PULSE, "A0_FIELD_IMPULSE", true));
        }
        if (wantEmptyOr(want, "A0_FIELD_STEP")) {
            arms.add(runFieldArm(
                    fieldDev, A0IdealSource.Command.U_STEP, "A0_FIELD_STEP", false));
        }
        if (wantEmptyOr(want, "A0_J0_BACTERIA")) {
            ArmResult r = runBacteria(device, peaks, A0IdealSource.Command.U_ZERO, "A0_J0_BACTERIA");
            arms.add(r);
            j0 = r.traj;
        }
        if (wantEmptyOr(want, "A0_J_PULSE_BACTERIA")) {
            if (j0 == null && want.isEmpty()) {
                throw new IllegalStateException("A0_J_PULSE_BACTERIA needs A0_J0_BACTERIA first A0_IDEAL_SOURCE NOT_FIG4B");
            }
            if (j0 == null) {
                System.out.println("  A0_J_PULSE_BACTERIA A0_IDEAL_SOURCE NOT_FIG4B running J=0 first for He-rise gate");
                ArmResult z = runBacteria(device, peaks, A0IdealSource.Command.U_ZERO, "A0_J0_BACTERIA");
                j0 = z.traj;
            }
            arms.add(runPulseBacteria(device, peaks, j0));
        }
        if (wantEmptyOr(want, "C1_OPEN_DILUTE") || want.isEmpty()) {
            arms.add(forbidOpenDilute());
        }

        if (arms.isEmpty()) {
            throw new IllegalArgumentException("no A0 arms selected: " + want);
        }

        boolean pass = true;
        for (ArmResult a : arms) {
            if ("C1_OPEN_DILUTE".equals(a.name)) {
                continue;
            }
            pass = pass && a.pass;
        }

        writeSummary(device, arms, pass);
        System.out.println(pass
                ? "A0=PASS A0_IDEAL_SOURCE NOT_FIG4B on C1c_FILLED_POCKET. "
                + "C1 open chip remains FAIL. A1 may start only on C1c_FILLED_POCKET, "
                + "still NOT_FIG4B, still not calibrated. Not Fig. 4b. Not C1_OPEN_DILUTE. W0 later."
                : "A0=FAIL A0_IDEAL_SOURCE NOT_FIG4B. Do not retune Object B / J_max / D1_spatial. "
                + "Do not start A1/A2/W0. Do not place the AC on C1_OPEN_DILUTE.");
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
    }

    private static void requireC1cFilledPass() throws IOException {
        String standing = Files.readString(C1C_STANDING, StandardCharsets.UTF_8);
        if (!standing.contains("**Status: PASS**")) {
            throw new IllegalStateException("C1c filled-pocket standing must remain PASS before A0");
        }
        if (!standing.contains("C1c_FILLED_POCKET") || !standing.contains("NOT_FIG4B")) {
            throw new IllegalStateException("C1c standing must name C1c_FILLED_POCKET NOT_FIG4B");
        }
    }

    private static void requireProtocolFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces") || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("A0 protocol is not frozen_before_traces");
        }
        if (!text.contains("A0_IDEAL_SOURCE") || !text.contains("NOT_FIG4B")) {
            throw new IllegalStateException("A0 protocol must freeze A0_IDEAL_SOURCE NOT_FIG4B");
        }
        if (!text.contains("C1c_FILLED_POCKET") || !json.contains("\"device\": \"C1c_FILLED_POCKET\"")) {
            throw new IllegalStateException("A0 protocol must freeze device C1c_FILLED_POCKET");
        }
        if (!json.contains("\"J_max\": 825.0") || !json.contains("\"ac_xyz_um\": [50.0, 50.0, 0.825]")) {
            throw new IllegalStateException("A0 protocol must freeze J_max=825 and AC xyz");
        }
        if (json.contains("\"TIME_ADJ\": true") || json.contains("\"narma\": true")) {
            throw new IllegalStateException("A0 protocol must freeze TIME_ADJ unused and NARMA off");
        }
        if (!json.contains("\"a0_on_open_dilute\": false")) {
            throw new IllegalStateException("A0 must not be declared on C1_OPEN_DILUTE");
        }
    }

    private static ArmResult runFieldArm(
            DaninoSIA0Field fieldDev,
            A0IdealSource.Command command,
            String name,
            boolean impulse) throws IOException {
        long t0 = System.nanoTime();
        System.out.printf("  %s A0_IDEAL_SOURCE NOT_FIG4B bacteria=OFF command=%s%n", name, command.name());
        DaninoSIA0Field.Trajectory a0 = fieldDev.run(
                command, true, A0IdealSource.T_END_FIELD, A0IdealSource.SAMPLE_DT_FIELD, A0IdealSource.RK_DT);
        DaninoSIA0Field.Trajectory n0 = fieldDev.run(
                command, false, A0IdealSource.T_END_FIELD, A0IdealSource.SAMPLE_DT_FIELD, A0IdealSource.RK_DT);
        writeFieldCsv(RESULTS.resolve("java_" + name + ".csv"), a0);
        writeFieldCsv(RESULTS.resolve("java_" + name + "_N0_ORACLE.csv"), n0);

        double heRel = maxHeRel(a0, n0);
        double cAbs = maxCentroidAbs(a0, n0);
        boolean massOk = a0.massReplayRel() <= LEDGER_REL && n0.massReplayRel() <= LEDGER_REL;
        boolean transportOk = a0.transportRel() <= LEDGER_REL && n0.transportRel() <= LEDGER_REL;
        boolean matchOk = heRel <= FIELD_HE_REL && cAbs <= CENTROID_ABS;
        boolean extraOk = impulse ? impulseExtra(a0) : stepExtra(a0);
        boolean armPass = massOk && transportOk && matchOk && extraOk;

        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  %s A0_IDEAL_SOURCE NOT_FIG4B cmd=%.6g src=%.6g mass_rel=%.3g transport_rel=%.3g "
                        + "He_vs_N0_rel=%.3g centroid_abs=%.3g extra=%s wall_s=%.1f gate=%s%n",
                name, a0.commandedMass, a0.ledger.getSourceAdded(), a0.massReplayRel(),
                a0.transportRel(), heRel, cAbs, extraOk ? "PASS" : "FAIL", sec,
                armPass ? "PASS" : "FAIL");
        return new ArmResult(
                name, command.name(), "OFF", "FIELD", Double.NaN, Double.NaN,
                a0.massReplayRel(), a0.transportRel(), heRel, cAbs,
                a0.commandedMass, a0.ledger.getSourceAdded(),
                Double.NaN, Double.NaN, Double.NaN,
                armPass, "bacteria OFF field vs N0 oracle A0_IDEAL_SOURCE NOT_FIG4B", null);
    }

    private static boolean impulseExtra(DaninoSIA0Field.Trajectory a0) {
        int iOn = firstIndexAtLeast(a0.t, A0IdealSource.T_OFF_PULSE);
        if (iOn < 0 || iOn >= a0.t.length - 1) {
            return false;
        }
        boolean falls = a0.heAc[a0.heAc.length - 1] < a0.heAc[iOn] - 1e-18;
        boolean farRises = false;
        for (int i = 1; i < a0.heFar.length; i++) {
            if (a0.heFar[i] > a0.heFar[0] + 1e-18) {
                farRises = true;
                break;
            }
        }
        return falls && farRises;
    }

    private static boolean stepExtra(DaninoSIA0Field.Trajectory a0) {
        boolean hotter = true;
        for (int i = 1; i < a0.t.length; i++) {
            if (a0.t[i] > A0IdealSource.T_OFF_STEP + 1e-12) {
                break;
            }
            if (!(a0.heAc[i] > a0.heFar[i])) {
                hotter = false;
                break;
            }
        }
        double maxAc = max(a0.heAc);
        double maxFar = max(a0.heFar);
        double tAc = timeToFraction(a0.t, a0.heAc, 0.10 * maxAc);
        double tFar = timeToFraction(a0.t, a0.heFar, 0.10 * maxFar);
        return hotter && tAc < tFar;
    }

    private static ArmResult runBacteria(
            DaninoSIC1cFilledPocket device,
            DaninoSIPeriodCheck.Protocol peaks,
            A0IdealSource.Command command,
            String name) throws IOException {
        long t0 = System.nanoTime();
        System.out.printf("  %s A0_IDEAL_SOURCE NOT_FIG4B bacteria=ON command=%s%n", name, command.name());
        A0IdealSource ac = new A0IdealSource(command);
        DaninoSIC1cFilledPocket.Probe probe = new DaninoSIC1cFilledPocket.Probe(
                A0IdealSource.I_AC, A0IdealSource.J_AC, A0IdealSource.K_AC,
                A0IdealSource.I_FAR, A0IdealSource.J_FAR, A0IdealSource.K_FAR);
        DaninoSIC1cFilledPocket.Trajectory traj = device.integrate(
                0.40, "AHL_KICK_005",
                A0IdealSource.T_END_BACTERIA, A0IdealSource.SAMPLE_DT_BACTERIA, A0IdealSource.RK_DT,
                ac::deposit, probe);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        if (traj.negativeState) {
            per = new DaninoSIPeriodCheck.Result(
                    "NEGATIVE_STATE", per.nPeaks, per.nPersistPeaks, Double.NaN,
                    per.peakTimes, per.relAmplitude, per.ampPersist, per.meanI, per.p5I, per.p95I);
        }
        double relT = Double.NaN;
        double relTC1c = Double.NaN;
        boolean periodOk = false;
        if ("OSC".equals(per.flag) && !Double.isNaN(per.period)) {
            relT = Math.abs(per.period - PY_T_040) / PY_T_040;
            relTC1c = Math.abs(per.period - C1C_T_040) / C1C_T_040;
            periodOk = relT <= PERIOD_REL_BAND;
        }
        double massRel = massReplayRel(ac.commandedMass(), traj.extra.getSourceAdded());
        double ahlRel = ahlResidualRel(traj);
        boolean massOk = massRel <= LEDGER_REL && ahlRel <= LEDGER_REL;
        boolean armPass = periodOk && massOk && !traj.negativeState
                && Math.abs(ac.commandedMass()) <= 1e-18;
        writeBacteriaCsv(RESULTS.resolve("java_" + name + ".csv"), traj, ac);
        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  %s A0_IDEAL_SOURCE NOT_FIG4B flag=%s T=%s relT_B=%s relT_C1c=%s "
                        + "mass_rel=%.3g ahl_rel=%.3g cmd=%.6g wall_s=%.1f gate=%s%n",
                name, per.flag, fmt(per.period), fmt(relT), fmt(relTC1c),
                massRel, ahlRel, ac.commandedMass(), sec,
                armPass ? "PASS" : "FAIL");
        ArmResult r = new ArmResult(
                name, command.name(), "ON", per.flag, per.period, relT,
                massRel, ahlRel, Double.NaN, Double.NaN,
                ac.commandedMass(), traj.extra.getSourceAdded(),
                mean(traj.luxIMean()), mean(traj.iAc), relTC1c,
                armPass, "J=0 Object B C1c class A0_IDEAL_SOURCE NOT_FIG4B", traj);
        return r;
    }

    private static ArmResult runPulseBacteria(
            DaninoSIC1cFilledPocket device,
            DaninoSIPeriodCheck.Protocol peaks,
            DaninoSIC1cFilledPocket.Trajectory j0) throws IOException {
        long t0 = System.nanoTime();
        String name = "A0_J_PULSE_BACTERIA";
        System.out.printf("  %s A0_IDEAL_SOURCE NOT_FIG4B bacteria=ON command=U_PULSE%n", name);
        A0IdealSource ac = new A0IdealSource(A0IdealSource.Command.U_PULSE);
        DaninoSIC1cFilledPocket.Probe probe = new DaninoSIC1cFilledPocket.Probe(
                A0IdealSource.I_AC, A0IdealSource.J_AC, A0IdealSource.K_AC,
                A0IdealSource.I_FAR, A0IdealSource.J_FAR, A0IdealSource.K_FAR);
        DaninoSIC1cFilledPocket.Trajectory traj = device.integrate(
                0.40, "AHL_KICK_005",
                A0IdealSource.T_END_BACTERIA, A0IdealSource.SAMPLE_DT_BACTERIA, A0IdealSource.RK_DT,
                ac::deposit, probe);
        DaninoSIPeriodCheck.Result per = DaninoSIPeriodCheck.extract(traj.t, traj.luxIMean(), peaks);
        int iOn = firstIndexAtLeast(traj.t, A0IdealSource.T_OFF_PULSE);
        boolean heRise = iOn >= 0 && traj.heAc[iOn] > j0.heAc[iOn];
        double massRel = massReplayRel(ac.commandedMass(), traj.extra.getSourceAdded());
        double ahlRel = ahlResidualRel(traj);
        boolean massOk = massRel <= LEDGER_REL && ahlRel <= LEDGER_REL
                && Math.abs(ac.commandedMass() - A0IdealSource.DELTA_M_KICK) <= 1e-9;
        boolean armPass = heRise && massOk && !traj.negativeState;
        writeBacteriaCsv(RESULTS.resolve("java_" + name + ".csv"), traj, ac);
        double sec = (System.nanoTime() - t0) / 1e9;
        System.out.printf(
                Locale.US,
                "  %s A0_IDEAL_SOURCE NOT_FIG4B flag=%s T=%s He_ac_rise=%s I_mean=%.6g I_ac=%.6g "
                        + "mass_rel=%.3g ahl_rel=%.3g cmd=%.6g wall_s=%.1f gate=%s%n",
                name, per.flag, fmt(per.period), heRise ? "YES" : "NO",
                mean(traj.luxIMean()), mean(traj.iAc),
                massRel, ahlRel, ac.commandedMass(), sec,
                armPass ? "PASS" : "FAIL");
        return new ArmResult(
                name, "U_PULSE", "ON", per.flag, per.period, Double.NaN,
                massRel, ahlRel, Double.NaN, Double.NaN,
                ac.commandedMass(), traj.extra.getSourceAdded(),
                mean(traj.luxIMean()), mean(traj.iAc), Double.NaN,
                armPass, "pulse vs J=0 reported I; He near AC rises A0_IDEAL_SOURCE NOT_FIG4B", traj);
    }

    private static ArmResult forbidOpenDilute() {
        System.out.println(
                "  C1_OPEN_DILUTE A0_IDEAL_SOURCE NOT_FIG4B forbidden; AC not placed; not an identity arm");
        return new ArmResult(
                "C1_OPEN_DILUTE", "FORBIDDEN", "NA", "FORBIDDEN", Double.NaN, Double.NaN,
                Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                true, "forbidden A0 device; not run A0_IDEAL_SOURCE NOT_FIG4B", null);
    }

    private static double ahlResidualRel(DaninoSIC1cFilledPocket.Trajectory traj) {
        double r = traj.ledger.residual(traj.remainingMass, traj.extra) + traj.extra.getSourceAdded();
        double mStar = Math.max(
                Math.max(Math.abs(traj.remainingMass), Math.abs(traj.ledger.synthSource())),
                Math.max(Math.abs(traj.extra.getSourceAdded()), 1e-15));
        return Math.abs(r) / mStar;
    }

    private static double massReplayRel(double cmd, double src) {
        double mStar = Math.max(Math.max(Math.abs(cmd), Math.abs(src)),
                Math.max(A0IdealSource.DELTA_M_KICK, 1e-15));
        return Math.abs(cmd - src) / mStar;
    }

    private static double maxHeRel(DaninoSIA0Field.Trajectory a, DaninoSIA0Field.Trajectory b) {
        double m = 0.0;
        for (int i = 0; i < a.t.length; i++) {
            m = Math.max(m, rel(a.heAc[i], b.heAc[i]));
            m = Math.max(m, rel(a.heFar[i], b.heFar[i]));
        }
        return m;
    }

    private static double maxCentroidAbs(DaninoSIA0Field.Trajectory a, DaninoSIA0Field.Trajectory b) {
        double m = 0.0;
        for (int i = 0; i < a.t.length; i++) {
            if (Double.isNaN(a.cx[i]) && Double.isNaN(b.cx[i])) {
                continue;
            }
            m = Math.max(m, Math.hypot(a.cx[i] - b.cx[i], a.cy[i] - b.cy[i]));
        }
        return m;
    }

    private static double rel(double a, double b) {
        return Math.abs(a - b) / Math.max(Math.abs(b), 1e-15);
    }

    private static int firstIndexAtLeast(double[] t, double value) {
        for (int i = 0; i < t.length; i++) {
            if (t[i] + 1e-12 >= value) {
                return i;
            }
        }
        return -1;
    }

    private static double max(double[] x) {
        double m = Double.NEGATIVE_INFINITY;
        for (double v : x) {
            m = Math.max(m, v);
        }
        return m;
    }

    private static double timeToFraction(double[] t, double[] y, double thresh) {
        for (int i = 0; i < t.length; i++) {
            if (y[i] >= thresh) {
                return t[i];
            }
        }
        return Double.POSITIVE_INFINITY;
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

    private static void writeFieldCsv(Path path, DaninoSIA0Field.Trajectory traj) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.println("t,He_ac,He_far,He_mean,cx,cy,cmd_mass,src_mass,NOT_FIG4B,A0_IDEAL_SOURCE");
            for (int i = 0; i < traj.t.length; i++) {
                w.printf(Locale.US,
                        "%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,NOT_FIG4B,A0_IDEAL_SOURCE%n",
                        traj.t[i], traj.heAc[i], traj.heFar[i], traj.heMean[i],
                        traj.cx[i], traj.cy[i], traj.cmdMass[i], traj.srcMass[i]);
            }
        }
    }

    private static void writeBacteriaCsv(
            Path path,
            DaninoSIC1cFilledPocket.Trajectory traj,
            A0IdealSource ac) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path, StandardCharsets.UTF_8))) {
            w.println("t,A_mean,I_mean,Hi_mean,He_pocket_mean,He_pocket_min,He_pocket_max,"
                    + "He_ac,He_far,I_ac,cmd_mass,NOT_FIG4B,A0_IDEAL_SOURCE");
            for (int i = 0; i < traj.t.length; i++) {
                w.printf(Locale.US,
                        "%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,NOT_FIG4B,A0_IDEAL_SOURCE%n",
                        traj.t[i], traj.Y[i][0], traj.Y[i][1], traj.Y[i][2],
                        traj.Y[i][3], traj.hePocketMin[i], traj.hePocketMax[i],
                        traj.heAc[i], traj.heFar[i], traj.iAc[i], ac.commandedMass());
            }
        }
    }

    private static void writeSummary(
            DaninoSIC1cFilledPocket device,
            List<ArmResult> arms,
            boolean pass) throws IOException {
        StringBuilder json = new StringBuilder();
        json.append("{\n");
        json.append("  \"gate\": \"A0\",\n");
        json.append("  \"status_label\": \"A0_IDEAL_SOURCE\",\n");
        json.append("  \"object\": \"DaninoSI_OccupiedDF\",\n");
        json.append("  \"A0\": \"").append(pass ? "PASS" : "FAIL").append("\",\n");
        json.append("  \"NOT_FIG4B\": true,\n");
        json.append("  \"device\": \"C1c_FILLED_POCKET\",\n");
        json.append("  \"C1\": \"FAIL\",\n");
        json.append("  \"C1b\": \"PASS_ISLAND_ONLY\",\n");
        json.append("  \"C1c\": \"PASS\",\n");
        json.append("  \"a0_on_open_dilute\": false,\n");
        json.append("  \"a1_may_start\": ").append(pass).append(",\n");
        json.append("  \"a1_calibrated\": false,\n");
        json.append("  \"a2_may_start\": false,\n");
        json.append("  \"w0_may_start\": false,\n");
        json.append("  \"J_max\": ").append(A0IdealSource.J_MAX).append(",\n");
        json.append("  \"ac_xyz_um\": [50.0, 50.0, 0.825],\n");
        json.append("  \"ac_voxel_ijk\": [25, 25, 0],\n");
        json.append("  \"payload\": \"extracellular_AHL\",\n");
        json.append("  \"TIME_ADJ\": false,\n");
        json.append("  \"n_cells\": ").append(device.nCells()).append(",\n");
        json.append("  \"D1_spatial\": ").append(device.d1Spatial()).append(",\n");
        json.append("  \"arms\": [\n");
        for (int i = 0; i < arms.size(); i++) {
            ArmResult a = arms.get(i);
            json.append("    {\"name\":\"").append(a.name).append("\",");
            json.append("\"flag\":\"").append(a.flag).append("\",");
            json.append("\"period\":").append(num(a.period)).append(",");
            json.append("\"period_rel_err\":").append(num(a.relPeriod)).append(",");
            json.append("\"mass_rel\":").append(num(a.massRel)).append(",");
            json.append("\"ahl_or_transport_rel\":").append(num(a.transportRel)).append(",");
            json.append("\"he_vs_n0_rel\":").append(num(a.heRel)).append(",");
            json.append("\"cmd_mass\":").append(num(a.cmd)).append(",");
            json.append("\"src_mass\":").append(num(a.src)).append(",");
            json.append("\"I_mean\":").append(num(a.iMean)).append(",");
            json.append("\"I_ac\":").append(num(a.iAc)).append(",");
            json.append("\"pass\":").append(a.pass).append(",");
            json.append("\"NOT_FIG4B\": true,\"A0_IDEAL_SOURCE\": true}");
            json.append(i + 1 < arms.size() ? ",\n" : "\n");
        }
        json.append("  ]\n}\n");
        Files.writeString(RESULTS.resolve("a0_summary.json"), json.toString(), StandardCharsets.UTF_8);
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
        final double transportRel;
        final double heRel;
        final double centroidAbs;
        final double cmd;
        final double src;
        final double iMean;
        final double iAc;
        final double relTC1c;
        final boolean pass;
        final String note;
        final DaninoSIC1cFilledPocket.Trajectory traj;

        ArmResult(String name, String command, String bacteria, String flag, double period, double relPeriod,
                  double massRel, double transportRel, double heRel, double centroidAbs,
                  double cmd, double src, double iMean, double iAc, double relTC1c,
                  boolean pass, String note, DaninoSIC1cFilledPocket.Trajectory traj) {
            this.name = name;
            this.command = command;
            this.bacteria = bacteria;
            this.flag = flag;
            this.period = period;
            this.relPeriod = relPeriod;
            this.massRel = massRel;
            this.transportRel = transportRel;
            this.heRel = heRel;
            this.centroidAbs = centroidAbs;
            this.cmd = cmd;
            this.src = src;
            this.iMean = iMean;
            this.iAc = iAc;
            this.relTC1c = relTC1c;
            this.pass = pass;
            this.note = note;
            this.traj = traj;
        }
    }
}
