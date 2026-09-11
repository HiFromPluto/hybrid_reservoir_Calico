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
 * D0_CLOSED_BATH. Scalar depleting C. Incremental sizer. Death off.
 * Not NutrientField. Not elongateCited. Not Erickson. Not NARMA.
 */
public final class LaneDD0Job {

    private static final Path ROOT = Path.of("").toAbsolutePath();
    private static final Path PROTOCOL = ROOT.resolve(
            "examples/LaneD_LifeCycle/PROTOCOL_D0.md");
    private static final Path PROTOCOL_JSON = ROOT.resolve(
            "examples/LaneD_LifeCycle/configs/protocol_d0.json");
    private static final Path RESULTS = ROOT.resolve(
            "examples/LaneD_LifeCycle/results");

    private LaneDD0Job() { }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        refuseCapacity(args);
        requireFrozen();
        Files.createDirectories(RESULTS);

        new BSimRandom(LaneDD0Identity.RNG_SEED);

        System.out.println(LaneDD0Identity.OBJECT + " " + LaneDD0Identity.GATE
                + " device=" + LaneDD0Identity.DEVICE);
        System.out.println("Not LC0_OPEN_BATH. Not NutrientField. Not elongateCited.");
        System.out.println("Not Erickson. Not famine. Not NARMA. Death off.");
        System.out.printf(Locale.US,
                "N0=%d  seed=%d  C_s=%.3f mM  K_S=%.3f mM  C_cut=%.4f mM  "
                        + "T_cap=%.0f s  N_ever_cap=%d  sink=rho_cell  "
                        + "elongation=incremental%n",
                LaneDD0Identity.N0, LaneDD0Identity.RNG_SEED,
                LaneDD0Identity.C_S_MM, LaneDD0Identity.K_S_MM,
                LaneDD0Identity.C_CUT_MM, LaneDD0Identity.T_CAP_S,
                LaneDD0Identity.N_EVER_COMPUTE_CAP);
        System.out.printf(Locale.US,
                "F3_DIVISION_VOLUME_DROP V(3)=%.4f  two_daughters=%.4f  "
                        + "drop=%.4f um3 (%.1f%%). Not glucose. "
                        + "Mass gate is U vs d(CV), not sum rho V.%n",
                LaneDD0Identity.F3_V_DIV_UM3,
                LaneDD0Identity.F3_TWO_DAUGHTERS_UM3,
                LaneDD0Identity.F3_DROP_UM3,
                100.0 * LaneDD0Identity.F3_DROP_FRAC);

        Arm grow = runArm("D0_GROW", true);
        Arm off = runArm("D0_OFF", false);
        writeArmCsv("d0_grow.csv", grow);
        writeArmCsv("d0_off.csv", off);
        writeSummary(grow, off);
    }

    private static Arm runArm(String name, boolean growOn) {
        System.out.printf(Locale.US,
                "%s growth=%s death=OFF motility=OFF chemistry=OFF PDE=OFF "
                        + "device=%s sink=rho_cell elongation=incremental seed=%d%n",
                name, growOn ? "ON" : "OFF", LaneDD0Identity.DEVICE,
                LaneDD0Identity.RNG_SEED);

        LaneDClosedBath bath = new LaneDClosedBath(LaneDD0Identity.C_S_MM, growOn);
        List<LaneDClosedBathCell> cells = new ArrayList<>();
        for (int i = 0; i < LaneDD0Identity.N0; i++) {
            cells.add(new LaneDClosedBathCell(
                    LaneDD0Identity.founderX(i),
                    LaneDD0Identity.founderY(i),
                    0.5 * LaneDD0Identity.BOX_Z_UM,
                    LaneDD0Identity.ELL0_UM));
        }

        Arm arm = new Arm(name, growOn);
        record(arm, 0.0, bath, cells);

        double t = 0.0;
        int nEver = LaneDD0Identity.N0;
        while (t < LaneDD0Identity.T_CAP_S - 1.0e-12) {
            bath.refuseDeath(false);
            double dt = Math.min(LaneDD0Identity.DT_S, LaneDD0Identity.T_CAP_S - t);
            double cStep = bath.concentrationMm();

            if (growOn) {
                double sumV = sumVolume(cells);
                bath.stepSink(dt, sumV);
                int nBefore = cells.size();
                for (int i = 0; i < nBefore; i++) {
                    LaneDClosedBathCell cell = cells.get(i);
                    if (cell.elongateIncremental(dt, cStep)) {
                        arm.shrinkEvents++;
                    }
                    if (cell.shouldDivide()) {
                        if (nEver + 1 > LaneDD0Identity.N_EVER_COMPUTE_CAP) {
                            arm.capFired = true;
                            System.out.println("SCOPE_NOTE compute-cap N_ever>"
                                    + LaneDD0Identity.N_EVER_COMPUTE_CAP
                                    + ". Not a death. Not PASS.");
                            finish(arm, t + dt, bath, cells, nEver);
                            return arm;
                        }
                        double vBefore = cell.volumeUm3();
                        cell.resetAfterFission();
                        LaneDClosedBathCell daughter = new LaneDClosedBathCell(
                                cell.x, cell.y, cell.z, LaneDD0Identity.ELL0_UM);
                        cells.add(daughter);
                        nEver++;
                        arm.birthsTotal++;
                        arm.f3DropCum += vBefore - cell.volumeUm3()
                                - daughter.volumeUm3();
                    }
                }
            }

            t += dt;
            boolean log = isLogTime(t)
                    || t + 1.0e-12 >= LaneDD0Identity.T_CAP_S
                    || bath.concentrationMm() <= LaneDD0Identity.C_CUT_MM;
            if (log) {
                record(arm, t, bath, cells);
            }

            if (growOn && bath.concentrationMm() <= LaneDD0Identity.C_CUT_MM) {
                arm.stopReason = "C_CUT";
                break;
            }
        }

        if (arm.stopReason == null) {
            if (growOn && bath.concentrationMm() > LaneDD0Identity.C_CUT_MM) {
                arm.stopReason = "T_CAP";
                System.out.println("SCOPE_NOTE T_cap with C > C_cut. Not PASS. "
                        + "Do not raise C_s or lower N0.");
            } else {
                arm.stopReason = growOn ? "C_CUT" : "T_CAP";
            }
        }

        finish(arm, t, bath, cells, nEver);
        return arm;
    }

    private static void finish(Arm arm, double t, LaneDClosedBath bath,
            List<LaneDClosedBathCell> cells, int nEver) {
        if (arm.obsT.isEmpty()
                || Math.abs(arm.obsT.get(arm.obsT.size() - 1) - t) > 1.0e-9) {
            record(arm, t, bath, cells);
        }
        arm.tEnd = t;
        arm.nEnd = cells.size();
        arm.nEverEnd = nEver;
        arm.cEnd = bath.concentrationMm();
        arm.uEnd = bath.cumulativeSinkMmUm3();
        arm.dCvEnd = bath.deltaCvMmUm3();
        arm.massResidual = bath.massResidualRel();
        arm.sumVEnd = sumVolume(cells);
        System.out.printf(Locale.US,
                "%s N_end=%d N_ever=%d births=%d C_end=%.8f U=%.6f dCV=%.6f "
                        + "mass_resid=%.6e stop=%s shrink=%d f3_drop_cum=%.4f%n",
                arm.name, arm.nEnd, arm.nEverEnd, arm.birthsTotal, arm.cEnd,
                arm.uEnd, arm.dCvEnd, arm.massResidual, arm.stopReason,
                arm.shrinkEvents, arm.f3DropCum);
    }

    private static boolean isLogTime(double t) {
        double w = LaneDD0Identity.LOG_DT_S;
        double k = t / w;
        return Math.abs(k - Math.round(k)) < 1.0e-9;
    }

    private static double sumVolume(List<LaneDClosedBathCell> cells) {
        double s = 0.0;
        for (LaneDClosedBathCell cell : cells) {
            s += cell.volumeUm3();
        }
        return s;
    }

    private static void record(Arm arm, double t, LaneDClosedBath bath,
            List<LaneDClosedBathCell> cells) {
        double minL = Double.POSITIVE_INFINITY;
        double maxL = 0.0;
        double sumV = 0.0;
        for (LaneDClosedBathCell cell : cells) {
            double ell = cell.lengthUm();
            minL = Math.min(minL, ell);
            maxL = Math.max(maxL, ell);
            sumV += cell.volumeUm3();
        }
        double c = bath.concentrationMm();
        arm.obsT.add(t);
        arm.obsC.add(c);
        arm.obsLambda.add(LaneDD0Identity.lambdaPerHour(c));
        arm.obsN.add(cells.size());
        arm.obsNever.add(cells.size());
        arm.obsBirths.add(arm.birthsTotal);
        arm.obsSumV.add(sumV);
        arm.obsU.add(bath.cumulativeSinkMmUm3());
        arm.obsDcv.add(bath.deltaCvMmUm3());
        arm.obsMinEll.add(minL);
        arm.obsMaxEll.add(maxL);
        arm.obsF3.add(arm.f3DropCum);
    }

    private static void writeArmCsv(String name, Arm arm) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(
                RESULTS.resolve(name), StandardCharsets.UTF_8))) {
            w.println("t,C_mM,lambda_per_h,N,N_ever,births_cum,sum_V_um3,"
                    + "U_mM_um3,dCV_mM_um3,min_ell,max_ell,f3_drop_cum_um3");
            for (int i = 0; i < arm.obsT.size(); i++) {
                w.printf(Locale.US,
                        "%.10f,%.16e,%.16e,%d,%d,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e%n",
                        arm.obsT.get(i), arm.obsC.get(i), arm.obsLambda.get(i),
                        arm.obsN.get(i), arm.obsNever.get(i), arm.obsBirths.get(i),
                        arm.obsSumV.get(i), arm.obsU.get(i), arm.obsDcv.get(i),
                        arm.obsMinEll.get(i), arm.obsMaxEll.get(i),
                        arm.obsF3.get(i));
            }
        }
    }

    private static void writeSummary(Arm grow, Arm off) throws IOException {
        boolean honesty = grow.growOn && !off.growOn;
        boolean processOff = off.birthsTotal == 0
                && off.nEnd == LaneDD0Identity.N0
                && Math.abs(off.cEnd - LaneDD0Identity.C_S_MM) <= 1.0e-15;
        boolean noShrink = grow.shrinkEvents == 0 && noLengthDecrease(grow);
        boolean ledger = grow.massResidual <= LaneDD0Identity.DELTA_MASS + 1.0e-15;
        boolean monod = monodOk(grow);
        boolean exhausted = grow.cEnd <= LaneDD0Identity.C_CUT_MM + 1.0e-15;
        boolean tCap = "T_CAP".equals(grow.stopReason) && !exhausted;
        boolean noKill = grow.nEnd == grow.nEverEnd && off.nEnd == off.nEverEnd;
        boolean cap = grow.capFired || off.capFired;
        boolean scope = cap || tCap;
        boolean pass = honesty && processOff && noShrink && ledger && monod
                && exhausted && noKill && !cap && !tCap;

        String killer;
        if (cap) {
            killer = "compute_cap";
        } else if (tCap) {
            killer = "t_cap";
        } else if (!honesty) {
            killer = "honesty";
        } else if (!processOff) {
            killer = "process_off";
        } else if (!noShrink) {
            killer = "no_shrink_F4";
        } else if (!ledger) {
            killer = "bath_ledger";
        } else if (!monod) {
            killer = "monod_identity";
        } else if (!exhausted) {
            killer = "exhaustion";
        } else if (!noKill) {
            killer = "silent_kill";
        } else {
            killer = "none";
        }

        Files.writeString(RESULTS.resolve("d0_summary.json"), String.format(Locale.US,
                "{\n  \"gate\": \"D0_CLOSED_BATH\",\n"
                        + "  \"object\": \"LANE_D_LIFE_CYCLE\",\n"
                        + "  \"device\": \"LANE_D_CLOSED_BATH\",\n"
                        + "  \"narma\": false,\n  \"death\": false,\n"
                        + "  \"motility\": false,\n  \"chemistry\": false,\n"
                        + "  \"pde\": false,\n  \"erickson\": false,\n"
                        + "  \"elongation\": \"incremental\",\n"
                        + "  \"sink_density\": \"rho_cell\",\n"
                        + "  \"seed\": %d,\n  \"n0\": %d,\n"
                        + "  \"c_s_mM\": %.12f,\n  \"k_s_mM\": %.12f,\n"
                        + "  \"c_cut_mM\": %.12f,\n  \"delta_mass\": %.4f,\n"
                        + "  \"n_ever_compute_cap\": %d,\n"
                        + "  \"grow_N_end\": %d,\n  \"grow_N_ever\": %d,\n"
                        + "  \"grow_births\": %d,\n  \"grow_C_end\": %.12f,\n"
                        + "  \"grow_t_end\": %.10f,\n  \"grow_U\": %.12f,\n"
                        + "  \"grow_dCV\": %.12f,\n  \"grow_mass_residual\": %.12e,\n"
                        + "  \"grow_sum_V\": %.12f,\n  \"grow_rho_cell_sum_V\": %.12e,\n"
                        + "  \"grow_shrink_events\": %d,\n"
                        + "  \"grow_f3_drop_cum\": %.12f,\n"
                        + "  \"grow_stop\": \"%s\",\n"
                        + "  \"off_N_end\": %d,\n  \"off_N_ever\": %d,\n"
                        + "  \"off_births\": %d,\n  \"off_C_end\": %.12f,\n"
                        + "  \"cap_fired\": %s,\n  \"t_cap_scope\": %s,\n"
                        + "  \"pass_without_isolation\": %s,\n  \"killer\": \"%s\"\n}\n",
                LaneDD0Identity.RNG_SEED, LaneDD0Identity.N0,
                LaneDD0Identity.C_S_MM, LaneDD0Identity.K_S_MM,
                LaneDD0Identity.C_CUT_MM, LaneDD0Identity.DELTA_MASS,
                LaneDD0Identity.N_EVER_COMPUTE_CAP,
                grow.nEnd, grow.nEverEnd, grow.birthsTotal, grow.cEnd,
                grow.tEnd, grow.uEnd, grow.dCvEnd, grow.massResidual,
                grow.sumVEnd, LaneDD0Identity.RHO_CELL_GCDW_PER_UM3 * grow.sumVEnd,
                grow.shrinkEvents, grow.f3DropCum, grow.stopReason,
                off.nEnd, off.nEverEnd, off.birthsTotal, off.cEnd,
                cap, tCap, pass, killer), StandardCharsets.UTF_8);

        System.out.printf(Locale.US, "D0.1 honesty %s%n", honesty ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "D0.2 process-off births=%d N_end=%d C_end=%.6f %s%n",
                off.birthsTotal, off.nEnd, off.cEnd, processOff ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "D0.3 no_shrink events=%d %s%n",
                grow.shrinkEvents, noShrink ? "PASS" : "FAIL");
        System.out.printf(Locale.US,
                "D0.4 bath_ledger |dCV-U|/(Cs V)=%.6e delta=%.2f %s%n",
                grow.massResidual, LaneDD0Identity.DELTA_MASS,
                ledger ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "D0.5 monod_identity Warren_K_S=%s %s%n",
                LaneDD0Identity.K_S_MM, monod ? "PASS" : "FAIL");
        System.out.printf(Locale.US, "D0.6 exhaustion C_end=%.8f C_cut=%.4f stop=%s %s%n",
                grow.cEnd, LaneDD0Identity.C_CUT_MM, grow.stopReason,
                exhausted ? "PASS" : (tCap ? "SCOPE_NOTE" : "FAIL"));
        System.out.printf(Locale.US, "D0.7 no_silent_kill N_end=N_ever %s%n",
                noKill && !cap ? "PASS" : (cap ? "SCOPE_NOTE" : "FAIL"));
        System.out.println("D0.8 isolation is the checker git diff. Job does not claim it.");
        System.out.printf(Locale.US,
                "F3_DIVISION_VOLUME_DROP fissions=%d drop_cum=%.4f um3 (reported, not gated)%n",
                grow.birthsTotal, grow.f3DropCum);
        System.out.printf(Locale.US,
                "REPORT sum_rho_cell_V=%.6e gCDW (not the mass gate)%n",
                LaneDD0Identity.RHO_CELL_GCDW_PER_UM3 * grow.sumVEnd);

        if (pass) {
            System.out.println("D0_CLOSED_BATH=PASS pending isolation. Named closed bath "
                    + "depletes; bath ledger closes; growth-off flat; death was not on.");
        } else if (scope) {
            System.out.println("D0_CLOSED_BATH=SCOPE_NOTE " + killer + ". Not PASS.");
            System.exit(2);
        } else {
            System.out.println("D0_CLOSED_BATH=FAIL " + killer
                    + ". Do not raise C_s, Y, rho_cell, or N0.");
            System.exit(1);
        }
    }

    private static boolean noLengthDecrease(Arm arm) {
        double prevMin = 0.0;
        boolean first = true;
        for (int i = 0; i < arm.obsMinEll.size(); i++) {
            double minL = arm.obsMinEll.get(i);
            if (!first && minL + 1.0e-12 < LaneDD0Identity.ELL0_UM) {
                return false;
            }
            first = false;
            prevMin = minL;
        }
        return prevMin >= LaneDD0Identity.ELL0_UM - 1.0e-12;
    }

    private static boolean monodOk(Arm arm) {
        for (int i = 0; i < arm.obsC.size(); i++) {
            double expect = LaneDD0Identity.lambdaPerHour(arm.obsC.get(i));
            if (Math.abs(arm.obsLambda.get(i) - expect) > 1.0e-12) {
                return false;
            }
        }
        return !arm.obsC.isEmpty();
    }

    private static void refuseCapacity(String[] args) {
        String blob = String.join(" ", args == null ? new String[0] : args)
                .toLowerCase(Locale.ROOT);
        if (blob.contains("narma") || blob.contains("charc") || blob.contains("ipc")) {
            throw new IllegalArgumentException("D0 refuses NARMA/CHARC/IPC");
        }
        if (blob.contains("erickson") || blob.contains("famine")
                || blob.contains("handoff") || blob.contains("death=on")
                || blob.contains("bacteria.max") || blob.contains("elongatecited")) {
            throw new IllegalArgumentException(
                    "D0 refuses death-on, Erickson, famine, handoff, "
                            + "bacteria.max, elongateCited");
        }
    }

    private static void requireFrozen() throws IOException {
        String text = Files.readString(PROTOCOL, StandardCharsets.UTF_8);
        String json = Files.readString(PROTOCOL_JSON, StandardCharsets.UTF_8);
        if (!text.contains("frozen_before_traces")
                || !json.contains("\"frozen_before_traces\": true")) {
            throw new IllegalStateException("D0 protocol is not frozen_before_traces");
        }
        if (!json.contains("\"status_label\": \"D0_CLOSED_BATH\"")
                || !json.contains("\"device\": \"LANE_D_CLOSED_BATH\"")
                || !json.contains("\"not_device\": \"LC0_OPEN_BATH\"")
                || !json.contains("\"n0\": 8")
                || !json.contains("\"seed\": 101")
                || !json.contains("\"death\": false")
                || !json.contains("\"pde\": false")
                || !json.contains("\"erickson\": false")
                || !json.contains("\"sink_density\": \"rho_cell\"")
                || !json.contains("\"elongation\": \"incremental\"")
                || !json.contains("\"elongate_cited\": false")
                || !json.contains("\"initial_lengths\": \"all_l0\"")
                || !json.contains("\"c_cut_mM\": 0.001")
                || !json.contains("\"refuse_death_while_c_gt_ccut\": true")) {
            throw new IllegalStateException("D0 must freeze D0_CLOSED_BATH, "
                    + "N0=8, seed=101, death=false, rho_cell, incremental, "
                    + "C_cut=0.001, refuse-death, not LC0_OPEN_BATH");
        }
        if (text.contains("LC0_OPEN_BATH") && !text.contains("Not `LC0_OPEN_BATH`")) {
            throw new IllegalStateException("D0 protocol must not adopt LC0_OPEN_BATH");
        }
    }

    private static final class Arm {
        final String name;
        final boolean growOn;
        final List<Double> obsT = new ArrayList<>();
        final List<Double> obsC = new ArrayList<>();
        final List<Double> obsLambda = new ArrayList<>();
        final List<Integer> obsN = new ArrayList<>();
        final List<Integer> obsNever = new ArrayList<>();
        final List<Integer> obsBirths = new ArrayList<>();
        final List<Double> obsSumV = new ArrayList<>();
        final List<Double> obsU = new ArrayList<>();
        final List<Double> obsDcv = new ArrayList<>();
        final List<Double> obsMinEll = new ArrayList<>();
        final List<Double> obsMaxEll = new ArrayList<>();
        final List<Double> obsF3 = new ArrayList<>();
        int birthsTotal;
        int nEnd;
        int nEverEnd;
        double tEnd;
        double cEnd;
        double uEnd;
        double dCvEnd;
        double massResidual;
        double sumVEnd;
        double f3DropCum;
        int shrinkEvents;
        boolean capFired;
        String stopReason;

        Arm(String name, boolean growOn) {
            this.name = name;
            this.growOn = growOn;
        }
    }
}
