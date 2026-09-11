package bsim.c1b;

import bsim.BSim;
import bsim.BSimRandom;
import bsim.BSimStepScheduler;
import bsim.c0.DaninoSIC0WellMixed;
import bsim.c1.DaninoSIC1AhlLedger;
import bsim.circuit.DaninoSIOccupiedDF;
import bsim.transport.BSimTransportField;
import bsim.transport.BSimTransportLedger;

/**
 * C1b island identity: C0-equivalent extracellular operator. NOT_FIG4B.
 *
 * He, mu, and Ve live on one fluid voxel of volume 12 µm³. Other voxels
 * are not constructed. D1_spatial=0. Not the open P0 chip.
 *
 * Traces use the C1 membrane / field path on that voxel (frozen-He RK4,
 * signed transferQuantity, field.decay). Algebra uses C0 vs Object B.
 */
public final class DaninoSIC1bIsland {

    public static final String OBJECT = DaninoSIOccupiedDF.OBJECT;
    public static final String NOT_FIG4B = "NOT_FIG4B";
    public static final double D1_SPATIAL = 0.0;
    public static final double V_BOUND_X = 2.0;
    public static final double V_BOUND_Y = 2.0;
    public static final double V_BOUND_Z = 3.0;

    private final BSimRandom unusedRng;
    private final DaninoSIC0WellMixed wellMixed;
    private final double boxVolume;

    public DaninoSIC1bIsland(BSimRandom rng) {
        this.unusedRng = rng;
        this.wellMixed = new DaninoSIC0WellMixed(rng);
        BSim sim = new BSim();
        sim.setBound(V_BOUND_X, V_BOUND_Y, V_BOUND_Z);
        sim.setSolid(true, true, true);
        boolean[][][] fluid = new boolean[1][1][1];
        fluid[0][0][0] = true;
        BSimTransportField field = new BSimTransportField(
                sim, new int[] {1, 1, 1}, D1_SPATIAL, 0.0, fluid);
        this.boxVolume = field.getBox()[0] * field.getBox()[1] * field.getBox()[2];
        if (Math.abs(boxVolume - wellMixed.vE()) > 1e-12) {
            throw new IllegalStateException(
                    "island voxel volume " + boxVolume + " != Ve " + wellMixed.vE() + " NOT_FIG4B");
        }
        if (Math.abs(wellMixed.d() - 0.5) > 1e-15) {
            throw new IllegalStateException("island d is " + wellMixed.d() + " not 0.5 NOT_FIG4B");
        }
        if (wellMixed.nCells() != 8) {
            throw new IllegalStateException("island N must be C0's 8 NOT_FIG4B");
        }
    }

    public BSimRandom unusedRng() {
        return unusedRng;
    }

    public DaninoSIC0WellMixed wellMixed() {
        return wellMixed;
    }

    public int nCells() {
        return wellMixed.nCells();
    }

    public double vCell() {
        return wellMixed.vCell();
    }

    public double vE() {
        return wellMixed.vE();
    }

    public double d() {
        return wellMixed.d();
    }

    public double d1Spatial() {
        return D1_SPATIAL;
    }

    public double boxVolume() {
        return boxVolume;
    }

    public int nFluidVoxels() {
        return 1;
    }

    /**
     * C1 field/membrane operator on C0 volumes. He(0) and mu hit the single
     * fluid voxel only. D1=0. NOT_FIG4B.
     */
    public Trajectory integrate(
            double mu,
            String ensemble,
            double tEnd,
            double sampleDt,
            double rkDt) {
        int n = nCells();
        double vCellVol = vCell();
        double dFrac = d();
        if (!(rkDt > 0.0) || Math.abs(tEnd / rkDt - Math.round(tEnd / rkDt)) > 1e-9) {
            throw new IllegalArgumentException("rkDt must divide tEnd NOT_FIG4B");
        }
        int nMicro = (int) Math.round(tEnd / rkDt);
        int nSample = (int) Math.round(tEnd / sampleDt);
        if (Math.abs(nSample * sampleDt - tEnd) > 1e-12) {
            throw new IllegalArgumentException("sampleDt must divide tEnd NOT_FIG4B");
        }
        int stride = (int) Math.round(sampleDt / rkDt);
        if (Math.abs(stride * rkDt - sampleDt) > 1e-12) {
            throw new IllegalArgumentException("rkDt must divide sampleDt NOT_FIG4B");
        }

        boolean kick = !"SI_BASAL_PERTURB".equals(ensemble);
        double a0 = 0.0;
        double i0 = kick ? 0.0 : 1.0;
        double hi0 = kick ? 0.05 : 0.0;
        double he0 = kick ? 0.05 : 0.0;
        double hist = hi0;

        double[] a = new double[n];
        double[] i = new double[n];
        double[] hi = new double[n];
        java.util.Arrays.fill(a, a0);
        java.util.Arrays.fill(i, i0);
        java.util.Arrays.fill(hi, hi0);

        BSim sim = new BSim();
        sim.setDt(rkDt);
        sim.setSimulationTime(tEnd);
        sim.setBound(V_BOUND_X, V_BOUND_Y, V_BOUND_Z);
        sim.setSolid(true, true, true);
        sim.setLeaky(false, false, false, false, false, false);
        if (unusedRng != null) {
            sim.setRandom(unusedRng);
        }

        boolean[][][] fluid = new boolean[1][1][1];
        fluid[0][0][0] = true;
        BSimTransportField field = new BSimTransportField(
                sim, new int[] {1, 1, 1}, D1_SPATIAL, mu, fluid);
        if (Math.abs(field.getBox()[0] * field.getBox()[1] * field.getBox()[2] - vE()) > 1e-12) {
            throw new IllegalStateException("integrate voxel volume != Ve NOT_FIG4B");
        }
        if (he0 > 0.0) {
            field.setConc(he0);
        }
        field.getTransportLedger().reset();

        DaninoSIC1AhlLedger ahl = new DaninoSIC1AhlLedger();
        ahl.reset(totalAhlMass(hi, field, vCellVol));

        int cap = (int) Math.round(DaninoSIOccupiedDF.TAU / rkDt) + 2;
        double[][] hiTape = new double[n][cap];
        for (int c = 0; c < n; c++) {
            hiTape[c][0] = hi[c];
        }

        double[] tOut = new double[nSample + 1];
        double[][] meanOut = new double[nSample + 1][4];
        recordSample(0, 0.0, a, i, hi, field.getConc(0, 0, 0), tOut, meanOut);

        double[] heSampled = new double[n];
        double[] dA1 = new double[n];
        double[] dI1 = new double[n];
        double[] dHi1 = new double[n];
        double[] dA2 = new double[n];
        double[] dI2 = new double[n];
        double[] dHi2 = new double[n];
        double[] dA3 = new double[n];
        double[] dI3 = new double[n];
        double[] dHi3 = new double[n];
        double[] dA4 = new double[n];
        double[] dI4 = new double[n];
        double[] dHi4 = new double[n];
        double[] aTmp = new double[n];
        double[] iTmp = new double[n];
        double[] hiTmp = new double[n];
        double[] hTau = new double[n];
        double[] memDeposit = new double[n];
        Rates r1 = new Rates();
        Rates r2 = new Rates();
        Rates r3 = new Rates();
        Rates r4 = new Rates();
        int[] micro = {0};
        double[] maxHiSpread = {0.0};

        BSimStepScheduler scheduler = new BSimStepScheduler(sim);
        scheduler.run(new BSimStepScheduler.Adapter() {
            @Override
            public void transport(BSimStepScheduler.Context context, double duration) {
                field.diffuse(duration);
                field.decay(duration);
            }

            @Override
            public void sampleState(BSimStepScheduler.Context context) {
                double he = field.getConc(0, 0, 0);
                java.util.Arrays.fill(heSampled, he);
            }

            @Override
            public void integrateModels(BSimStepScheduler.Context context) {
                int step = micro[0];
                double t = step * rkDt;
                fillTau(hiTape, cap, step, rkDt, t - DaninoSIOccupiedDF.TAU, hist, hi0, hTau);
                rhsFrozenHe(a, i, hi, heSampled, hTau, dFrac, vCellVol, dA1, dI1, dHi1, r1);
                axpy(a, i, hi, dA1, dI1, dHi1, 0.5 * rkDt, aTmp, iTmp, hiTmp);
                fillTau(hiTape, cap, step, rkDt, t + 0.5 * rkDt - DaninoSIOccupiedDF.TAU, hist, hi0, hTau);
                rhsFrozenHe(aTmp, iTmp, hiTmp, heSampled, hTau, dFrac, vCellVol, dA2, dI2, dHi2, r2);
                axpy(a, i, hi, dA2, dI2, dHi2, 0.5 * rkDt, aTmp, iTmp, hiTmp);
                rhsFrozenHe(aTmp, iTmp, hiTmp, heSampled, hTau, dFrac, vCellVol, dA3, dI3, dHi3, r3);
                axpy(a, i, hi, dA3, dI3, dHi3, rkDt, aTmp, iTmp, hiTmp);
                fillTau(hiTape, cap, step, rkDt, t + rkDt - DaninoSIOccupiedDF.TAU, hist, hi0, hTau);
                rhsFrozenHe(aTmp, iTmp, hiTmp, heSampled, hTau, dFrac, vCellVol, dA4, dI4, dHi4, r4);

                for (int c = 0; c < n; c++) {
                    double hiS1 = hi[c];
                    double hiS2 = hi[c] + 0.5 * rkDt * dHi1[c];
                    double hiS3 = hi[c] + 0.5 * rkDt * dHi2[c];
                    double hiS4 = hi[c] + rkDt * dHi3[c];
                    double memConc = (DaninoSIOccupiedDF.D_MEM / 6.0) * (
                            (heSampled[c] - hiS1)
                                    + 2.0 * (heSampled[c] - hiS2)
                                    + 2.0 * (heSampled[c] - hiS3)
                                    + (heSampled[c] - hiS4));
                    memDeposit[c] = rkDt * vCellVol * memConc;
                    a[c] += (rkDt / 6.0) * (dA1[c] + 2.0 * dA2[c] + 2.0 * dA3[c] + dA4[c]);
                    i[c] += (rkDt / 6.0) * (dI1[c] + 2.0 * dI2[c] + 2.0 * dI3[c] + dI4[c]);
                    hi[c] += (rkDt / 6.0) * (dHi1[c] + 2.0 * dHi2[c] + 2.0 * dHi3[c] + dHi4[c]);
                }
                ahl.addWeighted(
                        (rkDt / 6.0) * (r1.synth + 2.0 * r2.synth + 2.0 * r3.synth + r4.synth),
                        (rkDt / 6.0) * (r1.gammaH + 2.0 * r2.gammaH + 2.0 * r3.gammaH + r4.gammaH),
                        sum(memDeposit));
                for (int c = 0; c < n; c++) {
                    hiTape[c][(step + 1) % cap] = hi[c];
                }
                micro[0] = step + 1;
            }

            @Override
            public void depositFluxes(BSimStepScheduler.Context context) {
                for (int c = 0; c < n; c++) {
                    field.transferQuantity(0, 0, 0, -memDeposit[c]);
                }
            }

            @Override
            public void observe(BSimStepScheduler.Context context) {
                if (context.isInitialObservation()) {
                    return;
                }
                int completed = context.getCompletedUpdates();
                if (completed % stride == 0) {
                    int s = completed / stride;
                    recordSample(s, completed * rkDt, a, i, hi, field.getConc(0, 0, 0), tOut, meanOut);
                    maxHiSpread[0] = Math.max(maxHiSpread[0], hiSpread(hi));
                }
            }
        });

        if (micro[0] != nMicro) {
            throw new IllegalStateException("scheduler steps " + micro[0] + " != " + nMicro + " NOT_FIG4B");
        }
        boolean negative = false;
        for (double[] row : meanOut) {
            for (double v : row) {
                if (v < -1e-12) {
                    negative = true;
                    break;
                }
            }
        }
        double remaining = totalAhlMass(hi, field, vCellVol);
        return new Trajectory(
                tOut, meanOut, dFrac, mu, rkDt, ensemble, negative, ahl,
                remaining, field.getTransportLedger(), maxHiSpread[0], hiSpread(hi));
    }

    private static void rhsFrozenHe(
            double[] a,
            double[] i,
            double[] hi,
            double[] he,
            double[] hTau,
            double d,
            double vCell,
            double[] dA,
            double[] dI,
            double[] dHi,
            Rates rates) {
        double dens = 1.0 - Math.pow(d / DaninoSIOccupiedDF.D0, 4.0);
        double synth = 0.0;
        double gLoss = 0.0;
        for (int c = 0; c < a.length; c++) {
            double p = DaninoSIOccupiedDF.production(hTau[c]);
            double denomF = 1.0 + DaninoSIOccupiedDF.F * (a[c] + i[c]);
            dA[c] = DaninoSIOccupiedDF.C_A * dens * p - DaninoSIOccupiedDF.GAMMA_A * a[c] / denomF;
            dI[c] = DaninoSIOccupiedDF.C_I * dens * p - DaninoSIOccupiedDF.GAMMA_I * i[c] / denomF;
            double synthConc = DaninoSIOccupiedDF.B * i[c] / (1.0 + DaninoSIOccupiedDF.K * i[c]);
            double enzConc = (DaninoSIOccupiedDF.GAMMA_H * a[c] * hi[c])
                    / (1.0 + DaninoSIOccupiedDF.G * a[c]);
            double memConc = DaninoSIOccupiedDF.D_MEM * (he[c] - hi[c]);
            dHi[c] = synthConc - enzConc + memConc;
            synth += vCell * synthConc;
            gLoss += vCell * enzConc;
        }
        rates.synth = synth;
        rates.gammaH = gLoss;
    }

    private static double totalAhlMass(double[] hi, BSimTransportField field, double vCell) {
        double s = 0.0;
        for (double h : hi) {
            s += h * vCell;
        }
        return s + field.totalFluidQuantity();
    }

    private static void recordSample(
            int s,
            double t,
            double[] a,
            double[] i,
            double[] hi,
            double he,
            double[] tOut,
            double[][] meanOut) {
        tOut[s] = t;
        meanOut[s][0] = mean(a);
        meanOut[s][1] = mean(i);
        meanOut[s][2] = mean(hi);
        meanOut[s][3] = he;
    }

    private static void axpy(
            double[] a, double[] i, double[] hi,
            double[] dA, double[] dI, double[] dHi,
            double h,
            double[] aOut, double[] iOut, double[] hiOut) {
        for (int c = 0; c < a.length; c++) {
            aOut[c] = a[c] + h * dA[c];
            iOut[c] = i[c] + h * dI[c];
            hiOut[c] = hi[c] + h * dHi[c];
        }
    }

    private static void fillTau(
            double[][] tape,
            int cap,
            int lastFilled,
            double rkDt,
            double tDelay,
            double hiHistory,
            double hiAtZero,
            double[] hTau) {
        for (int c = 0; c < tape.length; c++) {
            hTau[c] = hiAt(tape[c], cap, lastFilled, rkDt, tDelay, hiHistory, hiAtZero);
        }
    }

    private static double hiAt(
            double[] tape,
            int cap,
            int lastFilled,
            double rkDt,
            double tDelay,
            double hiHistory,
            double hiAtZero) {
        if (tDelay < 0.0) {
            return hiHistory;
        }
        if (tDelay <= 1e-15) {
            return hiAtZero;
        }
        double idx = tDelay / rkDt;
        int i = (int) Math.floor(idx);
        if (i < 0) {
            return hiHistory;
        }
        if (i >= lastFilled) {
            return tape[lastFilled % cap];
        }
        double frac = idx - i;
        int i1 = i + 1;
        if (i1 > lastFilled) {
            return tape[i % cap];
        }
        return tape[i % cap] * (1.0 - frac) + tape[i1 % cap] * frac;
    }

    private static double mean(double[] x) {
        double s = 0.0;
        for (double v : x) {
            s += v;
        }
        return s / x.length;
    }

    private static double sum(double[] x) {
        double s = 0.0;
        for (double v : x) {
            s += v;
        }
        return s;
    }

    private static double hiSpread(double[] hi) {
        double lo = hi[0];
        double hiVal = hi[0];
        for (double v : hi) {
            lo = Math.min(lo, v);
            hiVal = Math.max(hiVal, v);
        }
        return hiVal - lo;
    }

    private static final class Rates {
        double synth;
        double gammaH;
    }

    public static final class Trajectory {
        public final double[] t;
        /** Mean (A, I, Hi) and island He. */
        public final double[][] Y;
        public final double d;
        public final double mu;
        public final double rkDt;
        public final String ensemble;
        public final boolean negativeState;
        public final LedgerView ledger;
        public final double remainingMass;
        public final double maxHiSpread;
        public final double finalHiSpread;
        private final DaninoSIC1AhlLedger ahl;
        private final BSimTransportLedger extra;

        Trajectory(
                double[] t,
                double[][] Y,
                double d,
                double mu,
                double rkDt,
                String ensemble,
                boolean negativeState,
                DaninoSIC1AhlLedger ahl,
                double remainingMass,
                BSimTransportLedger extra,
                double maxHiSpread,
                double finalHiSpread) {
            this.t = t;
            this.Y = Y;
            this.d = d;
            this.mu = mu;
            this.rkDt = rkDt;
            this.ensemble = ensemble;
            this.negativeState = negativeState;
            this.ahl = ahl;
            this.remainingMass = remainingMass;
            this.extra = extra;
            this.maxHiSpread = maxHiSpread;
            this.finalHiSpread = finalHiSpread;
            this.ledger = new LedgerView(ahl, extra);
        }

        public double[] luxIMean() {
            double[] out = new double[Y.length];
            for (int n = 0; n < Y.length; n++) {
                out[n] = Y[n][1];
            }
            return out;
        }

        public double residual() {
            return ahl.residual(remainingMass, extra);
        }

        public double relativeResidual() {
            return ahl.relativeResidual(remainingMass, extra);
        }
    }

    public static final class LedgerView {
        private final DaninoSIC1AhlLedger ahl;
        private final BSimTransportLedger extra;

        LedgerView(DaninoSIC1AhlLedger ahl, BSimTransportLedger extra) {
            this.ahl = ahl;
            this.extra = extra;
        }

        public double muLoss() {
            return extra.getDecayLoss();
        }

        public double membraneCancel() {
            return ahl.membraneCancel();
        }

        public double gammaHLoss() {
            return ahl.gammaHLoss();
        }

        public double synthSource() {
            return ahl.synthSource();
        }
    }
}
