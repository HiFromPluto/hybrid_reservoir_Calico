package bsim.c0;

import bsim.BSimRandom;
import bsim.BSimStepScheduler;
import bsim.circuit.DaninoSIOccupiedDF;

/**
 * Gate C0: N explicit Object B cells and one well-mixed He compartment.
 * NOT_FIG4B. Not a painted-N copy of {@link DaninoSIOccupiedDF}.
 *
 * Per-cell membrane flux is {@code D (He - Hi_i)}. He is updated from the
 * sum of those fluxes times {@code v_cell / V_e}, plus {@code -mu He}.
 * When all cells are identical and {@code d = N v_cell / (N v_cell + V_e)
 * = 0.5}, the system reduces to the Object B 4-DDE, including
 * {@code d/(1-d)}.
 *
 * Injected {@link BSimRandom} is unused. Mechanics OFF. He is not a PDE.
 */
public final class DaninoSIC0WellMixed {

    public static final String OBJECT = DaninoSIOccupiedDF.OBJECT;
    public static final String OBJECT_STATUS = DaninoSIOccupiedDF.OBJECT_STATUS;
    public static final String NOT_FIG4B = "NOT_FIG4B";

    public static final int N_CELLS = 8;
    public static final double V_CELL = 1.5;
    public static final double V_E = 12.0;

    private final int n;
    private final double vCell;
    private final double vE;
    private final double d;
    private final BSimRandom unusedRng;
    private final double bProd;
    private final double gammaH;

    public DaninoSIC0WellMixed(BSimRandom rng) {
        this(N_CELLS, V_CELL, V_E, DaninoSIOccupiedDF.B, DaninoSIOccupiedDF.GAMMA_H, rng);
    }

    public DaninoSIC0WellMixed(
            int n,
            double vCell,
            double vE,
            double bProd,
            double gammaH,
            BSimRandom rng) {
        if (n < 1) {
            throw new IllegalArgumentException("n must be >= 1");
        }
        if (!(vCell > 0.0) || !(vE > 0.0)) {
            throw new IllegalArgumentException("volumes must be positive");
        }
        this.n = n;
        this.vCell = vCell;
        this.vE = vE;
        this.d = (n * vCell) / (n * vCell + vE);
        this.bProd = bProd;
        this.gammaH = gammaH;
        this.unusedRng = rng;
    }

    public static DaninoSIC0WellMixed closedMembrane(BSimRandom rng) {
        return new DaninoSIC0WellMixed(N_CELLS, V_CELL, V_E, 0.0, 0.0, rng);
    }

    public BSimRandom unusedRng() {
        return unusedRng;
    }

    public int nCells() {
        return n;
    }

    public double vCell() {
        return vCell;
    }

    public double vE() {
        return vE;
    }

    public double d() {
        return d;
    }

    public double dOverOneMinusD() {
        return (n * vCell) / vE;
    }

    public double intracellularMass(double[] hi) {
        double s = 0.0;
        for (double h : hi) {
            s += h * vCell;
        }
        return s;
    }

    public double extracellularMass(double he) {
        return he * vE;
    }

    public double totalAhlMass(double[] hi, double he) {
        return intracellularMass(hi) + extracellularMass(he);
    }

    /**
     * Algebra: identical-cell He RHS must match Object B, including
     * {@code d/(1-d)}. NOT_FIG4B.
     */
    public void assertReducesToObjectB(double[] y4, double hTau, double mu) {
        if (Math.abs(d - 0.5) > 1e-15) {
            throw new IllegalStateException("C0 d from volumes is " + d + " not 0.5 NOT_FIG4B");
        }
        if (Math.abs(dOverOneMinusD() - (d / (1.0 - d))) > 1e-15) {
            throw new IllegalStateException("N v_cell / V_e != d/(1-d) NOT_FIG4B");
        }
        double[] bulk = new double[4];
        DaninoSIOccupiedDF.rhs(y4, hTau, d, mu, bulk);

        double[] a = fill(y4[0]);
        double[] i = fill(y4[1]);
        double[] hi = fill(y4[2]);
        double[] hTauArr = fill(hTau);
        double[] dA = new double[n];
        double[] dI = new double[n];
        double[] dHi = new double[n];
        double[] dHe = new double[1];
        Rates rates = new Rates();
        rhs(a, i, hi, y4[3], hTauArr, mu, dA, dI, dHi, dHe, rates);

        if (Math.abs(dA[0] - bulk[0]) > 1e-14
                || Math.abs(dI[0] - bulk[1]) > 1e-14
                || Math.abs(dHi[0] - bulk[2]) > 1e-14
                || Math.abs(dHe[0] - bulk[3]) > 1e-14) {
            throw new IllegalStateException(
                    "C0 identical RHS does not reduce to Object B NOT_FIG4B "
                            + "dA=" + dA[0] + " vs " + bulk[0]
                            + " dI=" + dI[0] + " vs " + bulk[1]
                            + " dHi=" + dHi[0] + " vs " + bulk[2]
                            + " dHe=" + dHe[0] + " vs " + bulk[3]);
        }
        if (Math.abs(rates.memIntoCells + (-rates.memIntoCells)) > 0.0) {
            throw new IllegalStateException("membrane bookkeeping broken NOT_FIG4B");
        }
    }

    public Trajectory integrate(
            double mu,
            double[][] cellY0,
            double[] hiHistory,
            double he0,
            double tEnd,
            double sampleDt,
            double rkDt) {
        if (cellY0.length != n || hiHistory.length != n) {
            throw new IllegalArgumentException("IVP length must equal N");
        }
        int nMicro = (int) Math.round(tEnd / rkDt);
        if (Math.abs(nMicro * rkDt - tEnd) > 1e-12) {
            throw new IllegalArgumentException("rkDt must divide tEnd");
        }
        int nSample = (int) Math.round(tEnd / sampleDt);
        if (Math.abs(nSample * sampleDt - tEnd) > 1e-12) {
            throw new IllegalArgumentException("sampleDt must divide tEnd");
        }
        int stride = (int) Math.round(sampleDt / rkDt);
        if (Math.abs(stride * rkDt - sampleDt) > 1e-12) {
            throw new IllegalArgumentException("rkDt must divide sampleDt");
        }

        double[] a = new double[n];
        double[] i = new double[n];
        double[] hi = new double[n];
        for (int c = 0; c < n; c++) {
            a[c] = cellY0[c][0];
            i[c] = cellY0[c][1];
            hi[c] = cellY0[c][2];
        }
        double[] he = {he0};
        double[][] hiTape = new double[n][nMicro + 1];
        for (int c = 0; c < n; c++) {
            hiTape[c][0] = hi[c];
        }

        DaninoSIC0AhlLedger ledger = new DaninoSIC0AhlLedger();
        ledger.reset(totalAhlMass(hi, he[0]));

        double[] tOut = new double[nSample + 1];
        double[][] meanOut = new double[nSample + 1][4];
        double[][] hiOut = new double[nSample + 1][n];
        double[][] iOut = new double[nSample + 1][n];
        recordSample(0, 0.0, a, i, hi, he[0], tOut, meanOut, hiOut, iOut);

        double[] dA1 = new double[n];
        double[] dI1 = new double[n];
        double[] dHi1 = new double[n];
        double[] dHe1 = new double[1];
        double[] dA2 = new double[n];
        double[] dI2 = new double[n];
        double[] dHi2 = new double[n];
        double[] dHe2 = new double[1];
        double[] dA3 = new double[n];
        double[] dI3 = new double[n];
        double[] dHi3 = new double[n];
        double[] dHe3 = new double[1];
        double[] dA4 = new double[n];
        double[] dI4 = new double[n];
        double[] dHi4 = new double[n];
        double[] dHe4 = new double[1];
        double[] aTmp = new double[n];
        double[] iTmp = new double[n];
        double[] hiTmp = new double[n];
        double[] hTau = new double[n];
        Rates r1 = new Rates();
        Rates r2 = new Rates();
        Rates r3 = new Rates();
        Rates r4 = new Rates();

        int[] micro = {0};
        double[] maxHiSpread = {0.0};
        BSimStepScheduler scheduler = new BSimStepScheduler(rkDt, tEnd);
        scheduler.run(new BSimStepScheduler.Adapter() {
            @Override
            public void integrateModels(BSimStepScheduler.Context context) {
                int step = micro[0];
                double t = step * rkDt;
                fillTau(hiTape, step, rkDt, t - DaninoSIOccupiedDF.TAU, hiHistory, cellY0, hTau);
                rhs(a, i, hi, he[0], hTau, mu, dA1, dI1, dHi1, dHe1, r1);

                axpy(a, i, hi, he[0], dA1, dI1, dHi1, dHe1[0], 0.5 * rkDt, aTmp, iTmp, hiTmp);
                double he2 = he[0] + 0.5 * rkDt * dHe1[0];
                fillTau(hiTape, step, rkDt, t + 0.5 * rkDt - DaninoSIOccupiedDF.TAU, hiHistory, cellY0, hTau);
                rhs(aTmp, iTmp, hiTmp, he2, hTau, mu, dA2, dI2, dHi2, dHe2, r2);

                axpy(a, i, hi, he[0], dA2, dI2, dHi2, dHe2[0], 0.5 * rkDt, aTmp, iTmp, hiTmp);
                double he3 = he[0] + 0.5 * rkDt * dHe2[0];
                rhs(aTmp, iTmp, hiTmp, he3, hTau, mu, dA3, dI3, dHi3, dHe3, r3);

                axpy(a, i, hi, he[0], dA3, dI3, dHi3, dHe3[0], rkDt, aTmp, iTmp, hiTmp);
                double he4 = he[0] + rkDt * dHe3[0];
                fillTau(hiTape, step, rkDt, t + rkDt - DaninoSIOccupiedDF.TAU, hiHistory, cellY0, hTau);
                rhs(aTmp, iTmp, hiTmp, he4, hTau, mu, dA4, dI4, dHi4, dHe4, r4);

                for (int c = 0; c < n; c++) {
                    a[c] += (rkDt / 6.0) * (dA1[c] + 2.0 * dA2[c] + 2.0 * dA3[c] + dA4[c]);
                    i[c] += (rkDt / 6.0) * (dI1[c] + 2.0 * dI2[c] + 2.0 * dI3[c] + dI4[c]);
                    hi[c] += (rkDt / 6.0) * (dHi1[c] + 2.0 * dHi2[c] + 2.0 * dHi3[c] + dHi4[c]);
                }
                he[0] += (rkDt / 6.0) * (dHe1[0] + 2.0 * dHe2[0] + 2.0 * dHe3[0] + dHe4[0]);
                ledger.addWeighted(
                        (rkDt / 6.0) * (r1.synth + 2.0 * r2.synth + 2.0 * r3.synth + r4.synth),
                        (rkDt / 6.0) * (r1.gammaH + 2.0 * r2.gammaH + 2.0 * r3.gammaH + r4.gammaH),
                        (rkDt / 6.0) * (r1.muLoss + 2.0 * r2.muLoss + 2.0 * r3.muLoss + r4.muLoss),
                        (rkDt / 6.0) * (r1.memIntoCells + 2.0 * r2.memIntoCells
                                + 2.0 * r3.memIntoCells + r4.memIntoCells));
                for (int c = 0; c < n; c++) {
                    hiTape[c][step + 1] = hi[c];
                }
                micro[0] = step + 1;
            }

            @Override
            public void observe(BSimStepScheduler.Context context) {
                if (context.isInitialObservation()) {
                    return;
                }
                int completed = context.getCompletedUpdates();
                if (completed % stride == 0) {
                    int s = completed / stride;
                    recordSample(s, completed * rkDt, a, i, hi, he[0], tOut, meanOut, hiOut, iOut);
                    maxHiSpread[0] = Math.max(maxHiSpread[0], hiSpread(hi));
                }
            }
        });

        if (micro[0] != nMicro) {
            throw new IllegalStateException("scheduler steps " + micro[0] + " != " + nMicro);
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
        double remaining = totalAhlMass(hi, he[0]);
        return new Trajectory(
                tOut,
                meanOut,
                hiOut,
                iOut,
                d,
                mu,
                rkDt,
                negative,
                ledger,
                remaining,
                maxHiSpread[0],
                hiSpread(hi));
    }

    private void rhs(
            double[] a,
            double[] i,
            double[] hi,
            double he,
            double[] hTau,
            double mu,
            double[] dA,
            double[] dI,
            double[] dHi,
            double[] dHe,
            Rates rates) {
        double dens = 1.0 - Math.pow(d / DaninoSIOccupiedDF.D0, 4.0);
        double memSum = 0.0;
        double synth = 0.0;
        double gLoss = 0.0;
        for (int c = 0; c < n; c++) {
            double p = DaninoSIOccupiedDF.production(hTau[c]);
            double denomF = 1.0 + DaninoSIOccupiedDF.F * (a[c] + i[c]);
            dA[c] = DaninoSIOccupiedDF.C_A * dens * p - DaninoSIOccupiedDF.GAMMA_A * a[c] / denomF;
            dI[c] = DaninoSIOccupiedDF.C_I * dens * p - DaninoSIOccupiedDF.GAMMA_I * i[c] / denomF;
            double synthConc = bProd * i[c] / (1.0 + DaninoSIOccupiedDF.K * i[c]);
            double enzConc = (gammaH * a[c] * hi[c]) / (1.0 + DaninoSIOccupiedDF.G * a[c]);
            double memConc = DaninoSIOccupiedDF.D_MEM * (he - hi[c]);
            dHi[c] = synthConc - enzConc + memConc;
            memSum += memConc;
            synth += vCell * synthConc;
            gLoss += vCell * enzConc;
        }
        dHe[0] = -(vCell / vE) * memSum - mu * he;
        rates.synth = synth;
        rates.gammaH = gLoss;
        rates.muLoss = mu * he * vE;
        rates.memIntoCells = vCell * memSum;
    }

    private double[] fill(double v) {
        double[] x = new double[n];
        java.util.Arrays.fill(x, v);
        return x;
    }

    private static void axpy(
            double[] a,
            double[] i,
            double[] hi,
            double he,
            double[] dA,
            double[] dI,
            double[] dHi,
            double dHe,
            double h,
            double[] aOut,
            double[] iOut,
            double[] hiOut) {
        for (int c = 0; c < a.length; c++) {
            aOut[c] = a[c] + h * dA[c];
            iOut[c] = i[c] + h * dI[c];
            hiOut[c] = hi[c] + h * dHi[c];
        }
    }

    private static void fillTau(
            double[][] tape,
            int lastFilled,
            double rkDt,
            double tDelay,
            double[] hiHistory,
            double[][] cellY0,
            double[] hTau) {
        for (int c = 0; c < tape.length; c++) {
            hTau[c] = hiAt(tape[c], lastFilled, rkDt, tDelay, hiHistory[c], cellY0[c][2]);
        }
    }

    private static double hiAt(
            double[] tape,
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
            return tape[lastFilled];
        }
        double frac = idx - i;
        int i1 = i + 1;
        if (i1 > lastFilled) {
            return tape[i];
        }
        return tape[i] * (1.0 - frac) + tape[i1] * frac;
    }

    private static void recordSample(
            int s,
            double t,
            double[] a,
            double[] i,
            double[] hi,
            double he,
            double[] tOut,
            double[][] meanOut,
            double[][] hiOut,
            double[][] iOut) {
        tOut[s] = t;
        meanOut[s][0] = mean(a);
        meanOut[s][1] = mean(i);
        meanOut[s][2] = mean(hi);
        meanOut[s][3] = he;
        System.arraycopy(hi, 0, hiOut[s], 0, hi.length);
        System.arraycopy(i, 0, iOut[s], 0, i.length);
    }

    private static double mean(double[] x) {
        double s = 0.0;
        for (double v : x) {
            s += v;
        }
        return s / x.length;
    }

    private static double hiSpread(double[] hi) {
        double lo = hi[0];
        double hiV = hi[0];
        for (double v : hi) {
            if (v < lo) {
                lo = v;
            }
            if (v > hiV) {
                hiV = v;
            }
        }
        return hiV - lo;
    }

    private static final class Rates {
        double synth;
        double gammaH;
        double muLoss;
        double memIntoCells;
    }

    public static final class Trajectory {
        public final double[] t;
        /** Mean (A, I, Hi) and shared He. */
        public final double[][] Y;
        public final double[][] Hi;
        public final double[][] I;
        public final double d;
        public final double mu;
        public final double rkDt;
        public final boolean negativeState;
        public final DaninoSIC0AhlLedger ledger;
        public final double remainingMass;
        public final double maxHiSpread;
        public final double finalHiSpread;

        Trajectory(
                double[] t,
                double[][] Y,
                double[][] Hi,
                double[][] I,
                double d,
                double mu,
                double rkDt,
                boolean negativeState,
                DaninoSIC0AhlLedger ledger,
                double remainingMass,
                double maxHiSpread,
                double finalHiSpread) {
            this.t = t;
            this.Y = Y;
            this.Hi = Hi;
            this.I = I;
            this.d = d;
            this.mu = mu;
            this.rkDt = rkDt;
            this.negativeState = negativeState;
            this.ledger = ledger;
            this.remainingMass = remainingMass;
            this.maxHiSpread = maxHiSpread;
            this.finalHiSpread = finalHiSpread;
        }

        public double[] luxIMean() {
            double[] out = new double[Y.length];
            for (int n = 0; n < Y.length; n++) {
                out[n] = Y[n][1];
            }
            return out;
        }

        public double[] he() {
            double[] out = new double[Y.length];
            for (int n = 0; n < Y.length; n++) {
                out[n] = Y[n][3];
            }
            return out;
        }

        public double residual() {
            return ledger.residual(remainingMass);
        }

        public double relativeResidual() {
            double mStar = DaninoSIC0AhlLedger.characteristicMass(
                    ledger.initialMass(), remainingMass, ledger.synthSource());
            return Math.abs(residual()) / mStar;
        }
    }
}
