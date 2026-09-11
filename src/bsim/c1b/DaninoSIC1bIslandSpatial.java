package bsim.c1b;

import bsim.BSim;
import bsim.BSimRandom;
import bsim.BSimStepScheduler;
import bsim.c1.DaninoSIC1AhlLedger;
import bsim.circuit.DaninoSIOccupiedDF;
import bsim.p0.DaninoPocketGeometry;
import bsim.transport.BSimTransportField;
import bsim.transport.BSimTransportLedger;

/**
 * C1b island extra: C1 packed patch voxels are the only fluid. NOT_FIG4B.
 * Not occupancy identity. D1_spatial=800 inside the island; solids outside.
 */
public final class DaninoSIC1bIslandSpatial {

    public static final String OBJECT = DaninoSIOccupiedDF.OBJECT;
    public static final String OBJECT_STATUS = DaninoSIOccupiedDF.OBJECT_STATUS;
    public static final String NOT_FIG4B = "NOT_FIG4B";
    public static final String MASK_SHA256 =
            "2b53bb0b872477ab28afab8b20da48e1c63daa834614caf55aab63043081f99b";

    public static final double V_CELL = 1.5;
    public static final int N_CELLS = 110;
    public static final int NX_ISLAND = 5;
    public static final int NY_ISLAND = 5;
    public static final int NZ_ISLAND = 1;
    public static final int PATCH_I0 = 0;
    public static final int PATCH_I1 = 4;
    public static final int PATCH_J0 = 0;
    public static final int PATCH_J1 = 4;
    public static final int PATCH_K = 0;
    public static final int N_PATCH_VOXELS = 25;
    public static final double V_VOXEL = DaninoPocketGeometry.DX_UM
            * DaninoPocketGeometry.DY_UM
            * DaninoPocketGeometry.DZ_UM;
    public static final double V_E_LOCAL = N_PATCH_VOXELS * V_VOXEL;
    public static final double D_LOCAL = (N_CELLS * V_CELL) / (N_CELLS * V_CELL + V_E_LOCAL);
    public static final double D1_SPATIAL = 800.0;
    public static final String D1_LABEL = "SI_TAKEN_SPATIAL";
    public static final double BOUND_X = NX_ISLAND * DaninoPocketGeometry.DX_UM;
    public static final double BOUND_Y = NY_ISLAND * DaninoPocketGeometry.DY_UM;
    public static final double BOUND_Z = NZ_ISLAND * DaninoPocketGeometry.DZ_UM;

    private final int n;
    private final double vCell;
    private final double d;
    private final double d1Spatial;
    private final BSimRandom unusedRng;
    private final boolean[][][] fluid;
    private final String maskSha;
    private final int[] vi;
    private final int[] vj;
    private final int[] vk;
    private final double[] xBsim;
    private final double[] yBsim;
    private final double[] zBsim;
    private final int[] patchI;
    private final int[] patchJ;

    public DaninoSIC1bIslandSpatial(BSimRandom rng) {
        this(N_CELLS, V_CELL, D1_SPATIAL, rng);
    }

    public DaninoSIC1bIslandSpatial(int n, double vCell, double d1Spatial, BSimRandom rng) {
        if (n != N_CELLS) {
            throw new IllegalArgumentException("C1b island extra N is frozen at " + N_CELLS + " NOT_FIG4B");
        }
        if (Math.abs(vCell - V_CELL) > 0.0) {
            throw new IllegalArgumentException("C1b v_cell is frozen at " + V_CELL + " NOT_FIG4B");
        }
        if (Math.abs(d1Spatial - D1_SPATIAL) > 0.0) {
            throw new IllegalArgumentException("C1b extra D1_spatial is frozen at " + D1_SPATIAL + " NOT_FIG4B");
        }
        this.n = n;
        this.vCell = vCell;
        this.d = D_LOCAL;
        this.d1Spatial = d1Spatial;
        this.unusedRng = rng;
        this.fluid = buildIslandMask();
        this.maskSha = "island_patch_25_voxels_NOT_FIG4B";
        int nFluid = countIslandFluid(fluid);
        if (nFluid != N_PATCH_VOXELS) {
            throw new IllegalStateException(
                    "island fluid voxels " + nFluid + " != " + N_PATCH_VOXELS + " NOT_FIG4B");
        }
        if (Math.abs(d - 0.5) > 1e-15) {
            throw new IllegalStateException("C1b island d is " + d + " not 0.5 NOT_FIG4B");
        }
        this.vi = new int[n];
        this.vj = new int[n];
        this.vk = new int[n];
        this.xBsim = new double[n];
        this.yBsim = new double[n];
        this.zBsim = new double[n];
        this.patchI = new int[N_PATCH_VOXELS];
        this.patchJ = new int[N_PATCH_VOXELS];
        placeHeldPack();
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

    public double vELocal() {
        return V_E_LOCAL;
    }

    public double d() {
        return d;
    }

    public double d1Spatial() {
        return d1Spatial;
    }

    public String maskSha256() {
        return maskSha;
    }

    public int nPatchVoxels() {
        return N_PATCH_VOXELS;
    }

    public Trajectory integrate(
            double mu,
            String ensemble,
            double tEnd,
            double sampleDt,
            double rkDt) {
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
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);
        sim.setLeaky(false, false, false, false, false, false);
        if (unusedRng != null) {
            sim.setRandom(unusedRng);
        }

        int[] boxes = {NX_ISLAND, NY_ISLAND, NZ_ISLAND};
        BSimTransportField field = new BSimTransportField(sim, boxes, d1Spatial, mu, fluid);
        if (he0 > 0.0) {
            field.setConc(he0);
        }
        field.getTransportLedger().reset();

        DaninoSIC1AhlLedger ledger = new DaninoSIC1AhlLedger();
        ledger.reset(totalAhlMass(hi, field));

        int cap = (int) Math.round(DaninoSIOccupiedDF.TAU / rkDt) + 2;
        double[][] hiTape = new double[n][cap];
        for (int c = 0; c < n; c++) {
            hiTape[c][0] = hi[c];
        }

        double[] tOut = new double[nSample + 1];
        double[][] meanOut = new double[nSample + 1][4];
        double[][] iOut = new double[nSample + 1][n];
        double[] hePatchMean = new double[nSample + 1];
        double[] hePatchMin = new double[nSample + 1];
        double[] hePatchMax = new double[nSample + 1];
        double[] heFluidMin = new double[nSample + 1];
        double[] heFluidMax = new double[nSample + 1];
        double[] heRange = {0.0, 0.0};
        recordSample(0, 0.0, a, i, hi, field, tOut, meanOut, iOut,
                hePatchMean, hePatchMin, hePatchMax, heFluidMin, heFluidMax, heRange);

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

        BSimStepScheduler scheduler = new BSimStepScheduler(sim);
        scheduler.run(new BSimStepScheduler.Adapter() {
            @Override
            public void transport(BSimStepScheduler.Context context, double duration) {
                field.diffuse(duration);
                field.decay(duration);
            }

            @Override
            public void sampleState(BSimStepScheduler.Context context) {
                for (int c = 0; c < n; c++) {
                    heSampled[c] = field.getConc(vi[c], vj[c], vk[c]);
                }
            }

            @Override
            public void integrateModels(BSimStepScheduler.Context context) {
                int step = micro[0];
                double t = step * rkDt;
                fillTau(hiTape, cap, step, rkDt, t - DaninoSIOccupiedDF.TAU, hist, hi0, hTau);
                rhsFrozenHe(a, i, hi, heSampled, hTau, dA1, dI1, dHi1, r1);
                axpy(a, i, hi, dA1, dI1, dHi1, 0.5 * rkDt, aTmp, iTmp, hiTmp);
                fillTau(hiTape, cap, step, rkDt, t + 0.5 * rkDt - DaninoSIOccupiedDF.TAU, hist, hi0, hTau);
                rhsFrozenHe(aTmp, iTmp, hiTmp, heSampled, hTau, dA2, dI2, dHi2, r2);
                axpy(a, i, hi, dA2, dI2, dHi2, 0.5 * rkDt, aTmp, iTmp, hiTmp);
                rhsFrozenHe(aTmp, iTmp, hiTmp, heSampled, hTau, dA3, dI3, dHi3, r3);
                axpy(a, i, hi, dA3, dI3, dHi3, rkDt, aTmp, iTmp, hiTmp);
                fillTau(hiTape, cap, step, rkDt, t + rkDt - DaninoSIOccupiedDF.TAU, hist, hi0, hTau);
                rhsFrozenHe(aTmp, iTmp, hiTmp, heSampled, hTau, dA4, dI4, dHi4, r4);

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
                    memDeposit[c] = rkDt * vCell * memConc;
                    a[c] += (rkDt / 6.0) * (dA1[c] + 2.0 * dA2[c] + 2.0 * dA3[c] + dA4[c]);
                    i[c] += (rkDt / 6.0) * (dI1[c] + 2.0 * dI2[c] + 2.0 * dI3[c] + dI4[c]);
                    hi[c] += (rkDt / 6.0) * (dHi1[c] + 2.0 * dHi2[c] + 2.0 * dHi3[c] + dHi4[c]);
                }
                ledger.addWeighted(
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
                    field.transferQuantity(vi[c], vj[c], vk[c], -memDeposit[c]);
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
                    recordSample(s, completed * rkDt, a, i, hi, field, tOut, meanOut, iOut,
                            hePatchMean, hePatchMin, hePatchMax, heFluidMin, heFluidMax, heRange);
                }
                int progress = (int) Math.round(100.0 / rkDt);
                if (progress > 0 && completed % progress == 0) {
                    System.err.printf(java.util.Locale.US,
                            "    t=%.0f NOT_FIG4B I_mean=%.6g He_fluid_range=%.3g%n",
                            completed * rkDt, mean(i), heRange[1]);
                    System.err.flush();
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
        double remaining = totalAhlMass(hi, field);
        return new Trajectory(
                tOut, meanOut, iOut,
                hePatchMean, hePatchMin, hePatchMax, heFluidMin, heFluidMax,
                d, mu, d1Spatial, rkDt, ensemble, negative, ledger,
                remaining, field.getTransportLedger(),
                heRange[0], heRange[1]);
    }

    private static boolean[][][] buildIslandMask() {
        boolean[][][] fluid = new boolean[NX_ISLAND][NY_ISLAND][NZ_ISLAND];
        for (int j = PATCH_J0; j <= PATCH_J1; j++) {
            for (int i = PATCH_I0; i <= PATCH_I1; i++) {
                fluid[i][j][PATCH_K] = true;
            }
        }
        return fluid;
    }

    private static int countIslandFluid(boolean[][][] fluid) {
        int n = 0;
        for (int i = 0; i < NX_ISLAND; i++) {
            for (int j = 0; j < NY_ISLAND; j++) {
                for (int k = 0; k < NZ_ISLAND; k++) {
                    if (fluid[i][j][k]) {
                        n++;
                    }
                }
            }
        }
        return n;
    }

    private void placeHeldPack() {
        int cell = 0;
        int nVox = 0;
        for (int j = PATCH_J0; j <= PATCH_J1; j++) {
            for (int i = PATCH_I0; i <= PATCH_I1; i++) {
                patchI[nVox] = i;
                patchJ[nVox] = j;
                int nHere = nVox < 10 ? 5 : 4;
                double cx = (i + 0.5) * DaninoPocketGeometry.DX_UM;
                double cy = (j + 0.5) * DaninoPocketGeometry.DY_UM;
                double cz = (PATCH_K + 0.5) * DaninoPocketGeometry.DZ_UM;
                double[][] off = nHere == 5
                        ? new double[][] {{0, 0}, {0.5, 0.5}, {0.5, -0.5}, {-0.5, 0.5}, {-0.5, -0.5}}
                        : new double[][] {{0.5, 0.5}, {0.5, -0.5}, {-0.5, 0.5}, {-0.5, -0.5}};
                for (int p = 0; p < nHere; p++) {
                    vi[cell] = i;
                    vj[cell] = j;
                    vk[cell] = PATCH_K;
                    xBsim[cell] = cx + off[p][0];
                    yBsim[cell] = cy + off[p][1];
                    zBsim[cell] = cz;
                    if ((int) (xBsim[cell] / DaninoPocketGeometry.DX_UM) != i
                            || (int) (yBsim[cell] / DaninoPocketGeometry.DY_UM) != j) {
                        throw new IllegalStateException(
                                "centre left its voxel NOT_FIG4B cell=" + cell);
                    }
                    cell++;
                }
                nVox++;
            }
        }
        if (cell != n || nVox != N_PATCH_VOXELS) {
            throw new IllegalStateException(
                    "held-pack placement N=" + cell + " nVox=" + nVox + " NOT_FIG4B");
        }
        for (int c = 0; c < n; c++) {
            if (!fluid[vi[c]][vj[c]][vk[c]]) {
                throw new IllegalStateException("packed centre in solid voxel NOT_FIG4B");
            }
        }
    }

    private void rhsFrozenHe(
            double[] a,
            double[] i,
            double[] hi,
            double[] he,
            double[] hTau,
            double[] dA,
            double[] dI,
            double[] dHi,
            Rates rates) {
        double dens = 1.0 - Math.pow(d / DaninoSIOccupiedDF.D0, 4.0);
        double synth = 0.0;
        double gLoss = 0.0;
        for (int c = 0; c < n; c++) {
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

    private double totalAhlMass(double[] hi, BSimTransportField field) {
        double s = 0.0;
        for (double h : hi) {
            s += h * vCell;
        }
        return s + field.totalFluidQuantity();
    }

    private void recordSample(
            int s,
            double t,
            double[] a,
            double[] i,
            double[] hi,
            BSimTransportField field,
            double[] tOut,
            double[][] meanOut,
            double[][] iOut,
            double[] hePatchMean,
            double[] hePatchMin,
            double[] hePatchMax,
            double[] heFluidMin,
            double[] heFluidMax,
            double[] heRange) {
        tOut[s] = t;
        meanOut[s][0] = mean(a);
        meanOut[s][1] = mean(i);
        meanOut[s][2] = mean(hi);
        System.arraycopy(i, 0, iOut[s], 0, n);
        double pMin = Double.POSITIVE_INFINITY;
        double pMax = Double.NEGATIVE_INFINITY;
        double pSum = 0.0;
        for (int v = 0; v < N_PATCH_VOXELS; v++) {
            double c = field.getConc(patchI[v], patchJ[v], PATCH_K);
            pMin = Math.min(pMin, c);
            pMax = Math.max(pMax, c);
            pSum += c;
        }
        hePatchMean[s] = pSum / N_PATCH_VOXELS;
        hePatchMin[s] = pMin;
        hePatchMax[s] = pMax;
        meanOut[s][3] = hePatchMean[s];
        double fMin = Double.POSITIVE_INFINITY;
        double fMax = Double.NEGATIVE_INFINITY;
        for (int ii = 0; ii < NX_ISLAND; ii++) {
            for (int jj = 0; jj < NY_ISLAND; jj++) {
                if (!fluid[ii][jj][PATCH_K]) {
                    continue;
                }
                double c = field.getConc(ii, jj, PATCH_K);
                fMin = Math.min(fMin, c);
                fMax = Math.max(fMax, c);
            }
        }
        heFluidMin[s] = fMin;
        heFluidMax[s] = fMax;
        heRange[0] = Math.max(heRange[0], pMax - pMin);
        heRange[1] = Math.max(heRange[1], fMax - fMin);
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

    private static final class Rates {
        double synth;
        double gammaH;
    }

    public static final class Trajectory {
        public final double[] t;
        /** Mean (A, I, Hi) and patch-mean He. */
        public final double[][] Y;
        public final double[][] I;
        public final double[] hePatchMean;
        public final double[] hePatchMin;
        public final double[] hePatchMax;
        public final double[] heFluidMin;
        public final double[] heFluidMax;
        public final double d;
        public final double mu;
        public final double d1Spatial;
        public final double rkDt;
        public final String ensemble;
        public final boolean negativeState;
        public final DaninoSIC1AhlLedger ledger;
        public final double remainingMass;
        public final BSimTransportLedger extra;
        public final double hePatchRangeMax;
        public final double heFluidRangeMax;

        Trajectory(
                double[] t,
                double[][] Y,
                double[][] I,
                double[] hePatchMean,
                double[] hePatchMin,
                double[] hePatchMax,
                double[] heFluidMin,
                double[] heFluidMax,
                double d,
                double mu,
                double d1Spatial,
                double rkDt,
                String ensemble,
                boolean negativeState,
                DaninoSIC1AhlLedger ledger,
                double remainingMass,
                BSimTransportLedger extra,
                double hePatchRangeMax,
                double heFluidRangeMax) {
            this.t = t;
            this.Y = Y;
            this.I = I;
            this.hePatchMean = hePatchMean;
            this.hePatchMin = hePatchMin;
            this.hePatchMax = hePatchMax;
            this.heFluidMin = heFluidMin;
            this.heFluidMax = heFluidMax;
            this.d = d;
            this.mu = mu;
            this.d1Spatial = d1Spatial;
            this.rkDt = rkDt;
            this.ensemble = ensemble;
            this.negativeState = negativeState;
            this.ledger = ledger;
            this.remainingMass = remainingMass;
            this.extra = extra;
            this.hePatchRangeMax = hePatchRangeMax;
            this.heFluidRangeMax = heFluidRangeMax;
        }

        public double[] luxIMean() {
            double[] out = new double[Y.length];
            for (int n = 0; n < Y.length; n++) {
                out[n] = Y[n][1];
            }
            return out;
        }

        public double residual() {
            return ledger.residual(remainingMass, extra);
        }

        public double relativeResidual() {
            return ledger.relativeResidual(remainingMass, extra);
        }
    }
}
