package bsim.c1c;

import bsim.BSim;
import bsim.BSimRandom;
import bsim.BSimStepScheduler;
import bsim.c1.DaninoSIC1AhlLedger;
import bsim.circuit.DaninoSIOccupiedDF;
import bsim.p0.DaninoPocketGeometry;
import bsim.transport.BSimTransportField;
import bsim.transport.BSimTransportLedger;

/**
 * C1c filled-pocket identity: cells throughout the 100×100×1.65 µm pocket.
 * Bus voxels are chemical solids. D1_spatial=800. NOT_FIG4B.
 *
 * One delay-DDE per voxel, occupancy-weighted, algebraically identical to
 * N=11000 cells that share a voxel He and start identical. Not a retune.
 */
public final class DaninoSIC1cFilledPocket {

    public static final String OBJECT = DaninoSIOccupiedDF.OBJECT;
    public static final String OBJECT_STATUS = DaninoSIOccupiedDF.OBJECT_STATUS;
    public static final String NOT_FIG4B = "NOT_FIG4B";

    public static final double V_CELL = 1.5;
    public static final int NX = 50;
    public static final int NY = 50;
    public static final int NZ = 1;
    public static final int N_VOXELS = NX * NY * NZ;
    public static final int N_FIVE = 1000;
    public static final int N_CELLS = N_FIVE * 5 + (N_VOXELS - N_FIVE) * 4;
    public static final double V_VOXEL = DaninoPocketGeometry.DX_UM
            * DaninoPocketGeometry.DY_UM
            * DaninoPocketGeometry.DZ_UM;
    public static final double V_E_POCKET = N_VOXELS * V_VOXEL;
    public static final double D_POCKET = (N_CELLS * V_CELL) / (N_CELLS * V_CELL + V_E_POCKET);
    public static final double D1_SPATIAL = 800.0;
    public static final String D1_LABEL = "SI_TAKEN_SPATIAL";
    public static final String OPEN_EDGE_BC = "no-flux";
    public static final double BOUND_X = NX * DaninoPocketGeometry.DX_UM;
    public static final double BOUND_Y = NY * DaninoPocketGeometry.DY_UM;
    public static final double BOUND_Z = NZ * DaninoPocketGeometry.DZ_UM;

    private final int nCells;
    private final int nVox;
    private final double vCell;
    private final double d;
    private final double d1Spatial;
    private final BSimRandom unusedRng;
    private final boolean[][][] fluid;
    private final int[] nHere;
    private final int[] vi;
    private final int[] vj;
    private final int[] vk;

    public DaninoSIC1cFilledPocket(BSimRandom rng) {
        this.unusedRng = rng;
        this.nCells = N_CELLS;
        this.nVox = N_VOXELS;
        this.vCell = V_CELL;
        this.d = D_POCKET;
        this.d1Spatial = D1_SPATIAL;
        if (nCells != 11000) {
            throw new IllegalStateException("C1c N must be 11000 NOT_FIG4B, got " + nCells);
        }
        if (Math.abs(V_E_POCKET - 16500.0) > 1e-12) {
            throw new IllegalStateException("C1c Ve_pocket must be 16500 NOT_FIG4B, got " + V_E_POCKET);
        }
        if (Math.abs(d - 0.5) > 1e-15) {
            throw new IllegalStateException("C1c pocket d is " + d + " not 0.5 NOT_FIG4B");
        }
        if (Math.abs(BOUND_X - DaninoPocketGeometry.POCKET_LX_UM) > 0.0
                || Math.abs(BOUND_Y - DaninoPocketGeometry.POCKET_LY_UM) > 0.0
                || Math.abs(BOUND_Z - DaninoPocketGeometry.POCKET_LZ_UM) > 0.0) {
            throw new IllegalStateException("C1c pocket bounds must match P0 pocket NOT_FIG4B");
        }
        this.fluid = buildPocketMask();
        this.nHere = new int[nVox];
        this.vi = new int[nVox];
        this.vj = new int[nVox];
        this.vk = new int[nVox];
        placeFilledPocket();
        BSim probe = new BSim();
        probe.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        probe.setSolid(true, true, true);
        BSimTransportField field = new BSimTransportField(
                probe, new int[] {NX, NY, NZ}, D1_SPATIAL, 0.0, fluid);
        double vox = field.getBox()[0] * field.getBox()[1] * field.getBox()[2];
        if (Math.abs(vox - V_VOXEL) > 1e-12) {
            throw new IllegalStateException("C1c voxel volume " + vox + " != " + V_VOXEL + " NOT_FIG4B");
        }
        if (Math.abs(field.totalFluidQuantity()) > 0.0) {
            throw new IllegalStateException("probe field must start empty NOT_FIG4B");
        }
        int nFluid = 0;
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                if (fluid[i][j][0]) {
                    nFluid++;
                }
            }
        }
        if (nFluid != N_VOXELS) {
            throw new IllegalStateException("C1c fluid voxels " + nFluid + " != " + N_VOXELS + " NOT_FIG4B");
        }
    }

    public BSimRandom unusedRng() {
        return unusedRng;
    }

    public int nCells() {
        return nCells;
    }

    public int nVoxels() {
        return nVox;
    }

    public double vCell() {
        return vCell;
    }

    public double vEPocket() {
        return V_E_POCKET;
    }

    public double d() {
        return d;
    }

    public double d1Spatial() {
        return d1Spatial;
    }

    public String openEdgeBc() {
        return OPEN_EDGE_BC;
    }

    /**
     * Optional scheduler flux-deposition hook. A0 uses this for the ideal
     * source. Null means C1c identity (no extra source). NOT_FIG4B.
     */
    public interface ExtraDeposit {
        void deposit(BSimTransportField field, double tStart, double dt);
    }

    /**
     * Interactive BSim-ticker handle. Same operator as {@link #integrate}, one
     * RK microstep per {@link LiveView#step()}. Viewer only; not a gate.
     * NOT_FIG4B.
     */
    public LiveView openLive(double mu, String ensemble, double rkDt, ExtraDeposit extra) {
        if (!(rkDt > 0.0)) {
            throw new IllegalArgumentException("rkDt must be positive NOT_FIG4B");
        }
        boolean kick = !"SI_BASAL_PERTURB".equals(ensemble);
        double a0 = 0.0;
        double i0 = kick ? 0.0 : 1.0;
        double hi0 = kick ? 0.05 : 0.0;
        double he0 = kick ? 0.05 : 0.0;
        double hist = hi0;

        double[] a = new double[nVox];
        double[] i = new double[nVox];
        double[] hi = new double[nVox];
        java.util.Arrays.fill(a, a0);
        java.util.Arrays.fill(i, i0);
        java.util.Arrays.fill(hi, hi0);

        BSim sim = new BSim();
        sim.setDt(rkDt);
        sim.setSimulationTime(1.0e9);
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);
        sim.setLeaky(false, false, false, false, false, false);
        if (unusedRng != null) {
            sim.setRandom(unusedRng);
        }

        BSimTransportField field = new BSimTransportField(
                sim, new int[] {NX, NY, NZ}, d1Spatial, mu, fluid);
        if (he0 > 0.0) {
            field.setConc(he0);
        }
        field.getTransportLedger().reset();

        int cap = (int) Math.round(DaninoSIOccupiedDF.TAU / rkDt) + 2;
        double[][] hiTape = new double[nVox][cap];
        for (int c = 0; c < nVox; c++) {
            hiTape[c][0] = hi[c];
        }
        return new LiveView(sim, field, extra, a, i, hi, hiTape, cap, rkDt, hist, hi0);
    }

    /** Live C1c state for the Processing window. Viewer only. NOT_FIG4B. */
    public final class LiveView {
        public final BSim sim;
        public final BSimTransportField field;
        public final double[] luxI;
        public final int[] voxelI;
        public final int[] voxelJ;
        public final int[] voxelK;
        private final ExtraDeposit extra;
        private final double[] a;
        private final double[] hi;
        private final double[][] hiTape;
        private final int cap;
        private final double rkDt;
        private final double hist;
        private final double hi0;
        private final double[] heSampled;
        private final double[] dA1;
        private final double[] dI1;
        private final double[] dHi1;
        private final double[] dA2;
        private final double[] dI2;
        private final double[] dHi2;
        private final double[] dA3;
        private final double[] dI3;
        private final double[] dHi3;
        private final double[] dA4;
        private final double[] dI4;
        private final double[] dHi4;
        private final double[] aTmp;
        private final double[] iTmp;
        private final double[] hiTmp;
        private final double[] hTau;
        private final double[] memDeposit;
        private final Rates r1 = new Rates();
        private final Rates r2 = new Rates();
        private final Rates r3 = new Rates();
        private final Rates r4 = new Rates();
        private int micro;

        private LiveView(
                BSim sim,
                BSimTransportField field,
                ExtraDeposit extra,
                double[] a,
                double[] i,
                double[] hi,
                double[][] hiTape,
                int cap,
                double rkDt,
                double hist,
                double hi0) {
            this.sim = sim;
            this.field = field;
            this.extra = extra;
            this.a = a;
            this.luxI = i;
            this.hi = hi;
            this.hiTape = hiTape;
            this.cap = cap;
            this.rkDt = rkDt;
            this.hist = hist;
            this.hi0 = hi0;
            this.voxelI = vi;
            this.voxelJ = vj;
            this.voxelK = vk;
            this.heSampled = new double[nVox];
            this.dA1 = new double[nVox];
            this.dI1 = new double[nVox];
            this.dHi1 = new double[nVox];
            this.dA2 = new double[nVox];
            this.dI2 = new double[nVox];
            this.dHi2 = new double[nVox];
            this.dA3 = new double[nVox];
            this.dI3 = new double[nVox];
            this.dHi3 = new double[nVox];
            this.dA4 = new double[nVox];
            this.dI4 = new double[nVox];
            this.dHi4 = new double[nVox];
            this.aTmp = new double[nVox];
            this.iTmp = new double[nVox];
            this.hiTmp = new double[nVox];
            this.hTau = new double[nVox];
            this.memDeposit = new double[nVox];
        }

        public double time() {
            return micro * rkDt;
        }

        public int nVoxels() {
            return nVox;
        }

        /** One C1c RK microstep. Same phase order as the identity scheduler. */
        public void step() {
            double half = 0.5 * rkDt;
            field.diffuse(half);
            field.decay(half);
            for (int c = 0; c < nVox; c++) {
                heSampled[c] = field.getConc(vi[c], vj[c], vk[c]);
            }
            int step = micro;
            double t = step * rkDt;
            fillTau(hiTape, cap, step, rkDt, t - DaninoSIOccupiedDF.TAU, hist, hi0, hTau);
            rhsFrozenHe(a, luxI, hi, heSampled, hTau, dA1, dI1, dHi1, r1);
            axpy(a, luxI, hi, dA1, dI1, dHi1, half, aTmp, iTmp, hiTmp);
            fillTau(hiTape, cap, step, rkDt, t + half - DaninoSIOccupiedDF.TAU, hist, hi0, hTau);
            rhsFrozenHe(aTmp, iTmp, hiTmp, heSampled, hTau, dA2, dI2, dHi2, r2);
            axpy(a, luxI, hi, dA2, dI2, dHi2, half, aTmp, iTmp, hiTmp);
            rhsFrozenHe(aTmp, iTmp, hiTmp, heSampled, hTau, dA3, dI3, dHi3, r3);
            axpy(a, luxI, hi, dA3, dI3, dHi3, rkDt, aTmp, iTmp, hiTmp);
            fillTau(hiTape, cap, step, rkDt, t + rkDt - DaninoSIOccupiedDF.TAU, hist, hi0, hTau);
            rhsFrozenHe(aTmp, iTmp, hiTmp, heSampled, hTau, dA4, dI4, dHi4, r4);
            for (int c = 0; c < nVox; c++) {
                double hiS1 = hi[c];
                double hiS2 = hi[c] + half * dHi1[c];
                double hiS3 = hi[c] + half * dHi2[c];
                double hiS4 = hi[c] + rkDt * dHi3[c];
                double memConc = (DaninoSIOccupiedDF.D_MEM / 6.0) * (
                        (heSampled[c] - hiS1)
                                + 2.0 * (heSampled[c] - hiS2)
                                + 2.0 * (heSampled[c] - hiS3)
                                + (heSampled[c] - hiS4));
                memDeposit[c] = rkDt * nHere[c] * vCell * memConc;
                a[c] += (rkDt / 6.0) * (dA1[c] + 2.0 * dA2[c] + 2.0 * dA3[c] + dA4[c]);
                luxI[c] += (rkDt / 6.0) * (dI1[c] + 2.0 * dI2[c] + 2.0 * dI3[c] + dI4[c]);
                hi[c] += (rkDt / 6.0) * (dHi1[c] + 2.0 * dHi2[c] + 2.0 * dHi3[c] + dHi4[c]);
            }
            for (int c = 0; c < nVox; c++) {
                field.transferQuantity(vi[c], vj[c], vk[c], -memDeposit[c]);
            }
            if (extra != null) {
                extra.deposit(field, t, rkDt);
            }
            field.diffuse(half);
            field.decay(half);
            for (int c = 0; c < nVox; c++) {
                hiTape[c][(step + 1) % cap] = hi[c];
            }
            micro = step + 1;
        }
    }

    /** Optional He/I probes at the A0 voxel and a far voxel. NOT_FIG4B. */
    public static final class Probe {
        public final int iAc;
        public final int jAc;
        public final int kAc;
        public final int iFar;
        public final int jFar;
        public final int kFar;

        public Probe(int iAc, int jAc, int kAc, int iFar, int jFar, int kFar) {
            this.iAc = iAc;
            this.jAc = jAc;
            this.kAc = kAc;
            this.iFar = iFar;
            this.jFar = jFar;
            this.kFar = kFar;
        }
    }

    /**
     * Filled-pocket field/membrane operator. He(0) and mu hit pocket fluid
     * only. Bus is not constructed. D1=800, no-flux faces. NOT_FIG4B.
     */
    public Trajectory integrate(
            double mu,
            String ensemble,
            double tEnd,
            double sampleDt,
            double rkDt) {
        return integrate(mu, ensemble, tEnd, sampleDt, rkDt, null, null);
    }

    /**
     * Same operator as {@link #integrate(double, String, double, double, double)}
     * with an optional extra source in {@code depositFluxes} and optional
     * voxel probes. J=0 / null extra must recover C1c identity. NOT_FIG4B.
     */
    public Trajectory integrate(
            double mu,
            String ensemble,
            double tEnd,
            double sampleDt,
            double rkDt,
            ExtraDeposit extra,
            Probe probe) {
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

        double[] a = new double[nVox];
        double[] i = new double[nVox];
        double[] hi = new double[nVox];
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

        BSimTransportField field = new BSimTransportField(
                sim, new int[] {NX, NY, NZ}, d1Spatial, mu, fluid);
        if (he0 > 0.0) {
            field.setConc(he0);
        }
        field.getTransportLedger().reset();

        DaninoSIC1AhlLedger ledger = new DaninoSIC1AhlLedger();
        ledger.reset(totalAhlMass(hi, field));

        int cap = (int) Math.round(DaninoSIOccupiedDF.TAU / rkDt) + 2;
        double[][] hiTape = new double[nVox][cap];
        for (int c = 0; c < nVox; c++) {
            hiTape[c][0] = hi[c];
        }

        double[] tOut = new double[nSample + 1];
        double[][] meanOut = new double[nSample + 1][4];
        double[] heMin = new double[nSample + 1];
        double[] heMax = new double[nSample + 1];
        double[] heRange = {0.0};
        double[] heAc = probe == null ? null : new double[nSample + 1];
        double[] heFar = probe == null ? null : new double[nSample + 1];
        double[] iAc = probe == null ? null : new double[nSample + 1];
        int acIndexFound = -1;
        if (probe != null) {
            for (int c = 0; c < nVox; c++) {
                if (vi[c] == probe.iAc && vj[c] == probe.jAc && vk[c] == probe.kAc) {
                    acIndexFound = c;
                    break;
                }
            }
            if (acIndexFound < 0) {
                throw new IllegalStateException("A0 AC voxel is not a C1c fluid voxel NOT_FIG4B");
            }
        }
        final int acIndex = acIndexFound;
        recordSample(0, 0.0, a, i, hi, field, tOut, meanOut, heMin, heMax, heRange,
                heAc, heFar, iAc, probe, acIndex);

        double[] heSampled = new double[nVox];
        double[] dA1 = new double[nVox];
        double[] dI1 = new double[nVox];
        double[] dHi1 = new double[nVox];
        double[] dA2 = new double[nVox];
        double[] dI2 = new double[nVox];
        double[] dHi2 = new double[nVox];
        double[] dA3 = new double[nVox];
        double[] dI3 = new double[nVox];
        double[] dHi3 = new double[nVox];
        double[] dA4 = new double[nVox];
        double[] dI4 = new double[nVox];
        double[] dHi4 = new double[nVox];
        double[] aTmp = new double[nVox];
        double[] iTmp = new double[nVox];
        double[] hiTmp = new double[nVox];
        double[] hTau = new double[nVox];
        double[] memDeposit = new double[nVox];
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
                for (int c = 0; c < nVox; c++) {
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

                for (int c = 0; c < nVox; c++) {
                    double hiS1 = hi[c];
                    double hiS2 = hi[c] + 0.5 * rkDt * dHi1[c];
                    double hiS3 = hi[c] + 0.5 * rkDt * dHi2[c];
                    double hiS4 = hi[c] + rkDt * dHi3[c];
                    double memConc = (DaninoSIOccupiedDF.D_MEM / 6.0) * (
                            (heSampled[c] - hiS1)
                                    + 2.0 * (heSampled[c] - hiS2)
                                    + 2.0 * (heSampled[c] - hiS3)
                                    + (heSampled[c] - hiS4));
                    memDeposit[c] = rkDt * nHere[c] * vCell * memConc;
                    a[c] += (rkDt / 6.0) * (dA1[c] + 2.0 * dA2[c] + 2.0 * dA3[c] + dA4[c]);
                    i[c] += (rkDt / 6.0) * (dI1[c] + 2.0 * dI2[c] + 2.0 * dI3[c] + dI4[c]);
                    hi[c] += (rkDt / 6.0) * (dHi1[c] + 2.0 * dHi2[c] + 2.0 * dHi3[c] + dHi4[c]);
                }
                ledger.addWeighted(
                        (rkDt / 6.0) * (r1.synth + 2.0 * r2.synth + 2.0 * r3.synth + r4.synth),
                        (rkDt / 6.0) * (r1.gammaH + 2.0 * r2.gammaH + 2.0 * r3.gammaH + r4.gammaH),
                        sum(memDeposit));
                for (int c = 0; c < nVox; c++) {
                    hiTape[c][(step + 1) % cap] = hi[c];
                }
                micro[0] = step + 1;
            }

            @Override
            public void depositFluxes(BSimStepScheduler.Context context) {
                for (int c = 0; c < nVox; c++) {
                    field.transferQuantity(vi[c], vj[c], vk[c], -memDeposit[c]);
                }
                if (extra != null) {
                    extra.deposit(field, context.getStartTime(), context.getDt());
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
                    recordSample(s, completed * rkDt, a, i, hi, field, tOut, meanOut, heMin, heMax, heRange,
                            heAc, heFar, iAc, probe, acIndex);
                }
                int progress = (int) Math.round(100.0 / rkDt);
                if (progress > 0 && completed % progress == 0) {
                    System.err.printf(java.util.Locale.US,
                            "    t=%.0f NOT_FIG4B I_mean=%.6g He_pocket_range=%.3g%n",
                            completed * rkDt, weightedMean(i), heRange[0]);
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
                tOut, meanOut, heMin, heMax,
                d, mu, d1Spatial, rkDt, ensemble, negative, ledger,
                remaining, field.getTransportLedger(), heRange[0],
                heAc, heFar, iAc);
    }

    private static boolean[][][] buildPocketMask() {
        boolean[][][] mask = new boolean[NX][NY][NZ];
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                mask[i][j][0] = true;
            }
        }
        return mask;
    }

    private void placeFilledPocket() {
        int vox = 0;
        int placed = 0;
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                nHere[vox] = vox < N_FIVE ? 5 : 4;
                vi[vox] = i;
                vj[vox] = j;
                vk[vox] = 0;
                placed += nHere[vox];
                vox++;
            }
        }
        if (vox != nVox || placed != nCells) {
            throw new IllegalStateException(
                    "filled-pocket placement vox=" + vox + " N=" + placed + " NOT_FIG4B");
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
        for (int c = 0; c < nVox; c++) {
            double p = DaninoSIOccupiedDF.production(hTau[c]);
            double denomF = 1.0 + DaninoSIOccupiedDF.F * (a[c] + i[c]);
            dA[c] = DaninoSIOccupiedDF.C_A * dens * p - DaninoSIOccupiedDF.GAMMA_A * a[c] / denomF;
            dI[c] = DaninoSIOccupiedDF.C_I * dens * p - DaninoSIOccupiedDF.GAMMA_I * i[c] / denomF;
            double synthConc = DaninoSIOccupiedDF.B * i[c] / (1.0 + DaninoSIOccupiedDF.K * i[c]);
            double enzConc = (DaninoSIOccupiedDF.GAMMA_H * a[c] * hi[c])
                    / (1.0 + DaninoSIOccupiedDF.G * a[c]);
            double memConc = DaninoSIOccupiedDF.D_MEM * (he[c] - hi[c]);
            dHi[c] = synthConc - enzConc + memConc;
            double mass = nHere[c] * vCell;
            synth += mass * synthConc;
            gLoss += mass * enzConc;
        }
        rates.synth = synth;
        rates.gammaH = gLoss;
    }

    private double totalAhlMass(double[] hi, BSimTransportField field) {
        double s = 0.0;
        for (int c = 0; c < nVox; c++) {
            s += hi[c] * vCell * nHere[c];
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
            double[] heMinOut,
            double[] heMaxOut,
            double[] heRange,
            double[] heAc,
            double[] heFar,
            double[] iAc,
            Probe probe,
            int acIndex) {
        tOut[s] = t;
        meanOut[s][0] = weightedMean(a);
        meanOut[s][1] = weightedMean(i);
        meanOut[s][2] = weightedMean(hi);
        double lo = Double.POSITIVE_INFINITY;
        double hiVal = Double.NEGATIVE_INFINITY;
        double sum = 0.0;
        for (int ii = 0; ii < NX; ii++) {
            for (int jj = 0; jj < NY; jj++) {
                double c = field.getConc(ii, jj, 0);
                lo = Math.min(lo, c);
                hiVal = Math.max(hiVal, c);
                sum += c;
            }
        }
        meanOut[s][3] = sum / nVox;
        heMinOut[s] = lo;
        heMaxOut[s] = hiVal;
        heRange[0] = Math.max(heRange[0], hiVal - lo);
        if (probe != null) {
            heAc[s] = field.getConc(probe.iAc, probe.jAc, probe.kAc);
            heFar[s] = field.getConc(probe.iFar, probe.jFar, probe.kFar);
            iAc[s] = i[acIndex];
        }
    }

    private double weightedMean(double[] x) {
        double s = 0.0;
        for (int c = 0; c < nVox; c++) {
            s += nHere[c] * x[c];
        }
        return s / nCells;
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
        if (tDelay < 0.0) {
            java.util.Arrays.fill(hTau, hiHistory);
            return;
        }
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
        /** Occupancy-weighted mean (A, I, Hi) and pocket-mean He. */
        public final double[][] Y;
        public final double[] hePocketMin;
        public final double[] hePocketMax;
        public final double d;
        public final double mu;
        public final double d1Spatial;
        public final double rkDt;
        public final String ensemble;
        public final boolean negativeState;
        public final DaninoSIC1AhlLedger ledger;
        public final double remainingMass;
        public final BSimTransportLedger extra;
        public final double hePocketRangeMax;
        /** He at the A0 voxel; null on C1c identity (no probe). */
        public final double[] heAc;
        /** He at the far voxel; null on C1c identity (no probe). */
        public final double[] heFar;
        /** Local I at the A0 voxel; null on C1c identity (no probe). */
        public final double[] iAc;

        Trajectory(
                double[] t,
                double[][] Y,
                double[] hePocketMin,
                double[] hePocketMax,
                double d,
                double mu,
                double d1Spatial,
                double rkDt,
                String ensemble,
                boolean negativeState,
                DaninoSIC1AhlLedger ledger,
                double remainingMass,
                BSimTransportLedger extra,
                double hePocketRangeMax) {
            this(t, Y, hePocketMin, hePocketMax, d, mu, d1Spatial, rkDt, ensemble,
                    negativeState, ledger, remainingMass, extra, hePocketRangeMax,
                    null, null, null);
        }

        Trajectory(
                double[] t,
                double[][] Y,
                double[] hePocketMin,
                double[] hePocketMax,
                double d,
                double mu,
                double d1Spatial,
                double rkDt,
                String ensemble,
                boolean negativeState,
                DaninoSIC1AhlLedger ledger,
                double remainingMass,
                BSimTransportLedger extra,
                double hePocketRangeMax,
                double[] heAc,
                double[] heFar,
                double[] iAc) {
            this.t = t;
            this.Y = Y;
            this.hePocketMin = hePocketMin;
            this.hePocketMax = hePocketMax;
            this.d = d;
            this.mu = mu;
            this.d1Spatial = d1Spatial;
            this.rkDt = rkDt;
            this.ensemble = ensemble;
            this.negativeState = negativeState;
            this.ledger = ledger;
            this.remainingMass = remainingMass;
            this.extra = extra;
            this.hePocketRangeMax = hePocketRangeMax;
            this.heAc = heAc;
            this.heFar = heFar;
            this.iAc = iAc;
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
