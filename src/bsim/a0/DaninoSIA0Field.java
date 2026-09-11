package bsim.a0;

import bsim.BSim;
import bsim.BSimRandom;
import bsim.BSimStepScheduler;
import bsim.c1c.DaninoSIC1cFilledPocket;
import bsim.transport.BSimTransportField;
import bsim.transport.BSimTransportLedger;

/**
 * A0 field operator on the C1c filled-pocket mask. Bacteria OFF.
 * A0_IDEAL_SOURCE. NOT_FIG4B. Not C1_OPEN_DILUTE.
 */
public final class DaninoSIA0Field {

    public static final String NOT_FIG4B = "NOT_FIG4B";
    public static final String STATUS = "A0_IDEAL_SOURCE";

    private final BSimRandom unusedRng;
    private final boolean[][][] fluid;

    public DaninoSIA0Field(BSimRandom rng) {
        this.unusedRng = rng;
        this.fluid = buildPocketMask();
    }

    public BSimRandom unusedRng() {
        return unusedRng;
    }

    /**
     * Bacteria-OFF pocket field with a prescribed source.
     *
     * @param acModule if true, deposit through {@link A0IdealSource}; if false,
     *                 N0 oracle uses the inline frozen J_max and voxel indices
     *                 and does not call {@link A0IdealSource}.
     */
    public Trajectory run(
            A0IdealSource.Command command,
            boolean acModule,
            double tEnd,
            double sampleDt,
            double rkDt) {
        if (!(rkDt > 0.0) || Math.abs(tEnd / rkDt - Math.round(tEnd / rkDt)) > 1e-9) {
            throw new IllegalArgumentException("rkDt must divide tEnd A0_IDEAL_SOURCE NOT_FIG4B");
        }
        int nMicro = (int) Math.round(tEnd / rkDt);
        int nSample = (int) Math.round(tEnd / sampleDt);
        if (Math.abs(nSample * sampleDt - tEnd) > 1e-12) {
            throw new IllegalArgumentException("sampleDt must divide tEnd A0_IDEAL_SOURCE NOT_FIG4B");
        }
        int stride = (int) Math.round(sampleDt / rkDt);
        if (Math.abs(stride * rkDt - sampleDt) > 1e-12) {
            throw new IllegalArgumentException("rkDt must divide sampleDt A0_IDEAL_SOURCE NOT_FIG4B");
        }

        A0IdealSource ac = acModule ? new A0IdealSource(command) : null;
        double jMax = A0IdealSource.J_MAX;
        double tOn = command.tOn;
        double tOff = command.tOff;
        boolean zero = command == A0IdealSource.Command.U_ZERO;

        BSim sim = new BSim();
        sim.setDt(rkDt);
        sim.setSimulationTime(tEnd);
        sim.setBound(
                DaninoSIC1cFilledPocket.BOUND_X,
                DaninoSIC1cFilledPocket.BOUND_Y,
                DaninoSIC1cFilledPocket.BOUND_Z);
        sim.setSolid(true, true, true);
        sim.setLeaky(false, false, false, false, false, false);
        if (unusedRng != null) {
            sim.setRandom(unusedRng);
        }

        BSimTransportField field = new BSimTransportField(
                sim,
                new int[] {DaninoSIC1cFilledPocket.NX, DaninoSIC1cFilledPocket.NY, DaninoSIC1cFilledPocket.NZ},
                DaninoSIC1cFilledPocket.D1_SPATIAL,
                0.40,
                fluid);
        field.getTransportLedger().reset();

        double[] tOut = new double[nSample + 1];
        double[] heAc = new double[nSample + 1];
        double[] heFar = new double[nSample + 1];
        double[] heMean = new double[nSample + 1];
        double[] cx = new double[nSample + 1];
        double[] cy = new double[nSample + 1];
        double[] cmdMass = new double[nSample + 1];
        double[] srcMass = new double[nSample + 1];
        int[] micro = {0};
        double[] inlineCommanded = {0.0};

        recordField(0, 0.0, field, ac, inlineCommanded[0], tOut, heAc, heFar, heMean, cx, cy, cmdMass, srcMass);

        BSimStepScheduler scheduler = new BSimStepScheduler(sim);
        scheduler.run(new BSimStepScheduler.Adapter() {
            @Override
            public void transport(BSimStepScheduler.Context context, double duration) {
                field.diffuse(duration);
                field.decay(duration);
            }

            @Override
            public void depositFluxes(BSimStepScheduler.Context context) {
                double tStart = context.getStartTime();
                double dt = context.getDt();
                if (acModule) {
                    ac.deposit(field, tStart, dt);
                } else if (!zero) {
                    double dM = jMax * A0IdealSource.overlap(tStart, dt, tOn, tOff);
                    if (dM > 0.0) {
                        field.addQuantity(A0IdealSource.I_AC, A0IdealSource.J_AC, A0IdealSource.K_AC, dM);
                        inlineCommanded[0] += dM;
                    }
                }
                micro[0]++;
            }

            @Override
            public void observe(BSimStepScheduler.Context context) {
                if (context.isInitialObservation()) {
                    return;
                }
                int completed = context.getCompletedUpdates();
                if (completed % stride == 0) {
                    int s = completed / stride;
                    recordField(s, completed * rkDt, field, ac, inlineCommanded[0],
                            tOut, heAc, heFar, heMean, cx, cy, cmdMass, srcMass);
                }
            }
        });

        if (micro[0] != nMicro) {
            throw new IllegalStateException(
                    "scheduler steps " + micro[0] + " != " + nMicro + " A0_IDEAL_SOURCE NOT_FIG4B");
        }
        double commanded = acModule ? ac.commandedMass() : inlineCommanded[0];
        BSimTransportLedger ledger = field.getTransportLedger();
        return new Trajectory(
                tOut, heAc, heFar, heMean, cx, cy, cmdMass, srcMass,
                commanded, ledger, field.totalFluidQuantity(), acModule, command.name());
    }

    private static boolean[][][] buildPocketMask() {
        boolean[][][] mask = new boolean[DaninoSIC1cFilledPocket.NX][DaninoSIC1cFilledPocket.NY][DaninoSIC1cFilledPocket.NZ];
        for (int i = 0; i < DaninoSIC1cFilledPocket.NX; i++) {
            for (int j = 0; j < DaninoSIC1cFilledPocket.NY; j++) {
                mask[i][j][0] = true;
            }
        }
        return mask;
    }

    private static void recordField(
            int s,
            double t,
            BSimTransportField field,
            A0IdealSource ac,
            double inlineCommanded,
            double[] tOut,
            double[] heAc,
            double[] heFar,
            double[] heMean,
            double[] cx,
            double[] cy,
            double[] cmdMass,
            double[] srcMass) {
        tOut[s] = t;
        heAc[s] = field.getConc(A0IdealSource.I_AC, A0IdealSource.J_AC, A0IdealSource.K_AC);
        heFar[s] = field.getConc(A0IdealSource.I_FAR, A0IdealSource.J_FAR, A0IdealSource.K_FAR);
        double mass = 0.0;
        double mx = 0.0;
        double my = 0.0;
        double sumC = 0.0;
        for (int i = 0; i < DaninoSIC1cFilledPocket.NX; i++) {
            for (int j = 0; j < DaninoSIC1cFilledPocket.NY; j++) {
                double c = field.getConc(i, j, 0);
                sumC += c;
                double q = c * DaninoSIC1cFilledPocket.V_VOXEL;
                mass += q;
                mx += q * ((i + 0.5) * DaninoSIC1cFilledPocket.BOUND_X / DaninoSIC1cFilledPocket.NX);
                my += q * ((j + 0.5) * DaninoSIC1cFilledPocket.BOUND_Y / DaninoSIC1cFilledPocket.NY);
            }
        }
        heMean[s] = sumC / DaninoSIC1cFilledPocket.N_VOXELS;
        if (mass > 0.0) {
            cx[s] = mx / mass;
            cy[s] = my / mass;
        } else {
            cx[s] = Double.NaN;
            cy[s] = Double.NaN;
        }
        cmdMass[s] = ac != null ? ac.commandedMass() : inlineCommanded;
        srcMass[s] = field.getTransportLedger().getSourceAdded();
    }

    public static final class Trajectory {
        public final double[] t;
        public final double[] heAc;
        public final double[] heFar;
        public final double[] heMean;
        public final double[] cx;
        public final double[] cy;
        public final double[] cmdMass;
        public final double[] srcMass;
        public final double commandedMass;
        public final BSimTransportLedger ledger;
        public final double remainingMass;
        public final boolean acModule;
        public final String command;

        Trajectory(
                double[] t,
                double[] heAc,
                double[] heFar,
                double[] heMean,
                double[] cx,
                double[] cy,
                double[] cmdMass,
                double[] srcMass,
                double commandedMass,
                BSimTransportLedger ledger,
                double remainingMass,
                boolean acModule,
                String command) {
            this.t = t;
            this.heAc = heAc;
            this.heFar = heFar;
            this.heMean = heMean;
            this.cx = cx;
            this.cy = cy;
            this.cmdMass = cmdMass;
            this.srcMass = srcMass;
            this.commandedMass = commandedMass;
            this.ledger = ledger;
            this.remainingMass = remainingMass;
            this.acModule = acModule;
            this.command = command;
        }

        public double massReplayRel() {
            double mStar = Math.max(
                    Math.max(Math.abs(commandedMass), Math.abs(ledger.getSourceAdded())),
                    Math.max(A0IdealSource.DELTA_M_KICK, 1e-15));
            return Math.abs(commandedMass - ledger.getSourceAdded()) / mStar;
        }

        public double transportRel() {
            double r = ledger.residual(0.0, remainingMass);
            double mStar = Math.max(
                    Math.max(Math.abs(ledger.getSourceAdded()), Math.abs(remainingMass)),
                    1e-15);
            return Math.abs(r) / mStar;
        }
    }
}
