package bsim.a1;

import bsim.BSim;
import bsim.BSimRandom;
import bsim.BSimStepScheduler;
import bsim.a0.A0IdealSource;
import bsim.c1c.DaninoSIC1cFilledPocket;
import bsim.transport.BSimTransportField;
import bsim.transport.BSimTransportLedger;

/**
 * A1 bacteria-OFF field operator on the C1c filled-pocket mask.
 * A1_BOUNDED_TRANSDUCER. HYPOTHETICAL_DESIGN_ENVELOPE. NOT_FIG4B.
 */
public final class DaninoSIA1Field {

    private final BSimRandom unusedRng;
    private final boolean[][][] fluid;

    public DaninoSIA1Field(BSimRandom rng) {
        this.unusedRng = rng;
        this.fluid = new boolean[DaninoSIC1cFilledPocket.NX][DaninoSIC1cFilledPocket.NY][DaninoSIC1cFilledPocket.NZ];
        for (int i = 0; i < DaninoSIC1cFilledPocket.NX; i++) {
            for (int j = 0; j < DaninoSIC1cFilledPocket.NY; j++) {
                fluid[i][j][0] = true;
            }
        }
    }

    public BSimRandom unusedRng() {
        return unusedRng;
    }

    public Trajectory run(
            A1BoundedTransducer.Command command,
            double tEnd,
            double sampleDt,
            double rkDt) {
        if (!(rkDt > 0.0) || Math.abs(tEnd / rkDt - Math.round(tEnd / rkDt)) > 1e-9) {
            throw new IllegalArgumentException("rkDt must divide tEnd A1_BOUNDED_TRANSDUCER NOT_FIG4B");
        }
        int nMicro = (int) Math.round(tEnd / rkDt);
        int nSample = (int) Math.round(tEnd / sampleDt);
        int stride = (int) Math.round(sampleDt / rkDt);
        A1BoundedTransducer ac = new A1BoundedTransducer(command, false);

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
        double[] jS = new double[nSample + 1];
        double[] jA0 = new double[nSample + 1];
        double[] x1 = new double[nSample + 1];
        double[] m = new double[nSample + 1];
        double[] heAc = new double[nSample + 1];
        double[] heFar = new double[nSample + 1];
        double[] released = new double[nSample + 1];
        double[] src = new double[nSample + 1];
        int[] micro = {0};

        record(0, 0.0, field, ac, command, tOut, jS, jA0, x1, m, heAc, heFar, released, src);

        BSimStepScheduler scheduler = new BSimStepScheduler(sim);
        scheduler.run(new BSimStepScheduler.Adapter() {
            @Override
            public void transport(BSimStepScheduler.Context context, double duration) {
                field.diffuse(duration);
                field.decay(duration);
            }

            @Override
            public void depositFluxes(BSimStepScheduler.Context context) {
                ac.deposit(field, context.getStartTime(), context.getDt());
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
                    record(s, completed * rkDt, field, ac, command,
                            tOut, jS, jA0, x1, m, heAc, heFar, released, src);
                }
            }
        });

        if (micro[0] != nMicro) {
            throw new IllegalStateException("scheduler steps " + micro[0] + " != " + nMicro
                    + " A1_BOUNDED_TRANSDUCER NOT_FIG4B");
        }
        return new Trajectory(
                tOut, jS, jA0, x1, m, heAc, heFar, released, src,
                ac, field.getTransportLedger(), field.totalFluidQuantity());
    }

    private static void record(
            int s,
            double t,
            BSimTransportField field,
            A1BoundedTransducer ac,
            A1BoundedTransducer.Command command,
            double[] tOut,
            double[] jS,
            double[] jA0,
            double[] x1,
            double[] m,
            double[] heAc,
            double[] heFar,
            double[] released,
            double[] src) {
        tOut[s] = t;
        jS[s] = ac.lastJ();
        jA0[s] = A1BoundedTransducer.a0J(command, t);
        x1[s] = ac.x1();
        m[s] = ac.payload();
        heAc[s] = field.getConc(A0IdealSource.I_AC, A0IdealSource.J_AC, A0IdealSource.K_AC);
        heFar[s] = field.getConc(A0IdealSource.I_FAR, A0IdealSource.J_FAR, A0IdealSource.K_FAR);
        released[s] = ac.released();
        src[s] = field.getTransportLedger().getSourceAdded();
    }

    public static final class Trajectory {
        public final double[] t;
        public final double[] jS;
        public final double[] jA0;
        public final double[] x1;
        public final double[] payload;
        public final double[] heAc;
        public final double[] heFar;
        public final double[] released;
        public final double[] src;
        public final A1BoundedTransducer ac;
        public final BSimTransportLedger ledger;
        public final double remainingMass;

        Trajectory(
                double[] t,
                double[] jS,
                double[] jA0,
                double[] x1,
                double[] payload,
                double[] heAc,
                double[] heFar,
                double[] released,
                double[] src,
                A1BoundedTransducer ac,
                BSimTransportLedger ledger,
                double remainingMass) {
            this.t = t;
            this.jS = jS;
            this.jA0 = jA0;
            this.x1 = x1;
            this.payload = payload;
            this.heAc = heAc;
            this.heFar = heFar;
            this.released = released;
            this.src = src;
            this.ac = ac;
            this.ledger = ledger;
            this.remainingMass = remainingMass;
        }

        public double massReplayRel() {
            return rel(ac.released(), ledger.getSourceAdded());
        }

        public double payloadRel() {
            return rel(ac.payloadDrop(), ac.released());
        }

        public double transportRel() {
            double r = ledger.residual(0.0, remainingMass);
            return Math.abs(r) / Math.max(Math.max(Math.abs(ledger.getSourceAdded()), Math.abs(remainingMass)), 1e-15);
        }

        public double a0ContrastRel() {
            double m = 0.0;
            for (int i = 0; i < t.length; i++) {
                m = Math.max(m, Math.abs(jS[i] - jA0[i]) / A1BoundedTransducer.J_MAX);
            }
            return m;
        }

        private static double rel(double a, double b) {
            double mStar = Math.max(Math.max(Math.abs(a), Math.abs(b)),
                    Math.max(A0IdealSource.DELTA_M_KICK, 1e-15));
            return Math.abs(a - b) / mStar;
        }
    }
}
