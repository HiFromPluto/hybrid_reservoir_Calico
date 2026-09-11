package bsim.laneb;

import bsim.BSim;
import bsim.BSimRandom;
import bsim.BSimStepScheduler;
import bsim.transport.BSimConservativeTransport;
import bsim.transport.BSimTransportField;
import bsim.transport.BSimTransportLedger;

import javax.vecmath.Vector3d;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Millimetre FLOW=0 dish on the N0 kernel. B0_MOTILITY_MEMORY.
 * Not Fig. 4b. Not C1c. Not HybridDish Hill.
 */
public final class LaneBDish {

    public enum Arm {
        MSD_MOTILE,
        MIX_BLOB,
        MIX_STRIPE,
        PAINT_MOTILE,
        PAINT_OFF,
        PAINT_FIELD,
        PAINT_SILENT,
        STRIPE_MOTILE,
        STRIPE_OFF,
        STRIPE_FIELD,
        STRIPE_SILENT,
        SUSTAIN_MOTILE,
        SUSTAIN_OFF,
        SUSTAIN_FIELD,
        SUSTAIN_SILENT,
        THERMAL
    }

    public static final class Result {
        public final String name;
        public final Arm arm;
        public final double duration;
        public final double[] time;
        public final double[] msdXy;
        public final double[] kappaCell;
        public final double[] kappaField;
        public final double[] autoCorr;
        public final double mu2d;
        public final double thermalD2d;
        public final double attResidualRel;
        public final double repResidualRel;
        public final double attMass;
        public final double repMass;
        public final double commandedAtt;
        public final double commandedRep;
        public final boolean ledgerPass;
        public final double[] meanAbsL;
        public final double pulseMeanAbsL;
        public final double[] deltaLEll;

        Result(String name, Arm arm, double duration, double[] time, double[] msdXy,
               double[] kappaCell, double[] kappaField, double[] autoCorr,
               double mu2d, double thermalD2d,
               double attResidualRel, double repResidualRel,
               double attMass, double repMass,
               double commandedAtt, double commandedRep, boolean ledgerPass,
               double[] meanAbsL, double pulseMeanAbsL, double[] deltaLEll) {
            this.name = name;
            this.arm = arm;
            this.duration = duration;
            this.time = time;
            this.msdXy = msdXy;
            this.kappaCell = kappaCell;
            this.kappaField = kappaField;
            this.autoCorr = autoCorr;
            this.mu2d = mu2d;
            this.thermalD2d = thermalD2d;
            this.attResidualRel = attResidualRel;
            this.repResidualRel = repResidualRel;
            this.attMass = attMass;
            this.repMass = repMass;
            this.commandedAtt = commandedAtt;
            this.commandedRep = commandedRep;
            this.ledgerPass = ledgerPass;
            this.meanAbsL = meanAbsL;
            this.pulseMeanAbsL = pulseMeanAbsL;
            this.deltaLEll = deltaLEll;
        }
    }

    private final BSimRandom rng;

    public LaneBDish(BSimRandom rng) {
        if (rng == null) throw new IllegalArgumentException("rng is required");
        this.rng = rng;
    }

    public Result run(Arm arm, double duration) {
        return run(arm, duration, LaneBIdentity.T_OFF);
    }

    public Result run(Arm arm, double duration, double tOff) {
        duration = LaneBIdentity.DT * Math.round(duration / LaneBIdentity.DT);
        BSim sim = new BSim();
        sim.setBound(LaneBIdentity.LX, LaneBIdentity.LY, LaneBIdentity.LZ);
        sim.setSolid(true, true, true);
        sim.setDt(LaneBIdentity.DT);
        sim.setSimulationTime(duration);
        sim.setTemperature(LaneBIdentity.T_K);
        sim.setVisc(LaneBIdentity.ETA);
        sim.setRandom(rng);

        boolean paint = isPaintArm(arm);
        boolean stripeArm = isStripeArm(arm);
        boolean sustainArm = isSustainArm(arm);
        boolean chemistry = paint || stripeArm || sustainArm;
        boolean cellsOn = arm != Arm.PAINT_FIELD && arm != Arm.STRIPE_FIELD
                && arm != Arm.SUSTAIN_FIELD;
        LaneBBacterium.Motility motility = motilityOf(arm);
        int nCells = nCellsOf(arm);
        boolean silent = arm == Arm.PAINT_SILENT || arm == Arm.STRIPE_SILENT
                || arm == Arm.SUSTAIN_SILENT;
        boolean stripeIc = stripeArm && !silent;
        boolean blob = arm == Arm.MIX_BLOB;
        boolean stripe = arm == Arm.MIX_STRIPE;

        int[] boxes = {LaneBIdentity.NX, LaneBIdentity.NY, LaneBIdentity.NZ};
        boolean[][][] fluid = BSimConservativeTransport.allFluid(boxes);
        BSimTransportField att = chemistry
                ? new BSimTransportField(sim, boxes, LaneBIdentity.D_ATT, LaneBIdentity.K_ATT, fluid)
                : null;
        BSimTransportField rep = chemistry
                ? new BSimTransportField(sim, boxes, LaneBIdentity.D_REP, LaneBIdentity.K_REP, fluid)
                : null;
        double initialAtt = 0.0;
        double initialRep = 0.0;
        if (stripeIc) {
            paintSquareWave(att, rep, LaneB2Identity.L_SLAB);
            initialAtt = att.totalFluidQuantity();
            initialRep = rep.totalFluidQuantity();
        }
        if (att != null) att.getTransportLedger().reset();
        if (rep != null) rep.getTransportLedger().reset();

        List<LaneBA0Source> sources = new ArrayList<LaneBA0Source>();
        if (paint) {
            double j = silent ? 0.0 : LaneBIdentity.J_MAX;
            double pulseOff = silent ? Double.NaN : tOff;
            for (int a = 0; a < LaneBIdentity.AC_COUNT; a++) {
                double[] xyz = LaneBIdentity.AC_XYZ[a];
                sources.add(new LaneBA0Source(
                        xyz[0], xyz[1], xyz[2], LaneBIdentity.AC_PAYLOAD[a],
                        j, 0.0, pulseOff));
            }
        }
        final LaneBSlabSource slab = (sustainArm && !silent)
                ? new LaneBSlabSource(LaneB3Identity.J_SLAB, 0.0, tOff)
                : null;

        LaneBBacterium[] cells = new LaneBBacterium[nCells];
        double[] x0 = new double[nCells];
        double[] y0 = new double[nCells];
        if (cellsOn) {
            for (int n = 0; n < nCells; n++) {
                Vector3d p = blob ? blobSpawn() : (stripe ? stripeSpawn() : interiorSpawn(arm));
                cells[n] = new LaneBBacterium(sim, p, rng, motility);
                x0[n] = cells[n].getPosition().x;
                y0[n] = cells[n].getPosition().y;
            }
        }

        int nSamp = (int) Math.round(duration / LaneBIdentity.SAMPLE_DT) + 1;
        double[] time = new double[nSamp];
        double[] msdXy = new double[nSamp];
        double[] kappaCell = new double[nSamp];
        double[] kappaField = new double[nSamp];
        double[] autoCorr = new double[nSamp];
        double[] meanAbsL = new double[nSamp];
        double[] deltaLEll = new double[nSamp];
        final double[][][] density0 = new double[1][][];
        final int[] samp = {0};
        final double pulseOff = tOff;
        final boolean stripeKappa = stripeArm || sustainArm;

        BSimStepScheduler scheduler = new BSimStepScheduler(sim);
        scheduler.addMechanicsEvent("move", LaneBIdentity.DT, context -> {
            for (LaneBBacterium cell : cells) {
                if (cell != null) cell.updatePosition();
            }
        });
        scheduler.run(new BSimStepScheduler.Adapter() {
            @Override
            public void transport(BSimStepScheduler.Context context, double half) {
                if (att != null) att.advance(half);
                if (rep != null) rep.advance(half);
            }

            @Override
            public void sampleState(BSimStepScheduler.Context context) {
                for (LaneBBacterium cell : cells) {
                    if (cell == null) continue;
                    double a = att == null ? 0.0 : att.getConc(cell.getPosition());
                    double r = rep == null ? 0.0 : rep.getConc(cell.getPosition());
                    cell.sampleLigand(a, r);
                }
            }

            @Override
            public void integrateModels(BSimStepScheduler.Context context) {
                for (LaneBBacterium cell : cells) {
                    if (cell != null) cell.action();
                }
            }

            @Override
            public void depositFluxes(BSimStepScheduler.Context context) {
                for (LaneBA0Source src : sources) {
                    BSimTransportField field = src.payload == LaneBA0Source.Payload.ATTRACTANT
                            ? att : rep;
                    src.deposit(field, context.getStartTime(), context.getDt());
                }
                if (slab != null) {
                    slab.deposit(att, rep, context.getStartTime(), context.getDt());
                }
            }

            @Override
            public void observe(BSimStepScheduler.Context context) {
                double t = context.getEndTime();
                int stride = (int) Math.round(LaneBIdentity.SAMPLE_DT / LaneBIdentity.DT);
                boolean take = context.isInitialObservation()
                        || context.getCompletedUpdates() % stride == 0;
                if (!take || samp[0] >= nSamp) {
                    return;
                }
                int s = samp[0]++;
                time[s] = t;
                double[][] dens = LaneBMaps.zeros();
                double msd = 0.0;
                for (int n = 0; n < cells.length; n++) {
                    if (cells[n] == null) continue;
                    Vector3d p = cells[n].getPosition();
                    LaneBMaps.addCell(dens, p.x, p.y);
                    double dx = p.x - x0[n];
                    double dy = p.y - y0[n];
                    msd += dx * dx + dy * dy;
                }
                if (nCells > 0) msd /= nCells;
                msdXy[s] = msd;
                if (stripeKappa) {
                    double[] xs = new double[nCells];
                    int nx = 0;
                    for (LaneBBacterium cell : cells) {
                        if (cell != null) xs[nx++] = cell.getPosition().x;
                    }
                    if (nx < xs.length) xs = java.util.Arrays.copyOf(xs, nx);
                    kappaCell[s] = nx > 0
                            ? LaneBMaps.cellStripeContrast(xs, LaneB2Identity.SLAB_UM) : Double.NaN;
                    kappaField[s] = (att != null && rep != null)
                            ? LaneBMaps.fieldStripeContrast(att, rep, LaneB2Identity.SLAB_UM)
                            : Double.NaN;
                    deltaLEll[s] = (att != null && rep != null)
                            ? LaneBMaps.medianAbsDeltaL(att, rep, LaneB2Identity.ELL_RUN_UM)
                            : Double.NaN;
                } else {
                    kappaCell[s] = nCells > 0 ? LaneBMaps.cellContrast(dens) : Double.NaN;
                    kappaField[s] = (att != null && rep != null)
                            ? LaneBMaps.fieldSceneContrast(att, rep) : Double.NaN;
                    deltaLEll[s] = Double.NaN;
                }
                meanAbsL[s] = (att != null && rep != null)
                        ? LaneBMaps.meanAbsLigand(att, rep) : Double.NaN;
                if (density0[0] == null) {
                    density0[0] = LaneBMaps.copy(dens);
                    autoCorr[s] = 1.0;
                } else {
                    autoCorr[s] = LaneBMaps.pearson(density0[0], dens);
                }
            }
        });

        int filled = samp[0];
        double mu = fitMu2d(time, msdXy, filled);
        double dTh = arm == Arm.THERMAL ? mu : Double.NaN;
        if (arm != Arm.MSD_MOTILE) {
            /* mu is still the MSD slope; only MSD_MOTILE uses it as swimming μ. */
        }
        double cmdA = slab != null
                ? slab.commandedAtt()
                : commanded(sources, LaneBA0Source.Payload.ATTRACTANT);
        double cmdR = slab != null
                ? slab.commandedRep()
                : commanded(sources, LaneBA0Source.Payload.REPELLENT);
        double attRes = relativeResidual(att, cmdA, initialAtt);
        double repRes = relativeResidual(rep, cmdR, initialRep);
        double attMass = att == null ? 0.0 : att.totalFluidQuantity();
        double repMass = rep == null ? 0.0 : rep.totalFluidQuantity();
        boolean ledgerOk = !chemistry
                || (Math.abs(attRes) < LaneBIdentity.LEDGER_REL
                && Math.abs(repRes) < LaneBIdentity.LEDGER_REL
                && attMass >= 0.0 && repMass >= 0.0
                && Double.isFinite(attMass) && Double.isFinite(repMass));
        double pulseL = pulseMean(time, meanAbsL, samp[0], pulseOff);
        return new Result(
                arm.name(), arm, duration, time, msdXy, kappaCell, kappaField, autoCorr,
                arm == Arm.MSD_MOTILE ? mu : Double.NaN,
                dTh, attRes, repRes, attMass, repMass, cmdA, cmdR, ledgerOk,
                meanAbsL, pulseL, deltaLEll);
    }

    private static boolean isPaintArm(Arm arm) {
        return arm == Arm.PAINT_MOTILE || arm == Arm.PAINT_OFF
                || arm == Arm.PAINT_FIELD || arm == Arm.PAINT_SILENT;
    }

    private static boolean isStripeArm(Arm arm) {
        return arm == Arm.STRIPE_MOTILE || arm == Arm.STRIPE_OFF
                || arm == Arm.STRIPE_FIELD || arm == Arm.STRIPE_SILENT;
    }

    private static boolean isSustainArm(Arm arm) {
        return arm == Arm.SUSTAIN_MOTILE || arm == Arm.SUSTAIN_OFF
                || arm == Arm.SUSTAIN_FIELD || arm == Arm.SUSTAIN_SILENT;
    }

    private static void paintSquareWave(BSimTransportField att, BSimTransportField rep, double lSlab) {
        for (int i = 0; i < LaneBIdentity.NX; i++) {
            double xc = (i + 0.5) * LaneBIdentity.DX;
            boolean even = LaneB2Identity.evenAttSlab(xc);
            for (int j = 0; j < LaneBIdentity.NY; j++) {
                if (even) att.setConc(i, j, 0, lSlab);
                else rep.setConc(i, j, 0, lSlab);
            }
        }
    }

    private static LaneBBacterium.Motility motilityOf(Arm arm) {
        if (arm == Arm.PAINT_OFF || arm == Arm.STRIPE_OFF || arm == Arm.SUSTAIN_OFF) {
            return LaneBBacterium.Motility.FROZEN;
        }
        if (arm == Arm.THERMAL) return LaneBBacterium.Motility.THERMAL;
        return LaneBBacterium.Motility.MOTILE;
    }

    private static int nCellsOf(Arm arm) {
        if (arm == Arm.PAINT_FIELD || arm == Arm.STRIPE_FIELD || arm == Arm.SUSTAIN_FIELD) {
            return 0;
        }
        if (arm == Arm.MSD_MOTILE || arm == Arm.THERMAL) return LaneBIdentity.N_MSD;
        return LaneBIdentity.N_DENSITY;
    }

    private Vector3d interiorSpawn(Arm arm) {
        double m = LaneBIdentity.SPAWN_MARGIN;
        double x = m + rng.nextDouble() * (LaneBIdentity.LX - 2.0 * m);
        double y = m + rng.nextDouble() * (LaneBIdentity.LY - 2.0 * m);
        if (arm == Arm.PAINT_MOTILE || arm == Arm.PAINT_OFF || arm == Arm.PAINT_SILENT
                || arm == Arm.STRIPE_MOTILE || arm == Arm.STRIPE_OFF || arm == Arm.STRIPE_SILENT
                || arm == Arm.SUSTAIN_MOTILE || arm == Arm.SUSTAIN_OFF || arm == Arm.SUSTAIN_SILENT) {
            x = rng.nextDouble() * LaneBIdentity.LX;
            y = rng.nextDouble() * LaneBIdentity.LY;
        }
        return new Vector3d(x, y, LaneBIdentity.LZ * 0.5);
    }

    private Vector3d blobSpawn() {
        double x;
        double y;
        do {
            x = 0.5 * LaneBIdentity.LX + LaneBIdentity.BLOB_SIGMA * rng.nextGaussian();
            y = 0.5 * LaneBIdentity.LY + LaneBIdentity.BLOB_SIGMA * rng.nextGaussian();
        } while (x < 1.0 || x > LaneBIdentity.LX - 1.0
                || y < 1.0 || y > LaneBIdentity.LY - 1.0);
        return new Vector3d(x, y, LaneBIdentity.LZ * 0.5);
    }

    private Vector3d stripeSpawn() {
        double x;
        double y;
        double high = LaneB1Identity.STRIPE_HIGH_UM;
        do {
            x = rng.nextDouble() * LaneBIdentity.LX;
            y = 1.0 + rng.nextDouble() * (LaneBIdentity.LY - 2.0);
        } while (x < 1.0 || x > LaneBIdentity.LX - 1.0
                || Math.floor(x / high) % 2.0 != 0.0);
        return new Vector3d(x, y, LaneBIdentity.LZ * 0.5);
    }

    static double pulseMean(double[] t, double[] y, int n, double tOff) {
        if (t == null || y == null || n <= 0 || !Double.isFinite(tOff)) {
            return Double.NaN;
        }
        double s = 0.0;
        int c = 0;
        for (int i = 0; i < n; i++) {
            if (t[i] <= 1e-12 || t[i] > tOff + 1e-12) continue;
            if (!Double.isFinite(y[i])) continue;
            s += y[i];
            c++;
        }
        return c == 0 ? Double.NaN : s / c;
    }

    static double fitMu2d(double[] t, double[] msd, int n) {
        double sT = 0.0;
        double sM = 0.0;
        double sTT = 0.0;
        double sTM = 0.0;
        int c = 0;
        for (int i = 0; i < n; i++) {
            if (t[i] < LaneBIdentity.MSD_FIT_T0 - 1e-12
                    || t[i] > LaneBIdentity.MSD_FIT_T1 + 1e-12) {
                continue;
            }
            if (!Double.isFinite(msd[i])) return Double.NaN;
            sT += t[i];
            sM += msd[i];
            sTT += t[i] * t[i];
            sTM += t[i] * msd[i];
            c++;
        }
        if (c < 3) return Double.NaN;
        double den = c * sTT - sT * sT;
        if (Math.abs(den) < 1e-18) return Double.NaN;
        double slope = (c * sTM - sT * sM) / den;
        return slope / 4.0;
    }

    static double tauMix(Result mix) {
        final double target = 1.0 / Math.E;
        for (int i = 1; i < mix.time.length; i++) {
            if (i > 1 && mix.time[i] <= mix.time[i - 1]) break;
            double c0 = mix.autoCorr[i - 1];
            double c1 = mix.autoCorr[i];
            if (!Double.isFinite(c0) || !Double.isFinite(c1)) continue;
            if (c0 >= target && c1 <= target) {
                double f = (c0 - target) / (c0 - c1 + 1e-18);
                return mix.time[i - 1] + f * (mix.time[i] - mix.time[i - 1]);
            }
        }
        return Double.NaN;
    }

    static double interpolate(double[] t, double[] y, double tq) {
        if (t == null || y == null || t.length == 0) return Double.NaN;
        if (tq <= t[0]) return y[0];
        for (int i = 1; i < t.length; i++) {
            if (t[i] < t[i - 1]) break;
            if (tq <= t[i]) {
                double f = (tq - t[i - 1]) / (t[i] - t[i - 1] + 1e-18);
                return y[i - 1] + f * (y[i] - y[i - 1]);
            }
        }
        return y[t.length - 1];
    }

    private static double commanded(List<LaneBA0Source> sources, LaneBA0Source.Payload payload) {
        double m = 0.0;
        for (LaneBA0Source s : sources) {
            if (s.payload == payload) m += s.commandedMass();
        }
        return m;
    }

    private static double relativeResidual(BSimTransportField field, double commanded, double initial) {
        if (field == null) return 0.0;
        BSimTransportLedger led = field.getTransportLedger();
        double remaining = field.totalFluidQuantity();
        double r = led.residual(initial, remaining);
        double scale = Math.max(Math.max(Math.max(commanded, remaining), initial), 1.0);
        return r / scale;
    }

    public static String line(Result r) {
        return String.format(Locale.US,
                "%s T=%.3g mu2d=%.6g Dth=%.6g attR=%.3e repR=%.3e cmdA=%.6g cmdR=%.6g ledger=%s",
                r.name, r.duration, r.mu2d, r.thermalD2d,
                r.attResidualRel, r.repResidualRel,
                r.commandedAtt, r.commandedRep, r.ledgerPass ? "PASS" : "FAIL");
    }
}
