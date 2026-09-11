package bsim.lanea;

import bsim.BSim;
import bsim.BSimRandom;
import bsim.BSimStepScheduler;
import bsim.particle.BSimParticle;
import bsim.transport.BSimTransportField;
import bsim.transport.BSimTransportLedger;

import javax.vecmath.Vector3d;
import java.util.ArrayList;
import java.util.List;

/**
 * N0 millimetre dish for Lane A occupancy. Immobilized Hill R,L.
 * Motility off on living arms. LANE_A_OCCUPIED_MILLIMETRE.
 */
public final class LaneAOccupiedDish {

    public enum Arm {
        CLOSED_CONSERVATIVE,
        DRIVEN,
        SILENT,
        BROWNIAN
    }

    public static final class Sample {
        public final double t;
        public final int window;
        public final int nBodies;
        public final double meanR;
        public final double meanL;
        public final double fracRgt05;
        public final double ahlUmMean;
        public final double ahlUmCenter;
        public final double remaining;
        public final double commanded;
        public final double sourceAdded;
        public final double decayLoss;
        public final double residual;
        public final double uCommand;
        public final double[] ahlReadout;
        public final double[] rReadout;
        public final double[] lReadout;

        Sample(double t, int window, int nBodies, double meanR, double meanL, double fracRgt05,
               double ahlUmMean, double ahlUmCenter, double remaining, double commanded,
               double sourceAdded, double decayLoss, double residual, double uCommand,
               double[] ahlReadout, double[] rReadout, double[] lReadout) {
            this.t = t;
            this.window = window;
            this.nBodies = nBodies;
            this.meanR = meanR;
            this.meanL = meanL;
            this.fracRgt05 = fracRgt05;
            this.ahlUmMean = ahlUmMean;
            this.ahlUmCenter = ahlUmCenter;
            this.remaining = remaining;
            this.commanded = commanded;
            this.sourceAdded = sourceAdded;
            this.decayLoss = decayLoss;
            this.residual = residual;
            this.uCommand = uCommand;
            this.ahlReadout = ahlReadout;
            this.rReadout = rReadout;
            this.lReadout = lReadout;
        }
    }

    public static final class Trajectory {
        public final Arm arm;
        public final List<Sample> samples;
        public final double commanded;
        public final double sourceAdded;
        public final double decayLoss;
        public final double remaining;
        public final double residual;
        public final double meanROccupancy;
        public final double meanLOccupancy;
        public final int occupancySamples;
        public final String occupancyFlag;

        Trajectory(Arm arm, List<Sample> samples, double commanded, double sourceAdded,
                   double decayLoss, double remaining, double residual,
                   double meanROccupancy, double meanLOccupancy, int occupancySamples,
                   String occupancyFlag) {
            this.arm = arm;
            this.samples = samples;
            this.commanded = commanded;
            this.sourceAdded = sourceAdded;
            this.decayLoss = decayLoss;
            this.remaining = remaining;
            this.residual = residual;
            this.meanROccupancy = meanROccupancy;
            this.meanLOccupancy = meanLOccupancy;
            this.occupancySamples = occupancySamples;
            this.occupancyFlag = occupancyFlag;
        }
    }

    private static final class Receiver {
        final Vector3d position;
        double r;
        double l;

        Receiver(Vector3d position) {
            this.position = position;
        }
    }

    public Trajectory run(Arm arm) {
        return run(arm, false, null);
    }

    public Trajectory run(Arm arm, boolean recordMaps) {
        return run(arm, recordMaps, null);
    }

    public Trajectory run(Arm arm, boolean recordMaps, double[] narmaU) {
        return run(arm, recordMaps, narmaU, LaneAConstants.RNG_SEED);
    }

    public Trajectory run(Arm arm, boolean recordMaps, double[] narmaU, long rngSeed) {
        boolean closed = arm == Arm.CLOSED_CONSERVATIVE;
        boolean brownian = arm == Arm.BROWNIAN;
        boolean living = arm == Arm.DRIVEN || arm == Arm.SILENT;
        boolean narma = narmaU != null;
        if (narma && (closed || brownian)) {
            throw new IllegalArgumentException(
                    "per-window extra is driven/silent only " + LaneAConstants.LABEL);
        }
        LaneAIdealSource.Command command;
        double tEnd;
        double decay;
        int occLo = LaneAConstants.OCCUPANCY_WINDOW_LO;
        int occHi = LaneAConstants.OCCUPANCY_WINDOW_HI;
        int nWin = narmaU != null ? narmaU.length : LaneAConstants.NUM_WINDOWS;
        if (closed) {
            command = LaneAIdealSource.Command.U_CLOSED_IMPULSE;
            tEnd = LaneAConstants.T_END_CLOSED;
            decay = 0.0;
        } else if (arm == Arm.SILENT) {
            command = LaneAIdealSource.Command.U_ZERO;
            tEnd = nWin * LaneAConstants.WINDOW_S;
            decay = LaneAConstants.K_AHL;
            if (narmaU != null && narmaU.length == LaneAConstants.KR_GR_WINDOWS) {
                occLo = LaneAConstants.KR_GR_OCC_LO;
                occHi = LaneAConstants.KR_GR_OCC_HI;
            } else if (narmaU != null && narmaU.length == LaneAConstants.WAVEFORM_WINDOWS) {
                occLo = LaneAConstants.WAVEFORM_OCC_LO;
                occHi = LaneAConstants.WAVEFORM_OCC_HI;
            } else if (narma) {
                occLo = LaneAConstants.NARMA10_OCC_LO;
                occHi = LaneAConstants.NARMA10_OCC_HI;
            }
        } else {
            command = narma ? LaneAIdealSource.Command.U_NARMA10 : LaneAIdealSource.Command.U_PULSE_TRAIN;
            tEnd = nWin * LaneAConstants.WINDOW_S;
            decay = LaneAConstants.K_AHL;
            if (narmaU != null && narmaU.length == LaneAConstants.KR_GR_WINDOWS) {
                occLo = LaneAConstants.KR_GR_OCC_LO;
                occHi = LaneAConstants.KR_GR_OCC_HI;
            } else if (narmaU != null && narmaU.length == LaneAConstants.WAVEFORM_WINDOWS) {
                occLo = LaneAConstants.WAVEFORM_OCC_LO;
                occHi = LaneAConstants.WAVEFORM_OCC_HI;
            } else if (narma) {
                occLo = LaneAConstants.NARMA10_OCC_LO;
                occHi = LaneAConstants.NARMA10_OCC_HI;
            }
        }

        LaneAIdealSource ac = narma && arm == Arm.DRIVEN
                ? new LaneAIdealSource(narmaU)
                : new LaneAIdealSource(command);
        BSim sim = new BSim();
        sim.setDt(LaneAConstants.DT);
        sim.setSimulationTime(tEnd);
        sim.setBound(LaneAConstants.BOUND_X, LaneAConstants.BOUND_Y, LaneAConstants.BOUND_Z);
        sim.setSolid(true, true, true);
        sim.setLeaky(false, false, false, false, false, false);
        sim.setRandomSeed(rngSeed);

        boolean[][][] fluid = allFluid();
        BSimTransportField field = new BSimTransportField(
                sim,
                new int[] {LaneAConstants.GRID_X, LaneAConstants.GRID_Y, LaneAConstants.GRID_Z},
                LaneAConstants.D_AHL,
                decay,
                fluid);
        int[] acBox = field.boxCoords(new Vector3d(
                LaneAConstants.AC_X, LaneAConstants.AC_Y, LaneAConstants.AC_Z));
        if (acBox[0] != LaneAConstants.I_AC || acBox[1] != LaneAConstants.J_AC
                || acBox[2] != LaneAConstants.K_AC) {
            throw new IllegalStateException(
                    "CENTER voxel mismatch " + acBox[0] + "," + acBox[1] + "," + acBox[2]
                            + " " + LaneAConstants.LABEL);
        }
        field.getTransportLedger().reset();

        BSimRandom rng = sim.getRandom();
        Receiver[] receivers = living ? placeReceivers(rng) : new Receiver[0];
        BSimParticle[] particles = brownian ? placeParticles(sim, rng) : new BSimParticle[0];

        List<Sample> samples = new ArrayList<Sample>();
        int[] steps = {0};
        int expectedSteps = BSimStepScheduler.exactSteps(tEnd, LaneAConstants.DT, "Lane A duration");

        BSimStepScheduler scheduler = new BSimStepScheduler(sim);
        if (brownian) {
            scheduler.addMechanicsEvent("brownian", LaneAConstants.DT, context -> {
                for (BSimParticle p : particles) {
                    p.action();
                    p.updatePosition();
                }
            });
        }
        scheduler.run(new BSimStepScheduler.Adapter() {
            @Override
            public void transport(BSimStepScheduler.Context context, double duration) {
                field.advance(duration);
            }

            @Override
            public void integrateModels(BSimStepScheduler.Context context) {
                if (!living) {
                    return;
                }
                double dt = context.getDt();
                for (Receiver cell : receivers) {
                    double cUm = field.getConc(cell.position) / LaneAConstants.MOLECULES_PER_UM3_PER_UM;
                    double cN = Math.pow(Math.max(0.0, cUm), LaneAConstants.HILL_N);
                    double kN = Math.pow(LaneAConstants.HILL_K_UM, LaneAConstants.HILL_N);
                    double target = cN / (kN + cN);
                    cell.r = target + (cell.r - target) * Math.exp(-dt / LaneAConstants.TAU_R_S);
                    cell.l = cell.r + (cell.l - cell.r) * Math.exp(-dt / LaneAConstants.TAU_L_S);
                }
            }

            @Override
            public void depositFluxes(BSimStepScheduler.Context context) {
                ac.deposit(field, context.getStartTime(), context.getDt());
                steps[0]++;
            }

            @Override
            public void observe(BSimStepScheduler.Context context) {
                if (context.isInitialObservation()) {
                    return;
                }
                double t = context.getEndTime();
                int completed = context.getCompletedUpdates();
                if (closed) {
                    if (isMultiple(t, 1.0) || t + 1e-12 >= tEnd) {
                        samples.add(record(t, 0, field, ac, receivers, particles, brownian, living, false));
                    }
                    return;
                }
                int stepsInWindow = BSimStepScheduler.exactSteps(
                        LaneAConstants.WINDOW_S, LaneAConstants.DT, "window");
                int sampleEvery = BSimStepScheduler.exactSteps(
                        LaneAConstants.SAMPLE_INTERVAL_S, LaneAConstants.DT, "sample");
                int stepInWindow = (completed - 1) % stepsInWindow;
                int window = (completed - 1) / stepsInWindow;
                boolean sample = stepInWindow % sampleEvery == 0
                        || stepInWindow == stepsInWindow - 1;
                if (sample) {
                    samples.add(record(t, window, field, ac, receivers, particles, brownian, living, recordMaps));
                    if (stepInWindow == stepsInWindow - 1 && (window % 10 == 9 || !narma)) {
                        System.out.printf(
                                "    %s %s window=%d t=%.2f mean_R=%.6g remaining=%.6g%n",
                                arm.name(),
                                narma ? LaneAConstants.NARMA10_LABEL : LaneAConstants.LABEL,
                                window, t,
                                samples.get(samples.size() - 1).meanR,
                                samples.get(samples.size() - 1).remaining);
                    }
                }
            }
        });

        if (steps[0] != expectedSteps) {
            throw new IllegalStateException(
                    "scheduler steps " + steps[0] + " != " + expectedSteps + " " + LaneAConstants.LABEL);
        }

        BSimTransportLedger ledger = field.getTransportLedger();
        double remaining = field.totalFluidQuantity();
        double residual = ledger.residual(0.0, remaining);
        double meanR = occupancyMean(samples, true, occLo, occHi);
        double meanL = occupancyMean(samples, false, occLo, occHi);
        int nOcc = occupancyCount(samples, occLo, occHi);
        String flag;
        if (brownian) {
            flag = "NO_RECEIVER";
        } else if (closed) {
            flag = "FIELD_ONLY";
        } else if (meanR >= LaneAConstants.OCCUPANCY_ALIVE_MIN_R) {
            flag = "ALIVE";
        } else {
            flag = "DEAD";
        }
        return new Trajectory(
                arm, samples, ac.commandedMass(), ledger.getSourceAdded(),
                ledger.getDecayLoss(), remaining, residual,
                meanR, meanL, nOcc, flag);
    }

    private static Receiver[] placeReceivers(BSimRandom rng) {
        Receiver[] cells = new Receiver[LaneAConstants.N_CELLS];
        for (int i = 0; i < cells.length; i++) {
            cells[i] = new Receiver(seedPosition(rng));
        }
        return cells;
    }

    private static BSimParticle[] placeParticles(BSim sim, BSimRandom rng) {
        BSimParticle[] particles = new BSimParticle[LaneAConstants.N_CELLS];
        for (int i = 0; i < particles.length; i++) {
            particles[i] = new BSimParticle(sim, seedPosition(rng), 1.0);
        }
        return particles;
    }

    private static Vector3d seedPosition(BSimRandom rng) {
        return new Vector3d(
                LaneAConstants.SEED_X0 + rng.nextDouble() * LaneAConstants.SEED_XSPAN,
                LaneAConstants.SEED_Y0 + rng.nextDouble() * LaneAConstants.SEED_YSPAN,
                LaneAConstants.SEED_Z);
    }

    private static boolean[][][] allFluid() {
        boolean[][][] fluid = new boolean[LaneAConstants.GRID_X][LaneAConstants.GRID_Y][LaneAConstants.GRID_Z];
        for (int i = 0; i < LaneAConstants.GRID_X; i++) {
            for (int j = 0; j < LaneAConstants.GRID_Y; j++) {
                fluid[i][j][0] = true;
            }
        }
        return fluid;
    }

    private static boolean inOccupancyWindow(int window, int lo, int hi) {
        return window >= lo && window <= hi;
    }

    private static boolean isMultiple(double t, double period) {
        double n = t / period;
        return Math.abs(n - Math.round(n)) * period <= 1e-9;
    }

    private static Sample record(
            double t,
            int window,
            BSimTransportField field,
            LaneAIdealSource ac,
            Receiver[] receivers,
            BSimParticle[] particles,
            boolean brownian,
            boolean living,
            boolean recordMaps) {
        BSimTransportLedger ledger = field.getTransportLedger();
        double remaining = field.totalFluidQuantity();
        double residual = ledger.residual(0.0, remaining);
        double ahlUmCenter = field.getConc(
                LaneAConstants.I_AC, LaneAConstants.J_AC, LaneAConstants.K_AC)
                / LaneAConstants.MOLECULES_PER_UM3_PER_UM;
        double ahlSum = 0.0;
        int nBoxes = LaneAConstants.GRID_X * LaneAConstants.GRID_Y * LaneAConstants.GRID_Z;
        for (int i = 0; i < LaneAConstants.GRID_X; i++) {
            for (int j = 0; j < LaneAConstants.GRID_Y; j++) {
                ahlSum += field.getConc(i, j, 0);
            }
        }
        double ahlUmMean = (ahlSum / nBoxes) / LaneAConstants.MOLECULES_PER_UM3_PER_UM;
        double meanR = 0.0;
        double meanL = 0.0;
        double frac = 0.0;
        int nBodies;
        if (brownian) {
            nBodies = particles.length;
        } else if (living) {
            nBodies = receivers.length;
            for (Receiver cell : receivers) {
                meanR += cell.r;
                meanL += cell.l;
                if (cell.r > 0.5) {
                    frac += 1.0;
                }
            }
            meanR /= nBodies;
            meanL /= nBodies;
            frac /= nBodies;
        } else {
            nBodies = 0;
        }
        double uCommand = ac.commandU(window, t);
        double[] ahlReadout = null;
        double[] rReadout = null;
        double[] lReadout = null;
        if (recordMaps) {
            ahlReadout = new double[LaneAConstants.READOUT_BINS];
            rReadout = new double[LaneAConstants.READOUT_BINS];
            lReadout = new double[LaneAConstants.READOUT_BINS];
            fillReadout(field, receivers, living, ahlReadout, rReadout, lReadout);
        }
        return new Sample(
                t, window, nBodies, meanR, meanL, frac,
                ahlUmMean, ahlUmCenter, remaining, ac.commandedMass(),
                ledger.getSourceAdded(), ledger.getDecayLoss(), residual,
                uCommand, ahlReadout, rReadout, lReadout);
    }

    private static void fillReadout(
            BSimTransportField field,
            Receiver[] receivers,
            boolean living,
            double[] ahlReadout,
            double[] rReadout,
            double[] lReadout) {
        double voxelX = LaneAConstants.BOUND_X / LaneAConstants.READOUT_X;
        double voxelY = LaneAConstants.BOUND_Y / LaneAConstants.READOUT_Y;
        double voxelZ = LaneAConstants.BOUND_Z / LaneAConstants.READOUT_Z;
        for (int x = 0; x < LaneAConstants.READOUT_X; x++) {
            for (int y = 0; y < LaneAConstants.READOUT_Y; y++) {
                int id = LaneAConstants.readoutIndex(x, y, 0);
                Vector3d centre = new Vector3d(
                        (x + 0.5) * voxelX, (y + 0.5) * voxelY, 0.5 * voxelZ);
                ahlReadout[id] = field.getConc(centre) / LaneAConstants.MOLECULES_PER_UM3_PER_UM;
            }
        }
        if (!living) {
            return;
        }
        double[] rSum = new double[LaneAConstants.READOUT_BINS];
        double[] lSum = new double[LaneAConstants.READOUT_BINS];
        int[] den = new int[LaneAConstants.READOUT_BINS];
        for (Receiver cell : receivers) {
            int xi = (int) (cell.position.x / voxelX);
            int yi = (int) (cell.position.y / voxelY);
            if (xi < 0) {
                xi = 0;
            } else if (xi >= LaneAConstants.READOUT_X) {
                xi = LaneAConstants.READOUT_X - 1;
            }
            if (yi < 0) {
                yi = 0;
            } else if (yi >= LaneAConstants.READOUT_Y) {
                yi = LaneAConstants.READOUT_Y - 1;
            }
            int id = LaneAConstants.readoutIndex(xi, yi, 0);
            den[id]++;
            rSum[id] += cell.r;
            lSum[id] += cell.l;
        }
        for (int i = 0; i < LaneAConstants.READOUT_BINS; i++) {
            if (den[i] > 0) {
                rReadout[i] = rSum[i] / den[i];
                lReadout[i] = lSum[i] / den[i];
            }
        }
    }

    private static double occupancyMean(List<Sample> samples, boolean rNotL, int lo, int hi) {
        double s = 0.0;
        int n = 0;
        for (Sample sample : samples) {
            if (inOccupancyWindow(sample.window, lo, hi)) {
                s += rNotL ? sample.meanR : sample.meanL;
                n++;
            }
        }
        return n == 0 ? Double.NaN : s / n;
    }

    private static int occupancyCount(List<Sample> samples, int lo, int hi) {
        int n = 0;
        for (Sample sample : samples) {
            if (inOccupancyWindow(sample.window, lo, hi)) {
                n++;
            }
        }
        return n;
    }
}
