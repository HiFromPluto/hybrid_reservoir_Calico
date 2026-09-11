package PocketNeckT3;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.draw.BSimP3DDrawer;
import bsim.particle.BSimBacterium;

import processing.core.PGraphics3D;

import javax.vecmath.Vector3d;
import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Properties;
import java.util.Random;
import java.util.Vector;

/**
 * PocketNeck-T3 living BSim occupancy. ENGINEERING port of T2 geometry
 * and clocks onto the HybridDish ticker / Hill receiver. Interior plastic
 * and bus AHL advection are wrapped here; src/BSimChemicalField is not
 * forked. No NARMA, no LuxI, no growth-on.
 */
public final class BSimPocketNeckT3 {
    enum Arm { FIELD_W100_FLUSH, FIELD_W20_L20, LIVE_W100_FLUSH, LIVE_W20_L20 }

    static final byte WALL = 0, POCKET = 1, NECK = 2, BUS = 3;

    static final double BOUND_X = 400.0, BOUND_Y = 200.0, BOUND_Z = 10.0;
    static final double DT_FROZEN = 0.02;
    static final int GRID_X = 80, GRID_Y = 40, GRID_Z = 1;
    static final double DX = BOUND_X / GRID_X, DY = BOUND_Y / GRID_Y, DZ = BOUND_Z / GRID_Z;

    static final double POCKET_X0 = 150.0, POCKET_Y0 = 0.0;
    static final double POCKET_LX = 100.0, POCKET_LY = 100.0;
    static final double BUS_X0 = 0.0, BUS_LX = 400.0, BUS_LY = 80.0;
    static final double FLUSH_BUS_Y0 = 100.0, NECKED_BUS_Y0 = 120.0;
    static final double NECK_Y0 = 100.0, NECK_LY = 20.0, NECK_W = 20.0;
    static final double SOURCE_X = 200.0, SOURCE_Y = 50.0, SOURCE_Z = 5.0;
    /** Preview only. GUV-class vesicle, 10 µm diameter. Not a payload or u→J model. */
    static final double AC_VISUAL_RADIUS_UM = 5.0;
    static final Color AC_VISUAL_COLOR = new Color(0, 180, 255);

    static final double AHL_DIFFUSIVITY = 159.0;
    static final double AHL_DECAY = 0.0033;
    static final double J_MAX = 636000.0;
    static final double RECEIVER_K_UM = 1.6;
    static final double RECEIVER_HILL_N = 2.0;
    static final double RECEIVER_TAU_S = 15.0;
    static final double TAU_L = 1500.0;
    static final double ALPHA_LUX = 1.0 / 1500.0;
    static final double DELTA_LUX = 1.0 / 1500.0;
    static final double MOLECULES_PER_UM3_PER_UM = 602.2;
    static final double BUS_FLOW = 3.0;
    static final double U_HOLD = 0.5;
    static final double T0_MEAN_R = 0.0318;
    static final double T0_MEAN_R_REL_TOL = 0.15;

    static double DT = DT_FROZEN;
    static double DURATION = 9000.0;
    static double SAMPLE_DT = 20.0;
    static boolean HEADLESS = true;
    static boolean CELLS = false;
    static boolean GROWTH_ON = false;
    static boolean CLAMP_OFF = true;
    static boolean SMOKE = false;
    static int INITIAL_POP = 0;
    static long RNG_SEED = 101;
    static Arm ARM = Arm.FIELD_W100_FLUSH;
    static String OUTPUT_DIR = "results/runs/FIELD_W100_FLUSH";
    static String RUN_LABEL = "field_w100_flush";
    static double ACID_CELL_PRODUCTION_RATE = 0.0;

    static final Vector<ReservoirBacterium> bacteria = new Vector<>();
    static final Vector<ReservoirBacterium> children = new Vector<>();
    static final Vector<ReservoirBacterium> removals = new Vector<>();
    static Random experimentRng;
    static int nextId;
    static ChipGeometry geom;

    static final class ChipGeometry {
        final boolean necked;
        final double busY0;
        final double neckX0, neckX1;
        final byte[][] region;
        final int[] pocketI, pocketJ, neckI, neckJ, busI, busJ;
        final int nPocket, nNeck, nBus, nActive;
        final int sourceI, sourceJ;
        final List<int[]> interfacePairs = new ArrayList<>();

        ChipGeometry(boolean necked) {
            this.necked = necked;
            this.busY0 = necked ? NECKED_BUS_Y0 : FLUSH_BUS_Y0;
            double mid = POCKET_X0 + 0.5 * POCKET_LX;
            this.neckX0 = mid - 0.5 * NECK_W;
            this.neckX1 = mid + 0.5 * NECK_W;
            region = new byte[GRID_X][GRID_Y];
            int[] pI = new int[GRID_X * GRID_Y], pJ = new int[GRID_X * GRID_Y];
            int[] nI = new int[GRID_X * GRID_Y], nJ = new int[GRID_X * GRID_Y];
            int[] bI = new int[GRID_X * GRID_Y], bJ = new int[GRID_X * GRID_Y];
            int np = 0, nn = 0, nb = 0;
            for (int i = 0; i < GRID_X; i++) {
                double x = (i + 0.5) * DX;
                for (int j = 0; j < GRID_Y; j++) {
                    double y = (j + 0.5) * DY;
                    boolean inP = x >= POCKET_X0 && x < POCKET_X0 + POCKET_LX
                            && y >= POCKET_Y0 && y < POCKET_Y0 + POCKET_LY;
                    boolean inB = x >= BUS_X0 && x < BUS_X0 + BUS_LX
                            && y >= busY0 && y < busY0 + BUS_LY;
                    boolean inN = false;
                    if (necked) {
                        inN = x >= neckX0 && x < neckX1
                                && y >= NECK_Y0 && y < NECK_Y0 + NECK_LY;
                    }
                    if (inP && inN)
                        throw new IllegalStateException("neck stole a pocket voxel at " + i + "," + j);
                    if (inP) {
                        region[i][j] = POCKET;
                        pI[np] = i;
                        pJ[np] = j;
                        np++;
                    } else if (inN) {
                        region[i][j] = NECK;
                        nI[nn] = i;
                        nJ[nn] = j;
                        nn++;
                    } else if (inB) {
                        region[i][j] = BUS;
                        bI[nb] = i;
                        bJ[nb] = j;
                        nb++;
                    } else {
                        region[i][j] = WALL;
                    }
                }
            }
            nPocket = np;
            nNeck = nn;
            nBus = nb;
            nActive = np + nn + nb;
            pocketI = trim(pI, np);
            pocketJ = trim(pJ, np);
            neckI = trim(nI, nn);
            neckJ = trim(nJ, nn);
            busI = trim(bI, nb);
            busJ = trim(bJ, nb);
            sourceI = (int) (SOURCE_X / DX);
            sourceJ = (int) (SOURCE_Y / DY);
            if (region[sourceI][sourceJ] != POCKET)
                throw new IllegalStateException("source is not in the pocket");
            byte downstream = necked ? NECK : BUS;
            for (int k = 0; k < nPocket; k++) {
                int i = pocketI[k], j = pocketJ[k];
                if (j + 1 < GRID_Y && region[i][j + 1] == downstream)
                    interfacePairs.add(new int[]{i, j, i, j + 1});
            }
            int expectedPocket = (int) Math.round(POCKET_LX / DX) * (int) Math.round(POCKET_LY / DY);
            if (nPocket != expectedPocket)
                throw new IllegalStateException("pocket voxels " + nPocket + " != " + expectedPocket);
            if (necked) {
                int expectedNeck = (int) Math.round(NECK_W / DX) * (int) Math.round(NECK_LY / DY);
                if (nNeck != expectedNeck)
                    throw new IllegalStateException("neck voxels " + nNeck + " != " + expectedNeck);
            } else if (nNeck != 0) {
                throw new IllegalStateException("flush arm has neck voxels");
            }
        }

        static int[] trim(int[] a, int n) {
            int[] out = new int[n];
            System.arraycopy(a, 0, out, 0, n);
            return out;
        }

        boolean active(int i, int j) {
            return region[i][j] != WALL;
        }

        boolean inPocket(double x, double y) {
            return x >= POCKET_X0 && x <= POCKET_X0 + POCKET_LX
                    && y >= POCKET_Y0 && y <= POCKET_Y0 + POCKET_LY;
        }

        boolean inNeck(double x, double y) {
            if (!necked) return false;
            return x >= neckX0 && x <= neckX1
                    && y >= NECK_Y0 && y <= NECK_Y0 + NECK_LY;
        }

        boolean inBus(double x, double y) {
            return y >= busY0;
        }

        boolean inGarage(double x, double y) {
            return inPocket(x, y) || inNeck(x, y);
        }
    }

    /**
     * Example-local wrap of BSimChemicalField. Engine diffuse/decay stay
     * in src/; walls and bus advection are applied here.
     */
    static final class ChipField extends BSimChemicalField {
        final ChipGeometry g;
        final double[][][] scratch;
        final double[][] rField, lField;

        ChipField(BSim sim, ChipGeometry g) {
            super(sim, new int[]{GRID_X, GRID_Y, GRID_Z}, AHL_DIFFUSIVITY, AHL_DECAY);
            this.g = g;
            scratch = new double[GRID_X][GRID_Y][GRID_Z];
            rField = new double[GRID_X][GRID_Y];
            lField = new double[GRID_X][GRID_Y];
        }

        void copyToScratch() {
            for (int i = 0; i < GRID_X; i++)
                for (int j = 0; j < GRID_Y; j++)
                    scratch[i][j][0] = quantity[i][j][0];
        }

        /**
         * Conservative FTCS on active–active faces only. Inactive neighbours
         * are no-flux (the T2 stencil), not Dirichlet sinks. z is one box
         * with solid lids, so this is the 2-D ENGINEERING port.
         */
        void maskedDiffuse() {
            copyToScratch();
            double kX = diffusivity * sim.getDt() / (box[0] * box[0]);
            double kY = diffusivity * sim.getDt() / (box[1] * box[1]);
            for (int i = 0; i < GRID_X - 1; i++) {
                for (int j = 0; j < GRID_Y; j++) {
                    if (!(g.active(i, j) && g.active(i + 1, j))) continue;
                    double qx = -kX * (scratch[i + 1][j][0] - scratch[i][j][0]);
                    quantity[i + 1][j][0] += qx;
                    quantity[i][j][0] -= qx;
                }
            }
            for (int i = 0; i < GRID_X; i++) {
                for (int j = 0; j < GRID_Y - 1; j++) {
                    if (!(g.active(i, j) && g.active(i, j + 1))) continue;
                    double qy = -kY * (scratch[i][j + 1][0] - scratch[i][j][0]);
                    quantity[i][j + 1][0] += qy;
                    quantity[i][j][0] -= qy;
                }
            }
        }

        void decayActive() {
            double keep = 1.0 - decayRate * sim.getDt();
            for (int n = 0; n < g.nPocket; n++)
                quantity[g.pocketI[n]][g.pocketJ[n]][0] *= keep;
            for (int n = 0; n < g.nNeck; n++)
                quantity[g.neckI[n]][g.neckJ[n]][0] *= keep;
            for (int n = 0; n < g.nBus; n++)
                quantity[g.busI[n]][g.busJ[n]][0] *= keep;
        }

        double activeMassBeforeDecay() {
            return totalActive();
        }

        void zeroWalls() {
            for (int i = 0; i < GRID_X; i++)
                for (int j = 0; j < GRID_Y; j++)
                    if (g.region[i][j] == WALL) quantity[i][j][0] = 0.0;
        }

        double wallMass() {
            double s = 0;
            for (int i = 0; i < GRID_X; i++)
                for (int j = 0; j < GRID_Y; j++)
                    if (g.region[i][j] == WALL) s += quantity[i][j][0];
            return s;
        }

        /**
         * First-order upwind vx=+3 µm/s in bus voxels only. Inlet C=0.
         * Conservative outlet at the last bus column. Pocket and neck v=0.
         */
        double advectBus() {
            copyToScratch();
            double cour = BUS_FLOW * sim.getDt() / box[0];
            double outlet = 0.0;
            for (int n = 0; n < g.nBus; n++) {
                int i = g.busI[n], j = g.busJ[n];
                double here = scratch[i][j][0];
                double west = 0.0;
                if (i > 0 && g.region[i - 1][j] == BUS) west = scratch[i - 1][j][0];
                quantity[i][j][0] = here + cour * (west - here);
                boolean eastBus = i + 1 < GRID_X && g.region[i + 1][j] == BUS;
                if (!eastBus) outlet += cour * here;
            }
            return outlet;
        }

        double totalActive() {
            double s = 0;
            for (int n = 0; n < g.nPocket; n++) s += quantity[g.pocketI[n]][g.pocketJ[n]][0];
            for (int n = 0; n < g.nNeck; n++) s += quantity[g.neckI[n]][g.neckJ[n]][0];
            for (int n = 0; n < g.nBus; n++) s += quantity[g.busI[n]][g.busJ[n]][0];
            return s;
        }

        double meanConc(int[] ii, int[] jj) {
            if (ii.length == 0) return Double.NaN;
            double s = 0;
            for (int n = 0; n < ii.length; n++) s += getConc(ii[n], jj[n], 0);
            return s / ii.length;
        }

        double meanRPocket() {
            double s = 0;
            for (int n = 0; n < g.nPocket; n++) s += rField[g.pocketI[n]][g.pocketJ[n]];
            return s / g.nPocket;
        }

        double meanLPocket() {
            double s = 0;
            for (int n = 0; n < g.nPocket; n++) s += lField[g.pocketI[n]][g.pocketJ[n]];
            return s / g.nPocket;
        }

        double minActiveConc() {
            double m = Double.POSITIVE_INFINITY;
            for (int n = 0; n < g.nPocket; n++) m = Math.min(m, getConc(g.pocketI[n], g.pocketJ[n], 0));
            for (int n = 0; n < g.nNeck; n++) m = Math.min(m, getConc(g.neckI[n], g.neckJ[n], 0));
            for (int n = 0; n < g.nBus; n++) m = Math.min(m, getConc(g.busI[n], g.busJ[n], 0));
            return m;
        }

        double interfaceFluxMolecules() {
            double rate = 0;
            double area = box[0] * box[2];
            for (int[] p : g.interfacePairs) {
                double cp = getConc(p[0], p[1], 0);
                double cd = getConc(p[2], p[3], 0);
                rate += diffusivity * (cp - cd) / box[1] * area;
            }
            return rate;
        }

        void advanceReporter() {
            double kn = Math.pow(RECEIVER_K_UM, RECEIVER_HILL_N);
            double decayR = Math.exp(-sim.getDt() / RECEIVER_TAU_S);
            for (int n = 0; n < g.nPocket; n++) {
                int i = g.pocketI[n], j = g.pocketJ[n];
                double cUm = getConc(i, j, 0) / MOLECULES_PER_UM3_PER_UM;
                if (cUm < 0) cUm = 0;
                double cn = Math.pow(cUm, RECEIVER_HILL_N);
                double target = cn / (kn + cn);
                rField[i][j] = target + (rField[i][j] - target) * decayR;
                lField[i][j] += (ALPHA_LUX * rField[i][j] - DELTA_LUX * lField[i][j]) * sim.getDt();
            }
        }

        double corner(int i, int j) {
            return getConc(i, j, 0) / MOLECULES_PER_UM3_PER_UM;
        }
    }

    public static final class ReservoirBacterium extends BSimBacterium {
        final int id;
        final ChipField ahlField;
        double receiver;
        double luminescence;
        boolean spillover;

        ReservoirBacterium(BSim sim, Vector3d position, ChipField ahlField) {
            super(sim, position);
            id = nextId++;
            this.ahlField = ahlField;
            this.receiver = 0.0;
            this.luminescence = 0.0;
            setSurfaceAreaGrowthRate(0.0);
        }

        double getResponse() { return receiver; }
        double getLuminescence() { return luminescence; }

        @Override
        public void action() {
            super.action();
            if (GROWTH_ON)
                throw new IllegalStateException("T3 growth must stay off");
            if (ACID_CELL_PRODUCTION_RATE != 0.0)
                throw new IllegalStateException("T3 acid cell production must stay off");
            if (!CLAMP_OFF)
                throw new IllegalStateException("T3 clamp death must stay off");
            double concentrationUm = ahlField.getConc(position) / MOLECULES_PER_UM3_PER_UM;
            if (concentrationUm < 0) concentrationUm = 0;
            double concentrationN = Math.pow(concentrationUm, RECEIVER_HILL_N);
            double target = concentrationN
                    / (Math.pow(RECEIVER_K_UM, RECEIVER_HILL_N) + concentrationN);
            receiver = target + (receiver - target) * Math.exp(-sim.getDt() / RECEIVER_TAU_S);
            luminescence += (ALPHA_LUX * receiver - DELTA_LUX * luminescence) * sim.getDt();
        }

        @Override
        public void updatePosition() {
            super.updatePosition();
            position.z = Math.max(0, Math.min(sim.getBound().z, position.z));
            if (geom.inBus(position.x, position.y)) {
                spillover = true;
                removals.add(this);
                return;
            }
            if (geom.inGarage(position.x, position.y)) return;
            bounceOffPlastic();
            if (geom.inBus(position.x, position.y)) {
                spillover = true;
                removals.add(this);
            }
        }

        /**
         * Mirror on pocket solid walls and neck walls. The pocket +y face is
         * solid except the door (flush: entire 100 µm; necked: W=20 mouth).
         * Do not pull roof-hits into the neck. Do not HybridDish-mirror the
         * outer 400 µm box.
         */
        void bounceOffPlastic() {
            double px0 = POCKET_X0, px1 = POCKET_X0 + POCKET_LX;
            double py0 = POCKET_Y0, py1 = POCKET_Y0 + POCKET_LY;
            double nx0 = geom.neckX0, nx1 = geom.neckX1;
            double ny0 = NECK_Y0, ny1 = NECK_Y0 + NECK_LY;
            boolean inNeckX = geom.necked && position.x >= nx0 && position.x <= nx1;

            if (position.y < py0) {
                position.y = 2 * py0 - position.y;
                position.x = clamp(position.x, px0, px1);
                position.y = clamp(position.y, py0, py1);
                return;
            }
            if (position.y <= py1) {
                if (position.x < px0) position.x = 2 * px0 - position.x;
                if (position.x > px1) position.x = 2 * px1 - position.x;
                position.x = clamp(position.x, px0, px1);
                position.y = clamp(position.y, py0, py1);
                return;
            }
            if (position.y < geom.busY0) {
                if (!inNeckX) {
                    position.y = 2 * py1 - position.y;
                    position.x = clamp(position.x, px0, px1);
                    position.y = clamp(position.y, py0, py1);
                    return;
                }
                if (position.x < nx0) position.x = 2 * nx0 - position.x;
                if (position.x > nx1) position.x = 2 * nx1 - position.x;
                if (position.y < ny0) position.y = 2 * ny0 - position.y;
                if (position.y > ny1) position.y = 2 * ny1 - position.y;
                position.x = clamp(position.x, nx0, nx1);
                position.y = clamp(position.y, ny0, ny1);
            }
        }

        static double clamp(double v, double lo, double hi) {
            return Math.max(lo, Math.min(hi, v));
        }

        @Override
        public void replicate() {
            throw new IllegalStateException("T3 growth is off; replicate() must not fire");
        }
    }

    public static void main(String[] args) {
        if (args.length > 0 && "--theory".equals(args[0])) {
            printTheory();
            return;
        }
        String configPath = args.length > 0 ? args[0] : "config/sim_config.properties";
        loadConfig(configPath);
        if (args.length > 1 && "preview".equals(args[1])) HEADLESS = false;
        if (args.length > 1 && "smoke".equals(args[1])) {
            SMOKE = true;
            DURATION = 50.0;
        }
        if (Math.abs(DT - DT_FROZEN) > 1e-12)
            throw new IllegalArgumentException("T3 freezes dt=0.02; got " + DT);
        if (GROWTH_ON)
            throw new IllegalArgumentException("T3 growth.on must be false");
        if (!CLAMP_OFF)
            throw new IllegalArgumentException("T3 clamp.off must be true");
        if (ACID_CELL_PRODUCTION_RATE != 0.0)
            throw new IllegalArgumentException("T3 acid cell production must be 0");
        if (SAMPLE_DT > 60.0 + 1e-12)
            throw new IllegalArgumentException("sample.dt must be <= 60 s");
        resetStaticState();
        experimentRng = new Random(RNG_SEED);
        seedParticleRng(RNG_SEED);

        boolean necked = ARM == Arm.FIELD_W20_L20 || ARM == Arm.LIVE_W20_L20;
        geom = new ChipGeometry(necked);
        printVoxelCheck(geom);

        double kFtcs = AHL_DIFFUSIVITY * DT / (DX * DX);
        double sixK = 6.0 * kFtcs;
        System.err.printf(Locale.US,
                "ENGINEERING dt=%.4f s  D dt/dx^2=%.4f  6D dt/dx^2=%.4f (need ~<1)%n",
                DT, kFtcs, sixK);
        if (sixK > 1.05)
            throw new IllegalStateException("FTCS stability gate failed: 6 D dt / dx^2 = " + sixK);

        final BSim sim = new BSim();
        sim.setDt(DT);
        sim.setSimulationTime(DURATION);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);

        final ChipField ahlField = new ChipField(sim, geom);
        if (CELLS) {
            for (int i = 0; i < INITIAL_POP; i++) {
                double x = POCKET_X0 + 3.0 + experimentRng.nextDouble() * (POCKET_LX - 6.0);
                double y = POCKET_Y0 + 3.0 + experimentRng.nextDouble() * (POCKET_LY - 6.0);
                ReservoirBacterium b = new ReservoirBacterium(
                        sim, new Vector3d(x, y, BOUND_Z / 2.0), ahlField);
                b.setRadius(1.0);
                b.setSurfaceAreaGrowthRate(0.0);
                b.setChildList(children);
                bacteria.add(b);
            }
        }

        final int nStart = bacteria.size();
        final int[] spillover = {0};
        final double[] injected = {0}, decayLoss = {0}, outletLoss = {0};
        final double[] wallSink = {0}, interfaceLoss = {0};
        final boolean[] cNonneg = {true};
        final List<Double> tSamples = new ArrayList<>();
        final List<Double> rSamples = new ArrayList<>();
        final int sampleEvery = Math.max(1, (int) Math.round(SAMPLE_DT / DT));
        final File outDir = new File(OUTPUT_DIR);
        if (!outDir.exists() && !outDir.mkdirs())
            throw new RuntimeException("Cannot create " + outDir);
        final BufferedWriter ts;
        try {
            ts = new BufferedWriter(new FileWriter(new File(outDir, "timeseries.csv")));
            ts.write("t_s,pocket_mean_AHL_uM,pocket_mean_R,pocket_mean_L,neck_mean_AHL_uM,"
                    + "bus_mean_AHL_uM,N,spillover_cum,cell_mean_R,cell_mean_L\n");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                ahlField.addQuantity(new Vector3d(SOURCE_X, SOURCE_Y, SOURCE_Z),
                        J_MAX * U_HOLD * DT);
                injected[0] += J_MAX * U_HOLD * DT;

                ahlField.maskedDiffuse();
                decayLoss[0] += AHL_DECAY * DT * ahlField.totalActive();
                ahlField.decayActive();
                double walls = ahlField.wallMass();
                wallSink[0] += walls;
                ahlField.zeroWalls();
                outletLoss[0] += ahlField.advectBus();
                interfaceLoss[0] += ahlField.interfaceFluxMolecules() * DT;
                ahlField.advanceReporter();
                if (ahlField.minActiveConc() < -1e-12) cNonneg[0] = false;

                if (CELLS) {
                    for (ReservoirBacterium b : bacteria) {
                        b.action();
                        b.updatePosition();
                    }
                    if (!children.isEmpty())
                        throw new IllegalStateException("T3 must not produce children");
                    for (ReservoirBacterium b : removals) {
                        if (b.spillover) spillover[0]++;
                    }
                    bacteria.removeAll(removals);
                    removals.clear();
                }

                int step = (int) sim.getTimestep();
                boolean last = step >= sim.timesteps(DURATION);
                if (step % sampleEvery == 0 || last) {
                    double t = sim.getTime();
                    double pC = ahlField.meanConc(geom.pocketI, geom.pocketJ) / MOLECULES_PER_UM3_PER_UM;
                    double pR = ahlField.meanRPocket();
                    double pL = ahlField.meanLPocket();
                    double nC = geom.nNeck == 0 ? Double.NaN
                            : ahlField.meanConc(geom.neckI, geom.neckJ) / MOLECULES_PER_UM3_PER_UM;
                    double bC = ahlField.meanConc(geom.busI, geom.busJ) / MOLECULES_PER_UM3_PER_UM;
                    double cellR = Double.NaN, cellL = Double.NaN;
                    if (CELLS && !bacteria.isEmpty()) {
                        double sr = 0, sl = 0;
                        for (ReservoirBacterium b : bacteria) {
                            sr += b.getResponse();
                            sl += b.getLuminescence();
                        }
                        cellR = sr / bacteria.size();
                        cellL = sl / bacteria.size();
                    }
                    tSamples.add(t);
                    rSamples.add(pR);
                    try {
                        ts.write(String.format(Locale.US,
                                "%.4f,%.9g,%.9g,%.9g,%.9g,%.9g,%d,%d,%.9g,%.9g%n",
                                t, pC, pR, pL, nC, bC, bacteria.size(), spillover[0],
                                cellR, cellL));
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                    if (step > 0 && step % (int) Math.round(500.0 / DT) == 0)
                        System.err.printf(Locale.US, "  %s t=%.0f s mean_R=%.5f N=%d%n",
                                ARM.name(), t, pR, bacteria.size());
                }
            }
        });

        if (!HEADLESS) {
            sim.setDrawer(new BSimP3DDrawer(sim, 800, 600) {
                @Override
                public void scene(PGraphics3D p3d) {
                    p3d.ortho(0, (float) BOUND_X, (float) BOUND_Y, 0, -1000, 10000);
                    p3d.camera((float) BOUND_X / 2f, (float) BOUND_Y / 2f, (float) BOUND_Y,
                            (float) BOUND_X / 2f, (float) BOUND_Y / 2f, 0f,
                            0f, 1f, 0f);
                    p3d.perspective((float) Math.PI / 2f,
                            (float) BOUND_X / (float) BOUND_Y, 0.1f, 10000f);
                    draw(ahlField, Color.ORANGE, (float) (255.0 / 400.0));
                    drawGarageOutline(p3d);
                    // Immobilized AC: visual only. Cyan, 10 µm diameter (GUV-class).
                    // Does not move, spill, grow, or change J. Bacteria spheres are
                    // drawn at 8 µm for visibility; true BSim radius is 1 µm.
                    sphere(new Vector3d(SOURCE_X, SOURCE_Y, SOURCE_Z),
                            AC_VISUAL_RADIUS_UM, AC_VISUAL_COLOR, 255);
                    for (ReservoirBacterium b : bacteria)
                        sphere(b.getPosition(), 8, Color.GREEN, 255);
                }

                void drawGarageOutline(PGraphics3D p3d) {
                    p3d.noFill();
                    p3d.stroke(80, 220, 255);
                    wire(p3d, POCKET_X0 + POCKET_LX / 2, POCKET_Y0 + POCKET_LY / 2, BOUND_Z / 2,
                            POCKET_LX, POCKET_LY, BOUND_Z);
                    if (geom.necked) {
                        p3d.stroke(255, 220, 80);
                        wire(p3d, (geom.neckX0 + geom.neckX1) / 2, NECK_Y0 + NECK_LY / 2, BOUND_Z / 2,
                                NECK_W, NECK_LY, BOUND_Z);
                    }
                    p3d.stroke(80, 255, 140);
                    wire(p3d, BUS_X0 + BUS_LX / 2, geom.busY0 + BUS_LY / 2, BOUND_Z / 2,
                            BUS_LX, BUS_LY, BOUND_Z);
                    p3d.noStroke();
                }

                void wire(PGraphics3D p3d, double cx, double cy, double cz,
                          double sx, double sy, double sz) {
                    p3d.pushMatrix();
                    p3d.translate((float) cx, (float) cy, (float) cz);
                    p3d.box((float) sx, (float) sy, (float) sz);
                    p3d.popMatrix();
                }
            });
        }

        PrintStream originalOut = System.out;
        try {
            if (HEADLESS) {
                System.setOut(new PrintStream(OutputStream.nullOutputStream()));
                sim.export();
            } else {
                System.err.println("Preview: close the window. This is not occupancy evidence.");
                sim.preview();
                return;
            }
        } finally {
            System.setOut(originalOut);
            try {
                ts.close();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

        writeSummary(ahlField, nStart, spillover[0], injected[0], decayLoss[0],
                outletLoss[0], wallSink[0], interfaceLoss[0], cNonneg[0],
                tSamples, rSamples);
        System.out.printf(Locale.US, "PocketNeckT3 done arm=%s dir=%s%n", ARM.name(), OUTPUT_DIR);
    }

    static void writeSummary(ChipField field, int nStart, int spillover,
                             double injected, double decay, double outlet, double wallSink,
                             double iface, boolean cNonneg,
                             List<Double> tSamples, List<Double> rSamples) {
        double remaining = field.totalActive();
        double residual = injected - remaining - decay - outlet - wallSink;
        double dominant = Math.max(Math.max(Math.abs(injected), Math.abs(decay)),
                Math.max(Math.abs(outlet), Math.abs(iface)));
        if (dominant < 1) dominant = 1;
        double vsDom = Math.abs(residual) / dominant;
        double pocketC = field.meanConc(geom.pocketI, geom.pocketJ) / MOLECULES_PER_UM3_PER_UM;
        double pocketR = field.meanRPocket();
        double pocketL = field.meanLPocket();
        double neckC = geom.nNeck == 0 ? Double.NaN
                : field.meanConc(geom.neckI, geom.neckJ) / MOLECULES_PER_UM3_PER_UM;
        double busC = field.meanConc(geom.busI, geom.busJ) / MOLECULES_PER_UM3_PER_UM;
        String occ = occupancyFlag(pocketR);
        boolean plateau = plateauOk(tSamples, rSamples, 0.10, 0.05);
        boolean t0Class = Math.abs(pocketR - T0_MEAN_R) / T0_MEAN_R <= T0_MEAN_R_REL_TOL;
        double cellR = Double.NaN, cellL = Double.NaN;
        if (CELLS && !bacteria.isEmpty()) {
            double sr = 0, sl = 0;
            for (ReservoirBacterium b : bacteria) {
                sr += b.getResponse();
                sl += b.getLuminescence();
            }
            cellR = sr / bacteria.size();
            cellL = sl / bacteria.size();
        }
        int i0 = minInt(geom.pocketI), i1 = maxInt(geom.pocketI);
        int j0 = minInt(geom.pocketJ), j1 = maxInt(geom.pocketJ);
        double cSW = field.corner(i0, j0), cSE = field.corner(i1, j0);
        double cNW = field.corner(i0, j1), cNE = field.corner(i1, j1);
        double interior = 0;
        int nInt = 0;
        for (int n = 0; n < geom.nPocket; n++) {
            int i = geom.pocketI[n], j = geom.pocketJ[n];
            if (i >= i0 + 2 && i <= i1 - 2 && j >= j0 + 2 && j <= j1 - 2) {
                interior += field.corner(i, j);
                nInt++;
            }
        }
        interior = nInt > 0 ? interior / nInt : Double.NaN;
        double closedEnd = 0.5 * (cSW + cSE);
        double ratio = interior > 0 ? closedEnd / interior : Double.NaN;
        double l2d = POCKET_LX * POCKET_LX / AHL_DIFFUSIVITY;
        double da = AHL_DECAY * l2d;
        double pe = BUS_FLOW * BUS_LY / AHL_DIFFUSIVITY;
        boolean massOk = vsDom <= 0.01;
        boolean wallsOk = Math.abs(wallSink) <= 1e-6 * Math.max(injected, 1.0);
        String transport = (massOk && cNonneg && wallsOk)
                ? "VALIDATED_FOR_SCREENING" : "MASS_BUDGET_FAIL";
        File dir = new File(OUTPUT_DIR);
        try (BufferedWriter w = new BufferedWriter(new FileWriter(new File(dir, "summary.csv")))) {
            w.write("Arm,Drive,Model,Label,Smoke,dx_um,dt_s,W_um,L_n_um,t_end_s,u_on,"
                    + "PocketMeanAHL_uM,PocketMeanR,PocketMeanL,Occupancy,PlateauR,"
                    + "NeckMeanAHL_uM,BusMeanAHL_uM,InjectedMass_molecules,RemainingMass_molecules,"
                    + "DecayLoss_molecules,BoundaryLoss_molecules,InterfaceLeak_molecules,"
                    + "WallSink_molecules,MassResidual_molecules,MassResidualVsDominant,"
                    + "DominantFlux_molecules,C_nonnegative,WallsZero,N_start,N_end,Spillover,"
                    + "CellMeanR,CellMeanL,CellsWriteAHL,FlushT0Class,C_SW,C_SE,C_NW,C_NE,"
                    + "C_interior,CornerInteriorRatio_closed_end,L2_over_D_s,Damkohler,Pe,"
                    + "TRANSPORT_MODEL_STATUS,n_pocket,n_neck,n_bus\n");
            w.write(String.format(Locale.US,
                    "%s,STEP_ON,bsim_java_chip,OCCUPANCY,%s,%.1f,%.4f,%.1f,%.1f,%.1f,%.2f,"
                            + "%.9g,%.9g,%.9g,%s,%s,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,"
                            + "%.9g,%s,%s,%d,%d,%d,%.9g,%.9g,%s,%s,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,"
                            + "%.9g,%.9g,%.9g,%s,%d,%d,%d%n",
                    ARM.name(), Boolean.toString(SMOKE), DX, DT,
                    ARM.name().contains("W20") ? 20.0 : 100.0,
                    ARM.name().contains("W20") ? 20.0 : 0.0,
                    DURATION, U_HOLD, pocketC, pocketR, pocketL, occ,
                    Boolean.toString(plateau), neckC, busC, injected, remaining, decay,
                    outlet, iface, wallSink, residual, vsDom, dominant,
                    Boolean.toString(cNonneg), Boolean.toString(wallsOk),
                    nStart, bacteria.size(), spillover, cellR, cellL,
                    "false", Boolean.toString(ARM.name().contains("W100") ? t0Class : false),
                    cSW, cSE, cNW, cNE, interior, ratio, l2d, da, pe, transport,
                    geom.nPocket, geom.nNeck, geom.nBus));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        System.err.printf(Locale.US,
                "%s occupancy=%s mean_R=%.6f C=%.4f uM residual/dom=%.3g N=%d->%d spill=%d%n",
                ARM.name(), occ, pocketR, pocketC, vsDom, nStart, bacteria.size(), spillover);
    }

    static String occupancyFlag(double meanR) {
        if (!Double.isFinite(meanR)) return "NA";
        if (meanR > 0.95) return "SATURATED";
        if (meanR >= 0.05) return "ALIVE";
        return "DEAD";
    }

    static boolean plateauOk(List<Double> t, List<Double> r, double lastFrac, double relTol) {
        if (t.isEmpty()) return false;
        double t1 = t.get(t.size() - 1);
        double t0 = t1 * (1.0 - lastFrac);
        double sum = 0;
        int n = 0;
        for (int i = 0; i < t.size(); i++) {
            if (t.get(i) >= t0) {
                sum += r.get(i);
                n++;
            }
        }
        if (n == 0) return false;
        double terminal = r.get(r.size() - 1);
        double ref = sum / n;
        if (Math.abs(ref) < 1e-15) return Math.abs(terminal) < relTol;
        return Math.abs(terminal - ref) / Math.abs(ref) <= relTol;
    }

    static int minInt(int[] a) {
        int m = a[0];
        for (int v : a) if (v < m) m = v;
        return m;
    }

    static int maxInt(int[] a) {
        int m = a[0];
        for (int v : a) if (v > m) m = v;
        return m;
    }

    static void printVoxelCheck(ChipGeometry g) {
        System.err.printf(Locale.US,
                "VOXEL %s: W/dx=%s Ln/dx=%s n_pocket=%d n_neck=%d n_bus=%d n_active=%d n_iface=%d%n",
                ARM.name(),
                g.necked ? String.format(Locale.US, "%.4g", NECK_W / DX) : "NA",
                g.necked ? String.format(Locale.US, "%.4g", NECK_LY / DY) : "NA",
                g.nPocket, g.nNeck, g.nBus, g.nActive, g.interfacePairs.size());
    }

    static void printTheory() {
        double v = POCKET_LX * POCKET_LY * BOUND_Z;
        double c1 = J_MAX / (AHL_DECAY * v * MOLECULES_PER_UM3_PER_UM);
        double c05 = 0.5 * c1;
        double h = (c05 * c05) / (RECEIVER_K_UM * RECEIVER_K_UM + c05 * c05);
        double kFtcs = AHL_DIFFUSIVITY * DT_FROZEN / (DX * DX);
        System.out.println("PocketNeck-T3 THEORY (well-mixed closed algebra + FTCS gate)");
        System.out.printf(Locale.US, "  J_max=%.3e molecules/s  V=%.0f um^3%n", J_MAX, v);
        System.out.printf(Locale.US,
                "  C_ss(u=1)=%.6f uM (target 3.2)  C_ss(u=0.5)=%.6f uM  H=%.4f%n",
                c1, c05, h);
        System.out.printf(Locale.US,
                "  dt=%.4f s  D dt/dx^2=%.4f  6D dt/dx^2=%.4f%n",
                DT_FROZEN, kFtcs, 6.0 * kFtcs);
        System.out.println("  Hypothesis: FIELD_W100_FLUSH DEAD T0-class mean_R~0.0318;");
        System.out.println("  FIELD_W20_L20 ALIVE T2-class mean_R~0.125; living follows field.");
        System.out.println("  Do not retune K or J_max. Occupancy voxels = pocket only.");
    }

    static void seedParticleRng(long seed) {
        try {
            java.lang.reflect.Field f = bsim.particle.BSimParticle.class.getDeclaredField("rng");
            f.setAccessible(true);
            Random rng = (Random) f.get(null);
            rng.setSeed(seed);
        } catch (ReflectiveOperationException e) {
            throw new RuntimeException(e);
        }
    }

    static void loadConfig(String path) {
        Properties p = loadProperties(new File(path));
        DT = getDouble(p, "dt", DT);
        DURATION = getDouble(p, "duration.s", DURATION);
        SAMPLE_DT = getDouble(p, "sample.dt.s", SAMPLE_DT);
        HEADLESS = Boolean.parseBoolean(p.getProperty("headless", String.valueOf(HEADLESS)));
        RNG_SEED = Long.parseLong(p.getProperty("rng.seed", String.valueOf(RNG_SEED)));
        INITIAL_POP = getInt(p, "initial.pop", INITIAL_POP);
        OUTPUT_DIR = p.getProperty("output.dir", OUTPUT_DIR);
        RUN_LABEL = p.getProperty("output.run.label", RUN_LABEL);
        GROWTH_ON = Boolean.parseBoolean(p.getProperty("growth.on", "false"));
        CLAMP_OFF = Boolean.parseBoolean(p.getProperty("clamp.off", "true"));
        ACID_CELL_PRODUCTION_RATE = getDouble(p, "acid.cell.production.rate", 0.0);
        String armName = p.getProperty("arm", ARM.name()).trim().toUpperCase(Locale.ROOT);
        ARM = Arm.valueOf(armName);
        CELLS = Boolean.parseBoolean(p.getProperty("cells",
                Boolean.toString(ARM == Arm.LIVE_W100_FLUSH || ARM == Arm.LIVE_W20_L20)));
        if (CELLS && INITIAL_POP <= 0) INITIAL_POP = 50;
        if (!CELLS) INITIAL_POP = 0;
    }

    static Properties loadProperties(File file) {
        Properties local = new Properties();
        try (FileInputStream in = new FileInputStream(file)) {
            local.load(in);
        } catch (IOException e) {
            throw new IllegalArgumentException("Cannot load config " + file, e);
        }
        String include = local.getProperty("config.include");
        if (include == null) return local;
        File parent = file.getAbsoluteFile().getParentFile();
        Properties merged = loadProperties(new File(parent, include));
        merged.putAll(local);
        return merged;
    }

    static int getInt(Properties p, String k, int d) {
        return Integer.parseInt(p.getProperty(k, String.valueOf(d)));
    }

    static double getDouble(Properties p, String k, double d) {
        return Double.parseDouble(p.getProperty(k, String.valueOf(d)));
    }

    static void resetStaticState() {
        bacteria.clear();
        children.clear();
        removals.clear();
        nextId = 0;
        geom = null;
    }
}
