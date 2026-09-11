package BSimDaninoD1;

import java.awt.Color;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Locale;
import java.util.Random;
import java.util.Vector;

import javax.vecmath.Vector3d;

import processing.core.PGraphics3D;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.ode.BSimOdeSolver;
import bsim.ode.BSimOdeSystem;
import bsim.particle.BSimBacterium;

/**
 * Isolated Danino D1 example: 4-ODE plus a correct µM ↔ molecules conversion,
 * open-loop well-mixed AHL bath. Not a reservoir window experiment.
 *
 * Stage 10 copies this ODE but couples y[1] (µM) to the field with 1e15 / 1e-15.
 * D1 replaces only that factor with 602 molecules/µm³ per µM. QS_* rates,
 * QS_KMLA, Hill n, and CELL_WALL_DIFF are copied from Stage 10 and not retuned.
 *
 * Plan Stage 3's FAIL of this ODE as a 300 s reservoir receiver stays archived.
 */
public class BSimDaninoD1 {

    static final double BOUND_X = 200.0;
    static final double BOUND_Y = 200.0;
    static final double BOUND_Z = 10.0;

    static final int GRID_X = 20;
    static final int GRID_Y = 20;
    static final int GRID_Z = 1;

    static final double AHL_DIFFUSIVITY = 159.0;
    static final double AHL_DECAY_RATE = 2.76e-3 / 60;

    static final double DT = 0.05;
    static final double SIM_TIME = 5400.0;
    static final int N_CELLS = 50;
    static final double GROWTH_RATE = 0.0;
    static final double FLOW_SPEED = 0.0;
    static final long RNG_SEED = 101L;

    static final double T_BATH_ON = 600.0;
    static final double T_BATH_OFF = 2400.0;
    static final double BATH_UM = 0.05;

    /** 1 µM = 602 molecules/µm³. Not Stage 10's 1e15. */
    static final double MOL_PER_UM3_PER_UM = 602.0;

    // QS kinetic parameters (from BSimReservoirStage10 / Stage 8, fixed).
    // All rates converted from min^-1 to s^-1 by dividing by 60.
    // Source: Danino et al. 2010 parameterization via BSimEntrainment_PIDCtrl.
    static final double TIME_ADJ = 60.0;

    static final double QS_DELTA1 = 0.8487 / TIME_ADJ;
    static final double QS_DELTA2 = 0.0234 / TIME_ADJ;
    static final double QS_G = 0.0412;
    static final double QS_KP2 = 9.0 / TIME_ADJ;
    static final double QS_KR1OFF = 6e-6 / TIME_ADJ;
    static final double QS_KR1ON = 5.99e-5 / TIME_ADJ;
    static final double QS_KCAT_AIIA = 2631.4 / TIME_ADJ;
    static final double QS_T_A = 0.00276 / TIME_ADJ;
    static final double QS_T_LA = 0.024 / TIME_ADJ;
    static final double QS_A0LI = 7.785e-6 / TIME_ADJ;
    static final double QS_A0AA = 6.183e-6 / TIME_ADJ;
    static final double QS_KPLI = 0.9 / TIME_ADJ;
    static final double QS_KPAA = 0.9 / TIME_ADJ;
    static final double QS_KMLA = 1e-2;
    static final double QS_KMAA = 1200.0;
    static final double QS_LTOT = 15.0;
    static final double QS_N = 2.0;

    // Cell-membrane AHL exchange rate (Kaplan & Greenberg 1985)
    // "conc. of AHL inside a cell and outside a cell [equilibrated] by 20 sec"
    static final double CELL_WALL_DIFF = 3.0 / TIME_ADJ;

    static final String EXPORT_DIR = "results/d1_seed101/";
    static final String CSV_NAME = "d1_timeseries.csv";

    static final Vector<DaninoBacterium> bacteria = new Vector<DaninoBacterium>();

    /**
     * Bacterium with the Stage 10 Danino 4-ODE. Open-loop: reads the bath,
     * does not write the field. Growth 0, no chemotaxis, no death.
     */
    static class DaninoBacterium extends BSimBacterium {
        double[] y;
        final BSimOdeSystem grn;
        final BSimChemicalField ahlField;
        double ahlExtUm;

        static {
            rng.setSeed(RNG_SEED);
        }

        public DaninoBacterium(BSim sim, Vector3d position, BSimChemicalField ahlField) {
            super(sim, position);
            this.ahlField = ahlField;
            this.y = new double[]{0.05, 0.05, 0.05, 0.05};
            this.grn = new QSGRN();
        }

        @Override
        public void action() {
            super.action();
            ahlExtUm = ahlField.getConc(position) / MOL_PER_UM3_PER_UM;
            y = BSimOdeSolver.rungeKutta45(grn, sim.getTime(), y, sim.getDt());
            for (int i = 0; i < y.length; i++) {
                if (y[i] < 0) y[i] = 0;
            }
            // Open loop: do not addQuantity from the membrane term.
        }

        @Override
        public void replicate() {
            // Growth rate is 0; replication is out of protocol.
        }

        class QSGRN implements BSimOdeSystem {
            @Override
            public double[] derivativeSystem(double t, double[] y) {
                double extraintradiff_uM = y[1] - ahlExtUm;
                double[] dy = new double[4];

                dy[0] = QS_A0LI
                        + QS_KPLI * (Math.pow(y[3], QS_N) / (Math.pow(QS_KMLA, QS_N) + Math.pow(y[3], QS_N)))
                        - (QS_DELTA1 * y[0]) / (QS_G * (y[0] + y[2]) + 1);

                dy[1] = QS_KP2 * y[0]
                        - QS_KR1ON * (QS_LTOT - y[3]) * y[1]
                        + QS_KR1OFF * y[3]
                        - (QS_KCAT_AIIA * y[2] * y[1]) / (QS_KMAA + y[1])
                        - QS_T_A * y[1]
                        - CELL_WALL_DIFF * extraintradiff_uM;

                dy[2] = QS_A0AA
                        + QS_KPAA * (Math.pow(y[3], QS_N) / (Math.pow(QS_KMLA, QS_N) + Math.pow(y[3], QS_N)))
                        - (QS_DELTA2 * y[2]) / (QS_G * (y[0] + y[2]) + 1);

                dy[3] = QS_KR1ON * (QS_LTOT - y[3]) * y[1]
                        - QS_KR1OFF * y[3]
                        - QS_T_LA * y[3];

                return dy;
            }

            @Override
            public int getNumEq() {
                return 4;
            }

            @Override
            public double[] getICs() {
                return new double[]{0.05, 0.05, 0.05, 0.05};
            }
        }
    }

    static double bathUm(double t) {
        return (t >= T_BATH_ON && t < T_BATH_OFF) ? BATH_UM : 0.0;
    }

    static double meanFieldUm(BSimChemicalField field) {
        double sum = 0.0;
        int n = 0;
        for (int i = 0; i < GRID_X; i++) {
            for (int j = 0; j < GRID_Y; j++) {
                for (int k = 0; k < GRID_Z; k++) {
                    sum += field.getConc(i, j, k);
                    n++;
                }
            }
        }
        return (sum / n) / MOL_PER_UM3_PER_UM;
    }

    public static void main(String[] args) {
        boolean preview = args.length > 0 && "preview".equals(args[0]);

        BSim sim = new BSim();
        sim.setDt(DT);
        sim.setSimulationTime(SIM_TIME);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);

        final BSimChemicalField ahlField = new BSimChemicalField(
                sim, new int[]{GRID_X, GRID_Y, GRID_Z}, AHL_DIFFUSIVITY, AHL_DECAY_RATE);

        Random placement = new Random(RNG_SEED);
        while (bacteria.size() < N_CELLS) {
            Vector3d p = new Vector3d(
                    10.0 + placement.nextDouble() * (BOUND_X - 20.0),
                    10.0 + placement.nextDouble() * (BOUND_Y - 20.0),
                    BOUND_Z / 2.0);
            DaninoBacterium b = new DaninoBacterium(sim, p, ahlField);
            b.setRadius(1.0);
            b.setSurfaceAreaGrowthRate(GROWTH_RATE);
            if (!b.intersection(bacteria)) bacteria.add(b);
        }

        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                ahlField.update();
                ahlField.setConc(bathUm(sim.getTime()) * MOL_PER_UM3_PER_UM);
                for (DaninoBacterium b : bacteria) {
                    b.action();
                    b.updatePosition();
                }
            }
        });

        sim.setDrawer(new BSimP3DDrawer(sim, 800, 600) {
            @Override
            public void scene(PGraphics3D p3d) {
                p3d.ortho(0, (float) BOUND_X,
                        (float) BOUND_Y, 0,
                        -1000, 10000);
                p3d.camera((float) BOUND_X / 2f, (float) BOUND_Y / 2f, (float) BOUND_Y,
                        (float) BOUND_X / 2f, (float) BOUND_Y / 2f, 0f,
                        0f, 1f, 0f);
                p3d.perspective((float) Math.PI / 2f,
                        (float) BOUND_X / (float) BOUND_Y,
                        0.1f, 10000f);

                draw(ahlField, Color.ORANGE, (float) (255.0 / (BATH_UM * MOL_PER_UM3_PER_UM)));
                for (DaninoBacterium b : bacteria) {
                    double luxIMax = QS_KPLI / QS_DELTA1;
                    double r = Math.min(1.0, Math.max(0.0, b.y[0] / luxIMax));
                    int red = (int) (220 * (1.0 - r)) + 30;
                    int green = (int) (220 * r) + 30;
                    sphere(b.getPosition(), 3.0, new Color(red, green, 30), 255);
                }
            }
        });

        if (preview) {
            System.out.println("Danino D1: preview mode");
            sim.preview();
            return;
        }

        BSimUtils.generateDirectoryPath(EXPORT_DIR);
        BSimLogger logger = new BSimLogger(sim, EXPORT_DIR + CSV_NAME) {
            @Override
            public void before() {
                super.before();
                write("t_s;AHL_ext_uM_mean;LuxI_mean;AHL_in_mean;AiiA_mean;LA_mean;N");
            }

            @Override
            public void during() {
                double meanLuxI = 0, meanAhlIn = 0, meanAiiA = 0, meanLA = 0;
                for (DaninoBacterium b : bacteria) {
                    meanLuxI += b.y[0];
                    meanAhlIn += b.y[1];
                    meanAiiA += b.y[2];
                    meanLA += b.y[3];
                }
                int n = bacteria.size();
                if (n > 0) {
                    meanLuxI /= n;
                    meanAhlIn /= n;
                    meanAiiA /= n;
                    meanLA /= n;
                }
                write(String.format(Locale.US, "%.0f;%.12f;%.12e;%.12e;%.12e;%.12e;%d",
                        (double) Math.round(sim.getTime()),
                        meanFieldUm(ahlField),
                        meanLuxI, meanAhlIn, meanAiiA, meanLA, n));
            }
        };
        logger.setDt(1.0);
        sim.addExporter(logger);

        System.out.println("Danino D1: exporting to " + EXPORT_DIR + CSV_NAME);
        System.out.println("  conversion = field/602 (not field*1e-15)");
        System.out.println("  bath ON " + T_BATH_ON + " <= t < " + T_BATH_OFF
                + " at " + BATH_UM + " uM");
        PrintStream originalOut = System.out;
        try {
            System.setOut(new PrintStream(OutputStream.nullOutputStream()));
            sim.export();
        } finally {
            System.setOut(originalOut);
        }
        System.out.println("Danino D1: export finished. Completeness is the CSV last row t=5400.");
    }
}
