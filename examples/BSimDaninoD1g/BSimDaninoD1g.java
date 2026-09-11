package BSimDaninoD1g;

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
 * Isolated Danino D1g citation box. Copies D1 (602 conversion, open-loop bath,
 * Stage 10 QS_*) and changes only zero ICs plus explicit linear dilution µ
 * on LuxI, AiiA, and LA. BSim growth stays off; N stays 50.
 *
 * D1 remains a closed FAIL. This ODE is not re-integrated into Plan Stage 6.
 */
public class BSimDaninoD1g {

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

    /** Explicit Danino dilution; ln 2 / 30 min doubling. Not BSim growth. */
    static final double MU = Math.log(2.0) / 1800.0;

    static final double[] ICS = {0.0, 0.0, 0.0, 0.0};

    // QS kinetic parameters (from D1 / BSimReservoirStage10, fixed).
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
    static final double CELL_WALL_DIFF = 3.0 / TIME_ADJ;

    static final String EXPORT_DIR = "results/d1g_seed101/";
    static final String CSV_NAME = "d1g_timeseries.csv";
    static final String PARAMS_NAME = "d1g_params.csv";

    static final Vector<DaninoBacterium> bacteria = new Vector<DaninoBacterium>();

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
            this.y = new double[]{ICS[0], ICS[1], ICS[2], ICS[3]};
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
                        - (QS_DELTA1 * y[0]) / (QS_G * (y[0] + y[2]) + 1)
                        - MU * y[0];

                dy[1] = QS_KP2 * y[0]
                        - QS_KR1ON * (QS_LTOT - y[3]) * y[1]
                        + QS_KR1OFF * y[3]
                        - (QS_KCAT_AIIA * y[2] * y[1]) / (QS_KMAA + y[1])
                        - QS_T_A * y[1]
                        - CELL_WALL_DIFF * extraintradiff_uM;

                dy[2] = QS_A0AA
                        + QS_KPAA * (Math.pow(y[3], QS_N) / (Math.pow(QS_KMLA, QS_N) + Math.pow(y[3], QS_N)))
                        - (QS_DELTA2 * y[2]) / (QS_G * (y[0] + y[2]) + 1)
                        - MU * y[2];

                dy[3] = QS_KR1ON * (QS_LTOT - y[3]) * y[1]
                        - QS_KR1OFF * y[3]
                        - QS_T_LA * y[3]
                        - MU * y[3];

                return dy;
            }

            @Override
            public int getNumEq() {
                return 4;
            }

            @Override
            public double[] getICs() {
                return new double[]{ICS[0], ICS[1], ICS[2], ICS[3]};
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
                    double r = Math.min(1.0, Math.max(0.0, b.y[0] / Math.max(luxIMax, 1e-12)));
                    int red = (int) (220 * (1.0 - r)) + 30;
                    int green = (int) (220 * r) + 30;
                    sphere(b.getPosition(), 3.0, new Color(red, green, 30), 255);
                }
            }
        });

        if (preview) {
            System.out.println("Danino D1g: preview mode");
            System.out.println("  mu = " + MU + " /s");
            System.out.println("  ICs = {0, 0, 0, 0}");
            sim.preview();
            return;
        }

        BSimUtils.generateDirectoryPath(EXPORT_DIR);
        BSimLogger params = new BSimLogger(sim, EXPORT_DIR + PARAMS_NAME) {
            @Override
            public void before() {
                super.before();
                write("parameter;value");
                write(String.format(Locale.US, "mu_per_s;%.16e", MU));
                write("mu_formula;Math.log(2.0)/1800.0");
                write("ics;0,0,0,0");
                write("d1_ics;0.05,0.05,0.05,0.05");
                write("conversion;field/602");
                write("bath_uM;" + BATH_UM);
                write("n_cells;" + N_CELLS);
                write("bsim_growth_rate;" + GROWTH_RATE);
                write("dilution_on;LuxI,AiiA,LA");
                write("dilution_off;AHL_in");
            }

            @Override
            public void during() { }
        };
        sim.addExporter(params);

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

        System.out.println("Danino D1g: exporting to " + EXPORT_DIR + CSV_NAME);
        System.out.println("  mu = " + MU + " /s  (ln2/1800)");
        System.out.println("  ICs = {0, 0, 0, 0}");
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
        System.out.println("Danino D1g: export finished. Completeness is the CSV last row t=5400.");
    }
}
