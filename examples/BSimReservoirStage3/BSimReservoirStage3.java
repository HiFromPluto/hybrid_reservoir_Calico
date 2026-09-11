package BSimReservoirStage3;

import java.awt.Color;

import javax.vecmath.Vector3d;

import processing.core.PGraphics3D;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;

/**
 * Stage 3: Single AC (artificial cell) source releasing into a chemical field.
 * Zero bacteria — purely validates the diffusion/decay profile.
 *
 * The AC continuously injects molecules into a BSimChemicalField at a
 * constant rate. The field diffuses and decays. Validation: the steady-state
 * concentration profile should decay with distance from the source, with
 * characteristic length scale L = sqrt(D/k).
 *
 * Chemical field API reused from BSimChemicalFieldTest / BSimQuorumOscillator.
 * The discrete-vesicle alternative (BSimVesiculation / BSimVesicleActivation)
 * plugs in at Stage 4 where triggered release is needed.
 */
public class BSimReservoirStage3 {

    // --- Domain (fixed across all stages) ---
    static final double BOUND_X = 1000.0;
    static final double BOUND_Y = 500.0;
    static final double BOUND_Z = 10.0;

    // --- Chemical field grid (from plan: 5 um/voxel) ---
    static final int GRID_X = 200;
    static final int GRID_Y = 100;
    static final int GRID_Z = 1;
    // Voxel size: 5 x 5 x 10 um

    // --- Chemical field parameters ---
    static final double DIFFUSIVITY  = 100.0;   // um^2/s
    static final double DECAY_RATE   = 0.01;    // 1/s
    static final double PROD_RATE    = 1e6;     // molecules/s

    // Characteristic diffusion length: sqrt(D/k) = sqrt(100/0.01) = 100 um
    static final double CHAR_LENGTH = Math.sqrt(DIFFUSIVITY / DECAY_RATE);

    // --- AC position (domain centre) ---
    static final Vector3d AC_POS = new Vector3d(BOUND_X / 2.0, BOUND_Y / 2.0, BOUND_Z / 2.0);

    public static void main(String[] args) {

        boolean exportData = true;
        if (args.length > 0 && args[0].equals("preview")) {
            exportData = false;
        }

        String timestamp = BSimUtils.timeStamp();
        String exportPath = "./results/" + timestamp + "/";

        // --- Simulation setup ---
        BSim sim = new BSim();
        // dt must satisfy stability: D*dt/dx^2 < 0.5
        // kX = 100*0.05/25 = 0.2 (safe)
        sim.setDt(0.05);
        sim.setSimulationTime(500);     // enough for steady state (~5x L^2/D = 500s)
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);

        // --- Chemical field ---
        final BSimChemicalField field = new BSimChemicalField(sim,
                new int[]{GRID_X, GRID_Y, GRID_Z}, DIFFUSIVITY, DECAY_RATE);

        // --- Ticker: inject molecules + update field ---
        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                // AC produces molecules at constant rate
                field.addQuantity(AC_POS, PROD_RATE * sim.getDt());
                // Diffuse + decay
                field.update();
            }
        });

        // --- Drawer ---
        sim.setDrawer(new BSimP3DDrawer(sim, 800, 600) {
            @Override
            public void scene(PGraphics3D p3d) {
                // Camera (same as Stage 1/2)
                p3d.ortho(0, (float) BOUND_X,
                          (float) BOUND_Y, 0,
                          -1000, 10000);
                p3d.camera((float) BOUND_X / 2f, (float) BOUND_Y / 2f, (float) BOUND_Y,
                           (float) BOUND_X / 2f, (float) BOUND_Y / 2f, 0f,
                           0f, 1f, 0f);
                p3d.perspective((float) Math.PI / 2f,
                        (float) BOUND_X / (float) BOUND_Y,
                        0.1f, 10000f);

                // Draw chemical field. Peak conc ~270 molecules/um^3.
                // alphaGrad scales conc → [0,255] alpha.
                // 255/270 ≈ 0.94 → peak is fully opaque; halve it for a gradient.
                draw(field, Color.CYAN, (float)(255.0 / 540.0));

                // Draw AC as a red sphere
                sphere(AC_POS, 15.0, Color.RED, 255);
            }
        });

        // --- Exporters ---
        if (exportData) {
            BSimUtils.generateDirectoryPath(exportPath);

            // Concentration profile along x-axis through centre (y=GRID_Y/2, z=0)
            BSimLogger profileLogger = new BSimLogger(sim, exportPath + "stage3_profile.csv") {
                @Override
                public void before() {
                    super.before();
                    // Header: time, then concentration at each x-voxel
                    StringBuilder header = new StringBuilder("time_s");
                    for (int i = 0; i < GRID_X; i++) {
                        double xMid = (i + 0.5) * (BOUND_X / GRID_X);
                        header.append(",x_").append(String.format("%.1f", xMid));
                    }
                    write(header.toString());
                }
                @Override
                public void during() {
                    int jCentre = GRID_Y / 2;
                    StringBuilder line = new StringBuilder(sim.getFormattedTime());
                    for (int i = 0; i < GRID_X; i++) {
                        line.append(",").append(field.getConc(i, jCentre, 0));
                    }
                    write(line.toString());
                }
            };
            profileLogger.setDt(10.0);  // every 10 seconds
            sim.addExporter(profileLogger);

            // Summary: time, peak concentration, total quantity, production-decay balance
            BSimLogger summaryLogger = new BSimLogger(sim, exportPath + "stage3_summary.csv") {
                @Override
                public void before() {
                    super.before();
                    write("time_s,peak_conc,total_quantity,total_produced");
                }
                double totalProduced = 0;
                @Override
                public void during() {
                    totalProduced += PROD_RATE * getDt();
                    // Peak = concentration at source voxel
                    int[] srcBox = field.boxCoords(AC_POS);
                    double peakConc = field.getConc(srcBox[0], srcBox[1], srcBox[2]);
                    write(sim.getFormattedTime()
                            + "," + peakConc
                            + "," + field.totalQuantity()
                            + "," + totalProduced);
                }
            };
            summaryLogger.setDt(10.0);
            sim.addExporter(summaryLogger);

            // Parameters
            BSimLogger paramLogger = new BSimLogger(sim, exportPath + "stage3_params.csv") {
                @Override
                public void before() {
                    super.before();
                    write("parameter,value");
                    write("diffusivity_um2_per_s," + DIFFUSIVITY);
                    write("decay_rate_per_s," + DECAY_RATE);
                    write("production_rate_mol_per_s," + PROD_RATE);
                    write("char_length_um," + CHAR_LENGTH);
                    write("grid," + GRID_X + "x" + GRID_Y + "x" + GRID_Z);
                    write("voxel_um," + (BOUND_X/GRID_X) + "x" + (BOUND_Y/GRID_Y) + "x" + (BOUND_Z/GRID_Z));
                    write("dt," + sim.getDt());
                    write("simTime_s," + sim.getSimulationTime());
                    write("domain," + BOUND_X + "x" + BOUND_Y + "x" + BOUND_Z);
                    write("kX," + (DIFFUSIVITY * sim.getDt() / Math.pow(BOUND_X/GRID_X, 2)));
                    write("kY," + (DIFFUSIVITY * sim.getDt() / Math.pow(BOUND_Y/GRID_Y, 2)));
                }
                @Override
                public void during() { }
            };
            sim.addExporter(paramLogger);

            System.out.println("Stage 3: exporting to " + exportPath);
            System.out.println("  Char. length = " + String.format("%.0f", CHAR_LENGTH) + " um");
            System.out.println("  Stability kX = " + (DIFFUSIVITY * sim.getDt() / Math.pow(BOUND_X/GRID_X, 2)));
            sim.export();

        } else {
            System.out.println("Stage 3: preview mode (GUI)");
            System.out.println("  Char. length = " + String.format("%.0f", CHAR_LENGTH) + " um");
            sim.preview();
        }
    }
}
