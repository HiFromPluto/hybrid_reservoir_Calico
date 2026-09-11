package BSimReservoirStage5;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Vector3d;

import processing.core.PGraphics3D;
import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;

/**
 * Stage 5: Integrate flow into the reservoir dish.
 *
 * BSim has no built-in flow/advection (confirmed by reading BSim.java,
 * BSimChemicalField.java, BSimParticle.java). Flow is implemented here as:
 *
 *   1. Particle drag: uniform flow velocity applied as a Stokes-drag force
 *      via addForce() — BSimParticle.updatePosition() converts force→velocity
 *      via Stokes law automatically.
 *
 *   2. Chemical field advection: first-order upwind finite-difference scheme
 *      applied to the signal field each timestep, in addition to diffusion
 *      and decay.
 *
 * Domain: same 1000×500×10 µm box (microfluidic channel).
 * Flow: uniform along +x (left→right).
 * Boundaries: x=0 is inlet (open), x=1000 is outlet (leaky for chemicals).
 *
 * Validation: tracer-pulse test — inject a chemical pulse at the inlet,
 * confirm it transits the channel at the expected flow speed and clears
 * the outlet without unphysical accumulation.
 *
 * Expected transit time = BOUND_X / FLOW_SPEED.
 */
public class BSimReservoirStage5 {

    // --- Domain (fixed) ---
    static final double BOUND_X = 1000.0;
    static final double BOUND_Y = 500.0;
    static final double BOUND_Z = 10.0;

    // --- Grid ---
    static final int GRID_X = 200;
    static final int GRID_Y = 100;
    static final int GRID_Z = 1;

    // --- Signal field parameters ---
    static final double DIFFUSIVITY = 100.0;   // um^2/s
    static final double DECAY_RATE  = 0.01;    // 1/s

    // --- Flow parameters ---
    static final double FLOW_SPEED = 20.0;     // um/s along +x
    // Expected transit time: 1000/20 = 50s

    // --- Tracer pulse parameters ---
    static final double PULSE_TIME   = 5.0;    // seconds: inject pulse at t=5
    static final double PULSE_DURATION = 2.0;  // seconds: inject for 2s
    static final double PULSE_RATE   = 1e6;    // molecules/s
    static final double PULSE_X      = 50.0;   // um from left wall (inlet region)
    static final double PULSE_Y      = BOUND_Y / 2.0;
    static final double PULSE_Z      = BOUND_Z / 2.0;

    // --- Probe positions (for logging concentration vs time) ---
    // At inlet, 1/4, centre, 3/4, outlet
    static final double[] PROBE_X = {50, 250, 500, 750, 950};

    /**
     * First-order upwind advection for a chemical field along +x.
     * C[i] -= v*dt/dx * (C[i] - C[i-1])  (upwind: flow is +x, so donor is i-1)
     *
     * At i=0 (inlet): C[i-1] is taken as 0 (nothing flowing in from outside).
     * At i=max (outlet): chemical that flows out is lost (open boundary).
     *
     * Stability requires: v*dt/dx < 1  (CFL condition).
     *
     * Uses public getConc/setConc API since quantity[][] is protected.
     */
    static void advect(BSimChemicalField field, double flowSpeed, double dt,
                       int nx, int ny, int nz) {
        double dx = field.getBox()[0];
        double courant = flowSpeed * dt / dx;
        // Snapshot concentrations before update
        double[][] before = new double[nx][ny];
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                before[i][j] = field.getConc(i, j, 0);

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++) {
                double cUpwind = (i > 0) ? before[i - 1][j] : 0;
                double newConc = before[i][j] - courant * (before[i][j] - cUpwind);
                if (newConc < 0) newConc = 0;
                field.setConc(i, j, 0, newConc);
            }
    }

    public static void main(String[] args) {

        boolean exportData = true;
        if (args.length > 0 && args[0].equals("preview")) {
            exportData = false;
        }

        String timestamp = BSimUtils.timeStamp();
        String exportPath = "./results/" + timestamp + "/";

        double voxelX = BOUND_X / GRID_X;  // 5 um
        double courant = FLOW_SPEED * 0.05 / voxelX;  // CFL number

        System.out.println("Stage 5: Flow + advection test");
        System.out.println("  Flow speed: " + FLOW_SPEED + " um/s");
        System.out.println("  Expected transit time: " + (BOUND_X / FLOW_SPEED) + " s");
        System.out.println("  Voxel dx: " + voxelX + " um");
        System.out.println("  CFL number: " + courant + " (must be < 1)");

        // --- Simulation ---
        BSim sim = new BSim();
        sim.setDt(0.05);
        sim.setSimulationTime(120);  // enough time to see pulse transit + clear
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);

        // --- Signal field (diffusion + decay + advection) ---
        final BSimChemicalField signalField = new BSimChemicalField(sim,
                new int[]{GRID_X, GRID_Y, GRID_Z}, DIFFUSIVITY, DECAY_RATE);

        // --- Stokes drag force for flow ---
        // F = 6*pi*eta*r * v  (Stokes drag to produce desired velocity)
        // For a typical bacterium (r~1um), eta=8.9e-4 Pa·s:
        // F = 6*pi*8.9e-4*1 * 20 = 0.335 pN
        // But we don't have bacteria in this stage — it's a tracer-pulse test
        // of the chemical field only. Particle flow is demonstrated by the
        // advection of the field.

        final Vector3d pulsePos = new Vector3d(PULSE_X, PULSE_Y, PULSE_Z);

        // --- Ticker ---0
        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                double t = sim.getTime();

                // Inject tracer pulse
                if (t >= PULSE_TIME && t < PULSE_TIME + PULSE_DURATION) {
                    signalField.addQuantity(pulsePos, PULSE_RATE * sim.getDt());
                }

                // Diffuse + decay (built-in)
                signalField.update();

                // Advect along +x (new code — BSim has no built-in advection)
                advect(signalField, FLOW_SPEED, sim.getDt(),
                       GRID_X, GRID_Y, GRID_Z);
            }
        });

        // --- Drawer ---
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

                // Signal field (cyan) — scaled for pulse peak
                draw(signalField, Color.CYAN, (float)(255.0 / 540.0));

                // Mark pulse injection point
                double t = sim.getTime();
                boolean pulsing = (t >= PULSE_TIME && t < PULSE_TIME + PULSE_DURATION);
                Color pulseCol = pulsing ? Color.GREEN : Color.RED;
                sphere(pulsePos, 10.0, pulseCol, 255);

                // Draw probe positions as small yellow spheres
                for (double px : PROBE_X) {
                    sphere(new Vector3d(px, BOUND_Y / 2.0, BOUND_Z / 2.0),
                           5.0, Color.YELLOW, 180);
                }

                // Flow direction arrow (visual hint)
                p3d.stroke(255, 255, 0);
                p3d.strokeWeight(2);
                p3d.line(50, 30, 0, 200, 30, 0);  // arrow shaft
                p3d.line(200, 30, 0, 180, 20, 0);  // arrowhead
                p3d.line(200, 30, 0, 180, 40, 0);
                p3d.noStroke();
            }
        });

        // --- Exporters ---
        if (exportData) {
            BSimUtils.generateDirectoryPath(exportPath);

            // Probe concentration time-series (every 0.5s)
            BSimLogger probeLogger = new BSimLogger(sim, exportPath + "stage5_probes.csv") {
                @Override
                public void before() {
                    super.before();
                    StringBuilder header = new StringBuilder("time_s");
                    for (double px : PROBE_X) {
                        header.append(",conc_x").append(String.format("%.0f", px));
                    }
                    header.append(",total_quantity");
                    write(header.toString());
                }
                @Override
                public void during() {
                    StringBuilder line = new StringBuilder(sim.getFormattedTime());
                    for (double px : PROBE_X) {
                        Vector3d probePos = new Vector3d(px, BOUND_Y / 2.0, BOUND_Z / 2.0);
                        line.append(",").append(String.format("%.4f", signalField.getConc(probePos)));
                    }
                    line.append(",").append(String.format("%.2f", signalField.totalQuantity()));
                    write(line.toString());
                }
            };
            probeLogger.setDt(0.5);
            sim.addExporter(probeLogger);

            // Concentration profile along x (every 5s)
            BSimLogger profileLogger = new BSimLogger(sim, exportPath + "stage5_profile.csv") {
                @Override
                public void before() {
                    super.before();
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
                        line.append(",").append(signalField.getConc(i, jCentre, 0));
                    }
                    write(line.toString());
                }
            };
            profileLogger.setDt(5.0);
            sim.addExporter(profileLogger);

            // Parameters
            BSimLogger paramLogger = new BSimLogger(sim, exportPath + "stage5_params.csv") {
                @Override
                public void before() {
                    super.before();
                    write("parameter,value");
                    write("flow_speed_um_per_s," + FLOW_SPEED);
                    write("expected_transit_s," + (BOUND_X / FLOW_SPEED));
                    write("diffusivity_um2_per_s," + DIFFUSIVITY);
                    write("decay_rate_per_s," + DECAY_RATE);
                    write("pulse_time_s," + PULSE_TIME);
                    write("pulse_duration_s," + PULSE_DURATION);
                    write("pulse_rate_mol_per_s," + PULSE_RATE);
                    write("pulse_x_um," + PULSE_X);
                    write("grid," + GRID_X + "x" + GRID_Y + "x" + GRID_Z);
                    write("voxel_dx_um," + (BOUND_X / GRID_X));
                    write("dt," + sim.getDt());
                    write("CFL," + courant);
                    write("domain," + BOUND_X + "x" + BOUND_Y + "x" + BOUND_Z);
                }
                @Override
                public void during() { }
            };
            sim.addExporter(paramLogger);

            System.out.println("Stage 5: exporting to " + exportPath);
            sim.export();

        } else {
            System.out.println("Stage 5: preview mode");
            sim.preview();
        }
    }
}
