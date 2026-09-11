package BSimReservoirStage4;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
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
 * Stage 4: ACs with per-cell independent sequence-driven poration.
 *
 * Each AC reads its own input file (input_sequence{id}.txt) with format:
 *   Line 1: bit_duration (seconds per bit for this AC)
 *   Line 2: space-separated 0/1 sequence
 *
 * This allows each AC to have unique timing and activation pattern,
 * enabling independent input channels for the reservoir.
 *
 * Two activation modes (global setting):
 *   CONTINUOUS — AC produces the entire time a bit is 1
 *   IMPULSE    — AC produces a brief burst (IMPULSE_DURATION) at the
 *                start of each bit=1 period
 */
public class BSimReservoirStage4 {

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
    static final double DECAY_RATE  = 0.01;    // 1/s  (tau = 100s)
    static final double PROD_RATE   = 1e6;     // molecules/s

    // --- Activation parameters ---
    static final double DOSE_THRESHOLD = 1.0;
    static final double DOSE_HIGH      = 10.0;
    static final double DOSE_RADIUS    = 30.0;  // um

    // --- Global mode ---
    static final double IMPULSE_DURATION = 1.0;  // seconds (impulse mode only)
    static final boolean IMPULSE_MODE   = false; // false = continuous, true = impulse

    // --- Sequence file pattern ---
    static final String SEQ_FILE_PATTERN = "input_sequence%d.txt";

    /**
     * An artificial cell (AC) with its own input sequence and bit duration.
     */
    static class AC {
        final int id;
        final Vector3d position;
        boolean activated = false;

        // Per-AC sequence parameters
        final int[] sequence;
        final double bitDuration;
        final double totalTime;

        AC(int id, Vector3d position, int[] sequence, double bitDuration) {
            this.id = id;
            this.position = position;
            this.sequence = sequence;
            this.bitDuration = bitDuration;
            this.totalTime = sequence.length * bitDuration;
        }

        /** Whether this AC should be producing at simulation time t. */
        boolean shouldProduce(double t) {
            int bitIndex = (int)(t / bitDuration);
            if (bitIndex >= sequence.length) bitIndex = sequence.length - 1;
            int bitValue = sequence[bitIndex];

            if (IMPULSE_MODE) {
                double tInBit = t - bitIndex * bitDuration;
                return (bitValue == 1) && (tInBit < IMPULSE_DURATION);
            } else {
                return (bitValue == 1);
            }
        }
    }

    /**
     * Read a per-AC sequence file.
     * Format: line 1 = bit_duration, line 2 = space-separated 0/1 values.
     * Returns [bitDuration, sequence...] packed — caller unpacks.
     */
    static Object[] readSequenceFile(String path, int acId) {
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            double bitDuration = Double.parseDouble(br.readLine().trim());
            String line = br.readLine().trim();
            String[] tokens = line.split("\\s+");
            int[] seq = new int[tokens.length];
            for (int i = 0; i < tokens.length; i++) {
                seq[i] = Integer.parseInt(tokens[i]);
            }
            return new Object[]{bitDuration, seq};
        } catch (IOException e) {
            System.err.println("Could not read " + path + ": " + e.getMessage());
            System.err.println("AC " + acId + ": using default sequence (5s, 1 0 1 0 1)");
            return new Object[]{5.0, new int[]{1, 0, 1, 0, 1}};
        }
    }

    public static void main(String[] args) {

        boolean exportData = true;
        if (args.length > 0 && args[0].equals("preview")) {
            exportData = false;
        }

        // --- Load per-AC sequences ---
        Vector3d[] positions = {
            new Vector3d(250, 250, 5),
            new Vector3d(500, 250, 5),
            new Vector3d(750, 250, 5)
        };

        final List<AC> acs = new ArrayList<>();
        double maxTotalTime = 0;

        for (int i = 0; i < positions.length; i++) {
            String seqFile = String.format(SEQ_FILE_PATTERN, i);
            Object[] result = readSequenceFile(seqFile, i);
            double bitDur = (Double) result[0];
            int[] seq = (int[]) result[1];

            AC ac = new AC(i, positions[i], seq, bitDur);
            acs.add(ac);

            if (ac.totalTime > maxTotalTime) maxTotalTime = ac.totalTime;

            // Print info
            System.out.print("AC " + i + " (" + bitDur + "s/bit, "
                    + seq.length + " bits, " + ac.totalTime + "s): ");
            for (int b : seq) System.out.print(b + " ");
            System.out.println();
        }

        System.out.println("Mode: " + (IMPULSE_MODE ? "IMPULSE (" + IMPULSE_DURATION + "s)" : "CONTINUOUS"));
        System.out.println("Simulation time: " + maxTotalTime + "s (longest AC sequence)");

        String timestamp = BSimUtils.timeStamp();
        String exportPath = "./results/" + timestamp + "/";

        // --- Simulation ---
        BSim sim = new BSim();
        sim.setDt(0.05);
        sim.setSimulationTime(maxTotalTime);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);

        // --- Signal field ---
        final BSimChemicalField signalField = new BSimChemicalField(sim,
                new int[]{GRID_X, GRID_Y, GRID_Z}, DIFFUSIVITY, DECAY_RATE);

        // --- Dose field ---
        final BSimChemicalField doseField = new BSimChemicalField(sim,
                new int[]{GRID_X, GRID_Y, GRID_Z}, 0, 0);

        final double voxelX = BOUND_X / GRID_X;
        final double voxelY = BOUND_Y / GRID_Y;

        // --- Ticker ---
        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                double t = sim.getTime();

                // Clear dose field
                for (int i = 0; i < GRID_X; i++)
                    for (int j = 0; j < GRID_Y; j++)
                        doseField.setConc(i, j, 0, 0);

                // Apply dose independently for each AC based on its own sequence
                for (AC ac : acs) {
                    if (t < ac.totalTime && ac.shouldProduce(t)) {
                        int cxMin = Math.max(0, (int)((ac.position.x - DOSE_RADIUS) / voxelX));
                        int cxMax = Math.min(GRID_X - 1, (int)((ac.position.x + DOSE_RADIUS) / voxelX));
                        int cyMin = Math.max(0, (int)((ac.position.y - DOSE_RADIUS) / voxelY));
                        int cyMax = Math.min(GRID_Y - 1, (int)((ac.position.y + DOSE_RADIUS) / voxelY));
                        for (int i = cxMin; i <= cxMax; i++)
                            for (int j = cyMin; j <= cyMax; j++)
                                doseField.setConc(i, j, 0, DOSE_HIGH);
                    }
                }

                // Each AC checks its local dose
                for (AC ac : acs) {
                    double dose = doseField.getConc(ac.position);
                    ac.activated = dose > DOSE_THRESHOLD;

                    if (ac.activated) {
                        signalField.addQuantity(ac.position, PROD_RATE * sim.getDt());
                    }
                }

                signalField.update();
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

                // Signal field (cyan)
                draw(signalField, Color.CYAN, (float)(255.0 / 540.0));

                // ACs: green if activated, red if not
                for (AC ac : acs) {
                    Color col = ac.activated ? Color.GREEN : Color.RED;
                    sphere(ac.position, 15.0, col, 255);
                }
            }
        });

        // --- Exporters ---
        if (exportData) {
            BSimUtils.generateDirectoryPath(exportPath);

            // State logger (every 0.5s)
            BSimLogger stateLogger = new BSimLogger(sim, exportPath + "stage4_response.csv") {
                @Override
                public void before() {
                    super.before();
                    StringBuilder header = new StringBuilder("time_s");
                    for (AC ac : acs) {
                        header.append(",ac").append(ac.id).append("_bit_index");
                        header.append(",ac").append(ac.id).append("_bit_value");
                        header.append(",ac").append(ac.id).append("_activated");
                        header.append(",ac").append(ac.id).append("_signal_conc");
                    }
                    write(header.toString());
                }
                @Override
                public void during() {
                    double t = sim.getTime();
                    StringBuilder line = new StringBuilder(sim.getFormattedTime());
                    for (AC ac : acs) {
                        int bi = (int)(t / ac.bitDuration);
                        if (bi >= ac.sequence.length) bi = ac.sequence.length - 1;
                        line.append(",").append(bi);
                        line.append(",").append(ac.sequence[bi]);
                        line.append(",").append(ac.activated ? 1 : 0);
                        line.append(",").append(String.format("%.2f", signalField.getConc(ac.position)));
                    }
                    write(line.toString());
                }
            };
            stateLogger.setDt(0.5);
            sim.addExporter(stateLogger);

            // Parameters
            BSimLogger paramLogger = new BSimLogger(sim, exportPath + "stage4_params.csv") {
                @Override
                public void before() {
                    super.before();
                    write("parameter,value");
                    write("num_acs," + acs.size());
                    for (AC ac : acs) {
                        StringBuilder seqStr = new StringBuilder();
                        for (int i = 0; i < ac.sequence.length; i++) {
                            if (i > 0) seqStr.append(" ");
                            seqStr.append(ac.sequence[i]);
                        }
                        write("ac" + ac.id + "_sequence," + seqStr);
                        write("ac" + ac.id + "_n_bits," + ac.sequence.length);
                        write("ac" + ac.id + "_bit_duration_s," + ac.bitDuration);
                        write("ac" + ac.id + "_total_time_s," + ac.totalTime);
                    }
                    write("impulse_mode," + IMPULSE_MODE);
                    write("impulse_duration_s," + IMPULSE_DURATION);
                    write("sim_total_time_s," + sim.getSimulationTime());
                    write("signal_D," + DIFFUSIVITY);
                    write("signal_decay," + DECAY_RATE);
                    write("signal_prod_rate," + PROD_RATE);
                    write("dose_threshold," + DOSE_THRESHOLD);
                    write("dt," + sim.getDt());
                }
                @Override
                public void during() { }
            };
            sim.addExporter(paramLogger);

            System.out.println("Stage 4: exporting to " + exportPath);
            sim.export();

        } else {
            System.out.println("Stage 4: preview mode");
            sim.preview();
        }
    }
}
