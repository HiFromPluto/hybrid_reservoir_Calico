package BSimReservoirStage10;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
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
 * Stage 10a: Repellent chemotaxis.
 *
 * Changes from Stage 8:
 *   1. NEW repellentField -- AC-driven, same D/decay/grid as signalField.
 *      ACs are typed as ATTRACTANT (produce into signalField) or REPELLENT
 *      (produce into repellentField). Type is set per-AC in the constructor.
 *
 *   2. ReservoirBacterium overrides pEndRun() to integrate both signals:
 *        - Moving up attractant gradient -> extend run (pEndRunUp, Berg 1972)
 *        - Moving up repellent gradient  -> shorten run (pEndRunElse)
 *        - Moving down repellent gradient -> extend run (pEndRunUp)
 *      This is a SIGN-FLIP of the attractant mechanism, labeled as a plumbing
 *      test per POST_STAGE6_ROADMAP.md Stage 10a. Real E. coli repellent taxis
 *      (leucine, nickel, low pH) has different receptor adaptation kinetics
 *      than attractant taxis. Literature-sourced repellent-specific parameters
 *      are deferred to a future grounding pass.
 *
 *   3. Repellent memory: separate memory array for repellent field gradient
 *      sensing, same short/long-term windows as attractant (Berg-cited 1s/3s).
 *
 * NOT changed from Stage 8:
 *   - 4-variable QS ODE (LuxI/AHL/AiiA/LA, Danino et al. 2010)
 *   - AHL field (bacterium-produced, D=159, decay=2.76e-3/60)
 *   - Growth rate 0.00698 um^2/s (30 min doubling)
 *   - Death: density-dependent stochastic (ungrounded, flagged for Stage 9)
 *   - Reflective Y boundaries, solid walls, flow=0
 */
public class BSimReservoirStage10Control {

    // =================================================================
    //  Domain (fixed across all stages)
    // =================================================================
    static final double BOUND_X = 1000.0;
    static final double BOUND_Y = 500.0;
    static final double BOUND_Z = 10.0;

    // =================================================================
    //  Grid
    // =================================================================
    static final int GRID_X = 200;
    static final int GRID_Y = 100;
    static final int GRID_Z = 1;

    // =================================================================
    //  Signal field parameters (AC-driven, for chemotaxis)
    //  Used for BOTH attractant and repellent fields (same molecule class)
    // =================================================================
    static final double DIFFUSIVITY = 100.0;   // um^2/s
    static final double DECAY_RATE  = 0.01;    // 1/s  (tau = 100s)
    static final double PROD_RATE   = 1e6;     // molecules/s per active AC

    // =================================================================
    //  AHL field parameters (Stage 8 -- bacterium-produced)
    // =================================================================
    static final double AHL_DIFFUSIVITY = 159.0;
    static final double AHL_DECAY_RATE  = 2.76e-3 / 60;

    // =================================================================
    //  QS kinetic parameters (from BSimEntrainment_PIDCtrl, fixed)
    // =================================================================
    static final double TIME_ADJ = 60.0;

    static final double QS_DELTA1   = 0.8487   / TIME_ADJ;
    static final double QS_DELTA2   = 0.0234   / TIME_ADJ;
    static final double QS_G        = 0.0412;
    static final double QS_KP2      = 9.0      / TIME_ADJ;
    static final double QS_KR1OFF   = 6e-6     / TIME_ADJ;
    static final double QS_KR1ON    = 5.99e-5  / TIME_ADJ;
    static final double QS_KCAT_AIIA = 2631.4  / TIME_ADJ;
    static final double QS_T_A      = 0.00276  / TIME_ADJ;
    static final double QS_T_LA     = 0.024    / TIME_ADJ;
    static final double QS_A0LI     = 7.785e-6 / TIME_ADJ;
    static final double QS_A0AA     = 6.183e-6 / TIME_ADJ;
    static final double QS_KPLI     = 0.9      / TIME_ADJ;
    static final double QS_KPAA     = 0.9      / TIME_ADJ;
    static final double QS_KMLA     = 1e-2;
    static final double QS_KMAA     = 1200.0;
    static final double QS_LTOT     = 15.0;
    static final double QS_N        = 2.0;

    static final double CELL_WALL_DIFF = 3.0 / TIME_ADJ;

    // =================================================================
    //  AC parameters (Stage 4)
    // =================================================================
    static final double DOSE_THRESHOLD = 1.0;
    static final double DOSE_HIGH      = 10.0;
    static final double DOSE_RADIUS    = 30.0;
    static final double IMPULSE_DURATION = 1.0;
    static final boolean IMPULSE_MODE   = false;
    static final String SEQ_FILE_PATTERN = "input_sequence%d.txt";

    // =================================================================
    //  AC signal types
    // =================================================================
    // CONTROL: AC1 is INERT (exists but produces nothing)
    // AC0 (x=250): attractant
    // AC1 (x=500): inert (no signal)
    // AC2 (x=750): attractant
    static final boolean[] AC_IS_REPELLENT = {false, false, false};
    static final boolean[] AC_IS_INERT     = {false, true, false};

    // =================================================================
    //  Flow parameters
    // =================================================================
    static final double FLOW_SPEED = 0.0;

    // =================================================================
    //  Growth / death parameters (Stage 7)
    // =================================================================
    static final double GROWTH_RATE       = 4.0 * Math.PI / 1800.0;
    static final double EXPECTED_T_GEN    = 4.0 * Math.PI / GROWTH_RATE;
    static final double DEATH_RATE_AT_CAP = 1.0 / EXPECTED_T_GEN;

    // =================================================================
    //  Population parameters
    // =================================================================
    static final int    INITIAL_POP       = 500;
    static final int    CARRYING_CAPACITY = 2000;

    // =================================================================
    //  AC class (from Stage 4, extended with signal type)
    // =================================================================
    static class AC {
        final int id;
        final Vector3d position;
        boolean activated = false;
        final int[] sequence;
        final double bitDuration;
        final double totalTime;
        final boolean repellent;  // true = produces repellent, false = attractant
        final boolean inert;     // true = produces nothing (control)

        AC(int id, Vector3d position, int[] sequence, double bitDuration, boolean repellent, boolean inert) {
            this.id = id;
            this.position = position;
            this.sequence = sequence;
            this.bitDuration = bitDuration;
            this.totalTime = sequence.length * bitDuration;
            this.repellent = repellent;
            this.inert = inert;
        }

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

    static Object[] readSequenceFile(String path, int acId) {
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            double bitDuration = Double.parseDouble(br.readLine().trim());
            String line = br.readLine().trim();
            String[] tokens = line.split("\\s+");
            int[] seq = new int[tokens.length];
            for (int i = 0; i < tokens.length; i++)
                seq[i] = Integer.parseInt(tokens[i]);
            return new Object[]{bitDuration, seq};
        } catch (IOException e) {
            System.err.println("Could not read " + path + " for AC " + acId + ", using default");
            return new Object[]{5.0, new int[]{1, 0, 1, 0, 1}};
        }
    }

    // =================================================================
    //  Bacterium with dual chemotaxis + QS + growth
    // =================================================================
    static int nextId = 0;
    static final Vector<ReservoirBacterium> bacteria = new Vector<>();
    static final Vector<ReservoirBacterium> children = new Vector<>();
    static final Vector<ReservoirBacterium> removals = new Vector<>();

    static class ReservoirBacterium extends BSimBacterium {
        final int id;
        double lastDivisionTime;
        double[] y;
        BSimOdeSystem grn;
        final BSimChemicalField signalField;
        final BSimChemicalField repellentField;
        final BSimChemicalField ahlField;
        double diffConc;

        // Repellent gradient memory (separate from BSimBacterium's attractant memory).
        // Same short/long-term windows as attractant (Berg 1972: 1s short, 3s long).
        // Sign-flip plumbing: uses identical algorithm, just on a different field.
        double[] repMemory;

        public ReservoirBacterium(BSim sim, Vector3d position, double birthTime,
                                  BSimChemicalField signalField,
                                  BSimChemicalField repellentField,
                                  BSimChemicalField ahlField) {
            super(sim, position);
            this.id = nextId++;
            this.lastDivisionTime = birthTime;
            this.y = new double[]{0.05, 0.05, 0.05, 0.05};
            this.signalField = signalField;
            this.repellentField = repellentField;
            this.ahlField = ahlField;
            this.grn = new QSGRN();
            // Attractant chemotaxis (Berg-cited, Stage 7)
            setGoal(signalField);
            // Allocate repellent memory (same length as attractant memory)
            int memLen = sim.timesteps(shortTermMemoryDuration + longTermMemoryDuration);
            repMemory = new double[memLen];
            double initConc = repellentField.getConc(position);
            for (int i = 0; i < repMemory.length; i++) repMemory[i] = initConc;
        }

        /**
         * Repellent gradient detection -- same algorithm as BSimBacterium's
         * movingUpGradient() (Schnitzer, Berg et al.), applied to repellentField.
         *
         * PLUMBING TEST (POST_STAGE6_ROADMAP.md Stage 10a): uses identical
         * sensitivity, memory windows, and comparison as attractant. Real E. coli
         * repellent taxis has different receptor adaptation kinetics -- deferred
         * to a future grounding pass with literature-sourced parameters.
         */
        public boolean movingUpRepellentGradient() {
            double shortTermCounter = 0, longTermCounter = 0;

            System.arraycopy(repMemory, 0, repMemory, 1, repMemory.length - 1);
            repMemory[0] = repellentField.getConc(position);

            for (int i = 0; i < repMemory.length; i++) {
                if (i < shortTermMemoryLength) {
                    shortTermCounter += repMemory[i];
                } else {
                    longTermCounter += repMemory[i];
                }
            }
            double shortTermMean = shortTermCounter / shortTermMemoryLength;
            double longTermMean  = longTermCounter / longTermMemoryLength;

            return shortTermMean - longTermMean > sensitivity;
        }

        /**
         * Integrated chemotaxis: attractant extends runs, repellent shortens them.
         *
         * Berg's cited values (BSimBacterium defaults):
         *   pEndRunUp   = 1/1.07 s^-1  (longer runs when moving up attractant)
         *   pEndRunElse = 1/0.86 s^-1  (shorter runs otherwise)
         *
         * Integration logic (sign-flip plumbing):
         *   Up attractant, not up repellent -> extend run  (pEndRunUp)
         *   Up repellent, not up attractant -> shorten run (pEndRunElse)
         *   Both or neither                 -> baseline    (pEndRunElse)
         *
         * This means moving DOWN a repellent gradient (fleeing) while not
         * near attractant still gives baseline tumbling -- the bacterium
         * doesn't actively "seek" low-repellent areas, it just avoids
         * climbing the gradient. This matches the simple sign-flip model.
         */
        @Override
        public double pEndRun() {
            boolean upAttractant = (goal != null) && movingUpGradient();
            boolean upRepellent  = movingUpRepellentGradient();

            if (upAttractant && !upRepellent) return pEndRunUp;
            // All other cases: baseline tumble rate
            return pEndRunElse;
        }

        @Override
        public void action() {
            super.action();

            addForce(new Vector3d(stokesCoefficient() * FLOW_SPEED, 0, 0));

            diffConc = y[1] * 1e15 - ahlField.getConc(position);

            y = BSimOdeSolver.rungeKutta45(grn, sim.getTime(), y, sim.getDt());

            for (int i = 0; i < y.length; i++)
                if (y[i] < 0) y[i] = 0;

            double cellVol = (4.0 / 3.0) * Math.PI * Math.pow(radius, 3);
            double amountToField = diffConc * CELL_WALL_DIFF * sim.getDt() * cellVol;
            ahlField.addQuantity(position, amountToField);

            double pDeath = sim.getDt() * DEATH_RATE_AT_CAP
                    * ((double) bacteria.size() / CARRYING_CAPACITY);
            if (Math.random() < pDeath) {
                removals.add(this);
            }
        }

        @Override
        public void updatePosition() {
            super.updatePosition();

            double boundY = sim.getBound().y;
            if (position.y < 0.0)         position.y = -position.y;
            else if (position.y > boundY) position.y = 2.0 * boundY - position.y;
            position.y = Math.max(0.0, Math.min(boundY, position.y));

            position.z = Math.max(0.0, Math.min(sim.getBound().z, position.z));
        }

        @SuppressWarnings("unchecked")
        @Override
        public void replicate() {
            double divisionTime = sim.getTime();
            setRadiusFromSurfaceArea(surfaceArea(replicationRadius) / 2);

            Vector3d childPos = new Vector3d(position);
            childPos.x += 2.0 * radius * (Math.random() - 0.5);
            childPos.y += 2.0 * radius * (Math.random() - 0.5);

            ReservoirBacterium child = new ReservoirBacterium(sim,
                    childPos, divisionTime, signalField, repellentField, ahlField);
            child.setRadius(radius);
            child.setSurfaceAreaGrowthRate(surfaceAreaGrowthRate);
            child.setChildList(childList);
            child.y = new double[]{ this.y[0], this.y[1], this.y[2], this.y[3] };

            childList.add(child);
            this.lastDivisionTime = divisionTime;
        }

        public double reporterNormalised() {
            double luxIMax = QS_KPLI / (QS_DELTA1);
            return Math.min(1.0, Math.max(0.0, y[0] / luxIMax));
        }

        public double extAHL_uM() {
            return ahlField.getConc(position) * 1e-15;
        }

        class QSGRN implements BSimOdeSystem {
            @Override
            public double[] derivativeSystem(double t, double[] y) {
                double extraintradiff_uM = diffConc * 1e-15;

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
            public int getNumEq() { return 4; }

            @Override
            public double[] getICs() { return new double[]{0.05, 0.05, 0.05, 0.05}; }
        }
    }

    // =================================================================
    //  Upwind advection (from Stage 5)
    // =================================================================
    static void advect(BSimChemicalField field, double flowSpeed, double dt,
                       int nx, int ny, int nz) {
        double dx = field.getBox()[0];
        double courant = flowSpeed * dt / dx;
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

    // =================================================================
    //  Main
    // =================================================================
    public static void main(String[] args) {

        boolean exportData = true;
        if (args.length > 0 && args[0].equals("preview")) {
            exportData = false;
        }

        // --- Load per-AC sequences ---
        Vector3d[] acPositions = {
                new Vector3d(250, 250, 5),
                new Vector3d(500, 250, 5),
                new Vector3d(750, 250, 5)
        };
        final List<AC> acs = new ArrayList<>();
        double maxACTime = 0;
        for (int i = 0; i < acPositions.length; i++) {
            String seqFile = String.format(SEQ_FILE_PATTERN, i);
            Object[] result = readSequenceFile(seqFile, i);
            double bitDur = (Double) result[0];
            int[] seq = (int[]) result[1];
            boolean repel = (i < AC_IS_REPELLENT.length) && AC_IS_REPELLENT[i];
            boolean inert = (i < AC_IS_INERT.length) && AC_IS_INERT[i];
            AC ac = new AC(i, acPositions[i], seq, bitDur, repel, inert);
            acs.add(ac);
            if (ac.totalTime > maxACTime) maxACTime = ac.totalTime;
            String typeLabel = inert ? "INERT" : (repel ? "REPELLENT" : "ATTRACTANT");
            System.out.print("AC " + i + " (" + typeLabel
                    + ", " + bitDur + "s/bit, " + seq.length + " bits, " + ac.totalTime + "s): ");
            for (int b : seq) System.out.print(b + " ");
            System.out.println();
        }

        final double simTime = Math.max(maxACTime + 30.0, 3600.0);

        String timestamp = BSimUtils.timeStamp();
        String exportPath = "./results/" + timestamp + "/";

        double voxelX = BOUND_X / GRID_X;
        double courant = FLOW_SPEED * 0.05 / voxelX;
        double kxAHL = AHL_DIFFUSIVITY * 0.05 / (voxelX * voxelX);

        System.out.println("Stage 10a CONTROL: AC1 inert (no repellent)");
        System.out.println("  Attractant field: D=" + DIFFUSIVITY + ", decay=" + DECAY_RATE);
        System.out.println("  Repellent field: D=" + DIFFUSIVITY + ", decay=" + DECAY_RATE);
        System.out.println("  AHL field: D=" + AHL_DIFFUSIVITY + ", decay=" + AHL_DECAY_RATE);
        System.out.println("    AHL kX=" + String.format("%.3f", kxAHL) + " (must be < 0.5)");
        System.out.println("  Flow: " + FLOW_SPEED + " um/s, CFL=" + courant);
        System.out.println("  Initial pop: " + INITIAL_POP);
        System.out.println("  Sim time: " + simTime + "s");

        // --- Simulation ---
        BSim sim = new BSim();
        sim.setDt(0.05);
        sim.setSimulationTime(simTime);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        sim.setSolid(true, true, true);

        // --- Attractant field (AC-driven, for chemotaxis toward) ---
        final BSimChemicalField signalField = new BSimChemicalField(sim,
                new int[]{GRID_X, GRID_Y, GRID_Z}, DIFFUSIVITY, DECAY_RATE);

        // --- Repellent field (AC-driven, for chemotaxis away from) ---
        final BSimChemicalField repellentField = new BSimChemicalField(sim,
                new int[]{GRID_X, GRID_Y, GRID_Z}, DIFFUSIVITY, DECAY_RATE);

        // --- AHL field (bacterium-produced, for QS) ---
        final BSimChemicalField ahlField = new BSimChemicalField(sim,
                new int[]{GRID_X, GRID_Y, GRID_Z}, AHL_DIFFUSIVITY, AHL_DECAY_RATE);

        // --- Dose field (for AC activation gating) ---
        final BSimChemicalField doseField = new BSimChemicalField(sim,
                new int[]{GRID_X, GRID_Y, GRID_Z}, 0, 0);

        final double voxelXval = BOUND_X / GRID_X;
        final double voxelYval = BOUND_Y / GRID_Y;

        // --- Seed bacteria ---
        for (int i = 0; i < INITIAL_POP; i++) {
            double bx = 300 + Math.random() * 400;
            double by = 150 + Math.random() * 200;
            double bz = BOUND_Z / 2.0;
            ReservoirBacterium b = new ReservoirBacterium(sim,
                    new Vector3d(bx, by, bz), 0.0, signalField, repellentField, ahlField);
            b.setRadius();
            b.setSurfaceAreaGrowthRate(GROWTH_RATE);
            b.setChildList(children);
            bacteria.add(b);
        }

        // --- Ticker ---
        sim.setTicker(new BSimTicker() {
            @Override
            public void tick() {
                double t = sim.getTime();

                // ---- AC activation (Stage 4) ----
                for (int i = 0; i < GRID_X; i++)
                    for (int j = 0; j < GRID_Y; j++)
                        doseField.setConc(i, j, 0, 0);

                for (AC ac : acs) {
                    if (t < ac.totalTime && ac.shouldProduce(t)) {
                        int cxMin = Math.max(0, (int)((ac.position.x - DOSE_RADIUS) / voxelXval));
                        int cxMax = Math.min(GRID_X - 1, (int)((ac.position.x + DOSE_RADIUS) / voxelXval));
                        int cyMin = Math.max(0, (int)((ac.position.y - DOSE_RADIUS) / voxelYval));
                        int cyMax = Math.min(GRID_Y - 1, (int)((ac.position.y + DOSE_RADIUS) / voxelYval));
                        for (int i = cxMin; i <= cxMax; i++)
                            for (int j = cyMin; j <= cyMax; j++)
                                doseField.setConc(i, j, 0, DOSE_HIGH);
                    }
                }

                for (AC ac : acs) {
                    double dose = doseField.getConc(ac.position);
                    ac.activated = dose > DOSE_THRESHOLD;
                    if (ac.activated && !ac.inert) {
                        // Route signal to attractant or repellent field based on AC type
                        BSimChemicalField target = ac.repellent ? repellentField : signalField;
                        target.addQuantity(ac.position, PROD_RATE * sim.getDt());
                    }
                }

                // ---- Fields: diffuse + decay + advect ----
                signalField.update();
                advect(signalField, FLOW_SPEED, sim.getDt(), GRID_X, GRID_Y, GRID_Z);

                repellentField.update();
                advect(repellentField, FLOW_SPEED, sim.getDt(), GRID_X, GRID_Y, GRID_Z);

                ahlField.update();
                advect(ahlField, FLOW_SPEED, sim.getDt(), GRID_X, GRID_Y, GRID_Z);

                // ---- Bacteria: action + position update ----
                for (ReservoirBacterium b : bacteria) {
                    b.action();
                    b.updatePosition();
                }

                // Remove dead / out-of-bounds
                for (ReservoirBacterium b : bacteria) {
                    if (b.getPosition().x > BOUND_X || b.getPosition().x < 0) {
                        if (!removals.contains(b)) removals.add(b);
                    }
                }
                bacteria.removeAll(removals);
                removals.clear();

                bacteria.addAll(children);
                children.clear();
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

                // Attractant field (cyan)
                draw(signalField, Color.CYAN, (float)(255.0 / 540.0));

                // Repellent field (magenta)
                draw(repellentField, new Color(255, 0, 255), (float)(255.0 / 540.0));

                // ACs: green=attractant active, red=repellent active, gray=inactive
                for (AC ac : acs) {
                    Color col;
                    if (!ac.activated) col = Color.DARK_GRAY;
                    else if (ac.repellent) col = Color.RED;
                    else col = Color.GREEN;
                    sphere(ac.position, 15.0, col, 255);
                }

                // Bacteria: brightness tracks LuxI (QS reporter)
                for (ReservoirBacterium b : bacteria) {
                    double r = b.reporterNormalised();
                    Color col = new Color(30, (int)(55 + 200 * r), 30);
                    sphere(b.getPosition(), 8.0, col, 255);
                }

                // Flow direction arrow (dimmed since flow=0)
                p3d.stroke(100, 100, 0);
                p3d.strokeWeight(1);
                p3d.line(50, 30, 0, 200, 30, 0);
                p3d.line(200, 30, 0, 180, 20, 0);
                p3d.line(200, 30, 0, 180, 40, 0);
                p3d.noStroke();
            }
        });

        // --- Exporters ---
        if (exportData) {
            BSimUtils.generateDirectoryPath(exportPath);

            // State logger (every 1s)
            BSimLogger stateLogger = new BSimLogger(sim, exportPath + "stage10_control_state.csv") {
                @Override
                public void before() {
                    super.before();
                    StringBuilder header = new StringBuilder(
                            "time_s,population,mean_luxI,mean_AHL_intra,mean_AiiA,mean_LA");
                    for (AC ac : acs) {
                        header.append(",ac").append(ac.id).append("_activated");
                        header.append(",ac").append(ac.id).append("_signal_conc");
                    }
                    write(header.toString());
                }
                @Override
                public void during() {
                    double meanLuxI = 0, meanAHL = 0, meanAiiA = 0, meanLA = 0;
                    for (ReservoirBacterium b : bacteria) {
                        meanLuxI += b.y[0];
                        meanAHL  += b.y[1];
                        meanAiiA += b.y[2];
                        meanLA   += b.y[3];
                    }
                    int n = bacteria.size();
                    if (n > 0) { meanLuxI /= n; meanAHL /= n; meanAiiA /= n; meanLA /= n; }

                    StringBuilder line = new StringBuilder();
                    line.append(sim.getFormattedTime());
                    line.append(",").append(n);
                    line.append(",").append(String.format("%.6f", meanLuxI));
                    line.append(",").append(String.format("%.6f", meanAHL));
                    line.append(",").append(String.format("%.6f", meanAiiA));
                    line.append(",").append(String.format("%.6f", meanLA));

                    for (AC ac : acs) {
                        line.append(",").append(ac.activated ? 1 : 0);
                        // Report concentration from the field this AC produces into
                        BSimChemicalField f = ac.repellent ? repellentField : signalField;
                        line.append(",").append(String.format("%.2f", f.getConc(ac.position)));
                    }
                    write(line.toString());
                }
            };
            stateLogger.setDt(1.0);
            sim.addExporter(stateLogger);

            // Per-bacterium snapshot (every 30s)
            BSimLogger bacLogger = new BSimLogger(sim, exportPath + "stage10_control_bacteria.csv") {
                @Override
                public void before() {
                    super.before();
                    write("time_s,bac_id,x,y,z,luxI,AHL_intra,AiiA,LA,ext_AHL_uM,"
                            + "attractant_conc,repellent_conc,dist_nearest_ac");
                }
                @Override
                public void during() {
                    String t = sim.getFormattedTime();
                    for (ReservoirBacterium b : bacteria) {
                        Vector3d p = b.getPosition();
                        double minDist = Double.MAX_VALUE;
                        for (AC ac : acs) {
                            double dx = p.x - ac.position.x;
                            double dy = p.y - ac.position.y;
                            double d = Math.sqrt(dx * dx + dy * dy);
                            if (d < minDist) minDist = d;
                        }
                        write(t + "," + b.id
                                + "," + String.format("%.1f", p.x)
                                + "," + String.format("%.1f", p.y)
                                + "," + String.format("%.1f", p.z)
                                + "," + String.format("%.6f", b.y[0])
                                + "," + String.format("%.6f", b.y[1])
                                + "," + String.format("%.6f", b.y[2])
                                + "," + String.format("%.6f", b.y[3])
                                + "," + String.format("%.6e", b.extAHL_uM())
                                + "," + String.format("%.2f", signalField.getConc(p))
                                + "," + String.format("%.2f", repellentField.getConc(p))
                                + "," + String.format("%.1f", minDist));
                    }
                }
            };
            bacLogger.setDt(30.0);
            sim.addExporter(bacLogger);

            // Parameters
            BSimLogger paramLogger = new BSimLogger(sim, exportPath + "stage10_control_params.csv") {
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
                        write("ac" + ac.id + "_bit_duration_s," + ac.bitDuration);
                        write("ac" + ac.id + "_type," + (ac.repellent ? "REPELLENT" : "ATTRACTANT"));
                    }
                    write("flow_speed_um_per_s," + FLOW_SPEED);
                    write("signal_diffusivity," + DIFFUSIVITY);
                    write("signal_decay_rate," + DECAY_RATE);
                    write("ahl_diffusivity," + AHL_DIFFUSIVITY);
                    write("ahl_decay_rate," + AHL_DECAY_RATE);
                    write("prod_rate," + PROD_RATE);
                    write("initial_pop," + INITIAL_POP);
                    write("carrying_capacity," + CARRYING_CAPACITY);
                    write("growth_rate," + GROWTH_RATE);
                    write("expected_T_gen," + EXPECTED_T_GEN);
                    write("CellWallDiff," + CELL_WALL_DIFF);
                    write("dt," + sim.getDt());
                    write("sim_time," + sim.getSimulationTime());
                    write("chemotaxis,dual: setGoal(signalField) + repellentField sign-flip");
                    write("repellent_note,PLUMBING TEST -- sign-flip of attractant; same sensitivity/memory/pEndRun values");
                }
                @Override
                public void during() { }
            };
            sim.addExporter(paramLogger);

            System.out.println("Stage 10a CONTROL: exporting to" + exportPath);
            sim.export();

        } else {
            System.out.println("Stage 10a CONTROL: preview mode");
            sim.preview();
        }
    }
}
