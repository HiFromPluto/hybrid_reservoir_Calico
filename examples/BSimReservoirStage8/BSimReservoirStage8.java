package BSimReservoirStage8;

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
 * Stage 8: Replace threshold QS with continuous AHL-LuxR kinetics.
 *
 * Changes from Stage 7:
 *   1. NEW ahlField -- bacterium-produced AHL (D=159 um^2/s, decay=2.76e-3/60 s^-1).
 *      This un-conflates the two roles of signalField flagged in Stage 7:
 *        - signalField remains AC-driven, used ONLY for chemotaxis (setGoal)
 *        - ahlField is bacterium-sourced, used for quorum sensing / gene expression
 *
 *   2. DELETED the generic 1-variable reporter ODE (alpha*I^n/(K^n+I^n) + alpha0 - gamma*y)
 *      carried since Stage 1. REPLACED with 4-variable LuxI/AHL/AiiA/LuxR-AHL ODE system
 *      from BSimEntrainment_PIDCtrl (Danino et al. 2010 parameterization):
 *        y[0] = LuxI   (autoinducer synthase, QS reporter)
 *        y[1] = AHL    (intracellular autoinducer)
 *        y[2] = AiiA   (autoinducer degradation enzyme)
 *        y[3] = LA     (LuxR-AHL complex, transcription factor)
 *      QS response: Kpli*(LA^n/(Kmla^n + LA^n)) -- continuous Hill, no threshold.
 *
 *   3. Membrane AHL exchange: intracellular AHL (y[1]) exchanges with extracellular
 *      ahlField via first-order kinetics (CellWallDiff = 1/20 s^-1,
 *      Kaplan & Greenberg 1985). Exchange uses cell volume for mass conservation.
 *
 * NOT changed:
 *   - signalField (AC-driven, molecule-count units, for chemotaxis)
 *   - Chemotaxis: setGoal(signalField) (Stage 7, Berg-cited)
 *   - Growth rate: 0.00698 um^2/s (30 min doubling, Stage 7)
 *   - Death: density-dependent stochastic (ungrounded, flagged for Stage 9)
 *   - sensitivity: 1 molecule/um^3 (~1.7 nM) -- only used by chemotaxis on
 *     signalField; ahlField uses uM units internally, no cross-comparison
 *   - Flow, ACs, advection (Stages 4-5)
 *
 * All QS parameters from BSimEntrainment_PIDCtrl (fixed, not random):
 *   Source: Danino et al. 2010, Nature 463:326-330
 *   Membrane exchange: Kaplan & Greenberg 1985, J. Bacteriol. 163(3):1210-1214
 *   (Roadmap note: verify these against the actual papers before trusting.)
 */
public class BSimReservoirStage8 {

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
    //  Signal field parameters (Stage 3 -- AC-driven, for chemotaxis)
    // =================================================================
    static final double DIFFUSIVITY = 100.0;   // um^2/s
    static final double DECAY_RATE  = 0.01;    // 1/s  (tau = 100s)
    static final double PROD_RATE   = 1e6;     // molecules/s per active AC

    // =================================================================
    //  AHL field parameters (NEW in Stage 8 -- bacterium-produced)
    //  Units: field quantity stored so that getConc()*1e-15 = uM
    //  (same convention as BSimEntrainment_PIDCtrl)
    // =================================================================
    static final double AHL_DIFFUSIVITY = 159.0;        // um^2/s (manuscript figure, in-code)
    static final double AHL_DECAY_RATE  = 2.76e-3 / 60; // s^-1 (t_half ~ 12.6 h, in-code)

    // =================================================================
    //  QS kinetic parameters (from BSimEntrainment_PIDCtrl, fixed)
    //  All rates converted from min^-1 to s^-1 by dividing by 60.
    //  Source: Danino et al. 2010 parameterization.
    // =================================================================
    static final double TIME_ADJ = 60.0;  // min -> s conversion

    static final double QS_DELTA1   = 0.8487   / TIME_ADJ;  // LuxI degradation (s^-1)
    static final double QS_DELTA2   = 0.0234   / TIME_ADJ;  // AiiA degradation (s^-1)
    static final double QS_G        = 0.0412;                // enzymatic saturation (dimensionless)
    static final double QS_KP2      = 9.0      / TIME_ADJ;  // AHL production by LuxI (s^-1)
    static final double QS_KR1OFF   = 6e-6     / TIME_ADJ;  // LuxR-AHL dissociation (s^-1)
    static final double QS_KR1ON    = 5.99e-5  / TIME_ADJ;  // LuxR-AHL association (uM^-1 s^-1)
    static final double QS_KCAT_AIIA = 2631.4  / TIME_ADJ;  // AiiA catalytic rate (s^-1)
    static final double QS_T_A      = 0.00276  / TIME_ADJ;  // AHL intracellular decay (s^-1)
    static final double QS_T_LA     = 0.024    / TIME_ADJ;  // LA complex decay (s^-1)
    static final double QS_A0LI     = 7.785e-6 / TIME_ADJ;  // basal LuxI production (uM s^-1)
    static final double QS_A0AA     = 6.183e-6 / TIME_ADJ;  // basal AiiA production (uM s^-1)
    static final double QS_KPLI     = 0.9      / TIME_ADJ;  // max LuxI production rate (uM s^-1)
    static final double QS_KPAA     = 0.9      / TIME_ADJ;  // max AiiA production rate (uM s^-1)
    static final double QS_KMLA     = 1e-2;                  // Hill constant for LA (uM)
    static final double QS_KMAA     = 1200.0;                // Michaelis constant for AiiA (uM)
    static final double QS_LTOT     = 15.0;                  // total LuxR (uM)
    static final double QS_N        = 2.0;                   // Hill coefficient

    // Cell-membrane AHL exchange rate (Kaplan & Greenberg 1985)
    // "conc. of AHL inside a cell and outside a cell [equilibrated] by 20 sec"
    static final double CELL_WALL_DIFF = 3.0 / TIME_ADJ;    // = 1/20 = 0.05 s^-1

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
    //  Flow parameters (Stage 5)
    // =================================================================
    static final double FLOW_SPEED = 0.0;     // um/s along +x

    // =================================================================
    //  Growth / death parameters (Stage 7)
    // =================================================================
    static final double GROWTH_RATE       = 4.0 * Math.PI / 1800.0;  // 0.00698 um^2/s, 30 min doubling
    static final double EXPECTED_T_GEN    = 4.0 * Math.PI / GROWTH_RATE;
    static final double DEATH_RATE_AT_CAP = 1.0 / EXPECTED_T_GEN;

    // =================================================================
    //  Population parameters — derived from Danino et al. 2010
    // =================================================================
    // Danino et al. (2010, Nature 463:326-330) operate their LuxI/LuxR
    // oscillator in micro-traps at confluent density (~10^10 cells/mL).
    // Canonical QS onset for LuxI/LuxR systems is 10^7-10^9 cells/mL
    // (Fuqua et al. 1994, J. Bacteriol. 176:269-275).
    //
    // Derivation for this domain:
    //   Domain volume = 1000 * 500 * 10 um^3 = 5e6 um^3 = 5e-6 mL
    //   At 10^8 cells/mL (conservative QS onset): 10^8 * 5e-6 = 500 cells
    //   At 10^9 cells/mL (strong QS):             10^9 * 5e-6 = 5000 cells
    //
    // INITIAL_POP = 500 targets the 10^8 cells/mL regime under full
    // validated Stage 5 flow (10 um/s). If advection-limited washout
    // prevents QS activation at this density, that is a real physical
    // finding about open-channel geometry, not a parameter to tune.
    //
    // Note: Danino et al.'s device physically traps cells (near-zero
    // flow inside trap, ~100x50x1 um chamber) so AHL accumulates
    // locally. Our straight flow-through channel may not sustain the
    // same AHL concentrations even at equivalent cell count.
    static final int    INITIAL_POP       = 500;
    static final int    CARRYING_CAPACITY = 2000;  // allow growth; 10^9 cells/mL ceiling

    // =================================================================
    //  AC class (from Stage 4)
    // =================================================================
    static class AC {
        final int id;
        final Vector3d position;
        boolean activated = false;
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
    //  Bacterium with QS gene expression + growth + division + death
    // =================================================================
    static int nextId = 0;
    static final Vector<ReservoirBacterium> bacteria = new Vector<>();
    static final Vector<ReservoirBacterium> children = new Vector<>();
    static final Vector<ReservoirBacterium> removals = new Vector<>();

    static class ReservoirBacterium extends BSimBacterium {
        final int id;
        double lastDivisionTime;
        // y[0]=LuxI, y[1]=AHL_intra, y[2]=AiiA, y[3]=LA (all in uM)
        double[] y;
        BSimOdeSystem grn;
        final BSimChemicalField signalField;
        final BSimChemicalField ahlField;
        // Intra-extra AHL concentration difference (field units/um^3)
        // Computed before ODE solve, read inside derivativeSystem
        double diffConc;

        public ReservoirBacterium(BSim sim, Vector3d position, double birthTime,
                                  BSimChemicalField signalField,
                                  BSimChemicalField ahlField) {
            super(sim, position);
            this.id = nextId++;
            this.lastDivisionTime = birthTime;
            this.y = new double[]{0.05, 0.05, 0.05, 0.05};  // ICs from BSimEntrainment_PIDCtrl
            this.signalField = signalField;
            this.ahlField = ahlField;
            this.grn = new QSGRN();
            // Stage 7: chemotaxis toward AC-driven signal field
            setGoal(signalField);
        }

        @Override
        public void action() {
            // Full BSimBacterium physics: Brownian + run/tumble + growth
            super.action();

            // Apply flow drag force
            addForce(new Vector3d(stokesCoefficient() * FLOW_SPEED, 0, 0));

            // Compute intra-extra AHL concentration difference BEFORE ODE solve.
            // diffConc = y[1]*1e15 - getConc(position), in field concentration units.
            // The ODE reads this via diffConc*1e-15 to get the difference in uM.
            // (Same pattern as BSimEntrainment_PIDCtrl)
            diffConc = y[1] * 1e15 - ahlField.getConc(position);

            // Solve QS ODE
            y = BSimOdeSolver.rungeKutta45(grn, sim.getTime(), y, sim.getDt());

            // Clamp negative values (numerical safety)
            for (int i = 0; i < y.length; i++)
                if (y[i] < 0) y[i] = 0;

            // Transfer AHL between cell and field.
            // Amount = diffConc * CellWallDiff * dt * V_cell
            // Uses cell volume (not voxel volume) for correct mass conservation
            // in our sparse-cell setup (see derivation in Stage 8 javadoc).
            double cellVol = (4.0 / 3.0) * Math.PI * Math.pow(radius, 3);
            double amountToField = diffConc * CELL_WALL_DIFF * sim.getDt() * cellVol;
            ahlField.addQuantity(position, amountToField);

            // Density-dependent stochastic death (ungrounded, Stage 9)
            double pDeath = sim.getDt() * DEATH_RATE_AT_CAP
                    * ((double) bacteria.size() / CARRYING_CAPACITY);
            if (Math.random() < pDeath) {
                removals.add(this);
            }
        }

        /**
         * Confine Y (and Z) to the domain. X is deliberately left alone here --
         * it's handled separately in the tick loop as an absorbing/removal
         * boundary (Stage 5's flow-through channel design), not touched by this
         * override.
         *
         * Checked before this existed: nothing was constraining Y at all --
         * real run data showed bacteria at y=971.5 with BOUND_Y=500, i.e.
         * unbounded drift, not a reflection bug -- there was no boundary
         * handling on this axis whatsoever.
         *
         * Mirror-reflection, not clamp -- ported from the original project
         * (reservoir_new/.../ReservoirBacterium_1.java, updatePosition()):
         * a clamp lets a bacterium whose run direction points outward sit
         * pinned at the wall for the duration of that run, producing an
         * artificial standing population at the domain edges. That was a
         * real, found-and-fixed bug there (STRATEGY_NEW_PAPER.md sec 2 round 3),
         * ported here rather than reinvented.
         *
         * Z is a plain clamp, not a reflection -- same reasoning as the
         * original project: BOUND_Z=10 with GRID_Z=1 means z never changes
         * voxel assignment, so reflecting it would be a no-op difference
         * from clamping. Added because z was equally unconstrained before
         * this, not because it's expected to matter.
         */
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
                    childPos, divisionTime, signalField, ahlField);
            child.setRadius(radius);
            child.setSurfaceAreaGrowthRate(surfaceAreaGrowthRate);
            child.setChildList(childList);
            // Inherit QS state (symmetric division)
            child.y = new double[]{ this.y[0], this.y[1], this.y[2], this.y[3] };

            childList.add(child);
            this.lastDivisionTime = divisionTime;
        }

        /** LuxI level as QS reporter, normalized 0-1 for visualization. */
        public double reporterNormalised() {
            // Scale by approximate max LuxI (Kpli/delta1 ~ 0.9/0.8487 ~ 1.06 uM)
            double luxIMax = QS_KPLI / (QS_DELTA1);
            return Math.min(1.0, Math.max(0.0, y[0] / luxIMax));
        }

        /** Extracellular AHL concentration at this bacterium's position, in uM. */
        public double extAHL_uM() {
            return ahlField.getConc(position) * 1e-15;
        }

        /**
         * 4-variable QS ODE system: LuxI / AHL / AiiA / LuxR-AHL.
         * From BSimEntrainment_PIDCtrl (Danino et al. 2010 parameterization).
         * All variables in uM, all rates in s^-1.
         * QS response: continuous Hill function of LA (y[3]), no threshold.
         */
        class QSGRN implements BSimOdeSystem {
            @Override
            public double[] derivativeSystem(double t, double[] y) {
                // Convert intra-extra difference to uM for ODE
                double extraintradiff_uM = diffConc * 1e-15;

                double[] dy = new double[4];

                // y[0] = LuxI: basal + Hill(LA) - enzymatic degradation
                dy[0] = QS_A0LI
                        + QS_KPLI * (Math.pow(y[3], QS_N) / (Math.pow(QS_KMLA, QS_N) + Math.pow(y[3], QS_N)))
                        - (QS_DELTA1 * y[0]) / (QS_G * (y[0] + y[2]) + 1);

                // y[1] = AHL_intra: production by LuxI - LuxR binding + LuxR-AHL dissociation
                //        - AiiA degradation - spontaneous decay - membrane exchange
                dy[1] = QS_KP2 * y[0]
                        - QS_KR1ON * (QS_LTOT - y[3]) * y[1]
                        + QS_KR1OFF * y[3]
                        - (QS_KCAT_AIIA * y[2] * y[1]) / (QS_KMAA + y[1])
                        - QS_T_A * y[1]
                        - CELL_WALL_DIFF * extraintradiff_uM;

                // y[2] = AiiA: basal + Hill(LA) - enzymatic degradation
                dy[2] = QS_A0AA
                        + QS_KPAA * (Math.pow(y[3], QS_N) / (Math.pow(QS_KMLA, QS_N) + Math.pow(y[3], QS_N)))
                        - (QS_DELTA2 * y[2]) / (QS_G * (y[0] + y[2]) + 1);

                // y[3] = LA (LuxR-AHL complex): association - dissociation - decay
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
            AC ac = new AC(i, acPositions[i], seq, bitDur);
            acs.add(ac);
            if (ac.totalTime > maxACTime) maxACTime = ac.totalTime;
            System.out.print("AC " + i + " (" + bitDur + "s/bit, "
                    + seq.length + " bits, " + ac.totalTime + "s): ");
            for (int b : seq) System.out.print(b + " ");
            System.out.println();
        }

        // Sim time: enough for QS dynamics to develop + growth validation
        final double simTime = Math.max(maxACTime + 30.0, 3600.0);

        String timestamp = BSimUtils.timeStamp();
        String exportPath = "./results/" + timestamp + "/";

        double voxelX = BOUND_X / GRID_X;
        double courant = FLOW_SPEED * 0.05 / voxelX;
        double courantAHL = FLOW_SPEED * 0.05 / voxelX;
        double kxAHL = AHL_DIFFUSIVITY * 0.05 / (voxelX * voxelX);

        System.out.println("Stage 8: Continuous QS kinetics (LuxI/AHL/AiiA/LA)");
        System.out.println("  Signal field: D=" + DIFFUSIVITY + ", decay=" + DECAY_RATE);
        System.out.println("  AHL field: D=" + AHL_DIFFUSIVITY + ", decay=" + AHL_DECAY_RATE);
        System.out.println("    AHL kX=" + String.format("%.3f", kxAHL) + " (must be < 0.5)");
        System.out.println("    AHL CFL=" + String.format("%.3f", courantAHL) + " (must be < 1)");
        System.out.println("  Flow: " + FLOW_SPEED + " um/s, CFL=" + courant);
        System.out.println("  Growth rate: " + String.format("%.5f", GROWTH_RATE) + " um^2/s");
        System.out.println("  Initial pop: " + INITIAL_POP + ", carrying capacity: " + CARRYING_CAPACITY);
        System.out.println("  Sim time: " + simTime + "s");

        // --- Simulation ---
        BSim sim = new BSim();
        sim.setDt(0.05);
        sim.setSimulationTime(simTime);
        sim.setTimeFormat("0.00");
        sim.setBound(BOUND_X, BOUND_Y, BOUND_Z);
        // BSim's own default is solid={false,false,false} on all axes, which for
        // BSimChemicalField's diffuse() means PERIODIC WRAP (checked directly in
        // BSimChemicalField.java: a non-solid axis maps the last box's "above"
        // neighbour to box 0, and vice versa). Left at the default, signalField
        // and ahlField were wrapping concentration from x=999 back to x=0 (and
        // y=499 back to y=0) on every diffusion step -- independent of, and in
        // addition to, the one-way advect() below. Setting all three solid
        // stops that: X stays effectively open in practice because advect()'s
        // own upwind scheme already carries mass out the right edge correctly
        // on its own (checked: cUpwind=0 at i=0, no special case needed at the
        // outflow end) -- solid=true here only removes the diffusion step's
        // extra, unwanted wraparound leak, it doesn't undo the flow-through
        // design. Z doesn't matter either way (GRID_Z=1), set for consistency.
        sim.setSolid(true, true, true);

        // --- Signal field (AC-driven, for chemotaxis) ---
        final BSimChemicalField signalField = new BSimChemicalField(sim,
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
                    new Vector3d(bx, by, bz), 0.0, signalField, ahlField);
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
                    if (ac.activated) {
                        signalField.addQuantity(ac.position, PROD_RATE * sim.getDt());
                    }
                }

                // ---- Fields: diffuse + decay + advect ----
                signalField.update();
                advect(signalField, FLOW_SPEED, sim.getDt(), GRID_X, GRID_Y, GRID_Z);

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

                // Signal field (cyan) -- AC-driven chemotaxis target
                draw(signalField, Color.CYAN, (float)(255.0 / 540.0));

                // ACs: green if activated, red if not
                for (AC ac : acs) {
                    Color col = ac.activated ? Color.GREEN : Color.RED;
                    sphere(ac.position, 15.0, col, 255);
                }

                // Bacteria: brightness tracks LuxI (QS reporter)
                for (ReservoirBacterium b : bacteria) {
                    double r = b.reporterNormalised();
                    Color col = new Color(30, (int)(55 + 200 * r), 30);
                    sphere(b.getPosition(), 8.0, col, 255);
                }

                // Flow direction arrow
                p3d.stroke(255, 255, 0);
                p3d.strokeWeight(2);
                p3d.line(50, 30, 0, 200, 30, 0);
                p3d.line(200, 30, 0, 180, 20, 0);
                p3d.line(200, 30, 0, 180, 40, 0);
                p3d.noStroke();
            }
        });

        // --- Exporters ---
        if (exportData) {
            BSimUtils.generateDirectoryPath(exportPath);

            // State logger (every 1s): population, mean QS state, AC status
            BSimLogger stateLogger = new BSimLogger(sim, exportPath + "stage8_state.csv") {
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
                        line.append(",").append(String.format("%.2f",
                                signalField.getConc(ac.position)));
                    }
                    write(line.toString());
                }
            };
            stateLogger.setDt(1.0);
            sim.addExporter(stateLogger);

            // Per-bacterium snapshot (every 30s)
            BSimLogger bacLogger = new BSimLogger(sim, exportPath + "stage8_bacteria.csv") {
                @Override
                public void before() {
                    super.before();
                    write("time_s,bac_id,x,y,z,luxI,AHL_intra,AiiA,LA,ext_AHL_uM,dist_nearest_ac");
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
                                + "," + String.format("%.1f", minDist));
                    }
                }
            };
            bacLogger.setDt(30.0);
            sim.addExporter(bacLogger);

            // Parameters
            BSimLogger paramLogger = new BSimLogger(sim, exportPath + "stage8_params.csv") {
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
                    write("QS_delta1," + QS_DELTA1);
                    write("QS_delta2," + QS_DELTA2);
                    write("QS_g," + QS_G);
                    write("QS_Kp2," + QS_KP2);
                    write("QS_Kr1ON," + QS_KR1ON);
                    write("QS_Kr1OFF," + QS_KR1OFF);
                    write("QS_KcatAiiA," + QS_KCAT_AIIA);
                    write("QS_t_A," + QS_T_A);
                    write("QS_t_LA," + QS_T_LA);
                    write("QS_a0LI," + QS_A0LI);
                    write("QS_a0AA," + QS_A0AA);
                    write("QS_Kpli," + QS_KPLI);
                    write("QS_KpaA," + QS_KPAA);
                    write("QS_Kmla," + QS_KMLA);
                    write("QS_KmaA," + QS_KMAA);
                    write("QS_Ltot," + QS_LTOT);
                    write("QS_n," + QS_N);
                    write("CellWallDiff," + CELL_WALL_DIFF);
                    write("dt," + sim.getDt());
                    write("sim_time," + sim.getSimulationTime());
                    write("chemotaxis,setGoal(signalField)");
                    write("sensitivity,1 (uncited; only used by chemotaxis on signalField)");
                }
                @Override
                public void during() { }
            };
            sim.addExporter(paramLogger);

            System.out.println("Stage 8: exporting to " + exportPath);
            sim.export();

        } else {
            System.out.println("Stage 8: preview mode");
            sim.preview();
        }
    }
}
