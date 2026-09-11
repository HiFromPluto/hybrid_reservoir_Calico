package BSimReservoirStage11;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.ode.BSimOdeSolver;
import bsim.ode.BSimOdeSystem;
import bsim.particle.BSimBacterium;

import javax.vecmath.Vector3d;
import java.util.Vector;

/**
 * Merged reservoir bacterium: Stage10's Danino QSGRN + dual chemotaxis
 * + homeostatic clamp + Hill toxicity, plus Monod growth coupling from
 * reservoir_new.
 *
 * Key differences from Stage10:
 *   - Unit conversion: 602 molecules/um^3 per uM (not 1e15)
 *   - Monod growth: effectiveGrowthRate = GROWTH_RATE * G/(K_G + G)
 *   - Glucose uptake (Michaelis-Menten)
 *   - VoxelAnalyzer birth/death event tracking
 */
public class ReservoirBacterium extends BSimBacterium {

    // ── QS kinetic parameters (Danino et al. 2010, via BSimEntrainment_PIDCtrl) ──
    static final double TIME_ADJ = 60.0;

    static final double QS_DELTA1    = 0.8487   / TIME_ADJ;
    static final double QS_DELTA2    = 0.0234   / TIME_ADJ;
    static final double QS_G         = 0.0412;
    static final double QS_KP2       = 9.0      / TIME_ADJ;
    static final double QS_KR1OFF    = 6e-6     / TIME_ADJ;
    static final double QS_KR1ON     = 5.99e-5  / TIME_ADJ;
    static final double QS_KCAT_AIIA = 2631.4   / TIME_ADJ;
    static final double QS_T_A       = 0.00276  / TIME_ADJ;
    static final double QS_T_LA      = 0.024    / TIME_ADJ;
    static final double QS_A0LI      = 7.785e-6 / TIME_ADJ;
    static final double QS_A0AA      = 6.183e-6 / TIME_ADJ;
    static final double QS_KPLI      = 0.9      / TIME_ADJ;
    static final double QS_KPAA      = 0.9      / TIME_ADJ;
    static final double QS_KMLA      = 1e-2;
    static final double QS_KMAA      = 1200.0;
    static final double QS_LTOT      = 15.0;
    static final double QS_N         = 2.0;

    // Cell-wall AHL diffusion rate (molecules/um^3 -> uM: * 1/602)
    static final double CELL_WALL_DIFF = 3.0 / TIME_ADJ;

    // Unit conversion: 1 uM = 602 molecules/um^3
    static final double MOL_PER_UM3_PER_UM = 602.0;

    // ── Growth parameters ──
    static double GROWTH_RATE = 4.0 * Math.PI / 1800.0; // um^2/s surface area growth

    // Monod growth half-saturation.
    // Senn et al. 1994 K_s ≈ 0.18 µM (glucose, E. coli) — use k.g=0.18 to decouple.
    // Default 5.0 matches KM_UPTAKE (supply-capped; mu proportional to v). Do not imply Senn covers 5.0.
    static double K_G = 5.0;

    // Glucose uptake half-saturation. TUNED (inherited from reservoir_new); not Senn 1994.
    static double V_MAX_UPTAKE = 1.0;
    static double KM_UPTAKE    = 5.0;

    // ── Homeostatic clamp (BSimLacOperon) — modeling convenience, NOT primary pop control ──
    // Observed plateau (~950) is set by Monod glucose limitation + uptake/replenishment,
    // not by this clamp (N << CARRYING_CAPACITY at equilibrium). Kept as soft upper bound only.
    static double T_REMOVAL         = 2250.0;
    static int    CARRYING_CAPACITY = 2000;

    // ── Toxicity from repellent (Hill dose-response, placeholder) ──
    static double TOX_EC50  = 150.0;
    static double TOX_K_MAX = 0.001;
    static double TOX_N     = 2.0;

    // ── Flow ──
    static double FLOW_SPEED = 0.0;

    // ── State ──
    private boolean markedForDeath = false;
    double[] y; // QS GRN: [LuxI, AHL_intra, AiiA, LuxR-AHL]
    private final BSimOdeSystem grn;
    double diffConc; // AHL intra-extra gradient (molecules/um^3)

    // Repellent gradient memory (Berg 1972, dual chemotaxis from Stage10)
    private double[] repMemory;

    // Dependencies
    final BSimChemicalField signalField;
    final BSimChemicalField repellentField;
    final BSimChemicalField ahlField;
    final BSimChemicalField glucoseField;
    final VoxelAnalyzer voxelAnalyzer;
    final Vector<ReservoirBacterium> childrenList;

    // Static reference for homeostatic clamp population count
    static Vector<ReservoirBacterium> allBacteria;

    public ReservoirBacterium(BSim sim, Vector3d position,
                              BSimChemicalField signalField,
                              BSimChemicalField repellentField,
                              BSimChemicalField ahlField,
                              BSimChemicalField glucoseField,
                              VoxelAnalyzer voxelAnalyzer,
                              Vector<ReservoirBacterium> childrenList) {
        super(sim, position);
        this.signalField    = signalField;
        this.repellentField = repellentField;
        this.ahlField       = ahlField;
        this.glucoseField   = glucoseField;
        this.voxelAnalyzer  = voxelAnalyzer;
        this.childrenList   = childrenList;
        this.y = new double[]{0.05, 0.05, 0.05, 0.05};
        this.grn = new QSGRN();

        // Attractant chemotaxis (setGoal)
        setGoal(signalField);

        // Repellent memory
        int memLen = sim.timesteps(shortTermMemoryDuration + longTermMemoryDuration);
        repMemory = new double[memLen];
        double initConc = repellentField.getConc(position);
        for (int i = 0; i < repMemory.length; i++) repMemory[i] = initConc;
    }

    public boolean isMarkedForDeath() { return markedForDeath; }
    public double[] getY() { return y; }

    public double reporterNormalised() {
        double luxIMax = QS_KPLI / QS_DELTA1;
        return Math.min(1.0, Math.max(0.0, y[0] / luxIMax));
    }

    // ── Repellent gradient (Stage10 dual chemotaxis) ──

    private double repellentGradientDelta() {
        double shortTermCounter = 0, longTermCounter = 0;
        System.arraycopy(repMemory, 0, repMemory, 1, repMemory.length - 1);
        repMemory[0] = repellentField.getConc(position);
        for (int i = 0; i < repMemory.length; i++) {
            if (i < shortTermMemoryLength)
                shortTermCounter += repMemory[i];
            else
                longTermCounter += repMemory[i];
        }
        return (shortTermCounter / shortTermMemoryLength)
             - (longTermCounter / longTermMemoryLength);
    }

    @Override
    public double pEndRun() {
        boolean upAttractant = (goal != null) && movingUpGradient();
        double repDelta = repellentGradientDelta();
        boolean upRepellent   = repDelta > sensitivity;
        boolean downRepellent = repDelta < -sensitivity;

        if (upAttractant && !upRepellent)    return pEndRunUp;
        if (upRepellent && !upAttractant)    return pEndRunElse;
        if (downRepellent && !upAttractant)  return pEndRunUp;
        return pEndRunElse;
    }

    @Override
    public void action() {
        super.action(); // run-and-tumble

        // Flow force
        addForce(new Vector3d(stokesCoefficient() * FLOW_SPEED, 0, 0));

        // AHL intra-extra gradient: y[1] is in uM, field is molecules/um^3
        // Correct conversion: y[1] * 602 gives molecules/um^3
        diffConc = y[1] * MOL_PER_UM3_PER_UM - ahlField.getConc(position);

        // QS GRN ODE step
        y = BSimOdeSolver.rungeKutta45(grn, sim.getTime(), y, sim.getDt());
        for (int i = 0; i < y.length; i++)
            if (y[i] < 0) y[i] = 0;

        // AHL field exchange
        double cellVol = (4.0 / 3.0) * Math.PI * Math.pow(radius, 3);
        double amountToField = diffConc * CELL_WALL_DIFF * sim.getDt() * cellVol;
        ahlField.addQuantity(position, amountToField);

        // Glucose uptake (Michaelis-Menten)
        double localGlucose = glucoseField.getConc(position);
        double desiredUptake = V_MAX_UPTAKE * localGlucose / (KM_UPTAKE + localGlucose);
        double actualUptake  = Math.min(localGlucose, desiredUptake * sim.getDt()) / sim.getDt();
        actualUptake = Math.max(0.0, actualUptake);
        glucoseField.addQuantity(position, -actualUptake * sim.getDt());

        // Homeostatic clamp (BSimLacOperon): P_removal = (dt/T) * 2^(-(1-N/K))
        double pRemoval = (sim.getDt() / T_REMOVAL)
                * Math.pow(2.0, -(1.0 - (double) allBacteria.size() / CARRYING_CAPACITY));
        if (Math.random() < pRemoval) {
            markedForDeath = true;
            return;
        }

        // Hill toxicity from repellent
        double repConc = repellentField.getConc(position);
        if (repConc > 0) {
            double kKill = TOX_K_MAX * Math.pow(repConc, TOX_N)
                    / (Math.pow(TOX_EC50, TOX_N) + Math.pow(repConc, TOX_N));
            if (Math.random() < kKill * sim.getDt()) {
                markedForDeath = true;
            }
        }
    }

    @Override
    public void updatePosition() {
        super.updatePosition();

        // Reflective Y boundaries (Stage10)
        double boundY = sim.getBound().y;
        if (position.y < 0.0)         position.y = -position.y;
        else if (position.y > boundY) position.y = 2.0 * boundY - position.y;
        position.y = Math.max(0.0, Math.min(boundY, position.y));
        position.z = Math.max(0.0, Math.min(sim.getBound().z, position.z));
    }

    @Override
    public void grow() {
        if (markedForDeath) return;

        double localGlucose = glucoseField.getConc(position);

        // Monod only: chemotaxis changes position, not biomass synthesis rate.
        // Input-correlated births require local glucose depletion shadows (see GLU_DIFF).
        double monodFactor = localGlucose / (K_G + localGlucose);

        double savedRate = surfaceAreaGrowthRate;
        surfaceAreaGrowthRate = savedRate * monodFactor;
        super.grow();
        surfaceAreaGrowthRate = savedRate;
    }

    @SuppressWarnings("unchecked")
    @Override
    public void replicate() {
        setRadiusFromSurfaceArea(surfaceArea(replicationRadius) / 2);

        Vector3d childPos = new Vector3d(position);
        childPos.x += 2.0 * radius * (Math.random() - 0.5);
        childPos.y += 2.0 * radius * (Math.random() - 0.5);

        ReservoirBacterium child = new ReservoirBacterium(sim,
                childPos, signalField, repellentField, ahlField, glucoseField,
                voxelAnalyzer, childrenList);
        child.setRadius(radius);
        child.setSurfaceAreaGrowthRate(GROWTH_RATE); // base rate, not Monod-modified
        child.setChildList(childList);
        child.y = new double[]{this.y[0], this.y[1], this.y[2], this.y[3]};

        childrenList.add(child);

        // Record birth event in voxel analyzer
        if (voxelAnalyzer != null) {
            voxelAnalyzer.recordBirth(childPos);
        }
    }

    // ── Danino QSGRN ODE (4-variable) ──
    class QSGRN implements BSimOdeSystem {
        @Override
        public double[] derivativeSystem(double t, double[] y) {
            // Convert AHL gradient to uM for the ODE
            double extraintradiff_uM = diffConc / MOL_PER_UM3_PER_UM;

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
