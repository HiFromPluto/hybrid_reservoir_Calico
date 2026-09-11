package BSimReservoirStage11;

import bsim.BSimChemicalField;

import javax.vecmath.Vector3d;
import java.util.List;

/**
 * Dual-grid voxel readout for the reservoir state vector.
 *
 * STATE grid (default 20x10x1 = 200 voxels): chemical concentrations + density + activation.
 * COUNT grid (default 4x2x1 = 8 bins): births and deaths per window (coarse bins for Poisson ≥ 5).
 *
 * Channels in voxels.csv:
 *   Att(200), Rep(200), Den(200), Act(200), Prot(200), AHL(200), Births(8), Deaths(8)
 */
public class VoxelAnalyzer {

    public static final class VoxelReadout {
        public final double[] attractant;
        public final double[] repellent;
        public final int[]    density;
        public final double[] activatedFraction;
        public final double[] protein;
        public final double[] ahl;
        public final int[]    births;
        public final int[]    deaths;

        public VoxelReadout(double[] attractant, double[] repellent,
                            int[] density, double[] activatedFraction,
                            double[] protein, double[] ahl,
                            int[] births, int[] deaths) {
            this.attractant        = attractant;
            this.repellent         = repellent;
            this.density           = density;
            this.activatedFraction = activatedFraction;
            this.protein           = protein;
            this.ahl               = ahl;
            this.births            = births;
            this.deaths            = deaths;
        }
    }

    private final int[] stateGrid;
    private final int[] countGrid;
    private final double boundX, boundY, boundZ;

    // Accumulated birth/death counts per COUNT bin, reset each window
    private int[] birthAccum;
    private int[] deathAccum;

    public VoxelAnalyzer(int[] stateGrid, int[] countGrid,
                         double boundX, double boundY, double boundZ) {
        this.stateGrid = stateGrid.clone();
        this.countGrid = countGrid.clone();
        this.boundX    = boundX;
        this.boundY    = boundY;
        this.boundZ    = boundZ;
        this.birthAccum = new int[countGrid[0] * countGrid[1] * countGrid[2]];
        this.deathAccum = new int[countGrid[0] * countGrid[1] * countGrid[2]];
    }

    public int[] snapshotBirths() { return birthAccum.clone(); }
    public int[] snapshotDeaths() { return deathAccum.clone(); }

    public int getStateVoxels() { return stateGrid[0] * stateGrid[1] * stateGrid[2]; }
    public int getCountBins()   { return countGrid[0] * countGrid[1] * countGrid[2]; }
    public int[] getStateGrid() { return stateGrid; }
    public int[] getCountGrid() { return countGrid; }

    /** Reset birth/death accumulators at the start of each window. */
    public void resetWindowCounters() {
        java.util.Arrays.fill(birthAccum, 0);
        java.util.Arrays.fill(deathAccum, 0);
    }

    /** Record a birth event at position. */
    public void recordBirth(Vector3d pos) {
        int bin = countBinFor(pos);
        if (bin >= 0) birthAccum[bin]++;
    }

    /** Record a death event at position. */
    public void recordDeath(Vector3d pos) {
        int bin = countBinFor(pos);
        if (bin >= 0) deathAccum[bin]++;
    }

    private int countBinFor(Vector3d pos) {
        double vx = boundX / countGrid[0];
        double vy = boundY / countGrid[1];
        double vz = boundZ / countGrid[2];
        int xi = clamp((int)(pos.x / vx), 0, countGrid[0] - 1);
        int yi = clamp((int)(pos.y / vy), 0, countGrid[1] - 1);
        int zi = clamp((int)(pos.z / vz), 0, countGrid[2] - 1);
        return xi * countGrid[1] * countGrid[2] + yi * countGrid[2] + zi;
    }

    /**
     * Compute all readout channels for the current state.
     * Bacterium type is ReservoirBacterium (Stage11's merged bacterium class).
     */
    public VoxelReadout analyze(BSimChemicalField attractantField,
                                BSimChemicalField repellentField,
                                BSimChemicalField ahlField,
                                List<ReservoirBacterium> bacteria) {

        int total = getStateVoxels();
        double[] attractant   = new double[total];
        double[] repellent    = new double[total];
        double[] ahlConc      = new double[total];
        int[]    density      = new int[total];
        double[] actSum       = new double[total];
        double[] protSum      = new double[total];

        double voxelX = boundX / stateGrid[0];
        double voxelY = boundY / stateGrid[1];
        double voxelZ = boundZ / stateGrid[2];

        for (int xi = 0; xi < stateGrid[0]; xi++)
            for (int yi = 0; yi < stateGrid[1]; yi++)
                for (int zi = 0; zi < stateGrid[2]; zi++) {
                    int id = xi * stateGrid[1] * stateGrid[2] + yi * stateGrid[2] + zi;
                    Vector3d centre = new Vector3d(
                            (xi + 0.5) * voxelX,
                            (yi + 0.5) * voxelY,
                            (zi + 0.5) * voxelZ);
                    attractant[id] = attractantField.getConc(centre);
                    repellent[id]  = repellentField.getConc(centre);
                    ahlConc[id]    = ahlField.getConc(centre);
                }

        for (ReservoirBacterium b : bacteria) {
            Vector3d pos = b.getPosition();
            int xi = clamp((int)(pos.x / voxelX), 0, stateGrid[0] - 1);
            int yi = clamp((int)(pos.y / voxelY), 0, stateGrid[1] - 1);
            int zi = clamp((int)(pos.z / voxelZ), 0, stateGrid[2] - 1);
            int id = xi * stateGrid[1] * stateGrid[2] + yi * stateGrid[2] + zi;
            density[id]++;
            actSum[id]  += b.reporterNormalised();
            protSum[id] += b.getY()[0]; // LuxI as protein readout
        }

        double[] actFrac = new double[total];
        double[] prot    = new double[total];
        for (int i = 0; i < total; i++) {
            if (density[i] > 0) {
                actFrac[i] = actSum[i] / density[i];
                prot[i]    = protSum[i] / density[i];
            }
        }

        return new VoxelReadout(attractant, repellent, density, actFrac, prot, ahlConc,
                                birthAccum.clone(), deathAccum.clone());
    }

    private static int clamp(int val, int min, int max) {
        return Math.max(min, Math.min(max, val));
    }
}
