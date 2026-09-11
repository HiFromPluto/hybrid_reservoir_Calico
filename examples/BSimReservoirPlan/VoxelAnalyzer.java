package BSimReservoirPlan;

import bsim.BSimChemicalField;

import javax.vecmath.Vector3d;
import java.util.Arrays;
import java.util.List;

/**
 * I/O-only dual-grid analyzer. The state and count grids are explicit and
 * independent of the chemical-field grid.
 */
public final class VoxelAnalyzer {
    public static final class VoxelReadout {
        public final double[] attractant;
        public final double[] repellent;
        public final double[] ahl;
        public final int[] density;
        public final double[] meanLuxI;
        public final int[] births;
        public final int[] deaths;

        VoxelReadout(double[] attractant, double[] repellent, double[] ahl,
                     int[] density, double[] meanLuxI, int[] births, int[] deaths) {
            this.attractant = attractant;
            this.repellent = repellent;
            this.ahl = ahl;
            this.density = density;
            this.meanLuxI = meanLuxI;
            this.births = births;
            this.deaths = deaths;
        }
    }

    private final int[] stateGrid;
    private final int[] countGrid;
    private final double boundX;
    private final double boundY;
    private final double boundZ;
    private final int[] birthAccum;
    private final int[] deathAccum;

    public VoxelAnalyzer(int[] stateGrid, int[] countGrid,
                         double boundX, double boundY, double boundZ) {
        this.stateGrid = stateGrid.clone();
        this.countGrid = countGrid.clone();
        this.boundX = boundX;
        this.boundY = boundY;
        this.boundZ = boundZ;
        this.birthAccum = new int[getCountBins()];
        this.deathAccum = new int[getCountBins()];
    }

    public int getStateVoxels() {
        return stateGrid[0] * stateGrid[1] * stateGrid[2];
    }

    public int getCountBins() {
        return countGrid[0] * countGrid[1] * countGrid[2];
    }

    public int[] getStateGrid() {
        return stateGrid.clone();
    }

    public int[] getCountGrid() {
        return countGrid.clone();
    }

    public void resetWindowCounters() {
        Arrays.fill(birthAccum, 0);
        Arrays.fill(deathAccum, 0);
    }

    public int[] snapshotBirths() {
        return birthAccum.clone();
    }

    public int[] snapshotDeaths() {
        return deathAccum.clone();
    }

    public void recordBirth(Vector3d position) {
        birthAccum[countBinFor(position)]++;
    }

    public void recordDeath(Vector3d position) {
        deathAccum[countBinFor(position)]++;
    }

    public VoxelReadout analyze(BSimChemicalField attractantField,
                                BSimChemicalField repellentField,
                                BSimChemicalField ahlField,
                                List<BSimReservoirPlan.ReservoirBacterium> bacteria) {
        int total = getStateVoxels();
        double[] attractant = new double[total];
        double[] repellent = new double[total];
        double[] ahl = new double[total];
        int[] density = new int[total];
        double[] luxISum = new double[total];

        double voxelX = boundX / stateGrid[0];
        double voxelY = boundY / stateGrid[1];
        double voxelZ = boundZ / stateGrid[2];

        for (int x = 0; x < stateGrid[0]; x++) {
            for (int y = 0; y < stateGrid[1]; y++) {
                for (int z = 0; z < stateGrid[2]; z++) {
                    int id = index(x, y, z, stateGrid);
                    Vector3d centre = new Vector3d(
                            (x + 0.5) * voxelX,
                            (y + 0.5) * voxelY,
                            (z + 0.5) * voxelZ);
                    attractant[id] = attractantField.getConc(centre);
                    repellent[id] = repellentField.getConc(centre);
                    ahl[id] = ahlField.getConc(centre);
                }
            }
        }

        for (BSimReservoirPlan.ReservoirBacterium bacterium : bacteria) {
            Vector3d p = bacterium.getPosition();
            int x = clamp((int) (p.x / voxelX), 0, stateGrid[0] - 1);
            int y = clamp((int) (p.y / voxelY), 0, stateGrid[1] - 1);
            int z = clamp((int) (p.z / voxelZ), 0, stateGrid[2] - 1);
            int id = index(x, y, z, stateGrid);
            density[id]++;
            luxISum[id] += bacterium.getLuxI();
        }

        double[] meanLuxI = new double[total];
        for (int i = 0; i < total; i++) {
            if (density[i] > 0) meanLuxI[i] = luxISum[i] / density[i];
        }

        return new VoxelReadout(attractant, repellent, ahl, density, meanLuxI,
                birthAccum.clone(), deathAccum.clone());
    }

    private int countBinFor(Vector3d p) {
        int x = clamp((int) (p.x / (boundX / countGrid[0])), 0, countGrid[0] - 1);
        int y = clamp((int) (p.y / (boundY / countGrid[1])), 0, countGrid[1] - 1);
        int z = clamp((int) (p.z / (boundZ / countGrid[2])), 0, countGrid[2] - 1);
        return index(x, y, z, countGrid);
    }

    private static int index(int x, int y, int z, int[] grid) {
        return x * grid[1] * grid[2] + y * grid[2] + z;
    }

    private static int clamp(int value, int min, int max) {
        return Math.max(min, Math.min(max, value));
    }
}
