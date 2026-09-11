package BSimReservoirPlanStage4;

import bsim.BSimChemicalField;

import javax.vecmath.Vector3d;
import java.util.Arrays;
import java.util.List;

/** I/O-only analyzer; readout grids remain independent of the chemical grid. */
public final class VoxelAnalyzer {
    public static final class VoxelReadout {
        public final double[] attractant, repellent, ahl, acidMm, pH, meanResponse;
        public final double[] fractionQAboveHalf;
        public final int[] density, births, deaths;

        VoxelReadout(double[] attractant, double[] repellent, double[] ahl,
                     double[] acidMm, double[] pH, int[] density,
                     double[] meanResponse, double[] fractionQAboveHalf,
                     int[] births, int[] deaths) {
            this.attractant = attractant;
            this.repellent = repellent;
            this.ahl = ahl;
            this.acidMm = acidMm;
            this.pH = pH;
            this.density = density;
            this.meanResponse = meanResponse;
            this.fractionQAboveHalf = fractionQAboveHalf;
            this.births = births;
            this.deaths = deaths;
        }
    }

    private final int[] stateGrid, countGrid;
    private final double boundX, boundY, boundZ;
    private final int[] birthAccum, deathAccum;

    public VoxelAnalyzer(int[] stateGrid, int[] countGrid,
                         double boundX, double boundY, double boundZ) {
        this.stateGrid = stateGrid.clone();
        this.countGrid = countGrid.clone();
        this.boundX = boundX;
        this.boundY = boundY;
        this.boundZ = boundZ;
        birthAccum = new int[getCountBins()];
        deathAccum = new int[getCountBins()];
    }

    public int getStateVoxels() {
        return stateGrid[0] * stateGrid[1] * stateGrid[2];
    }

    public int getCountBins() {
        return countGrid[0] * countGrid[1] * countGrid[2];
    }

    public void resetWindowCounters() {
        Arrays.fill(birthAccum, 0);
        Arrays.fill(deathAccum, 0);
    }

    public void recordBirth(Vector3d p) { birthAccum[countBinFor(p)]++; }
    public void recordDeath(Vector3d p) { deathAccum[countBinFor(p)]++; }

    public VoxelReadout analyze(BSimChemicalField attractantField,
                                BSimChemicalField repellentField,
                                BSimChemicalField ahlField,
                                BSimChemicalField acidField,
                                List<BSimReservoirPlanStage4.ReservoirBacterium> bacteria) {
        int total = getStateVoxels();
        double[] attractant = new double[total], repellent = new double[total];
        double[] ahl = new double[total], acidMm = new double[total];
        double[] pH = new double[total], responseSum = new double[total];
        double[] qAbove = new double[total];
        int[] density = new int[total];
        double voxelX = boundX / stateGrid[0];
        double voxelY = boundY / stateGrid[1];
        double voxelZ = boundZ / stateGrid[2];

        for (int x = 0; x < stateGrid[0]; x++) {
            for (int y = 0; y < stateGrid[1]; y++) {
                for (int z = 0; z < stateGrid[2]; z++) {
                    int id = index(x, y, z, stateGrid);
                    Vector3d centre = new Vector3d((x + .5) * voxelX,
                            (y + .5) * voxelY, (z + .5) * voxelZ);
                    attractant[id] = attractantField.getConc(centre);
                    repellent[id] = repellentField.getConc(centre);
                    ahl[id] = ahlField.getConc(centre);
                    acidMm[id] = BSimReservoirPlanStage4.acidMmFromConc(
                            acidField.getConc(centre));
                    pH[id] = BSimReservoirPlanStage4.pHFromAcidMm(acidMm[id]);
                }
            }
        }

        for (BSimReservoirPlanStage4.ReservoirBacterium bacterium : bacteria) {
            Vector3d p = bacterium.getPosition();
            int id = index(clamp((int) (p.x / voxelX), 0, stateGrid[0] - 1),
                    clamp((int) (p.y / voxelY), 0, stateGrid[1] - 1),
                    clamp((int) (p.z / voxelZ), 0, stateGrid[2] - 1), stateGrid);
            density[id]++;
            responseSum[id] += bacterium.getResponse();
            if (bacterium.getQ() > .5) qAbove[id]++;
        }

        double[] meanResponse = new double[total];
        double[] fractionQAboveHalf = new double[total];
        for (int i = 0; i < total; i++) {
            if (density[i] > 0) {
                meanResponse[i] = responseSum[i] / density[i];
                fractionQAboveHalf[i] = qAbove[i] / density[i];
            }
        }
        return new VoxelReadout(attractant, repellent, ahl, acidMm, pH, density,
                meanResponse, fractionQAboveHalf, birthAccum.clone(), deathAccum.clone());
    }

    private int countBinFor(Vector3d p) {
        return index(clamp((int) (p.x / (boundX / countGrid[0])), 0, countGrid[0] - 1),
                clamp((int) (p.y / (boundY / countGrid[1])), 0, countGrid[1] - 1),
                clamp((int) (p.z / (boundZ / countGrid[2])), 0, countGrid[2] - 1),
                countGrid);
    }

    private static int index(int x, int y, int z, int[] grid) {
        return x * grid[1] * grid[2] + y * grid[2] + z;
    }

    private static int clamp(int value, int min, int max) {
        return Math.max(min, Math.min(max, value));
    }
}
