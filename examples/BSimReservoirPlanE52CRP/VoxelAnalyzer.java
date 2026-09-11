package BSimReservoirPlanE52CRP;

import bsim.BSimChemicalField;
import bsim.particle.BSimParticle;

import javax.vecmath.Vector3d;
import java.util.Arrays;
import java.util.List;

/** I/O-only analyzer; readout grids remain independent of the chemical grid. */
public final class VoxelAnalyzer {
    public static final class VoxelReadout {
        public final double[] ahlUm, pH, meanResponse, fractionQAboveHalf;
        public final double[] meanLum, sumLum;
        public final int[] density, births, clampDeaths, acidDeaths, oobDeaths;

        VoxelReadout(double[] ahlUm, double[] pH, int[] density,
                     double[] meanResponse, double[] fractionQAboveHalf,
                     double[] meanLum, double[] sumLum, int[] births,
                     int[] clampDeaths, int[] acidDeaths, int[] oobDeaths) {
            this.ahlUm = ahlUm;
            this.pH = pH;
            this.density = density;
            this.meanResponse = meanResponse;
            this.fractionQAboveHalf = fractionQAboveHalf;
            this.meanLum = meanLum;
            this.sumLum = sumLum;
            this.births = births;
            this.clampDeaths = clampDeaths;
            this.acidDeaths = acidDeaths;
            this.oobDeaths = oobDeaths;
        }
    }

    private final int[] stateGrid, countGrid;
    private final double boundX, boundY, boundZ;
    private final int[] birthAccum, clampDeathAccum, acidDeathAccum, oobDeathAccum;

    public VoxelAnalyzer(int[] stateGrid, int[] countGrid,
                         double boundX, double boundY, double boundZ) {
        this.stateGrid = stateGrid.clone();
        this.countGrid = countGrid.clone();
        this.boundX = boundX;
        this.boundY = boundY;
        this.boundZ = boundZ;
        birthAccum = new int[getCountBins()];
        clampDeathAccum = new int[getCountBins()];
        acidDeathAccum = new int[getCountBins()];
        oobDeathAccum = new int[getCountBins()];
    }

    public int getStateVoxels() {
        return stateGrid[0] * stateGrid[1] * stateGrid[2];
    }

    public int getCountBins() {
        return countGrid[0] * countGrid[1] * countGrid[2];
    }

    public void resetWindowCounters() {
        Arrays.fill(birthAccum, 0);
        Arrays.fill(clampDeathAccum, 0);
        Arrays.fill(acidDeathAccum, 0);
        Arrays.fill(oobDeathAccum, 0);
    }

    public void recordBirth(Vector3d p) { birthAccum[countBinFor(p)]++; }

    public void recordDeath(Vector3d p, BSimReservoirPlanE52CRP.DeathCause cause) {
        int bin = countBinFor(p);
        if (cause == BSimReservoirPlanE52CRP.DeathCause.CLAMP) clampDeathAccum[bin]++;
        else if (cause == BSimReservoirPlanE52CRP.DeathCause.ACID) acidDeathAccum[bin]++;
        else if (cause == BSimReservoirPlanE52CRP.DeathCause.OOB) oobDeathAccum[bin]++;
    }

    public VoxelReadout analyze(BSimChemicalField ahlField,
                                BSimChemicalField acidField,
                                List<BSimReservoirPlanE52CRP.ReservoirBacterium> bacteria) {
        int total = getStateVoxels();
        double[] ahlUm = new double[total], pH = new double[total];
        double[] responseSum = new double[total], lumSum = new double[total];
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
                    ahlUm[id] = ahlField.getConc(centre)
                            / BSimReservoirPlanE52CRP.MOLECULES_PER_UM3_PER_UM;
                    pH[id] = BSimReservoirPlanE52CRP.pHFromAcidMm(
                            BSimReservoirPlanE52CRP.acidMmFromConc(acidField.getConc(centre)));
                }
            }
        }

        for (BSimReservoirPlanE52CRP.ReservoirBacterium bacterium : bacteria) {
            Vector3d p = bacterium.getPosition();
            int id = index(clamp((int) (p.x / voxelX), 0, stateGrid[0] - 1),
                    clamp((int) (p.y / voxelY), 0, stateGrid[1] - 1),
                    clamp((int) (p.z / voxelZ), 0, stateGrid[2] - 1), stateGrid);
            density[id]++;
            responseSum[id] += bacterium.getResponse();
            lumSum[id] += bacterium.getLuminescence();
            if (bacterium.getQ() > .5) qAbove[id]++;
        }

        double[] meanResponse = new double[total];
        double[] fractionQAboveHalf = new double[total];
        double[] meanLum = new double[total];
        for (int i = 0; i < total; i++) {
            if (density[i] > 0) {
                meanResponse[i] = responseSum[i] / density[i];
                fractionQAboveHalf[i] = qAbove[i] / density[i];
                meanLum[i] = lumSum[i] / density[i];
            }
        }
        return new VoxelReadout(ahlUm, pH, density, meanResponse, fractionQAboveHalf,
                meanLum, lumSum, birthAccum.clone(), clampDeathAccum.clone(),
                acidDeathAccum.clone(), oobDeathAccum.clone());
    }

    /** Spatial density only. Used by the Brownian null; no chemical or cell-state features. */
    public int[] density(List<? extends BSimParticle> particles) {
        int total = getStateVoxels();
        int[] density = new int[total];
        double voxelX = boundX / stateGrid[0];
        double voxelY = boundY / stateGrid[1];
        double voxelZ = boundZ / stateGrid[2];
        for (BSimParticle particle : particles) {
            Vector3d p = particle.getPosition();
            int id = index(clamp((int) (p.x / voxelX), 0, stateGrid[0] - 1),
                    clamp((int) (p.y / voxelY), 0, stateGrid[1] - 1),
                    clamp((int) (p.z / voxelZ), 0, stateGrid[2] - 1), stateGrid);
            density[id]++;
        }
        return density;
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
