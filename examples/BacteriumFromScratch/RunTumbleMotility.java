package BacteriumFromScratch;

import java.util.Random;

/**
 * Unbiased 3-D E. coli run-and-tumble renewal process.
 *
 * <p>Quantitative parameters: Kurzthaler et al. 2024,
 * DOI 10.1103/PhysRevLett.132.038302. Units are um and s.
 */
public final class RunTumbleMotility {

    public static final double MEAN_SPEED_UM_S = 16.0; // TAKEN
    public static final double SPEED_SD_UM_S = 5.78;   // TAKEN
    public static final double MEAN_RUN_S = 2.39;      // TAKEN
    public static final double MEAN_TUMBLE_S = 0.38;   // TAKEN

    private RunTumbleMotility() {}

    public static double exponential(Random rng, double mean) {
        return -mean * Math.log(1.0 - rng.nextDouble());
    }

    /** Sample the measured Gaussian speed distribution, rejecting v<=0. */
    public static double speed(Random rng) {
        double value;
        do {
            value = MEAN_SPEED_UM_S + SPEED_SD_UM_S * rng.nextGaussian();
        } while (value <= 0.0);
        return value;
    }

    /** Uniform direction on the unit sphere. */
    public static void direction(Random rng, double[] out) {
        double z = 2.0 * rng.nextDouble() - 1.0;
        double phi = 2.0 * Math.PI * rng.nextDouble();
        double radial = Math.sqrt(1.0 - z * z);
        out[0] = radial * Math.cos(phi);
        out[1] = radial * Math.sin(phi);
        out[2] = z;
    }

    /** Published renewal-process long-time active diffusivity, um^2/s. */
    public static double theoreticalDiffusivity() {
        double v2 = MEAN_SPEED_UM_S * MEAN_SPEED_UM_S
                + SPEED_SD_UM_S * SPEED_SD_UM_S;
        return v2 * MEAN_RUN_S * MEAN_RUN_S
                / (3.0 * (MEAN_RUN_S + MEAN_TUMBLE_S));
    }

    public static double runFraction() {
        return MEAN_RUN_S / (MEAN_RUN_S + MEAN_TUMBLE_S);
    }
}
