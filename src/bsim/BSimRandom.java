package bsim;

import java.util.Random;

/**
 * Single injectable pseudo-random source for a simulation.
 * Instances are intentionally stateful and should be shared by all stochastic
 * components in one deterministic run.
 */
public final class BSimRandom {
    private final Random random;

    /** Legacy unseeded construction. New reproducible jobs should pass a seed. */
    public BSimRandom() {
        random = new Random();
    }

    public BSimRandom(long seed) {
        random = new Random(seed);
    }

    public double nextDouble() {
        return random.nextDouble();
    }

    public double nextGaussian() {
        return random.nextGaussian();
    }

    public int nextInt(int bound) {
        return random.nextInt(bound);
    }

    public long nextLong() {
        return random.nextLong();
    }

    /**
     * Shared stream for legacy {@code java.util.Random} fields (e.g. Scratch
     * {@code EcoliRodCell.divisionRng}). Draws consume this source; do not
     * reseed the returned object independently.
     */
    public Random asJavaRandom() {
        return random;
    }
}
