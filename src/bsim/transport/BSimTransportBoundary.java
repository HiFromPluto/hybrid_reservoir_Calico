package bsim.transport;

import bsim.BSim;

/**
 * Chemical boundary conditions in BSim's documented face order:
 * x-high, x-low, y-high, y-low, z-high, z-low.
 */
public final class BSimTransportBoundary {
    public static final int X_HIGH = 0;
    public static final int X_LOW = 1;
    public static final int Y_HIGH = 2;
    public static final int Y_LOW = 3;
    public static final int Z_HIGH = 4;
    public static final int Z_LOW = 5;

    public enum Type {
        NO_FLUX,
        PERIODIC,
        OPEN,
        LEAKY
    }

    private final Type[] type = new Type[6];
    private final double[] exteriorConcentration = new double[6];
    private final double[] leakRate = new double[6];

    public BSimTransportBoundary() {
        for (int face = 0; face < 6; face++) {
            type[face] = Type.NO_FLUX;
        }
    }

    public static BSimTransportBoundary fromSimulation(BSim sim) {
        BSimTransportBoundary boundary = new BSimTransportBoundary();
        boolean[] solid = sim.getSolid();
        boolean[] leaky = sim.getLeaky();
        double[] rates = sim.getLeakyRate();
        for (int axis = 0; axis < 3; axis++) {
            int high = 2 * axis;
            int low = high + 1;
            if (!solid[axis]) {
                boundary.setPeriodic(high);
                boundary.setPeriodic(low);
            } else {
                if (leaky[high]) boundary.setLeaky(high, rates[high]);
                if (leaky[low]) boundary.setLeaky(low, rates[low]);
            }
        }
        return boundary;
    }

    public BSimTransportBoundary setNoFlux(int face) {
        checkFace(face);
        type[face] = Type.NO_FLUX;
        exteriorConcentration[face] = 0.0;
        leakRate[face] = 0.0;
        return this;
    }

    public BSimTransportBoundary setPeriodic(int face) {
        checkFace(face);
        type[face] = Type.PERIODIC;
        exteriorConcentration[face] = 0.0;
        leakRate[face] = 0.0;
        return this;
    }

    public BSimTransportBoundary setOpen(int face, double concentration) {
        checkFace(face);
        requireNonnegativeFinite(concentration, "exterior concentration");
        type[face] = Type.OPEN;
        exteriorConcentration[face] = concentration;
        leakRate[face] = 0.0;
        return this;
    }

    /**
     * Set the legacy BSim Robin-like loss coefficient in micrometres squared
     * per second. The exterior concentration is zero.
     */
    public BSimTransportBoundary setLeaky(int face, double rate) {
        checkFace(face);
        requireNonnegativeFinite(rate, "leak rate");
        type[face] = Type.LEAKY;
        leakRate[face] = rate;
        exteriorConcentration[face] = 0.0;
        return this;
    }

    public Type getType(int face) {
        checkFace(face);
        return type[face];
    }

    public double getExteriorConcentration(int face) {
        checkFace(face);
        return exteriorConcentration[face];
    }

    public double getLeakRate(int face) {
        checkFace(face);
        return leakRate[face];
    }

    public void validatePeriodicPairs() {
        for (int axis = 0; axis < 3; axis++) {
            int high = 2 * axis;
            int low = high + 1;
            boolean highPeriodic = type[high] == Type.PERIODIC;
            boolean lowPeriodic = type[low] == Type.PERIODIC;
            if (highPeriodic != lowPeriodic) {
                throw new IllegalArgumentException(
                        "periodic boundaries must be paired on axis " + axis);
            }
        }
    }

    private static void checkFace(int face) {
        if (face < 0 || face >= 6) {
            throw new IllegalArgumentException("boundary face must be in [0,5]: " + face);
        }
    }

    private static void requireNonnegativeFinite(double value, String name) {
        if (!Double.isFinite(value) || value < 0.0) {
            throw new IllegalArgumentException(name + " must be finite and nonnegative: " + value);
        }
    }
}
