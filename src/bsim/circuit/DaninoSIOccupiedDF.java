package bsim.circuit;

import bsim.BSimRandom;

/**
 * Object B: {@code DaninoSI_OccupiedDF}. <b>NOT_FIG4B</b>.
 *
 * Occupied SI delay-DDE at the TAKEN table, {@code d=0.5}, bulk {@code D1=0},
 * SI Hill {@code P}, delay only in {@code Hi}. This is the numerical organism
 * for D1/C0/C1. It is not a twin of Danino 2010 Fig. 4b.
 *
 * Object A {@code Danino2010_Fig4b_bulk_twin} remains FAIL / FAIL_NO_IDENTITY
 * and is not this class.
 *
 * Injected {@link BSimRandom} is accepted for the N0 contract and is unused:
 * the circuit is deterministic.
 */
public final class DaninoSIOccupiedDF {

    public static final String OBJECT = "DaninoSI_OccupiedDF";
    public static final String OBJECT_STATUS = "OCCUPIED_NOT_FIG4B";

    public static final double C_A = 1.0;
    public static final double C_I = 4.0;
    public static final double DELTA = 0.001;
    public static final double ALPHA = 2500.0;
    public static final double TAU = 10.0;
    public static final double K = 1.0;
    public static final double K1 = 0.1;
    public static final double B = 0.06;
    public static final double GAMMA_A = 15.0;
    public static final double GAMMA_I = 24.0;
    public static final double GAMMA_H = 0.01;
    public static final double F = 0.3;
    public static final double G = 0.01;
    public static final double D0 = 0.88;
    public static final double D_MEM = 2.5;
    public static final double D_BULK = 0.5;
    public static final double D1_BULK = 0.0;

    private final BSimRandom unusedRng;

    public DaninoSIOccupiedDF() {
        this(null);
    }

    /** RNG is injected and never called. */
    public DaninoSIOccupiedDF(BSimRandom rng) {
        this.unusedRng = rng;
    }

    public BSimRandom unusedRng() {
        return unusedRng;
    }

    /** SI Hill: {@code P = (delta + alpha H^2) / (1 + k1 H^2)}. */
    public static double production(double hTau) {
        double h2 = hTau * hTau;
        return (DELTA + ALPHA * h2) / (1.0 + K1 * h2);
    }

    public static void rhs(double[] y, double hTau, double d, double mu, double[] out) {
        double a = y[0];
        double i = y[1];
        double hI = y[2];
        double hE = y[3];
        double p = production(hTau);
        double dens = 1.0 - Math.pow(d / D0, 4.0);
        double denomF = 1.0 + F * (a + i);
        out[0] = C_A * dens * p - GAMMA_A * a / denomF;
        out[1] = C_I * dens * p - GAMMA_I * i / denomF;
        out[2] = B * i / (1.0 + K * i) - (GAMMA_H * a * hI) / (1.0 + G * a) + D_MEM * (hE - hI);
        out[3] = -(d / (1.0 - d)) * D_MEM * (hE - hI) - mu * hE;
    }

    /**
     * Integrate with RK4 and a linear {@code Hi} delay tape.
     * {@code y0} is {@code (A, I, Hi, He)}. History for {@code t < 0} is {@code hiHistory}.
     */
    public Trajectory integrate(
            double d,
            double mu,
            double[] y0,
            double hiHistory,
            double tEnd,
            double sampleDt,
            double rkDt) {
        if (y0.length != 4) {
            throw new IllegalArgumentException("y0 must be (A, I, Hi, He)");
        }
        int nMicro = (int) Math.round(tEnd / rkDt);
        if (Math.abs(nMicro * rkDt - tEnd) > 1e-12) {
            throw new IllegalArgumentException("rkDt must divide tEnd");
        }
        int nSample = (int) Math.round(tEnd / sampleDt);
        if (Math.abs(nSample * sampleDt - tEnd) > 1e-12) {
            throw new IllegalArgumentException("sampleDt must divide tEnd");
        }
        int stride = (int) Math.round(sampleDt / rkDt);
        if (Math.abs(stride * rkDt - sampleDt) > 1e-12) {
            throw new IllegalArgumentException("rkDt must divide sampleDt");
        }

        double[] y = y0.clone();
        double[] k1 = new double[4];
        double[] k2 = new double[4];
        double[] k3 = new double[4];
        double[] k4 = new double[4];
        double[] yTmp = new double[4];
        double[] hiTape = new double[nMicro + 1];
        hiTape[0] = y[2];

        double[] tOut = new double[nSample + 1];
        double[][] yOut = new double[nSample + 1][4];
        tOut[0] = 0.0;
        System.arraycopy(y, 0, yOut[0], 0, 4);

        for (int n = 0; n < nMicro; n++) {
            double t = n * rkDt;
            rhs(y, hiAt(hiTape, n, rkDt, t - TAU, hiHistory, y0[2]), d, mu, k1);
            axpy(y, k1, 0.5 * rkDt, yTmp);
            rhs(yTmp, hiAt(hiTape, n, rkDt, t + 0.5 * rkDt - TAU, hiHistory, y0[2]), d, mu, k2);
            axpy(y, k2, 0.5 * rkDt, yTmp);
            rhs(yTmp, hiAt(hiTape, n, rkDt, t + 0.5 * rkDt - TAU, hiHistory, y0[2]), d, mu, k3);
            axpy(y, k3, rkDt, yTmp);
            rhs(yTmp, hiAt(hiTape, n, rkDt, t + rkDt - TAU, hiHistory, y0[2]), d, mu, k4);
            for (int j = 0; j < 4; j++) {
                y[j] += (rkDt / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
            }
            hiTape[n + 1] = y[2];
            if ((n + 1) % stride == 0) {
                int s = (n + 1) / stride;
                tOut[s] = (n + 1) * rkDt;
                System.arraycopy(y, 0, yOut[s], 0, 4);
            }
        }

        boolean negative = false;
        for (double[] row : yOut) {
            for (double v : row) {
                if (v < -1e-12) {
                    negative = true;
                    break;
                }
            }
        }
        return new Trajectory(tOut, yOut, d, mu, y0.clone(), hiHistory, rkDt, negative);
    }

    private static void axpy(double[] y, double[] k, double h, double[] out) {
        for (int j = 0; j < 4; j++) {
            out[j] = y[j] + h * k[j];
        }
    }

    private static double hiAt(
            double[] tape,
            int lastFilled,
            double rkDt,
            double tDelay,
            double hiHistory,
            double hiAtZero) {
        if (tDelay < 0.0) {
            return hiHistory;
        }
        if (tDelay <= 1e-15) {
            return hiAtZero;
        }
        double idx = tDelay / rkDt;
        int i = (int) Math.floor(idx);
        if (i < 0) {
            return hiHistory;
        }
        if (i >= lastFilled) {
            i = lastFilled;
            return tape[i];
        }
        double frac = idx - i;
        int i1 = i + 1;
        if (i1 > lastFilled) {
            return tape[i];
        }
        return tape[i] * (1.0 - frac) + tape[i1] * frac;
    }

    public static final class Trajectory {
        public final double[] t;
        public final double[][] Y;
        public final double d;
        public final double mu;
        public final double[] y0;
        public final double hiHistory;
        public final double rkDt;
        public final boolean negativeState;

        Trajectory(
                double[] t,
                double[][] Y,
                double d,
                double mu,
                double[] y0,
                double hiHistory,
                double rkDt,
                boolean negativeState) {
            this.t = t;
            this.Y = Y;
            this.d = d;
            this.mu = mu;
            this.y0 = y0;
            this.hiHistory = hiHistory;
            this.rkDt = rkDt;
            this.negativeState = negativeState;
        }

        public double[] luxI() {
            double[] i = new double[Y.length];
            for (int n = 0; n < Y.length; n++) {
                i[n] = Y[n][1];
            }
            return i;
        }
    }
}
