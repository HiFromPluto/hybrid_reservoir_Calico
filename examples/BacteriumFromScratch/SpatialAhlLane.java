package BacteriumFromScratch;

import java.util.Arrays;

/**
 * Dilanji et al. 2012 one-dimensional AHL transport problem.
 *
 * <p>DOI 10.1021/ja211593q. Units: mm, h, nM. This is a standalone
 * millimetre-scale transport validation, not the 60 um Job 3c chamber.
 */
public final class SpatialAhlLane {

    public static final double LENGTH_MM = 32.0;          // TAKEN
    public static final double LOAD_LENGTH_MM = 2.0;     // TAKEN
    public static final double D_MM2_PER_H = 1.98;       // TAKEN
    public static final double C_INFINITY_NM = 4.0;      // TAKEN
    public static final double HALF_ACTIVATION_NM = 1.5; // TAKEN
    public static final double CFL = 0.40;               // ENGINEERING

    private final double dx;
    private final double baseDt;
    private final double[] c;
    private final double[] next;
    private final double initialMass;
    private double time;
    private double maxMassRelativeError;
    private double minimum = Double.POSITIVE_INFINITY;

    public SpatialAhlLane(double dxMm) {
        int n = (int) Math.round(LENGTH_MM / dxMm);
        if (n <= 2 || Math.abs(n * dxMm - LENGTH_MM) > 1e-12) {
            throw new IllegalArgumentException("dx must tile the 32 mm lane: " + dxMm);
        }
        this.dx = dxMm;
        this.baseDt = CFL * dx * dx / D_MM2_PER_H;
        this.c = new double[n];
        this.next = new double[n];
        double loaded = C_INFINITY_NM * LENGTH_MM / LOAD_LENGTH_MM;
        for (int i = 0; i < n; i++) {
            double x = (i + 0.5) * dx;
            c[i] = x < LOAD_LENGTH_MM ? loaded : 0.0;
        }
        this.initialMass = mass();
        updateDiagnostics();
    }

    public double dxMm() {
        return dx;
    }

    public double timeHours() {
        return time;
    }

    public double baseDtHours() {
        return baseDt;
    }

    public double initialMassNmMm() {
        return initialMass;
    }

    public double massNmMm() {
        return mass();
    }

    public double maxMassRelativeError() {
        return maxMassRelativeError;
    }

    public double minimumNm() {
        return minimum;
    }

    public double[] copyProfile() {
        return Arrays.copyOf(c, c.length);
    }

    public double cellCentreMm(int index) {
        return (index + 0.5) * dx;
    }

    public double sampleNearest(double xMm) {
        int i = Math.max(0, Math.min(c.length - 1, (int) Math.floor(xMm / dx)));
        return c[i];
    }

    /** Advance exactly to target time using a conservative finite-volume step. */
    public void advanceTo(double targetHours) {
        if (targetHours < time) {
            throw new IllegalArgumentException("cannot integrate backwards");
        }
        while (time < targetHours - 1e-14) {
            double dt = Math.min(baseDt, targetHours - time);
            step(dt);
        }
    }

    private void step(double dt) {
        double q = D_MM2_PER_H * dt / (dx * dx);
        // No-flux finite-volume boundaries: boundary face flux is exactly zero.
        next[0] = c[0] + q * (c[1] - c[0]);
        for (int i = 1; i < c.length - 1; i++) {
            next[i] = c[i] + q * (c[i - 1] - 2.0 * c[i] + c[i + 1]);
        }
        int last = c.length - 1;
        next[last] = c[last] + q * (c[last - 1] - c[last]);
        System.arraycopy(next, 0, c, 0, c.length);
        time += dt;
        updateDiagnostics();
    }

    private double mass() {
        double sum = 0.0;
        for (double value : c) {
            sum += value;
        }
        return sum * dx;
    }

    private void updateDiagnostics() {
        maxMassRelativeError = Math.max(maxMassRelativeError,
                Math.abs(mass() - initialMass) / initialMass);
        for (double value : c) {
            minimum = Math.min(minimum, value);
        }
    }
}
