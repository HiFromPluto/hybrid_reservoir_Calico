package BacteriumFromScratch;

import java.io.File;
import java.io.PrintStream;
import java.util.Locale;

/** Job 5b — standalone Dilanji millimetre-scale spatial AHL transport. */
public final class Job5bSims {

    static final String DIR = "results/job5b_spatial_ahl/";
    static final double T_PROFILE_H = 10.0;
    static final double T_END_H = 30.0;
    static final double X_GATE_MM = 10.0;

    private Job5bSims() {}

    private static final class Result {
        double dx;
        double dt;
        double arrival;
        double massError;
        double minimum;
        double finalMean;
        double[] profile10;
    }

    private static PrintStream stream(String name) {
        try {
            return new PrintStream(DIR + name, "UTF-8");
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private static Result run(double dx) {
        SpatialAhlLane lane = new SpatialAhlLane(dx);
        Result result = new Result();
        result.dx = dx;
        result.dt = lane.baseDtHours();
        result.arrival = Double.NaN;
        double previousT = lane.timeHours();
        double previousC = lane.sampleNearest(X_GATE_MM);
        boolean profileWritten = false;

        while (lane.timeHours() < T_END_H - 1e-13) {
            double target = Math.min(T_END_H,
                    lane.timeHours() + lane.baseDtHours());
            if (!profileWritten && target > T_PROFILE_H) {
                target = T_PROFILE_H;
            }
            lane.advanceTo(target);
            double currentC = lane.sampleNearest(X_GATE_MM);
            if (Double.isNaN(result.arrival)
                    && previousC < SpatialAhlLane.HALF_ACTIVATION_NM
                    && currentC >= SpatialAhlLane.HALF_ACTIVATION_NM) {
                double fraction = (SpatialAhlLane.HALF_ACTIVATION_NM - previousC)
                        / (currentC - previousC);
                result.arrival = previousT
                        + fraction * (lane.timeHours() - previousT);
            }
            if (!profileWritten
                    && Math.abs(lane.timeHours() - T_PROFILE_H) < 1e-10) {
                result.profile10 = lane.copyProfile();
                profileWritten = true;
            }
            previousT = lane.timeHours();
            previousC = currentC;
        }

        result.massError = lane.maxMassRelativeError();
        result.minimum = lane.minimumNm();
        result.finalMean = lane.massNmMm() / SpatialAhlLane.LENGTH_MM;

        PrintStream profile = stream(String.format(Locale.US,
                "profile_dx%03d_t10h.csv", (int) Math.round(dx * 1000.0)));
        profile.println("x_mm;C_nM");
        for (int i = 0; i < result.profile10.length; i++) {
            profile.printf(Locale.US, "%.6f;%.12e%n",
                    (i + 0.5) * dx, result.profile10[i]);
        }
        profile.close();
        return result;
    }

    public static void main(String[] args) {
        new File(DIR).mkdirs();
        System.out.println("Job 5b spatial AHL transport (standalone)");
        System.out.println("  Dilanji lane: L=32 mm, load=2 mm, D=1.98 mm^2/h");
        System.out.println("  TAKEN C_inf=4 nM, half-activation=1.5 nM");
        System.out.println("  not the 60 um Job 3c dish; no Weber parameter changes");

        double[] dxs = {0.20, 0.10, 0.05};
        Result[] results = new Result[dxs.length];
        for (int i = 0; i < dxs.length; i++) {
            results[i] = run(dxs[i]);
        }

        PrintStream summary = stream("job5b_summary.csv");
        summary.println("dx_mm;dt_h;arrival_x10_h;mass_rel_err;min_nM;final_mean_nM");
        for (Result result : results) {
            summary.printf(Locale.US, "%.4f;%.12e;%.12f;%.6e;%.6e;%.12f%n",
                    result.dx, result.dt, result.arrival, result.massError,
                    result.minimum, result.finalMean);
            System.out.printf(Locale.US,
                    "  dx=%.2f mm arrival@10mm=%.4f h mass_err=%.2e min=%.2e%n",
                    result.dx, result.arrival, result.massError, result.minimum);
        }
        summary.close();

        PrintStream scales = stream("geometry_scales.csv");
        double dUm2S = SpatialAhlLane.D_MM2_PER_H * 1.0e6 / 3600.0;
        double t10h = 10.0 * 10.0
                / (4.0 * SpatialAhlLane.D_MM2_PER_H);
        double t60s = 0.060 * 0.060
                / (4.0 * SpatialAhlLane.D_MM2_PER_H) * 3600.0;
        scales.println("D_mm2_h;D_um2_s;t_diff_10mm_h;t_diff_60um_s");
        scales.printf(Locale.US, "%.12f;%.12f;%.12f;%.12f%n",
                SpatialAhlLane.D_MM2_PER_H, dUm2S, t10h, t60s);
        scales.close();
        System.out.println("wrote " + DIR);
    }
}
