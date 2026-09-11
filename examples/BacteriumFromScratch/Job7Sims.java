package BacteriumFromScratch;

import java.io.File;
import java.io.PrintStream;
import java.util.Locale;
import java.util.Random;

/** Job 7 — seeded ensemble validation of unbiased E. coli motility. */
public final class Job7Sims {

    static final String DIR = "results/job7_motility/";
    static final int N_TRAJECTORIES = 20_000; // ENGINEERING
    static final long SEED = 101L;             // ENGINEERING
    static final double[] OBS_S = {0.1, 0.5, 1.0, 10.0, 100.0, 1000.0};

    private Job7Sims() {}

    private static final class Moments {
        long n;
        double sum;
        double sumSq;

        void add(double x) {
            n++;
            sum += x;
            sumSq += x * x;
        }

        double mean() {
            return sum / n;
        }

        double sd() {
            double variance = (sumSq - sum * sum / n) / (n - 1.0);
            return Math.sqrt(Math.max(0.0, variance));
        }
    }

    private static PrintStream stream(String name) {
        try {
            return new PrintStream(DIR + name, "UTF-8");
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    public static void main(String[] args) {
        new File(DIR).mkdirs();
        Random rng = new Random(SEED);

        double[] sumR2 = new double[OBS_S.length];
        double[] sumX = new double[OBS_S.length];
        double[] sumY = new double[OBS_S.length];
        double[] sumZ = new double[OBS_S.length];
        Moments runs = new Moments();
        Moments tumbles = new Moments();
        Moments speeds = new Moments();
        Moments runLengths = new Moments();
        double[] direction = new double[3];

        for (int particle = 0; particle < N_TRAJECTORIES; particle++) {
            double x = 0.0, y = 0.0, z = 0.0, time = 0.0;
            boolean running = true;
            double remaining = RunTumbleMotility.exponential(
                    rng, RunTumbleMotility.MEAN_RUN_S);
            double speed = RunTumbleMotility.speed(rng);
            RunTumbleMotility.direction(rng, direction);
            runs.add(remaining);
            speeds.add(speed);
            runLengths.add(remaining * speed);

            int observation = 0;
            while (observation < OBS_S.length) {
                double target = OBS_S[observation];
                double segment = Math.min(remaining, target - time);
                if (running) {
                    x += speed * direction[0] * segment;
                    y += speed * direction[1] * segment;
                    z += speed * direction[2] * segment;
                }
                time += segment;
                remaining -= segment;

                if (Math.abs(time - target) < 1e-10) {
                    double r2 = x * x + y * y + z * z;
                    sumR2[observation] += r2;
                    sumX[observation] += x;
                    sumY[observation] += y;
                    sumZ[observation] += z;
                    observation++;
                }

                if (remaining <= 1e-12) {
                    running = !running;
                    if (running) {
                        remaining = RunTumbleMotility.exponential(
                                rng, RunTumbleMotility.MEAN_RUN_S);
                        speed = RunTumbleMotility.speed(rng);
                        RunTumbleMotility.direction(rng, direction);
                        runs.add(remaining);
                        speeds.add(speed);
                        runLengths.add(remaining * speed);
                    } else {
                        remaining = RunTumbleMotility.exponential(
                                rng, RunTumbleMotility.MEAN_TUMBLE_S);
                        tumbles.add(remaining);
                    }
                }
            }
        }

        PrintStream msd = stream("msd.csv");
        msd.println("t_s;msd_um2;mean_x_um;mean_y_um;mean_z_um");
        for (int i = 0; i < OBS_S.length; i++) {
            msd.printf(Locale.US, "%.6f;%.12e;%.12e;%.12e;%.12e%n",
                    OBS_S[i], sumR2[i] / N_TRAJECTORIES,
                    sumX[i] / N_TRAJECTORIES,
                    sumY[i] / N_TRAJECTORIES,
                    sumZ[i] / N_TRAJECTORIES);
        }
        msd.close();

        double finalMsd = sumR2[OBS_S.length - 1] / N_TRAJECTORIES;
        double simulatedD = finalMsd / (6.0 * OBS_S[OBS_S.length - 1]);
        PrintStream summary = stream("job7_summary.csv");
        summary.println("metric;value;reference;units");
        summary.printf(Locale.US, "mean_run;%.12f;%.12f;s%n",
                runs.mean(), RunTumbleMotility.MEAN_RUN_S);
        summary.printf(Locale.US, "mean_tumble;%.12f;%.12f;s%n",
                tumbles.mean(), RunTumbleMotility.MEAN_TUMBLE_S);
        summary.printf(Locale.US, "mean_speed;%.12f;%.12f;um/s%n",
                speeds.mean(), RunTumbleMotility.MEAN_SPEED_UM_S);
        summary.printf(Locale.US, "speed_sd;%.12f;%.12f;um/s%n",
                speeds.sd(), RunTumbleMotility.SPEED_SD_UM_S);
        summary.printf(Locale.US, "run_fraction;%.12f;0.86;1%n",
                runs.mean() / (runs.mean() + tumbles.mean()));
        summary.printf(Locale.US, "mean_run_length;%.12f;%.12f;um%n",
                runLengths.mean(),
                RunTumbleMotility.MEAN_SPEED_UM_S
                        * RunTumbleMotility.MEAN_RUN_S);
        summary.printf(Locale.US, "D_theory;%.12f;198;um2/s%n",
                RunTumbleMotility.theoreticalDiffusivity());
        summary.printf(Locale.US, "D_simulated;%.12f;185;um2/s%n", simulatedD);
        summary.printf(Locale.US, "run_events;%d;0;count%n", runs.n);
        summary.printf(Locale.US, "tumble_events;%d;0;count%n", tumbles.n);
        summary.printf(Locale.US, "trajectories;%d;%d;count%n",
                N_TRAJECTORIES, N_TRAJECTORIES);
        summary.printf(Locale.US, "seed;%d;%d;1%n", SEED, SEED);
        summary.close();

        System.out.println("Job 7 run-and-tumble motility (standalone)");
        System.out.printf(Locale.US,
                "  sampled run=%.4f s tumble=%.4f s speed=%.4f +/- %.4f um/s%n",
                runs.mean(), tumbles.mean(), speeds.mean(), speeds.sd());
        System.out.printf(Locale.US,
                "  D_theory=%.3f D_simulated=%.3f um^2/s, N=%d seed=%d%n",
                RunTumbleMotility.theoreticalDiffusivity(), simulatedD,
                N_TRAJECTORIES, SEED);
        System.out.println("  basal motility only; no chemotaxis or chassis coupling");
        System.out.println("wrote " + DIR);
    }
}
