package BacteriumFromScratch;

import java.io.File;
import java.io.PrintStream;
import java.util.Locale;

/**
 * Job 6 — standalone carbon-starvation viability.
 *
 * <p>No BSim and no chassis imports. Death remains OFF in Jobs 2/3/3b/3c:
 * Biselli's law applies after transfer to zero carbon and does not supply a
 * cited interpolation for the nutrient-replete Job 3c dish.
 */
public final class Job6Sims {

    static final String DIR = "results/job6_starvation/";
    static final double DT_DAY = 0.01;       // ENGINEERING
    static final double T_END_DAY = 10.0;    // ENGINEERING
    static final double VIABILITY_MU = 0.70; // TAKEN comparison condition
    static final double OPT_GRID = 1.0e-4;   // ENGINEERING, h^-1

    private Job6Sims() {}

    private static PrintStream stream(String name) {
        try {
            return new PrintStream(DIR + name, "UTF-8");
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private static double numericalOptimum(double feastHours,
            double starvationDays) {
        double bestMu = 0.0;
        double bestFitness = Double.NEGATIVE_INFINITY;
        for (double mu = 0.0; mu <= 3.0 + 0.5 * OPT_GRID; mu += OPT_GRID) {
            double f = StarvationViability.feastFamineFitness(
                    mu, feastHours, starvationDays);
            if (f > bestFitness) {
                bestFitness = f;
                bestMu = mu;
            }
        }
        return bestMu;
    }

    private static void writeViability() {
        PrintStream out = stream("viability_timeseries.csv");
        out.println("t_day;mu_pre_h-1;gamma_day-1;N_rk4;N_exact;abs_err");
        double n = 1.0;
        int steps = (int) Math.round(T_END_DAY / DT_DAY);
        for (int i = 0; i <= steps; i++) {
            double t = i * DT_DAY;
            double exact = StarvationViability.exactRelativeViability(
                    t, VIABILITY_MU);
            out.printf(Locale.US, "%.4f;%.4f;%.12f;%.15e;%.15e;%.6e%n",
                    t, VIABILITY_MU,
                    StarvationViability.deathRatePerDay(VIABILITY_MU),
                    n, exact, Math.abs(n - exact));
            if (i < steps) {
                n = StarvationViability.rk4ViabilityStep(
                        n, DT_DAY, VIABILITY_MU);
            }
        }
        out.close();
    }

    private static void writeSummary() {
        PrintStream out = stream("job6_summary.csv");
        out.println("section;case;arg1;arg2;value;reference");

        double[][] ev1 = {
            {0.10, 0.24, 0.02},
            {0.30, 0.29, 0.03},
            {0.50, 0.33, 0.03},
            {0.70, 0.40, 0.04},
        };
        for (double[] row : ev1) {
            out.printf(Locale.US,
                    "units;EV1_WT_chemostat;%.4f;%.4f;%.12f;%.4f%n",
                    row[0], row[2],
                    StarvationViability.deathRatePerDay(row[0]), row[1]);
        }

        out.printf(Locale.US, "lag;UV_50_50;1.0;0;%.12f;2.3%n",
                StarvationViability.uvKilledLagDays(1.0));
        out.printf(Locale.US, "lag;UV_30_70;%.12f;0;%.12f;1.2%n",
                30.0 / 70.0,
                StarvationViability.uvKilledLagDays(30.0 / 70.0));
        out.printf(Locale.US, "lag;maintenance;%.4f;%.4f;%.12f;0.825%n",
                StarvationViability.SCHINK_E0_FMOL_PER_CFU,
                StarvationViability.SCHINK_BETA_FMOL_PER_CFU_DAY,
                StarvationViability.maintenanceLagDays());

        double[][] cycles = {
            {3.0, 6.0},
            {6.0, 6.0},
            {3.0, 12.0},
        };
        for (double[] cycle : cycles) {
            double closed = StarvationViability.optimumMuPerHour(
                    cycle[0], cycle[1]);
            double numeric = numericalOptimum(cycle[0], cycle[1]);
            out.printf(Locale.US,
                    "optimum;cycle;%.4f;%.4f;%.12f;%.12f%n",
                    cycle[0], cycle[1], numeric, closed);
        }

        double correctedMu = StarvationViability.optimumMuPerHour(3.0, 6.0);
        double oldMu = StarvationViability.optimumMuPerHour(3.0, 3.0);
        out.printf(Locale.US, "erratum;corrected_3h_6d;3;6;%.12f;%.12f%n",
                correctedMu,
                StarvationViability.deathRatePerDay(correctedMu));
        out.printf(Locale.US, "erratum;old_3h_3d;3;3;%.12f;%.12f%n",
                oldMu, StarvationViability.deathRatePerDay(oldMu));

        out.printf(Locale.US, "mechanism;alpha_beta_over_gamma;%.4f;%.4f;%.12f;1.14%n",
                StarvationViability.SCHINK_BETA_FMOL_PER_CFU_DAY,
                StarvationViability.SCHINK_GAMMA_PER_DAY,
                StarvationViability.recycledNutrientPerDeathFmol());
        out.printf(Locale.US, "mechanism;wrong_alphaG_substitution;%.4f;%.4f;%.12f;0.43%n",
                StarvationViability.SCHINK_BETA_FMOL_PER_CFU_DAY, 0.18,
                StarvationViability.SCHINK_BETA_FMOL_PER_CFU_DAY / 0.18);
        out.close();
    }

    public static void main(String[] args) {
        new File(DIR).mkdirs();
        System.out.println("Job 6 starvation viability (standalone)");
        System.out.println("  TAKEN gamma(mu) = 0.21 day^-1 exp[(1.0 h) mu]");
        System.out.println("  TAKEN Schink gamma=0.43 day^-1 beta=0.49 fmol/CFU/day");
        System.out.println("  ENGINEERING dt=" + DT_DAY + " day horizon="
                + T_END_DAY + " day optimum_grid=" + OPT_GRID + " h^-1");
        System.out.println("  death remains OFF in Jobs 2/3/3b/3c");
        writeViability();
        writeSummary();
        System.out.println("wrote " + DIR);
    }
}
