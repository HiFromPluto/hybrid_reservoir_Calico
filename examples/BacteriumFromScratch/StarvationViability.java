package BacteriumFromScratch;

/**
 * Biselli–Schink–Gerland carbon-starvation viability laws.
 *
 * <p>Biselli, Schink & Gerland 2020, Molecular Systems Biology 16:e9478,
 * DOI 10.15252/msb.20209478, with mandatory erratum
 * DOI 10.15252/msb.202010070. Supply/demand mechanism from Schink et al.
 * 2019, Cell Systems 9:64–73.e3, DOI 10.1016/j.cels.2019.06.003.
 *
 * <p><b>Units:</b> pre-starvation growth rate in h^-1, starvation time in
 * days, death rate in day^-1, nutrient amounts in fmol/CFU. This standalone
 * module is not a death switch for the nutrient-replete Job 3c colony.
 */
public final class StarvationViability {

    /** TAKEN: central rounded death-law prefactor, day^-1 (not hour^-1). */
    public static final double GAMMA0_PER_DAY = 0.21;

    /** TAKEN: central rounded growth-death slope, hours. */
    public static final double DEATH_SLOPE_H = 1.0;

    /** TAKEN: Schink WT glycerol death rate, day^-1. */
    public static final double SCHINK_GAMMA_PER_DAY = 0.43;

    /** TAKEN: Schink WT glycerol maintenance, fmol CFU^-1 day^-1. */
    public static final double SCHINK_BETA_FMOL_PER_CFU_DAY = 0.49;

    /** TAKEN: Schink glycerol pulse, fmol CFU^-1. */
    public static final double SCHINK_E0_FMOL_PER_CFU = 0.40;

    private StarvationViability() {}

    /** Biselli central growth-death relation, day^-1. */
    public static double deathRatePerDay(double preStarvationMuPerHour) {
        return GAMMA0_PER_DAY
                * Math.exp(DEATH_SLOPE_H * preStarvationMuPerHour);
    }

    /** Exact relative viability under constant starvation death rate. */
    public static double exactRelativeViability(double timeDays,
            double preStarvationMuPerHour) {
        return Math.exp(-deathRatePerDay(preStarvationMuPerHour) * timeDays);
    }

    /** Schink Eq. (6): lag after adding UV-killed cells, days. */
    public static double uvKilledLagDays(double uvKilledToViableRatio) {
        return uvKilledToViableRatio / SCHINK_GAMMA_PER_DAY;
    }

    /** Schink maintenance-limited lag T=E0/beta, days. */
    public static double maintenanceLagDays() {
        return SCHINK_E0_FMOL_PER_CFU / SCHINK_BETA_FMOL_PER_CFU_DAY;
    }

    /** Recycled nutrient per death from gamma=beta/alpha, fmol/CFU. */
    public static double recycledNutrientPerDeathFmol() {
        return SCHINK_BETA_FMOL_PER_CFU_DAY / SCHINK_GAMMA_PER_DAY;
    }

    /**
     * Feast-famine log fitness. Feast time is hours; starvation time is days.
     */
    public static double feastFamineFitness(double muPerHour,
            double feastHours, double starvationDays) {
        return muPerHour * feastHours
                - deathRatePerDay(muPerHour) * starvationDays;
    }

    /** Biselli closed-form optimum, h^-1. */
    public static double optimumMuPerHour(double feastHours,
            double starvationDays) {
        double ratio = feastHours
                / (DEATH_SLOPE_H * GAMMA0_PER_DAY * starvationDays);
        return Math.log(ratio) / DEATH_SLOPE_H;
    }

    /** RK4 step for dN/dt=-gamma*N, with t in days. */
    public static double rk4ViabilityStep(double n, double dtDays,
            double preStarvationMuPerHour) {
        double gamma = deathRatePerDay(preStarvationMuPerHour);
        double k1 = -gamma * n;
        double k2 = -gamma * (n + 0.5 * dtDays * k1);
        double k3 = -gamma * (n + 0.5 * dtDays * k2);
        double k4 = -gamma * (n + dtDays * k3);
        return n + dtDays * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
}
