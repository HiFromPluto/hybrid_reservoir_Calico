package BSimReservoirStage11;

import javax.vecmath.Vector3d;

/**
 * Configuration data object for an ArtificialCell.
 * Ported from reservoir_new (AC_OurReservoir_Scenario.ACConfig).
 */
public final class ACConfig {

    public enum SecretionType {
        ATTRACTANT,
        REPELLENT,
        AHL,
        MIXED,
        /** Injects into glucoseField. Stage 7 AC5 nutrient channel; Stage 2 birth-input fix. */
        NUTRIENT
    }

    public final Vector3d position;
    public final SecretionType secretionType;

    public double CS_max_mean = 2e5;
    public double CS_max_cv   = 0.3;
    public double gammaSym_mean = 1e4;
    public double gammaSym_cv   = 0.25;
    public double Km_mean = 1e5;
    public double Km_cv   = 0.2;

    public double  s_threshold      = 0.5;
    public double  P_leak           = 0.001;
    public double  k_refill         = 0.1;
    public boolean refillingEnabled = true;
    public boolean leakageEnabled   = true;
    public boolean heterogeneityEnabled = true;

    private ACConfig(Vector3d position, SecretionType secretionType) {
        this.position      = new Vector3d(position);
        this.secretionType = secretionType;
    }

    public static Builder builder(Vector3d position, SecretionType secretionType) {
        return new Builder(position, secretionType);
    }

    public static final class Builder {
        private final ACConfig cfg;

        private Builder(Vector3d position, SecretionType secretionType) {
            cfg = new ACConfig(position, secretionType);
        }

        public Builder csMax(double mean, double cv) {
            cfg.CS_max_mean = mean; cfg.CS_max_cv = cv; return this;
        }
        public Builder gammaSym(double mean, double cv) {
            cfg.gammaSym_mean = mean; cfg.gammaSym_cv = cv; return this;
        }
        public Builder km(double mean, double cv) {
            cfg.Km_mean = mean; cfg.Km_cv = cv; return this;
        }
        public Builder sThreshold(double s) { cfg.s_threshold = s; return this; }
        public Builder pLeak(double p) { cfg.P_leak = p; return this; }
        public Builder kRefill(double k) { cfg.k_refill = k; return this; }
        public Builder refillingEnabled(boolean e) { cfg.refillingEnabled = e; return this; }
        public Builder leakageEnabled(boolean e) { cfg.leakageEnabled = e; return this; }
        public Builder heterogeneityEnabled(boolean e) { cfg.heterogeneityEnabled = e; return this; }

        public ACConfig build() { return cfg; }
    }
}
