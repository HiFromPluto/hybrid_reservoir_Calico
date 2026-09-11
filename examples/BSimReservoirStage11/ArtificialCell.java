package BSimReservoirStage11;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.ode.BSimOdeSolver;
import bsim.particle.BSimBacterium;

import javax.vecmath.Vector3d;
import java.util.Random;

/**
 * Artificial Cell -- fixed-position signal source with ODE-driven vesicle release.
 * Ported from reservoir_new (AC_OurReservoir_Scenario.ArtificialCell).
 *
 * Release: i_R = gammaSym * (CS_in / (CS_in + Km)) * gate(s)
 * Gate:    gate(s) = 1 / (1 + exp(-GATE_SLOPE * (s - s_threshold)))
 */
public class ArtificialCell extends BSimBacterium {

    public static double GATE_SLOPE = 10.0;

    public double forcedInput = 0.0;

    private final ACConfig cfg;
    private final BSimChemicalField attractantField;
    private final BSimChemicalField repellentField;
    private final BSimChemicalField ahlField;
    private final BSimChemicalField glucoseField;

    private final ACInternalDynamics odeSystem;
    private double[] internalState;

    private final double CS_max;
    private final double gammaSym;
    private final double Km;
    private double CS_in;

    public ArtificialCell(BSim sim, ACConfig cfg,
                          BSimChemicalField attractantField,
                          BSimChemicalField repellentField,
                          BSimChemicalField ahlField,
                          BSimChemicalField glucoseField) {
        super(sim, new Vector3d(cfg.position));
        setRadius(3.0);
        this.cfg = cfg;
        this.attractantField = attractantField;
        this.repellentField  = repellentField;
        this.ahlField        = ahlField;
        this.glucoseField    = glucoseField;
        this.odeSystem       = new ACInternalDynamics();
        this.internalState   = odeSystem.getICs();

        if (cfg.heterogeneityEnabled) {
            Random rand = new Random();
            CS_max   = sampleLognormal(cfg.CS_max_mean,   cfg.CS_max_cv,   rand);
            gammaSym = sampleLognormal(cfg.gammaSym_mean, cfg.gammaSym_cv, rand);
            Km       = sampleLognormal(cfg.Km_mean,       cfg.Km_cv,       rand);
        } else {
            CS_max   = cfg.CS_max_mean;
            gammaSym = cfg.gammaSym_mean;
            Km       = cfg.Km_mean;
        }
        CS_in = CS_max;
    }

    @Override
    public void updatePosition() { }

    @Override
    public void action() {
        double dt = sim.getDt();

        odeSystem.externalInput = Math.max(0.0, Math.min(1.0, forcedInput));
        internalState = BSimOdeSolver.rungeKutta45(
                odeSystem, sim.getTime(), internalState, dt);
        if (internalState[0] < 0) internalState[0] = 0;
        if (internalState[1] < 0) internalState[1] = 0;

        double s = internalState[1];
        double gate = 1.0 / (1.0 + Math.exp(-GATE_SLOPE * (s - cfg.s_threshold)));

        double iR = 0.0;
        if (CS_in > 0.0) {
            iR = gammaSym * (CS_in / (CS_in + Km)) * gate;
        }

        double iLeak = 0.0;
        if (cfg.leakageEnabled && CS_in > 0.0) {
            double localConc  = localConcentration();
            double C_in_norm  = CS_in / CS_max;
            double C_out_norm = localConc / CS_max;
            double gradient   = C_in_norm - C_out_norm;
            if (gradient > 0.0) {
                double area = 4.0 * Math.PI * radius * radius;
                iLeak = cfg.P_leak * area * gradient * CS_max;
            }
        }

        double released = (iR + iLeak) * dt;
        CS_in = Math.max(0.0, CS_in - released);

        if (released > 0.0) inject(released);

        if (cfg.refillingEnabled) {
            double refill = cfg.k_refill * (CS_max - CS_in);
            CS_in = Math.min(CS_max, CS_in + refill * dt);
        }
    }

    private double localConcentration() {
        switch (cfg.secretionType) {
            case ATTRACTANT: return attractantField.getConc(position);
            case REPELLENT:  return repellentField.getConc(position);
            case AHL:        return ahlField.getConc(position);
            case NUTRIENT:   return glucoseField != null ? glucoseField.getConc(position) : 0.0;
            case MIXED:      return (attractantField.getConc(position)
                                   + ahlField.getConc(position)) * 0.5;
            default:         return 0.0;
        }
    }

    private void inject(double amount) {
        switch (cfg.secretionType) {
            case ATTRACTANT:
                attractantField.addQuantity(position, amount); break;
            case REPELLENT:
                repellentField.addQuantity(position, amount); break;
            case AHL:
                ahlField.addQuantity(position, amount); break;
            case NUTRIENT:
                if (glucoseField != null) glucoseField.addQuantity(position, amount);
                break;
            case MIXED:
                attractantField.addQuantity(position, amount * 0.5);
                ahlField.addQuantity(position, amount * 0.5); break;
        }
    }

    private static double sampleLognormal(double mean, double cv, Random rand) {
        double sigma = Math.sqrt(Math.log(1.0 + cv * cv));
        double mu    = Math.log(mean) - 0.5 * sigma * sigma;
        return Math.exp(mu + sigma * rand.nextGaussian());
    }

    public double getSecretionReadiness() { return internalState[1]; }
    public double getStoreFraction() { return CS_max > 0 ? CS_in / CS_max : 0.0; }
    public ACConfig.SecretionType getType() { return cfg.secretionType; }
}
