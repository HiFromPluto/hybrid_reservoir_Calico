package bsim.laneb;

import bsim.BSim;
import bsim.BSimRandom;
import bsim.ode.BSimOdeSolver;
import bsim.particle.BSimBacterium;

import javax.vecmath.Vector3d;

/**
 * Run-tumble E. coli with Barkai–Leibler att−rep chemotaxis.
 * Growth, death, AHL, QS, lookahead, and Brownian (motile arm) are OFF.
 * B0_MOTILITY_MEMORY. Not C1c. Not Fig. 4b.
 */
public final class LaneBBacterium extends BSimBacterium {

    public enum Motility { MOTILE, FROZEN, THERMAL }

    private final LaneBChemotaxisDynamics chemotaxis;
    private double[] chemotaxisState;
    private final Motility motility;
    private double localAtt;
    private double localRep;

    public LaneBBacterium(BSim sim, Vector3d position, BSimRandom rng, Motility motility) {
        super(sim, new Vector3d(position));
        setRadius(LaneBIdentity.RADIUS);
        setSurfaceAreaGrowthRate(0.0);
        setForceMagnitude(LaneBIdentity.FORCE_PN);
        pEndRunElse(LaneBIdentity.P_END_RUN_ELSE);
        pEndRunUp(LaneBIdentity.P_END_RUN_ELSE);
        pEndTumble(LaneBIdentity.P_END_TUMBLE);
        this.motility = motility;
        this.chemotaxis = new LaneBChemotaxisDynamics(rng);
        this.chemotaxisState = chemotaxis.getICs();
        this.localAtt = 0.0;
        this.localRep = 0.0;
    }

    public Motility motility() {
        return motility;
    }

    public double cheyP() {
        return chemotaxisState[2];
    }

    public void sampleLigand(double att, double rep) {
        this.localAtt = att;
        this.localRep = rep;
    }

    @Override
    public void grow() { }

    @Override
    public void replicate() { }

    @Override
    public double pEndRun() {
        double y = chemotaxisState[2] / LaneBIdentity.Y0;
        double factor = Math.max(0.2, Math.min(3.0, y));
        return LaneBIdentity.P_END_RUN_ELSE * factor;
    }

    @Override
    public void action() {
        if (motility == Motility.FROZEN) {
            return;
        }
        if (motility == Motility.THERMAL) {
            brownianForce();
            return;
        }
        chemotaxis.ligandConcentration = localAtt - localRep;
        chemotaxisState = BSimOdeSolver.rungeKutta45(
                chemotaxis, sim.getTime(), chemotaxisState, sim.getDt());
        chemotaxisState[0] = Math.max(0.0, Math.min(5.0, chemotaxisState[0]));
        chemotaxisState[1] = Math.max(0.0, Math.min(1.0, chemotaxisState[1]));
        chemotaxisState[2] = Math.max(0.0, chemotaxisState[2]);

        switch (motionState) {
            case RUNNING:
                if (sim.getRandom().nextDouble() < pEndRun() * sim.getDt()) {
                    motionState = MotionState.TUMBLING;
                }
                break;
            case TUMBLING:
                if (sim.getRandom().nextDouble() < pEndTumble() * sim.getDt()) {
                    bsim.BSimUtils.rotatePerp(direction, tumbleAngle(), sim.getRandom());
                    motionState = MotionState.RUNNING;
                }
                break;
            default:
                break;
        }
        if (motionState == MotionState.RUNNING) {
            rotationalDiffusion();
            flagellarForce();
        }
    }

    @Override
    public void updatePosition() {
        if (motility == Motility.FROZEN) {
            return;
        }
        super.updatePosition();
    }
}
