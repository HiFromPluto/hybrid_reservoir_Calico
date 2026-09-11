package bsim;

import java.util.ArrayList;
import java.util.List;

/**
 * Opt-in deterministic scheduler for new simulations. Legacy BSim.export()
 * keeps its historical inclusive-loop behavior.
 */
public final class BSimStepScheduler {
    public interface Phases {
        void observe(Context context);
        void applyBoundary(Context context);
        void transport(Context context, double duration);
        void sampleState(Context context);
        void integrateModels(Context context);
        void depositFluxes(Context context);
    }

    public abstract static class Adapter implements Phases {
        @Override public void observe(Context context) { }
        @Override public void applyBoundary(Context context) { }
        @Override public void transport(Context context, double duration) { }
        @Override public void sampleState(Context context) { }
        @Override public void integrateModels(Context context) { }
        @Override public void depositFluxes(Context context) { }
    }

    public interface MechanicsEvent {
        void run(Context context);
    }

    public static final class Context {
        private final int completedUpdates;
        private final double startTime;
        private final double endTime;
        private final double dt;
        private final boolean initialObservation;

        private Context(int completedUpdates, double startTime, double endTime,
                        double dt, boolean initialObservation) {
            this.completedUpdates = completedUpdates;
            this.startTime = startTime;
            this.endTime = endTime;
            this.dt = dt;
            this.initialObservation = initialObservation;
        }

        public int getCompletedUpdates() { return completedUpdates; }
        public double getStartTime() { return startTime; }
        public double getMidpointTime() { return startTime + 0.5 * dt; }
        public double getEndTime() { return endTime; }
        public double getDt() { return dt; }
        public boolean isInitialObservation() { return initialObservation; }
    }

    private static final class ScheduledMechanics {
        private final String name;
        private final int stride;
        private final MechanicsEvent event;

        private ScheduledMechanics(String name, int stride, MechanicsEvent event) {
            this.name = name;
            this.stride = stride;
            this.event = event;
        }
    }

    private final double dt;
    private final double duration;
    private final int updateCount;
    private final List<ScheduledMechanics> mechanics = new ArrayList<ScheduledMechanics>();
    private BSim simulation;

    public BSimStepScheduler(double dt, double duration) {
        requirePositiveFinite(dt, "dt");
        requireNonnegativeFinite(duration, "duration");
        this.dt = dt;
        this.duration = duration;
        this.updateCount = exactSteps(duration, dt, "simulation duration");
    }

    public BSimStepScheduler(BSim sim) {
        this(sim.getDt(), sim.getSimulationTime());
        this.simulation = sim;
    }

    public void addMechanicsEvent(String name, double period, MechanicsEvent event) {
        if (name == null || name.trim().isEmpty())
            throw new IllegalArgumentException("mechanics event name is required");
        if (event == null) throw new IllegalArgumentException("mechanics event is required");
        int stride = exactSteps(period, dt, "mechanics period " + name);
        if (stride <= 0)
            throw new IllegalArgumentException("mechanics period must be at least one base step");
        mechanics.add(new ScheduledMechanics(name, stride, event));
    }

    /**
     * Run an explicit initial observation followed by exactly duration/dt
     * updates. Observations after updates carry the resulting end timestamp.
     */
    public void run(Phases phases) {
        if (phases == null) throw new IllegalArgumentException("phase handler is required");
        setSimulationClock(0.0);
        phases.observe(new Context(0, 0.0, 0.0, dt, true));
        for (int zeroBased = 0; zeroBased < updateCount; zeroBased++) {
            int completed = zeroBased + 1;
            double start = zeroBased * dt;
            double end = completed * dt;
            Context context = new Context(completed, start, end, dt, false);
            setSimulationClock(zeroBased);
            phases.applyBoundary(context);
            phases.transport(context, 0.5 * dt);
            setSimulationClock(zeroBased + 0.5);
            phases.sampleState(context);
            phases.integrateModels(context);
            phases.depositFluxes(context);
            phases.transport(context, 0.5 * dt);
            setSimulationClock(completed);
            for (ScheduledMechanics scheduled : mechanics)
                if (completed % scheduled.stride == 0)
                    scheduled.event.run(context);
            phases.observe(context);
        }
    }

    public double getDt() { return dt; }
    public double getDuration() { return duration; }
    public int getUpdateCount() { return updateCount; }

    private void setSimulationClock(double timestep) {
        if (simulation != null) simulation.setSchedulerTimestep(timestep);
    }

    public static int exactSteps(double duration, double dt, String name) {
        requireNonnegativeFinite(duration, name);
        requirePositiveFinite(dt, "dt");
        double ratio = duration / dt;
        long rounded = Math.round(ratio);
        double tolerance = 1e-10 * Math.max(1.0, Math.abs(duration));
        if (rounded > Integer.MAX_VALUE
                || Math.abs(rounded * dt - duration) > tolerance) {
            throw new IllegalArgumentException(
                    name + " must be an integer number of base steps: "
                            + duration + " / " + dt);
        }
        return (int) rounded;
    }

    private static void requirePositiveFinite(double value, String name) {
        if (!Double.isFinite(value) || value <= 0.0)
            throw new IllegalArgumentException(name + " must be finite and positive: " + value);
    }

    private static void requireNonnegativeFinite(double value, String name) {
        if (!Double.isFinite(value) || value < 0.0)
            throw new IllegalArgumentException(name + " must be finite and nonnegative: " + value);
    }
}
