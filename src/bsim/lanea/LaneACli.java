package bsim.lanea;

import java.nio.file.Path;
import java.util.Locale;

/**
 * Placement-seed CLI for the named seed-replicate extra.
 * Allowed seeds are frozen: 101, 202, 303. Not a retune of J_max.
 */
public final class LaneACli {

    public final long seed;
    public final boolean drivenOnly;

    private LaneACli(long seed, boolean drivenOnly) {
        this.seed = seed;
        this.drivenOnly = drivenOnly;
    }

    public static LaneACli parse(String[] args) {
        long seed = LaneAConstants.RNG_SEED;
        boolean drivenOnly = false;
        for (int i = 0; i < args.length; i++) {
            String a = args[i];
            if ("--driven-only".equals(a)) {
                drivenOnly = true;
                continue;
            }
            if (a.startsWith("--seed=")) {
                seed = Long.parseLong(a.substring(7).trim());
                continue;
            }
            if ("--seed".equals(a)) {
                if (i + 1 >= args.length) {
                    throw new IllegalArgumentException("--seed needs a value");
                }
                seed = Long.parseLong(args[++i].trim());
                continue;
            }
        }
        if (seed != 101L && seed != 202L && seed != 303L) {
            throw new IllegalArgumentException(
                    "placement seed must be 101, 202, or 303 (frozen list); got " + seed);
        }
        return new LaneACli(seed, drivenOnly);
    }

    public Path outDir(Path defaultResults) {
        if (seed == LaneAConstants.RNG_SEED) {
            return defaultResults;
        }
        return defaultResults.resolve("seed_" + seed);
    }

    public void announce(String label) {
        System.out.printf(Locale.US, "  placement_seed=%d driven_only=%s %s%n",
                seed, drivenOnly, label);
    }
}
