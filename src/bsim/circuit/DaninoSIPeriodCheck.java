package bsim.circuit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

/**
 * NARMA-blind peak-to-peak period of intracellular LuxI {@code I}.
 * Rules copied from D0b {@code period_check.py}. Object B, NOT_FIG4B.
 */
public final class DaninoSIPeriodCheck {

    private DaninoSIPeriodCheck() { }

    public static void refuseNarma(String[] argv) {
        String blob = String.join(" ", argv == null ? new String[0] : argv).toLowerCase(Locale.ROOT);
        if (blob.contains("narma") || blob.contains("mackey")) {
            throw new IllegalArgumentException("D1 period checker cannot see a NARMA target.");
        }
    }

    public static Result extract(double[] t, double[] luxI, Protocol p) {
        int start = 0;
        while (start < t.length && t[start] < p.tDiscard - 1e-12) {
            start++;
        }
        int nW = t.length - start;
        if (nW < 8) {
            return Result.noPeriod(Double.NaN);
        }
        double[] yW = Arrays.copyOfRange(luxI, start, t.length);
        double[] tW = Arrays.copyOfRange(t, start, t.length);
        double p5 = percentile(yW, 5.0);
        double p95 = percentile(yW, 95.0);
        double amp = p95 - p5;
        double median = percentile(yW, 50.0);
        double relAmp = median > 0.0 ? amp / median : 0.0;
        if (!(amp > 0.0) || Double.isNaN(amp)) {
            return Result.noPeriod(mean(yW));
        }
        double prominence = p.peakRelProminence * amp;
        int distance = Math.max(1, (int) Math.ceil(p.peakMinSpacing / p.sampleDt));
        int[] idx = findPeaks(yW, prominence, distance);
        List<Double> times = new ArrayList<>();
        for (int j : idx) {
            times.add(tW[j]);
        }
        double tLo = tW[0];
        double tHi = tW[tW.length - 1];
        double persistCut = tLo + (1.0 - p.persistLastFraction) * (tHi - tLo);
        int nPersist = 0;
        for (double tau : times) {
            if (tau >= persistCut) {
                nPersist++;
            }
        }
        int nAmp = Math.max(8, (int) (p.ampLastFraction * yW.length));
        double ampFirst = percentile(Arrays.copyOfRange(yW, 0, nAmp), 95.0)
                - percentile(Arrays.copyOfRange(yW, 0, nAmp), 5.0);
        double ampLast = percentile(Arrays.copyOfRange(yW, yW.length - nAmp, yW.length), 95.0)
                - percentile(Arrays.copyOfRange(yW, yW.length - nAmp, yW.length), 5.0);
        double persist = ampFirst > 0.0 ? ampLast / ampFirst : 0.0;
        boolean osc = times.size() >= p.minPeaks
                && nPersist >= 1
                && persist >= p.ampPersistRatio
                && relAmp >= p.relAmplitudeMin;
        double period = Double.NaN;
        String flag = "NO_PERIOD";
        if (osc) {
            flag = "OSC";
            double sum = 0.0;
            for (int i = 1; i < times.size(); i++) {
                sum += times.get(i) - times.get(i - 1);
            }
            period = sum / (times.size() - 1);
        }
        return new Result(flag, times.size(), nPersist, period, times, relAmp, persist, mean(yW), p5, p95);
    }

    static int[] findPeaks(double[] y, double prominence, int distance) {
        List<Integer> cand = new ArrayList<>();
        for (int i = 1; i < y.length - 1; i++) {
            if (y[i] >= y[i - 1] && y[i] > y[i + 1]) {
                cand.add(i);
            }
        }
        List<Integer> kept = new ArrayList<>();
        for (int i : cand) {
            if (peakProminence(y, i) >= prominence - 1e-15) {
                kept.add(i);
            }
        }
        if (kept.isEmpty() || distance <= 1) {
            return toInt(kept);
        }
        Integer[] order = kept.toArray(new Integer[0]);
        Arrays.sort(order, (a, b) -> Double.compare(y[b], y[a]));
        boolean[] taken = new boolean[y.length];
        List<Integer> out = new ArrayList<>();
        for (int i : order) {
            boolean ok = true;
            int lo = Math.max(0, i - distance + 1);
            int hi = Math.min(y.length - 1, i + distance - 1);
            for (int j = lo; j <= hi; j++) {
                if (taken[j]) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                taken[i] = true;
                out.add(i);
            }
        }
        out.sort(Integer::compareTo);
        return toInt(out);
    }

    /** Prominence: height above the higher of the left and right flanking minima. */
    static double peakProminence(double[] y, int i) {
        double leftMin = y[i];
        for (int j = i - 1; j >= 0; j--) {
            if (y[j] < leftMin) {
                leftMin = y[j];
            }
            if (y[j] > y[i]) {
                break;
            }
        }
        double rightMin = y[i];
        for (int j = i + 1; j < y.length; j++) {
            if (y[j] < rightMin) {
                rightMin = y[j];
            }
            if (y[j] > y[i]) {
                break;
            }
        }
        return y[i] - Math.max(leftMin, rightMin);
    }

    static double percentile(double[] a, double q) {
        double[] s = a.clone();
        Arrays.sort(s);
        if (s.length == 1) {
            return s[0];
        }
        double pos = (q / 100.0) * (s.length - 1);
        int lo = (int) Math.floor(pos);
        int hi = (int) Math.ceil(pos);
        if (lo == hi) {
            return s[lo];
        }
        return s[lo] + (pos - lo) * (s[hi] - s[lo]);
    }

    static double mean(double[] a) {
        double s = 0.0;
        for (double v : a) {
            s += v;
        }
        return s / a.length;
    }

    private static int[] toInt(List<Integer> xs) {
        int[] out = new int[xs.size()];
        for (int i = 0; i < xs.size(); i++) {
            out[i] = xs.get(i);
        }
        return out;
    }

    public static final class Protocol {
        public final double tDiscard;
        public final double sampleDt;
        public final double peakMinSpacing;
        public final double peakRelProminence;
        public final int minPeaks;
        public final double persistLastFraction;
        public final double ampLastFraction;
        public final double ampPersistRatio;
        public final double relAmplitudeMin;

        public Protocol(
                double tDiscard,
                double sampleDt,
                double peakMinSpacing,
                double peakRelProminence,
                int minPeaks,
                double persistLastFraction,
                double ampLastFraction,
                double ampPersistRatio,
                double relAmplitudeMin) {
            this.tDiscard = tDiscard;
            this.sampleDt = sampleDt;
            this.peakMinSpacing = peakMinSpacing;
            this.peakRelProminence = peakRelProminence;
            this.minPeaks = minPeaks;
            this.persistLastFraction = persistLastFraction;
            this.ampLastFraction = ampLastFraction;
            this.ampPersistRatio = ampPersistRatio;
            this.relAmplitudeMin = relAmplitudeMin;
        }

        public static Protocol d0b() {
            return new Protocol(180.0, 0.5, 25.0, 0.20, 4, 0.30, 0.40, 0.40, 0.25);
        }
    }

    public static final class Result {
        public final String flag;
        public final int nPeaks;
        public final int nPersistPeaks;
        public final double period;
        public final List<Double> peakTimes;
        public final double relAmplitude;
        public final double ampPersist;
        public final double meanI;
        public final double p5I;
        public final double p95I;

        public Result(
                String flag,
                int nPeaks,
                int nPersistPeaks,
                double period,
                List<Double> peakTimes,
                double relAmplitude,
                double ampPersist,
                double meanI,
                double p5I,
                double p95I) {
            this.flag = flag;
            this.nPeaks = nPeaks;
            this.nPersistPeaks = nPersistPeaks;
            this.period = period;
            this.peakTimes = peakTimes;
            this.relAmplitude = relAmplitude;
            this.ampPersist = ampPersist;
            this.meanI = meanI;
            this.p5I = p5I;
            this.p95I = p95I;
        }

        public static Result noPeriod(double meanI) {
            return new Result("NO_PERIOD", 0, 0, Double.NaN, List.of(), Double.NaN, Double.NaN, meanI, Double.NaN, Double.NaN);
        }
    }
}
