package BacteriumFromScratch;

import java.io.PrintStream;
import java.util.Locale;

/**
 * Job 5 — Weber &amp; Buceta 2013 LuxI/LuxR quorum sensing, well-mixed.
 *
 * <p>Weber M, Buceta J (2013), BMC Systems Biology 7:6,
 * DOI 10.1186/1752-0509-7-6. Eq. (1) reaction set, Table 1 parameters.
 *
 * <p><b>Units are Weber's: nM and min.</b> Not the chassis's mM/um/s. There is
 * no conversion anywhere here because there is no coupling.
 *
 * <p><b>STANDALONE.</b> No chassis imports, no BSim. Job 5 must not touch
 * EcoliRodCell, ValdezHertzian or NutrientField, and does not use lambda_S,
 * K_S or Monod growth.
 *
 * <p><b>SPECIES CROSSING — ENGINEERING.</b> Weber models Vibrio fischeri; the
 * chassis is E. coli. Declared in PROTOCOL.md; every figure from this module
 * carries the label.
 *
 * <p><b>tau = 45 min is WEBER's</b> (RM/succinate, 30 C), not the chassis
 * T_div = 43.25 min. Do not substitute one for the other.
 *
 * <p><b>D = 10 /min is a transmembrane transport rate, not a diffusivity.</b>
 *
 * <p><b>Stiffness.</b> The ON branch is stiff; the OFF branch is not. In the
 * high state C2 reaches ~5400 nM, so the promoter mode
 * {@code Klux- + KFLUX*(C2+G)} reaches ~282 /min against RK4's real-axis
 * stability limit of 2.785/h. Fixed-step RK4 would need dt ~ 1e-3 min, i.e.
 * ~3e9 steps for the full sweep. This module uses ADAPTIVE step-doubling RK4
 * with Richardson correction. A quasi-steady-state reduction of the fast
 * binding equilibria would remove the stiffness but would CHANGE WEBER'S
 * MODEL, and is banned by PROTOCOL.md — this is a port.
 */
public final class Job5Sims {

    private Job5Sims() {}

    static final String DIR = "results/job5_seed101/";

    // ---- Weber Table 1. TAKEN unless noted. nM, min. ----
    static final double KD1 = 100.0;      // LuxR-A dissociation, nM
    static final double K1_MINUS = 10.0;  // LuxR-A unbinding, /min
    static final double KD2 = 20.0;       // dimerisation, nM
    static final double K2_MINUS = 1.0;   // dimer dissociation, /min
    static final double KDLUX = 200.0;    // complex->promoter, nM
    static final double KLUX_MINUS = 10.0;
    static final double B_BURST = 20.0;
    static final double KR = 200.0 / B_BURST;   // luxR transcription, /min
    static final double KI = 50.0 / B_BURST;    // luxI transcription, /min
    static final double D_MR = 0.347;     // luxR mRNA degradation, /min
    static final double D_MI = 0.347;
    static final double PR = B_BURST * D_MR;    // translation, /min
    static final double PI_ = B_BURST * D_MI;
    static final double ALPHA_R = 0.001;
    static final double ALPHA_I = 0.01;
    static final double D_A = 0.001;      // autoinducer degradation, /min
    static final double D_C2 = 0.002;
    static final double D_C = 0.002;
    static final double D_R = 0.002;
    static final double D_I = 0.01;
    static final double D_MEMBRANE = 10.0;      // transport RATE, /min
    static final double TAU_MIN = 45.0;         // Weber doubling time
    static final double KA_LUX02 = 0.04;        // A synthesis by LuxI, /min
    static final double KA_LUX01 = 0.0;         // lux01 lacks luxI

    static final double V0_UM3 = 1.5;           // cell volume, um^3
    static final double VTOT_UL = 2.0e-4;       // culture volume, uL
    static final int N_CELLS = 100;

    // ---- derived ----
    static final double GAMMA = Math.log(2.0) / TAU_MIN;          // /min
    static final double KF1 = K1_MINUS / KD1;                     // /nM/min
    static final double KF2 = K2_MINUS / KD2;
    static final double KFLUX = KLUX_MINUS / KDLUX;
    /** V_c,tot and V_ext in uL; 1 um^3 = 1e-9 uL. */
    static final double VC_TOT_UL = N_CELLS * V0_UM3 * 1.0e-9;
    static final double V_EXT_UL = VTOT_UL - VC_TOT_UL;
    static final double R_VOL = VC_TOT_UL / V_EXT_UL;

    /**
     * Total DNA concentration. DERIVED, not in Weber Table 1: one operon copy
     * per cell of volume V_0, c = 1/(N_A V_0). Stated explicitly because the
     * paper does not give it and the model needs it.
     */
    static final double DNA_TOT_NM =
            1.0 / (6.02214076e23 * V0_UM3 * 1.0e-15) * 1.0e9;

    // ---- integrator controls: ENGINEERING, frozen by the Gate 7 study ----
    /** Relative error tolerance. This is the parameter Gate 7 halves. */
    static double rtol = 1.0e-7;
    static final double ATOL = 1.0e-10;   // nM; below this a species is numerically zero
    static final double H_INIT = 1.0e-3;
    static final double H_MIN = 1.0e-9;
    static final double H_MAX = 1.0;

    // ---- state layout ----
    static final int MR = 0, MI = 1, R = 2, I = 3, A = 4,
                     C1 = 5, C2 = 6, G = 7, GC = 8, AE = 9;
    static final int NSPEC = 10;
    static final String[] NAMES = {
        "mRNA_luxR", "mRNA_luxI", "LuxR", "LuxI_GFP", "A",
        "LuxR_A", "LuxR_A_2", "DNA", "DNA_C2", "A_ext"
    };

    static double[] uninduced() {
        double[] y = new double[NSPEC];
        y[G] = DNA_TOT_NM;   // all DNA free, everything else zero
        return y;
    }

    /** Weber Eq. (1) as ODEs. cAStar is the exogenous control parameter (nM). */
    static void rhs(double[] y, double cAStar, double kA, double[] out) {
        double mr = y[MR], mi = y[MI], r = y[R], i = y[I], a = y[A];
        double c1 = y[C1], c2 = y[C2], g = y[G], gc = y[GC], ae = y[AE];

        double bind1 = KF1 * r * a - K1_MINUS * c1;          // luxR + A <-> luxR.A
        double dim = KF2 * c1 * c1 - K2_MINUS * c2;          // 2(luxR.A) <-> (luxR.A)2
        double prom = KFLUX * c2 * g - KLUX_MINUS * gc;      // (luxR.A)2 + DNA <-> DNA.(luxR.A)2
        double transport = D_MEMBRANE * (a - ae);            // A <-> A_ext

        out[MR] = ALPHA_R * KR * g + KR * gc - D_MR * mr - GAMMA * mr;
        out[MI] = ALPHA_I * KI * g + KI * gc - D_MI * mi - GAMMA * mi;
        out[R] = PR * mr - bind1 - D_R * r - GAMMA * r;
        out[I] = PI_ * mi - D_I * i - GAMMA * i;
        out[A] = kA * i - bind1 - transport - D_A * a - GAMMA * a;
        out[C1] = bind1 - 2.0 * dim - D_C * c1 - GAMMA * c1;
        out[C2] = dim - prom - D_C2 * c2 - GAMMA * c2;
        // DNA duplication (+gamma*(g+gc)) exactly cancels dilution (-gamma*g),
        // so total DNA is conserved -- Gate 5 checks this invariant.
        out[G] = -prom + GAMMA * gc;
        out[GC] = prom - GAMMA * gc;
        // A_ext takes NO growth dilution (Weber). Influx/efflux at gamma.
        out[AE] = R_VOL * transport - D_A * ae
                + GAMMA * (cAStar * VTOT_UL / V_EXT_UL - ae);
    }

    // scratch buffers
    private static final double[] KA_ = new double[NSPEC], KB = new double[NSPEC],
            KC = new double[NSPEC], KD_ = new double[NSPEC], TMP = new double[NSPEC],
            YBIG = new double[NSPEC], YSML = new double[NSPEC], YMID = new double[NSPEC];

    /** One classical RK4 step from y into out. */
    static void rk4(double[] y, double cAStar, double kA, double h, double[] out) {
        rhs(y, cAStar, kA, KA_);
        for (int s = 0; s < NSPEC; s++) TMP[s] = y[s] + 0.5 * h * KA_[s];
        rhs(TMP, cAStar, kA, KB);
        for (int s = 0; s < NSPEC; s++) TMP[s] = y[s] + 0.5 * h * KB[s];
        rhs(TMP, cAStar, kA, KC);
        for (int s = 0; s < NSPEC; s++) TMP[s] = y[s] + h * KC[s];
        rhs(TMP, cAStar, kA, KD_);
        for (int s = 0; s < NSPEC; s++) {
            out[s] = y[s] + h / 6.0 * (KA_[s] + 2.0 * KB[s] + 2.0 * KC[s] + KD_[s]);
        }
    }

    /** Integration statistics, so the cost of the tolerance choice stays visible. */
    static final class Stats {
        long accepted;
        long rejected;
        double hMin = Double.MAX_VALUE;
        double hMax;
        double worstNegative;
        double maxDnaDrift;
    }

    /**
     * Adaptive step-doubling RK4: one step of h against two of h/2, keeping the
     * h/2 result with Richardson correction. Order 4 error control, order 5
     * accepted value.
     */
    static void integrate(double[] y, double cAStar, double kA, double tEnd, Stats st) {
        double t = 0.0;
        double h = H_INIT;
        while (t < tEnd) {
            if (t + h > tEnd) h = tEnd - t;
            rk4(y, cAStar, kA, h, YBIG);
            rk4(y, cAStar, kA, 0.5 * h, YMID);
            rk4(YMID, cAStar, kA, 0.5 * h, YSML);

            double err = 0.0;
            for (int s = 0; s < NSPEC; s++) {
                double scale = ATOL + rtol * Math.max(Math.abs(YSML[s]), Math.abs(y[s]));
                err = Math.max(err, Math.abs(YSML[s] - YBIG[s]) / scale);
            }

            if (err <= 1.0 || h <= H_MIN) {
                for (int s = 0; s < NSPEC; s++) {
                    y[s] = YSML[s] + (YSML[s] - YBIG[s]) / 15.0;
                    if (!Double.isFinite(y[s])) {
                        throw new ArithmeticException("non-finite " + NAMES[s]
                                + " at cA*=" + cAStar + " t=" + t);
                    }
                    if (y[s] < st.worstNegative) st.worstNegative = y[s];
                }
                st.maxDnaDrift = Math.max(st.maxDnaDrift,
                        Math.abs((y[G] + y[GC]) - DNA_TOT_NM) / DNA_TOT_NM);
                t += h;
                st.accepted++;
                st.hMin = Math.min(st.hMin, h);
                st.hMax = Math.max(st.hMax, h);
            } else {
                st.rejected++;
            }

            double factor = (err > 0.0) ? 0.9 * Math.pow(err, -0.2) : 5.0;
            factor = Math.max(0.2, Math.min(5.0, factor));
            h = Math.max(H_MIN, Math.min(H_MAX, h * factor));
        }
    }

    /** Sweep grid: fine where the switch lives, coarse above. */
    static double[] sweepGrid() {
        java.util.ArrayList<Double> v = new java.util.ArrayList<Double>();
        for (double c = 0.0; c <= 25.0 + 1e-9; c += 0.25) v.add(c);
        for (double c = 26.0; c <= 100.0 + 1e-9; c += 2.0) v.add(c);
        double[] out = new double[v.size()];
        for (int k = 0; k < out.length; k++) out[k] = v.get(k);
        return out;
    }

    static final double SS_MINUTES = 6000.0;   // 100 h, Weber's steady state
    static final double BISTABLE_RATIO = 5.0;  // down/up GFP ratio marking bistability

    /** Continuation sweep up then down. Returns {lo, hi, gfpLow, gfpHigh}. */
    static double[] hysteresis(double kA, String tag, boolean write, Stats st) {
        double[] grid = sweepGrid();
        double[] up = new double[grid.length];
        double[] down = new double[grid.length];

        double[] y = uninduced();
        for (int k = 0; k < grid.length; k++) {
            integrate(y, grid[k], kA, SS_MINUTES, st);
            up[k] = y[I];
        }
        for (int k = grid.length - 1; k >= 0; k--) {
            integrate(y, grid[k], kA, SS_MINUTES, st);
            down[k] = y[I];
        }

        double lo = Double.NaN, hi = Double.NaN;
        for (int k = 0; k < grid.length; k++) {
            if (up[k] > 0.0 && down[k] / up[k] > BISTABLE_RATIO) {
                if (Double.isNaN(lo)) lo = grid[k];
                hi = grid[k];
            }
        }

        if (write) {
            try {
                new java.io.File(DIR).mkdirs();
                PrintStream ps = new PrintStream(DIR + "hysteresis_" + tag + ".csv", "UTF-8");
                ps.println("cA_star_nM;GFP_up_nM;GFP_down_nM;bistable;rtol;kA");
                for (int k = 0; k < grid.length; k++) {
                    boolean bi = up[k] > 0.0 && down[k] / up[k] > BISTABLE_RATIO;
                    ps.printf(Locale.US, "%.4f;%.9e;%.9e;%d;%.3e;%.4f%n",
                            grid[k], up[k], down[k], bi ? 1 : 0, rtol, kA);
                }
                ps.close();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }
        return new double[] {lo, hi, up[0], down[grid.length - 1]};
    }

    /** Time from uninduced to 95% of the 100 h GFP value at cA* = 100 nM. */
    static double settling(double kA, Stats st) {
        double[] ref = uninduced();
        integrate(ref, 100.0, kA, SS_MINUTES, st);
        double target = 0.95 * ref[I];

        double[] y = uninduced();
        double block = 0.5;                       // output resolution, min
        double hit = Double.NaN;
        try {
            new java.io.File(DIR).mkdirs();
            PrintStream ps = new PrintStream(DIR + "settling.csv", "UTF-8");
            ps.println("t_min;GFP_nM;target_nM;DNA_tot_nM;rtol");
            for (double t = 0.0; t < SS_MINUTES; t += block) {
                integrate(y, 100.0, kA, block, st);
                if (Double.isNaN(hit) && y[I] >= target) hit = t + block;
                if (t < 720.0 || Math.abs(t / 60.0 - Math.round(t / 60.0)) < block / 2.0) {
                    ps.printf(Locale.US, "%.4f;%.9e;%.9e;%.9e;%.3e%n",
                            t + block, y[I], target, y[G] + y[GC], rtol);
                }
            }
            ps.close();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        return hit;
    }

    /**
     * Increasing branch computed two ways, to decide what Weber's "steady state
     * concentration (induction time 100 hours)" actually specifies.
     *
     * <p>CONTINUATION carries the previous steady state forward, so it finds the
     * true saddle-node. FRESH restarts from the uninduced state at every c_A*
     * and holds 100 h; near a saddle-node, critical slowing down means the
     * switch may not complete within that window, so the apparent jump lands at
     * HIGHER c_A*. Which one Weber used is a reading of the paper, not a
     * parameter to tune.
     */
    static void protocolComparison(double kA, String tag) {
        double[] grid = sweepGrid();
        double[] cont = new double[grid.length];
        double[] fresh = new double[grid.length];
        Stats st = new Stats();

        double[] y = uninduced();
        for (int k = 0; k < grid.length; k++) {
            integrate(y, grid[k], kA, SS_MINUTES, st);
            cont[k] = y[I];
        }
        for (int k = 0; k < grid.length; k++) {
            double[] z = uninduced();
            integrate(z, grid[k], kA, SS_MINUTES, st);
            fresh[k] = z[I];
        }

        double ref = Math.min(cont[0], fresh[0]);
        double jc = Double.NaN, jf = Double.NaN;
        for (int k = 0; k < grid.length; k++) {
            if (Double.isNaN(jc) && cont[k] > 10.0 * ref) jc = grid[k];
            if (Double.isNaN(jf) && fresh[k] > 10.0 * ref) jf = grid[k];
        }
        System.out.printf(Locale.US,
                "  %-6s up-branch jump:  continuation %6.2f nM   fresh-IC %6.2f nM"
                + "   (Weber Fig. 3: ~15)%n", tag, jc, jf);

        try {
            new java.io.File(DIR).mkdirs();
            PrintStream ps = new PrintStream(DIR + "protocol_" + tag + ".csv", "UTF-8");
            ps.println("cA_star_nM;GFP_continuation_nM;GFP_freshIC_nM;rtol");
            for (int k = 0; k < grid.length; k++) {
                ps.printf(Locale.US, "%.4f;%.9e;%.9e;%.3e%n",
                        grid[k], cont[k], fresh[k], rtol);
            }
            ps.close();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    static void banner() {
        System.out.println("Job 5 - Weber & Buceta 2013 well-mixed LuxI/LuxR");
        System.out.println("  SPECIES CROSSING: V. fischeri lux on an E. coli chassis (ENGINEERING)");
        System.out.printf(Locale.US, "  tau = %.1f min (WEBER's, not chassis T_div)%n", TAU_MIN);
        System.out.printf(Locale.US, "  D = %.1f /min is a TRANSPORT RATE, not a diffusivity%n", D_MEMBRANE);
        System.out.printf(Locale.US, "  DNA_tot = %.4f nM (DERIVED: 1 copy per V_0 = %.1f um^3)%n",
                DNA_TOT_NM, V0_UM3);
        System.out.printf(Locale.US, "  V_c,tot = %.4e uL  V_ext = %.4e uL  r = %.6e%n",
                VC_TOT_UL, V_EXT_UL, R_VOL);
        System.out.println("  adaptive step-doubling RK4 + Richardson (NO QSSA reduction)");
    }

    public static void main(String[] args) {
        banner();
        boolean study = args.length > 0 && "tolstudy".equalsIgnoreCase(args[0]);
        if (args.length > 0 && "protocol".equalsIgnoreCase(args[0])) {
            System.out.println();
            System.out.println("Increasing-branch protocol comparison (down-branch unaffected)");
            protocolComparison(KA_LUX01, "lux01");
            protocolComparison(KA_LUX02, "lux02");
            return;
        }

        if (study) {
            System.out.println();
            System.out.println("GATE 7 tolerance study - bounds must move < 0.1 nM on halving");
            System.out.printf(Locale.US, "  %10s %10s %10s %10s %10s %12s %11s%n",
                    "rtol", "lux01_lo", "lux01_hi", "lux02_lo", "lux02_hi", "steps", "h_min");
            for (double tol : new double[] {1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8}) {
                rtol = tol;
                Stats st = new Stats();
                double[] a = hysteresis(KA_LUX01, "lux01", false, st);
                double[] b = hysteresis(KA_LUX02, "lux02", false, st);
                System.out.printf(Locale.US, "  %10.1e %10.2f %10.2f %10.2f %10.2f %12d %11.2e%n",
                        tol, a[0], a[1], b[0], b[1], st.accepted, st.hMin);
            }
            return;
        }

        System.out.printf(Locale.US, "%n  rtol = %.1e (ENGINEERING, Gate 7)%n", rtol);
        Stats sa = new Stats();
        double[] a = hysteresis(KA_LUX01, "lux01", true, sa);
        System.out.printf(Locale.US,
                "  lux01: bistable [%.2f, %.2f] nM   GFP low %.4e high %.4e   ratio %.1f%n",
                a[0], a[1], a[2], a[3], a[3] / a[2]);
        Stats sb = new Stats();
        double[] b = hysteresis(KA_LUX02, "lux02", true, sb);
        System.out.printf(Locale.US,
                "  lux02: bistable [%.2f, %.2f] nM   GFP low %.4e high %.4e   ratio %.1f%n",
                b[0], b[1], b[2], b[3], b[3] / b[2]);
        Stats sc = new Stats();
        double ts = settling(KA_LUX02, sc);
        System.out.printf(Locale.US, "  settling to 95%% at cA*=100 nM: %.1f min%n", ts);

        double worst = Math.min(sa.worstNegative, Math.min(sb.worstNegative, sc.worstNegative));
        double drift = Math.max(sa.maxDnaDrift, Math.max(sb.maxDnaDrift, sc.maxDnaDrift));
        System.out.printf(Locale.US,
                "  worst negative %.3e   max DNA drift %.3e   steps %d (%d rejected)%n",
                worst, drift, sa.accepted + sb.accepted + sc.accepted,
                sa.rejected + sb.rejected + sc.rejected);

        try {
            new java.io.File(DIR).mkdirs();
            PrintStream ps = new PrintStream(DIR + "summary.csv", "UTF-8");
            ps.println("construct;lo_nM;hi_nM;gfp_low_nM;gfp_high_nM;ratio;"
                    + "worst_negative;dna_drift_rel;settling_min;rtol;h_min;h_max");
            ps.printf(Locale.US, "lux01;%.4f;%.4f;%.9e;%.9e;%.4f;%.3e;%.3e;;%.3e;%.3e;%.3e%n",
                    a[0], a[1], a[2], a[3], a[3] / a[2], sa.worstNegative, sa.maxDnaDrift,
                    rtol, sa.hMin, sa.hMax);
            ps.printf(Locale.US, "lux02;%.4f;%.4f;%.9e;%.9e;%.4f;%.3e;%.3e;%.4f;%.3e;%.3e;%.3e%n",
                    b[0], b[1], b[2], b[3], b[3] / b[2], sb.worstNegative, sb.maxDnaDrift,
                    ts, rtol, sb.hMin, sb.hMax);
            ps.close();
            System.out.println("  wrote " + DIR + "summary.csv");
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}
