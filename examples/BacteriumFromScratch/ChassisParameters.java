package BacteriumFromScratch;

/**
 * Frozen chassis numbers. Every value is TAKEN from a BUILD paper or
 * labelled ENGINEERING with a named objective that is not "look alive".
 * Job 2 growth numbers are locked. Job 3 adds Hertzian packing only.
 */
public final class ChassisParameters {

    private ChassisParameters() {}

    /** Warren 2019 glucose minimal medium, batch culture. TAKEN. h^{-1}. */
    public static final double LAMBDA_S_PER_HOUR = 1.0;

    /** Same ν as λ_S. Valdez eq. (5). TAKEN. s^{-1}. */
    public static final double NU_PER_SECOND = LAMBDA_S_PER_HOUR / 3600.0;

    /** Warren K_S for glucose = 0.02 mM = 20 µM. TAKEN. mM. */
    public static final double KS_MM = 0.02;

    /** Warren saturating substrate C_s. TAKEN. mM. */
    public static final double NUTRIENT_BATH_MM = 0.5;

    /** Warren cell diameter w_0. TAKEN. µm. */
    public static final double W0_UM = 1.0;

    /** Hemisphere radius = w_0 / 2. TAKEN. µm. */
    public static final double RADIUS_UM = W0_UM / 2.0;

    /** Warren dividing cylindrical length. TAKEN. µm. */
    public static final double LDIV_UM = 3.0;

    /** Warren ℓ_div = 2 ℓ_0 + w_0. TAKEN. µm. */
    public static final double L0_UM = (LDIV_UM - W0_UM) / 2.0;

    /** Integrator step. ENGINEERING: growth timescale is hours. s. */
    public static final double DT_S = 1.0;

    /** CSV log step. ENGINEERING. s. */
    public static final double LOG_DT_S = 10.0;

    /** Isolated-cell export length. ENGINEERING. s. */
    public static final double SIM_TIME_S = 9000.0;

    /** Recorded seed. ODE is deterministic; unused by growth. ENGINEERING. */
    public static final long RNG_SEED = 101L;

    // ---- Job 3 Hertzian (Valdez 2025 eqs. 6–10; Warren Table S4) ----

    /**
     * Cell–cell Hertzian stiffness k_cc. TAKEN.
     * RodCellVSPhage Constants.cpp: "Mechanical parameters (Table S4 Warren2019.pdf)"
     * k_cc = 3.0e4. Do not substitute T7Plaque.txt (k_cc = 1e6) and do not
     * raise this after seeing overlap.
     */
    public static final double K_CC = 3.0e4;

    /**
     * Cell–agar / cell–wall stiffness k_ac. TAKEN. Same table, k_wc = 3.0e4.
     */
    public static final double K_AC = 3.0e4;

    /**
     * Dynamic friction μ_cc. TAKEN. Warren Fig. 9 standard set μ_cc = 0.1
     * and RodCellVSPhage cell_mu = 1.0e-1.
     */
    public static final double MU_CC = 0.1;

    /**
     * Dynamic friction μ_ac. TAKEN. Warren Fig. 9 μ_ca = 0.8
     * and RodCellVSPhage wall_mu = 8.0e-1.
     */
    public static final double MU_AC = 0.8;

    /**
     * Tangential dissipation γ_{cc,τ}. TAKEN. Warren Fig. 9
     * γ_cc,t = 10000 μm^{-1} h^{-1}; RodCellVSPhage gamma_t = 1.0e4.
     */
    public static final double GAMMA_T_PER_UM_PER_H = 1.0e4;

    /**
     * Normal dissipation γ_{cc,n}. TAKEN. RodCellVSPhage gamma_n = 5.0e2
     * citing Table S4.
     */
    public static final double GAMMA_N_PER_UM_PER_H = 5.0e2;

    /**
     * Water viscosity μ_liq. TAKEN. RodCellVSPhage viscosity = 1.00160e-3
     * (Pa s); Valdez eq. (12) ~ 1 mPa s.
     */
    public static final double MU_LIQ_PA_S = 1.00160e-3;

    /**
     * RodCellVSPhage Integrate.cpp: eta = 1e6 * viscosity before Tirado
     * drag (μm, hours). TAKEN from the BUILD integrator that implements
     * Valdez. Without it, water-scale η makes 1 s growth steps explode;
     * Job 3 keeps this factor and does not raise k_cc.
     */
    public static final double RODCELL_ETA_UNIT_FACTOR = 1.0e6;

    /**
     * Max pole displacement per mechanical substep. ENGINEERING CFL so
     * that water-scale Tirado mobility stays stable at Job 2 Δt = 1 s.
     */
    public static final double MECH_CFL_UM = 0.02;

    /** Cap on CFL substeps per growth step. ENGINEERING. */
    public static final int MECH_SUBSTEP_CAP = 20000;

    /**
     * Two-body seed overlap δ_cc(0). ENGINEERING: known overlap > 0 so
     * the Hertzian must move the capsules.
     */
    public static final double TWO_BODY_DELTA0_UM = 0.20;

    /** Two-body box. ENGINEERING: walls far from the pair. μm. */
    public static final double TWO_BODY_BOUND_X = 40.0;
    public static final double TWO_BODY_BOUND_Y = 40.0;
    public static final double TWO_BODY_BOUND_Z = 10.0;

    /** Two-body export. ENGINEERING: ~2 h at the paper's μm/h mobility. s. */
    public static final double TWO_BODY_TIME_S = 7200.0;

    /**
     * Allowed leftover overlap after two-body relaxation (surface gap
     * ≥ −this). ENGINEERING residual, not a retuned spring.
     */
    public static final double TWO_BODY_GAP_RESIDUAL_UM = 0.02;

    /**
     * Minimum centre displacement that counts as "positions moved".
     * ENGINEERING. μm.
     */
    public static final double TWO_BODY_MIN_DISP_UM = 0.05;

    /**
     * Growing-cluster closed box. ENGINEERING, named: growth must push.
     * μm.
     */
    public static final double CLUSTER_BOUND_X = 8.0;
    public static final double CLUSTER_BOUND_Y = 6.0;
    public static final double CLUSTER_BOUND_Z = 4.0;

    /**
     * Cluster export. ENGINEERING: two Job-2 division times plus margin.
     * 2 × T_div ≈ 5190 s. s.
     */
    public static final double CLUSTER_TIME_S = 6000.0;

    /**
     * Max pairwise overlap δ_cc during the cluster. ENGINEERING cap.
     * Frozen before looking at overlap; not a raised k_cc.
     */
    public static final double CLUSTER_OVERLAP_CAP_UM = 0.25;

    /**
     * Initial pairwise overlap of the two seed rods. ENGINEERING.
     * Known δ_cc > 0 so the cluster starts in contact.
     */
    public static final double CLUSTER_SEED_DELTA_UM = 0.15;

    /** Isolated packing check. ENGINEERING. s. */
    public static final double ISOLATED_PACK_TIME_S = 100.0;

    /**
     * Packing-force floor for an isolated cell. ENGINEERING: numerical
     * zero. Native force units of Valdez (8).
     */
    public static final double ISOLATED_F_PACK_MAX = 1.0e-8;

    // ---- Job 3b nutrient field (Valdez 16–17; Warren A1.2.2) ----

    /**
     * Colony diffusivity D_C. TAKEN. Warren Table S2 / RodCellVSPhage
     * DiffColony scaling (90 µm²/s in SI). DOI 10.7554/eLife.41093.
     */
    public static final double D_C_UM2_PER_S = 90.0;

    /**
     * Agar diffusivity D_A. TAKEN. Warren Table S2 (600 µm²/s in SI).
     */
    public static final double D_A_UM2_PER_S = 600.0;

    /**
     * Glucose yield Y. TAKEN. Warren Table S2 (0.5 gCDW/g glucose).
     */
    public static final double YIELD_GCDW_PER_GGLUCOSE = 0.5;

    /**
     * Single-cell dry mass density ρ_cell. TAKEN. Warren Table S3
     * (0.137×10⁻¹² gCDW/µm³).
     */
    public static final double RHO_CELL_GCDW_PER_UM3 = 0.137e-12;

    /**
     * Occupied-voxel mass density ρ₀ ≈ 0.68 ρ_cell. TAKEN. Warren
     * colony-density approximation; constant per occupied voxel v1.
     */
    public static final double RHO_0_GCDW_PER_UM3 = 0.68 * RHO_CELL_GCDW_PER_UM3;

    /**
     * g glucose per µm³ per mM (180 g/mol; 1 mM = 0.18 g/L; 1 L = 1e15 µm³).
     * TAKEN unit conversion — 1.8e-16 g/(µm³·mM).
     */
    public static final double MM_TO_G_PER_UM3 = 180.0 / 1.0e15 * 1.0e-3;

    /** Steady-state flux balance tolerance for Job 3b gate. ENGINEERING. */
    public static final double NUTRIENT_FLUX_TOL_REL = 0.005;

    /**
     * Consecutive GS sweeps with maxDelta &lt; tol required before accepting
     * quasi-steady. ENGINEERING: avoids single-sweep false convergence.
     */
    public static final int NUTRIENT_GS_STABLE_SWEEPS = 2;

    /** Depletion flux gate ignores rows with t &lt; this. ENGINEERING. s. */
    public static final double DEPLETION_FLUX_WARMUP_S = 500.0;

    /**
     * Division symmetry breaking. ENGINEERING.
     *
     * <p>Mechanism cited: Melke et al. 2010 Methods -- "At each cell division
     * we introduce some randomness in order to break the axial symmetry of the
     * system, giving two daughter cells with slightly different sizes and
     * imperfect alignment" -- citing Cho et al. 2007 (PLoS Biol 5:e302).
     * Melke's model is non-dimensional (see CITATION_DOSSIER erratum), so the
     * MECHANISM is TAKEN but these MAGNITUDES are ENGINEERING and must be
     * justified by a morphology sweep, not assumed.
     *
     * <p>Without this, division is exactly collinear and a colony seeded from
     * one founder in an open chamber grows as a 1-D filament: R ~ N instead of
     * R ~ sqrt(N), and lateral Hertzian contact never occurs.
     */
    public static final double DIVISION_ANGLE_SD_RAD = 0.05;

    /** Fractional daughter length asymmetry at division. ENGINEERING. */
    public static final double DIVISION_LENGTH_ASYM = 0.05;

    /** Nutrient grid spacing. ENGINEERING: resolve CLUSTER_BOX. µm. */
    public static final double NUTRIENT_DX_UM = 0.5;

    /** Agar padding layers at i=j=k=0. ENGINEERING. voxels. */
    public static final int NUTRIENT_AGAR_PAD = 2;

    /** Quasi-steady GS tolerance. ENGINEERING. mM. */
    public static final double NUTRIENT_CONV_TOL_MM = 1.0e-9;

    /** Max GS iterations per tick. ENGINEERING. */
    public static final int NUTRIENT_MAX_ITER = 8000;

    /** Job 3b uniform regression box (same as Job 2). ENGINEERING. µm. */
    public static final double UNIFORM_REG_BOUND_X = 100.0;
    public static final double UNIFORM_REG_BOUND_Y = 100.0;
    public static final double UNIFORM_REG_BOUND_Z = 10.0;

    /** Job 3b depletion / coupled export. ENGINEERING. s. */
    public static final double DEPLETION_TIME_S = 6000.0;

    public static double monodFactor(double nutrientMm) {
        double n = Math.max(0.0, nutrientMm);
        return n / (n + KS_MM);
    }

    /**
     * Valdez eq. (5): σ := ν log2(ℓ_div / ℓ_0). TAKEN. s^{-1}.
     */
    public static double sigmaPerSecond() {
        return NU_PER_SECOND * (Math.log(LDIV_UM / L0_UM) / Math.log(2.0));
    }

    /**
     * Time for ℓ_0 → ℓ_div at constant bath N. Valdez (4)–(5). s.
     */
    public static double analyticDivisionTimeS(double nutrientMm) {
        return Math.log(2.0) / (NU_PER_SECOND * monodFactor(nutrientMm));
    }

    public static double analyticLengthUm(double timeSinceBirthS, double nutrientMm) {
        double alpha = sigmaPerSecond() * monodFactor(nutrientMm);
        return L0_UM * Math.exp(alpha * timeSinceBirthS);
    }

    /** Mechanical dt in hours (RodCellVSPhage / Valdez native time). */
    public static double dtHours(double dtS) {
        return dtS / 3600.0;
    }

    public static double etaScale() {
        return RODCELL_ETA_UNIT_FACTOR * MU_LIQ_PA_S;
    }

    /**
     * Warren A1.2.2 Monod sink in colony voxels: (λ_S/Y)ρ₀ N/(N+K_S), mM/s.
     * Valdez (16) identifies ν=λ_S, κ=K_S; Warren supplies Y and ρ₀.
     */
    public static double consumptionSinkMmPerS(double nMm, double rhoOccupied) {
        if (rhoOccupied <= 0.0) {
            return 0.0;
        }
        double n = Math.max(0.0, nMm);
        double lambdaSPerS = LAMBDA_S_PER_HOUR / 3600.0;
        double monod = monodFactor(n);
        double massSink = (lambdaSPerS / YIELD_GCDW_PER_GGLUCOSE) * rhoOccupied * monod;
        return massSink / MM_TO_G_PER_UM3;
    }

    public static void printLedger() {
        System.out.println("CHASSIS PARAMETER LEDGER");
        System.out.println("  TAKEN     lambda_S = " + LAMBDA_S_PER_HOUR
                + " /h   Warren 2019 glucose MM  DOI 10.7554/eLife.41093");
        System.out.println("  TAKEN     nu       = " + NU_PER_SECOND
                + " /s   = lambda_S / 3600");
        System.out.println("  TAKEN     K_S      = " + KS_MM
                + " mM    Warren (20 uM glucose)");
        System.out.println("  TAKEN     N_bath   = " + NUTRIENT_BATH_MM
                + " mM    Warren C_s");
        System.out.println("  TAKEN     w0       = " + W0_UM + " um");
        System.out.println("  TAKEN     radius   = " + RADIUS_UM + " um");
        System.out.println("  TAKEN     L_div    = " + LDIV_UM + " um");
        System.out.println("  TAKEN     L0       = " + L0_UM
                + " um    from L_div = 2 L0 + w0");
        System.out.println("  TAKEN     sigma    = " + sigmaPerSecond()
                + " /s   Valdez eq 5");
        System.out.println("  TAKEN     Monod    = " + monodFactor(NUTRIENT_BATH_MM)
                + "      N/(N+K_S) at C_s");
        System.out.println("  TAKEN     T_div    = " + analyticDivisionTimeS(NUTRIENT_BATH_MM)
                + " s    ln2 / (nu * Monod)");
        System.out.println("  TAKEN     k_cc     = " + K_CC
                + "      Warren Table S4 via RodCellVSPhage Constants.cpp");
        System.out.println("  TAKEN     k_ac     = " + K_AC
                + "      Warren Table S4 k_wc");
        System.out.println("  TAKEN     mu_cc    = " + MU_CC
                + "      Warren Fig 9 / Table S4");
        System.out.println("  TAKEN     mu_ac    = " + MU_AC
                + "      Warren Fig 9 mu_ca");
        System.out.println("  TAKEN     gamma_t  = " + GAMMA_T_PER_UM_PER_H
                + " /um/h  Warren Fig 9 gamma_cc,t");
        System.out.println("  TAKEN     gamma_n  = " + GAMMA_N_PER_UM_PER_H
                + " /um/h  RodCellVSPhage Table S4");
        System.out.println("  TAKEN     mu_liq   = " + MU_LIQ_PA_S
                + " Pa s  water; Valdez eq 12");
        System.out.println("  TAKEN     eta_1e6  = " + RODCELL_ETA_UNIT_FACTOR
                + "      RodCellVSPhage Integrate.cpp eta = 1e6 * viscosity");
        System.out.println("  ENGINEERING cfl    = " + MECH_CFL_UM
                + " um/substep  quasi-static contact solve");
        System.out.println("  TAKEN     Hertz    = Valdez (8)-(10) elastic; (6)-(7) overdamped");
        System.out.println("  ENGINEERING dt     = " + DT_S
                + " s     hour-scale growth");
        System.out.println("  ENGINEERING log_dt = " + LOG_DT_S + " s");
        System.out.println("  ENGINEERING T_end  = " + SIM_TIME_S + " s   Job 2");
        System.out.println("  ENGINEERING seed   = " + RNG_SEED
                + "      recorded; growth has no RNG");
        System.out.println("  ENGINEERING divide = L >= L_div  (not Valdez stochastic P)");
        System.out.println("  ENGINEERING L_ran  = 0  (Warren allows a small fluctuation)");
        System.out.println("  ENGINEERING box_2b = "
                + TWO_BODY_BOUND_X + "x" + TWO_BODY_BOUND_Y + "x"
                + TWO_BODY_BOUND_Z + " um");
        System.out.println("  ENGINEERING dcc0   = " + TWO_BODY_DELTA0_UM + " um");
        System.out.println("  ENGINEERING T_2b   = " + TWO_BODY_TIME_S + " s");
        System.out.println("  ENGINEERING gap_res= " + TWO_BODY_GAP_RESIDUAL_UM + " um");
        System.out.println("  ENGINEERING box_cl = CLUSTER_BOX "
                + CLUSTER_BOUND_X + "x" + CLUSTER_BOUND_Y + "x"
                + CLUSTER_BOUND_Z + " um  closed");
        System.out.println("  ENGINEERING T_cl   = " + CLUSTER_TIME_S + " s");
        System.out.println("  ENGINEERING dcc_cap= " + CLUSTER_OVERLAP_CAP_UM + " um");
        System.out.println("  ENGINEERING dcc_seed_cl= " + CLUSTER_SEED_DELTA_UM + " um");
        System.out.println("  ENGINEERING T_iso  = " + ISOLATED_PACK_TIME_S + " s");
        System.out.println("  TAKEN     D_C      = " + D_C_UM2_PER_S
                + " um^2/s  Warren Table S2 colony");
        System.out.println("  TAKEN     D_A      = " + D_A_UM2_PER_S
                + " um^2/s  Warren Table S2 agar");
        System.out.println("  TAKEN     Y        = " + YIELD_GCDW_PER_GGLUCOSE
                + " gCDW/g  Warren Table S2");
        System.out.println("  TAKEN     rho_cell = " + RHO_CELL_GCDW_PER_UM3
                + " gCDW/um^3  Warren Table S3");
        System.out.println("  TAKEN     rho_0    = " + RHO_0_GCDW_PER_UM3
                + " gCDW/um^3  Warren 0.68 rho_cell (v1 constant/voxel)");
        System.out.println("  TAKEN     sink     = (lambda_S/Y) rho_0 N/(N+K_S)  Warren A1.2.2");
        System.out.println("  TAKEN     MM_conv  = " + MM_TO_G_PER_UM3
                + " g/(um^3 mM)  glucose unit conversion");
        System.out.println("  ENGINEERING BC_n   = low faces i=j=k=0 Dirichlet C_s;"
                + " high faces Neumann; corner agar shell (not Warren z<0 slab)");
        System.out.println("  ENGINEERING flux_tol= " + NUTRIENT_FLUX_TOL_REL
                + "      steady-state inflow vs consumption gate");
        System.out.println("  ENGINEERING gs_stable= " + NUTRIENT_GS_STABLE_SWEEPS
                + " sweeps + flux residual  quasi-steady acceptance");
        System.out.println("  ENGINEERING flux_wu = " + DEPLETION_FLUX_WARMUP_S
                + " s     depletion gate warm-up");
        System.out.println("  TAKEN     N(x) solve= quasi-steady GS each tick  Warren A1.2.7");
        System.out.println("  ENGINEERING dx_n   = " + NUTRIENT_DX_UM + " um  nutrient grid");
        System.out.println("  ENGINEERING agar_pad= " + NUTRIENT_AGAR_PAD + " voxels");
        System.out.println("  ENGINEERING occ_vol= capsule voxels at rho_0  Warren v1 footprint");
        System.out.println("  ENGINEERING T_dep  = " + DEPLETION_TIME_S + " s   Job 3b depletion");
        System.out.println("  OFF phage, QS, chemotaxis ODE, metabolism ODE, Brownian, F_s, P(L)");
        System.out.println("  OFF RelaxationMover / BSimCapsuleBacterium.k_cell spring");
        System.out.println("  OFF uncited k_ov from other examples");
    }
}
