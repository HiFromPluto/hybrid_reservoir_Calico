package DaninoPocketP0;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;

import javax.vecmath.Vector3d;

import BacteriumFromScratch.ChassisParameters;
import BacteriumFromScratch.EcoliRodCell;
import BacteriumFromScratch.ValdezHertzian;

import bsim.BSim;
import bsim.BSimStepScheduler;
import bsim.BSimUtils;
import bsim.p0.DaninoPocketGeometry;

/**
 * Gate P0 — packed Danino-class pocket, chemistry OFF, N0 scheduler.
 *
 * <p>Not C0, not C1, not Fig. 4b, not Object A, not Object B. NOT_FIG4B.
 * Motility OFF. Death clamp OFF. Not a silent edit of ChassisPocket I0c.
 */
public final class DaninoPocketP0 {

    private DaninoPocketP0() {}

    static final double BX = DaninoPocketGeometry.POCKET_LX_UM;
    static final double BY = DaninoPocketGeometry.POCKET_LY_UM;
    static final double BZ = DaninoPocketGeometry.POCKET_LZ_UM;
    static final double T_END = 25950.0;
    static final double T_SMOKE = 15570.0;
    static final int SMOKE_N = 32;
    static final double SPILL_TOL_UM = 1.0e-3;
    static final double FILAMENT_R_RDISC = 2.0;
    static final double OFFPLANE_EXCURSION_UM = 0.40;
    static final double FOUNDER_X = 50.0;
    static final double FOUNDER_Y = 80.0;
    static final long SEED = ChassisParameters.RNG_SEED;
    static final int LOG_STRIDE = 10;

    private static final ValdezHertzian.WallSpec WALLS = ValdezHertzian.WallSpec.openYPlus();

    private static Vector3d centreOfMass(List<EcoliRodCell> cells) {
        Vector3d c = new Vector3d();
        for (EcoliRodCell cell : cells) {
            c.add(cell.centre());
        }
        if (!cells.isEmpty()) {
            c.scale(1.0 / cells.size());
        }
        return c;
    }

    private static double outOfPlane(List<EcoliRodCell> cells) {
        int out = 0;
        for (EcoliRodCell cell : cells) {
            if (Math.abs(cell.centre().z - BZ / 2.0) > OFFPLANE_EXCURSION_UM) {
                out++;
            }
        }
        return cells.isEmpty() ? 0.0 : out / (double) cells.size();
    }

    private static double maxZExcursion(List<EcoliRodCell> cells) {
        double m = 0.0;
        for (EcoliRodCell cell : cells) {
            m = Math.max(m, Math.abs(cell.centre().z - BZ / 2.0));
        }
        return m;
    }

    private static double discRadius(int n) {
        double footprint = ChassisParameters.W0_UM
                * (0.5 * (ChassisParameters.L0_UM + ChassisParameters.LDIV_UM)
                   + 2.0 * ChassisParameters.RADIUS_UM);
        return Math.sqrt(n * footprint / 0.85 / Math.PI);
    }

    private static double colonyRadius(List<EcoliRodCell> cells, Vector3d com) {
        double r = 0.0;
        for (EcoliRodCell cell : cells) {
            r = Math.max(r, Math.hypot(cell.x1.x - com.x, cell.x1.y - com.y));
            r = Math.max(r, Math.hypot(cell.x2.x - com.x, cell.x2.y - com.y));
        }
        return r;
    }

    private static double maxPoleY(List<EcoliRodCell> cells) {
        double y = Double.NEGATIVE_INFINITY;
        for (EcoliRodCell cell : cells) {
            y = Math.max(y, cell.x1.y);
            y = Math.max(y, cell.x2.y);
        }
        return cells.isEmpty() ? 0.0 : y;
    }

    static int doorContact(List<EcoliRodCell> cells) {
        double thresh = BY - ChassisParameters.RADIUS_UM;
        for (EcoliRodCell cell : cells) {
            if (cell.x1.y >= thresh || cell.x2.y >= thresh) {
                return 1;
            }
        }
        return 0;
    }

    static boolean centreLeftThroughYPlus(EcoliRodCell cell) {
        return cell.centre().y > BY + SPILL_TOL_UM;
    }

    static boolean centreWallLeak(EcoliRodCell cell) {
        Vector3d c = cell.centre();
        return c.x < -SPILL_TOL_UM || c.x > BX + SPILL_TOL_UM
                || c.y < -SPILL_TOL_UM
                || c.z < -SPILL_TOL_UM || c.z > BZ + SPILL_TOL_UM;
    }

    static int removeSpilled(List<EcoliRodCell> cells, int[] wallLeak) {
        int spilled = 0;
        Iterator<EcoliRodCell> it = cells.iterator();
        while (it.hasNext()) {
            EcoliRodCell cell = it.next();
            if (centreWallLeak(cell)) {
                wallLeak[0]++;
                it.remove();
            } else if (centreLeftThroughYPlus(cell)) {
                spilled++;
                it.remove();
            }
        }
        return spilled;
    }

    private static void reportMorphology(List<EcoliRodCell> cells, Vector3d bound,
            double t, int spillCum, int wallLeak) {
        Vector3d com = centreOfMass(cells);
        double R = colonyRadius(cells, com);
        int n = cells.size();
        double rDisc = discRadius(n);
        ValdezHertzian.PackingStats pst = ValdezHertzian.stats(cells, bound, WALLS);
        System.out.printf(Locale.US,
                "  t=%7.0f  N=%4d  R=%6.2f  R/R_disc=%5.2f  offplane=%.3f  z_exc=%.3f  "
                        + "d_min=%.3f  dcc=%.4f  spill_cum=%d  door=%d  y_max=%.2f  wall_leak=%d%n",
                t, n, R, rDisc > 0.0 ? R / rDisc : 0.0, outOfPlane(cells),
                maxZExcursion(cells), pst.minCentreDist, pst.maxDeltaCc, spillCum,
                doorContact(cells), maxPoleY(cells), wallLeak);
    }

    private static void writeUtf8(Path path, String text) throws IOException {
        Files.createDirectories(path.getParent());
        Files.write(path, text.getBytes(StandardCharsets.UTF_8));
    }

    static void run(String mode, Path exportDir) throws IOException {
        boolean smoke = mode.startsWith("smoke");
        double tEnd = smoke ? T_SMOKE : T_END;
        BSimUtils.generateDirectoryPath(exportDir.toString());

        boolean[][][] fluid = DaninoPocketGeometry.buildFluidMask();
        writeUtf8(exportDir.resolve("geometry_mask.json"),
                DaninoPocketGeometry.toJson(fluid));

        String freeze = String.format(Locale.US,
                "{\n"
                + "  \"frozen_before_traces\": true,\n"
                + "  \"gate\": \"P0\",\n"
                + "  \"chemistry\": \"OFF\",\n"
                + "  \"motility\": \"OFF\",\n"
                + "  \"object_b\": \"OFF\",\n"
                + "  \"not_fig4b\": true,\n"
                + "  \"object_a\": \"FAIL\",\n"
                + "  \"Lx_um\": %.2f,\n"
                + "  \"Ly_um\": %.2f,\n"
                + "  \"Lz_um\": %.4f,\n"
                + "  \"open_edge\": \"+y\",\n"
                + "  \"bus\": \"particle_sink_only\",\n"
                + "  \"scheduler\": \"BSimStepScheduler\",\n"
                + "  \"dt_s\": %.1f,\n"
                + "  \"T_s\": %.1f,\n"
                + "  \"N_steps\": %d,\n"
                + "  \"log_dt_s\": %.1f,\n"
                + "  \"mechanics_period_s\": %.1f,\n"
                + "  \"seed\": %d,\n"
                + "  \"founder_count\": 1,\n"
                + "  \"founder_um\": [%.2f, %.2f, %.4f],\n"
                + "  \"death_clamp\": \"OFF\",\n"
                + "  \"height_class\": \"Danino_SI_bulk_trap_1.65um\"\n"
                + "}\n",
                BX, BY, BZ, ChassisParameters.DT_S, tEnd,
                (int) Math.round(tEnd / ChassisParameters.DT_S),
                ChassisParameters.LOG_DT_S, ChassisParameters.DT_S, SEED,
                FOUNDER_X, FOUNDER_Y, BZ / 2.0);
        writeUtf8(exportDir.resolve("protocol_freeze.json"), freeze);

        BSim sim = new BSim();
        sim.setDt(ChassisParameters.DT_S);
        sim.setSimulationTime(tEnd);
        sim.setTimeFormat("0.00");
        sim.setBound(BX, BY, BZ);
        sim.setSolid(true, true, true);
        sim.setRandomSeed(SEED);

        double halfL = ChassisParameters.L0_UM / 2.0;
        final List<EcoliRodCell> cells = new ArrayList<EcoliRodCell>();
        EcoliRodCell founder = new EcoliRodCell(sim,
                new Vector3d(FOUNDER_X - halfL, FOUNDER_Y, BZ / 2.0),
                new Vector3d(FOUNDER_X + halfL, FOUNDER_Y, BZ / 2.0));
        founder.divisionMode = EcoliRodCell.DivisionMode.SYMMETRY_BROKEN;
        founder.divisionRng = sim.getRandom().asJavaRandom();
        cells.add(founder);

        final Vector3d bound = sim.getBound();
        final double bathMm = ChassisParameters.NUTRIENT_BATH_MM;
        final int[] nextReport = {2};
        final boolean[] filament = {false};
        final int[] spillCum = {0};
        final int[] spillSinceLog = {0};
        final int[] wallLeak = {0};
        final int[] nAtLastLog = {1};
        final boolean[] mechanicsDone = {false};

        Path csvPath = exportDir.resolve("colony_timeseries.csv");
        final BufferedWriter csv = Files.newBufferedWriter(csvPath, StandardCharsets.UTF_8);

        System.out.println("P0 packed device -> " + csvPath.toString().replace('\\', '/'));
        System.out.println("  chemistry=OFF  motility=OFF  Object B=OFF  NOT_FIG4B");
        System.out.println("  Lx=" + BX + "  Ly=" + BY + "  Lz=" + BZ
                + " um  open=+y  bus=particle_sink_only");
        System.out.println("  scheduler=BSimStepScheduler  seed=" + SEED
                + "  founder=(" + FOUNDER_X + ", " + FOUNDER_Y + ", " + (BZ / 2.0) + ")");
        System.out.println("  height=Danino_SI_bulk_trap_1.65um  not HybridDish 10 um"
                + "  not ChassisPocket 1 um");
        System.out.println("  division=SYMMETRY_BROKEN  Hertzian Job 3  death=OFF"
                + "  AHL/Hill/NARMA/AC/LuxI=OFF");
        System.out.println("  mode=" + (smoke ? "SMOKE N=" + SMOKE_N + " T=" + (int) tEnd
                : "FULL T=" + (int) tEnd) + " s  N_steps="
                + (int) Math.round(tEnd / ChassisParameters.DT_S));

        BSimStepScheduler scheduler = new BSimStepScheduler(sim);
        scheduler.addMechanicsEvent("growth-pack-divide-spill", sim.getDt(), context -> {
            if (smoke && mechanicsDone[0]) {
                return;
            }
            for (EcoliRodCell cell : cells) {
                cell.elongateCited(sim.getDt(), bathMm);
            }
            ValdezHertzian.relaxContacts(cells, bound, WALLS);
            spillSinceLog[0] += removeSpilled(cells, wallLeak);
            List<EcoliRodCell> newborns = new ArrayList<EcoliRodCell>();
            for (EcoliRodCell cell : cells) {
                if (cell.shouldDivide()) {
                    newborns.add(cell.divideCited());
                }
            }
            cells.addAll(newborns);
            if (!newborns.isEmpty()) {
                ValdezHertzian.relaxContacts(cells, bound, WALLS);
                spillSinceLog[0] += removeSpilled(cells, wallLeak);
            }
            if (cells.size() >= nextReport[0]) {
                reportMorphology(cells, bound, context.getEndTime(),
                        spillCum[0] + spillSinceLog[0], wallLeak[0]);
                double rDisc = discRadius(cells.size());
                double R = colonyRadius(cells, centreOfMass(cells));
                double rr = rDisc > 0.0 ? R / rDisc : 0.0;
                if (rr > FILAMENT_R_RDISC) {
                    filament[0] = true;
                    System.out.println("FILAMENT: R/R_disc=" + rr
                            + "  check SYMMETRY_BROKEN + injected BSimRandom seed 101.");
                }
                nextReport[0] *= 2;
            }
            if (smoke && cells.size() >= SMOKE_N) {
                mechanicsDone[0] = true;
            }
        });

        try {
            csv.write("t_s;seed;N;R_um;R_disc_um;R_Rdisc;offplane;z_exc_max_um;"
                    + "delta_cc_max_um;d_centers_min_um;spill_tick;spill_cum;"
                    + "N_ever;door_contact;y_max_um;wall_leak;n_drop_unexplained");
            csv.newLine();

            scheduler.run(new BSimStepScheduler.Adapter() {
                @Override
                public void observe(BSimStepScheduler.Context context) {
                    boolean log = context.isInitialObservation()
                            || context.getCompletedUpdates() % LOG_STRIDE == 0;
                    if (!log) {
                        return;
                    }
                    spillCum[0] += spillSinceLog[0];
                    int tickSpill = spillSinceLog[0];
                    spillSinceLog[0] = 0;
                    int n = cells.size();
                    int unexplained = Math.max(0, nAtLastLog[0] - tickSpill - n);
                    nAtLastLog[0] = n;
                    ValdezHertzian.PackingStats pst = ValdezHertzian.stats(cells, bound, WALLS);
                    Vector3d com = centreOfMass(cells);
                    double R = colonyRadius(cells, com);
                    double rDisc = discRadius(n);
                    double t = context.isInitialObservation() ? 0.0 : context.getEndTime();
                    try {
                        csv.write(String.format(Locale.US,
                                "%.1f;%d;%d;%.6f;%.6f;%.6f;%.6f;%.6f;%.9e;%.9f;%d;%d;%d;%d;%.6f;%d;%d",
                                t, SEED, n, R, rDisc, rDisc > 0.0 ? R / rDisc : 0.0,
                                outOfPlane(cells), maxZExcursion(cells),
                                pst.maxDeltaCc, pst.minCentreDist, tickSpill, spillCum[0],
                                n + spillCum[0], doorContact(cells), maxPoleY(cells),
                                wallLeak[0], unexplained));
                        csv.newLine();
                        csv.flush();
                    } catch (IOException e) {
                        throw new IllegalStateException("P0 CSV write failed", e);
                    }
                }
            });
        } finally {
            csv.close();
        }

        spillCum[0] += spillSinceLog[0];
        System.out.println("P0 finished  t=" + (int) tEnd + "  N=" + cells.size()
                + "  spill_cum=" + spillCum[0] + "  door=" + doorContact(cells)
                + "  wall_leak=" + wallLeak[0]);
        reportMorphology(cells, bound, tEnd, spillCum[0], wallLeak[0]);

        String summary = String.format(Locale.US,
                "{\n"
                + "  \"gate\": \"P0\",\n"
                + "  \"chemistry\": \"OFF\",\n"
                + "  \"motility\": \"OFF\",\n"
                + "  \"object_b\": \"OFF\",\n"
                + "  \"not_fig4b\": true,\n"
                + "  \"scheduler\": \"BSimStepScheduler\",\n"
                + "  \"seed\": %d,\n"
                + "  \"Lx_um\": %.2f, \"Ly_um\": %.2f, \"Lz_um\": %.4f,\n"
                + "  \"N_end\": %d,\n"
                + "  \"N_ever\": %d,\n"
                + "  \"spill_cum\": %d,\n"
                + "  \"wall_leak\": %d,\n"
                + "  \"door_contact\": %d,\n"
                + "  \"filament\": %s,\n"
                + "  \"mode\": \"%s\"\n"
                + "}\n",
                SEED, BX, BY, BZ, cells.size(), cells.size() + spillCum[0],
                spillCum[0], wallLeak[0], doorContact(cells),
                filament[0] ? "true" : "false", smoke ? "smoke" : "full");
        writeUtf8(exportDir.resolve("p0_summary.json"), summary);

        if (filament[0]) {
            throw new IllegalStateException("P0 filament at N=" + cells.size());
        }
        if (smoke) {
            Vector3d com = centreOfMass(cells);
            double R = colonyRadius(cells, com);
            double rr = R / discRadius(Math.max(1, cells.size()));
            double off = outOfPlane(cells);
            if (cells.size() < SMOKE_N) {
                throw new IllegalStateException("P0 smoke FAIL  N=" + cells.size()
                        + " < " + SMOKE_N);
            }
            if (off > 0.05 || rr > FILAMENT_R_RDISC) {
                throw new IllegalStateException("P0 smoke FAIL  offplane=" + off
                        + "  R/R_disc=" + rr);
            }
            System.out.printf(Locale.US,
                    "P0 smoke OK  N=%d  offplane=%.3f  R/R_disc=%.2f  spill_cum=%d  y_max=%.2f%n",
                    cells.size(), off, rr, spillCum[0], maxPoleY(cells));
        }
    }

    public static void main(String[] args) throws IOException {
        String mode = (args.length > 0 && args[0] != null && !args[0].trim().isEmpty())
                ? args[0].trim() : "full";
        Path exportDir;
        if ("smokeA".equalsIgnoreCase(mode) || "smoke".equalsIgnoreCase(mode)) {
            exportDir = Paths.get("results", "p0_smoke_a");
            mode = "smokeA";
        } else if ("smokeB".equalsIgnoreCase(mode)) {
            exportDir = Paths.get("results", "p0_smoke_b");
            mode = "smokeB";
        } else {
            exportDir = Paths.get("results", "p0_seed101");
            mode = "full";
        }
        ChassisParameters.printLedger();
        run(mode, exportDir);
    }
}
