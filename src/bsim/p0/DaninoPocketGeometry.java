package bsim.p0;

import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Locale;

/**
 * Frozen Danino-class pocket + horizontal bus mask for P0 / later N0 transport.
 *
 * <p>P0 does not advance a chemical field. The mask is declared so C1 can
 * reuse the same voxels. Height is the SI bulk trap 1.65 µm, not HybridDish
 * 10 µm and not ChassisPocket 1 µm.
 */
public final class DaninoPocketGeometry {

    public static final double POCKET_LX_UM = 100.0;
    public static final double POCKET_LY_UM = 100.0;
    public static final double POCKET_LZ_UM = 1.65;

    public static final double BUS_LENGTH_UM = 400.0;
    public static final double BUS_WIDTH_UM = 80.0;

    public static final double DX_UM = 2.0;
    public static final double DY_UM = 2.0;
    public static final double DZ_UM = POCKET_LZ_UM;

    public static final double MASK_X0_UM = (POCKET_LX_UM - BUS_LENGTH_UM) / 2.0;
    public static final double MASK_Y0_UM = 0.0;
    public static final double MASK_Z0_UM = 0.0;

    public static final int NX = 200;
    public static final int NY = 90;
    public static final int NZ = 1;

    public static final String OPEN_EDGE = "+y";
    public static final String BUS_ROLE = "particle_sink_only";

    private DaninoPocketGeometry() {}

    public static boolean[][][] buildFluidMask() {
        boolean[][][] fluid = new boolean[NX][NY][NZ];
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                for (int k = 0; k < NZ; k++) {
                    fluid[i][j][k] = isFluidVoxel(i, j, k);
                }
            }
        }
        return fluid;
    }

    public static boolean isFluidVoxel(int i, int j, int k) {
        double cx = MASK_X0_UM + (i + 0.5) * DX_UM;
        double cy = MASK_Y0_UM + (j + 0.5) * DY_UM;
        double cz = MASK_Z0_UM + (k + 0.5) * DZ_UM;
        boolean inPocket = cx >= 0.0 && cx <= POCKET_LX_UM
                && cy >= 0.0 && cy <= POCKET_LY_UM
                && cz >= 0.0 && cz <= POCKET_LZ_UM;
        boolean inBus = cx >= MASK_X0_UM && cx <= MASK_X0_UM + BUS_LENGTH_UM
                && cy > POCKET_LY_UM && cy <= POCKET_LY_UM + BUS_WIDTH_UM
                && cz >= 0.0 && cz <= POCKET_LZ_UM;
        return inPocket || inBus;
    }

    public static int countFluid(boolean[][][] fluid) {
        int n = 0;
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                for (int k = 0; k < NZ; k++) {
                    if (fluid[i][j][k]) {
                        n++;
                    }
                }
            }
        }
        return n;
    }

    public static int countPocketFluid(boolean[][][] fluid) {
        int n = 0;
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                for (int k = 0; k < NZ; k++) {
                    if (!fluid[i][j][k]) {
                        continue;
                    }
                    double cy = MASK_Y0_UM + (j + 0.5) * DY_UM;
                    if (cy <= POCKET_LY_UM) {
                        n++;
                    }
                }
            }
        }
        return n;
    }

    public static String maskSha256(boolean[][][] fluid) {
        byte[] packed = new byte[NX * NY * NZ];
        int p = 0;
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                for (int k = 0; k < NZ; k++) {
                    packed[p++] = (byte) (fluid[i][j][k] ? 1 : 0);
                }
            }
        }
        try {
            MessageDigest md = MessageDigest.getInstance("SHA-256");
            byte[] digest = md.digest(packed);
            StringBuilder hex = new StringBuilder(digest.length * 2);
            for (byte b : digest) {
                hex.append(String.format(Locale.US, "%02x", b));
            }
            return hex.toString();
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 required", e);
        }
    }

    public static String toJson(boolean[][][] fluid) {
        int nFluid = countFluid(fluid);
        int nPocket = countPocketFluid(fluid);
        StringBuilder sb = new StringBuilder();
        sb.append("{\n");
        sb.append("  \"chemistry\": \"OFF\",\n");
        sb.append("  \"object_b\": \"OFF\",\n");
        sb.append("  \"not_fig4b\": true,\n");
        sb.append(String.format(Locale.US, "  \"pocket_um\": [%.2f, %.2f, %.4f],%n",
                POCKET_LX_UM, POCKET_LY_UM, POCKET_LZ_UM));
        sb.append("  \"open_edge\": \"").append(OPEN_EDGE).append("\",\n");
        sb.append("  \"bus_role\": \"").append(BUS_ROLE).append("\",\n");
        sb.append(String.format(Locale.US, "  \"bus_um\": [%.2f, %.2f, %.4f],%n",
                BUS_LENGTH_UM, BUS_WIDTH_UM, POCKET_LZ_UM));
        sb.append(String.format(Locale.US, "  \"mask_origin_um\": [%.2f, %.2f, %.2f],%n",
                MASK_X0_UM, MASK_Y0_UM, MASK_Z0_UM));
        sb.append(String.format(Locale.US, "  \"dx_um\": %.2f, \"dy_um\": %.2f, \"dz_um\": %.4f,%n",
                DX_UM, DY_UM, DZ_UM));
        sb.append(String.format(Locale.US, "  \"boxes\": [%d, %d, %d],%n", NX, NY, NZ));
        sb.append(String.format(Locale.US, "  \"n_fluid\": %d,%n", nFluid));
        sb.append(String.format(Locale.US, "  \"n_pocket_fluid\": %d,%n", nPocket));
        sb.append(String.format(Locale.US, "  \"n_bus_fluid\": %d,%n", nFluid - nPocket));
        sb.append("  \"mask_sha256\": \"").append(maskSha256(fluid)).append("\",\n");
        sb.append("  \"height_class\": \"Danino_SI_bulk_trap_1.65um\",\n");
        sb.append("  \"not_hybriddish_10um\": true,\n");
        sb.append("  \"not_chassispocket_1um\": true\n");
        sb.append("}\n");
        return sb.toString();
    }

    public static byte[] jsonBytes(boolean[][][] fluid) {
        return toJson(fluid).getBytes(StandardCharsets.UTF_8);
    }
}
