package BacteriumFromScratch;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import javax.vecmath.Vector3d;

import bsim.BSim;

/**
 * Job 3c step 3 diagnostic: is the consumed volume dx-dependent?
 *
 * <p>markCapsule dilates the capsule by +0.5*dx before testing voxel centres.
 * That dilation is a fixed FRACTION OF THE GRID, not of the cell, so the
 * marked volume -- and hence total Monod consumption -- may scale with dx.
 * If so, the source term itself changes under grid refinement and no
 * refinement study can converge.
 */
public final class Job3cOcc {

    private Job3cOcc() {}

    public static void main(String[] args) {
        double bx = 20.0, by = 20.0, bz = 4.0;
        BSim s = new BSim();
        s.setDt(ChassisParameters.DT_S);
        s.setSimulationTime(10.0);
        s.setBound(bx, by, bz);
        s.setSolid(true, true, true);

        double r = ChassisParameters.RADIUS_UM;
        double L = ChassisParameters.L0_UM;
        double halfL = L / 2.0;
        List<EcoliRodCell> one = new ArrayList<EcoliRodCell>();
        one.add(new EcoliRodCell(s,
                new Vector3d(bx / 2.0 - halfL, by / 2.0, bz / 2.0),
                new Vector3d(bx / 2.0 + halfL, by / 2.0, bz / 2.0)));

        // true sphero-cylinder volume: cylinder + sphere
        double vTrue = Math.PI * r * r * L + (4.0 / 3.0) * Math.PI * r * r * r;
        System.out.printf(Locale.US,
                "single cell  L = %.1f um, r = %.1f um   true capsule volume = %.4f um^3%n%n",
                L, r, vTrue);
        System.out.printf(Locale.US, "  %7s %14s %10s %14s %10s%n",
                "dx", "legacy_um3", "ratio", "conserv_um3", "ratio");

        for (double dx : new double[] {2.0, 1.0, 0.5, 0.25, 0.125}) {
            int pad = Math.max(1, (int) Math.round(1.0 / dx));
            NutrientField f = new NutrientField(bx, by, bz, false, true, dx, pad,
                    NutrientField.BcMode.DISH_LATERAL);
            f.markOccupancy(one);
            double vLegacy = f.markedMassUm3();
            NutrientField g = new NutrientField(bx, by, bz, false, true, dx, pad,
                    NutrientField.BcMode.DISH_LATERAL);
            g.setOccupancyMode(NutrientField.OccupancyMode.CONSERVATIVE_MASS);
            g.markOccupancy(one);
            double vCons = g.markedMassUm3();
            System.out.printf(Locale.US, "  %7.3f %14.4f %10.3f %14.4f %10.3f%n",
                    dx, vLegacy, vLegacy / vTrue, vCons, vCons / vTrue);
        }
        System.out.println();
        System.out.println("  ratio == 1.0 would mean consumption is grid independent.");
    }
}
