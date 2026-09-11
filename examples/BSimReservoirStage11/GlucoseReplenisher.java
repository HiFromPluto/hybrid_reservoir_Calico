package BSimReservoirStage11;

import bsim.BSimChemicalField;

/**
 * Model A glucose replenishment: dG/dt = k_supply * (G0 - G) per voxel.
 * Ported from reservoir_new (AC_OurReservoir_Scenario.GlucoseReplenisher).
 */
public class GlucoseReplenisher {

    private final double G0;
    private final double k_supply;

    public GlucoseReplenisher(double G0, double k_supply) {
        this.G0       = G0;
        this.k_supply = k_supply;
    }

    public void apply(BSimChemicalField glucoseField, double dt) {
        int[]    boxes  = glucoseField.getBoxes();
        double[] box    = glucoseField.getBox();
        double   volume = box[0] * box[1] * box[2];

        for (int i = 0; i < boxes[0]; i++)
            for (int j = 0; j < boxes[1]; j++)
                for (int k = 0; k < boxes[2]; k++) {
                    double G  = glucoseField.getConc(i, j, k);
                    double dG = k_supply * (G0 - G);
                    glucoseField.addQuantity(i, j, k, dG * volume * dt);
                }
    }
}
