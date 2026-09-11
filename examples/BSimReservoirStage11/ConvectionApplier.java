package BSimReservoirStage11;

import bsim.BSimChemicalField;
import javax.vecmath.Vector3d;

/**
 * Upwind advection for chemical fields.
 * Ported from reservoir_new (AC_OurReservoir_Scenario.ConvectionApplier).
 */
public class ConvectionApplier {

    public void apply(BSimChemicalField field, Vector3d velocity, double dt,
                      double cInlet, boolean applyInletBC) {
        int[]    boxes   = field.getBoxes();
        double[] boxSize = field.getBox();

        double[][][] conc = new double[boxes[0]][boxes[1]][boxes[2]];
        for (int i = 0; i < boxes[0]; i++)
            for (int j = 0; j < boxes[1]; j++)
                for (int k = 0; k < boxes[2]; k++)
                    conc[i][j][k] = field.getConc(i, j, k);

        for (int i = 0; i < boxes[0]; i++) {
            for (int j = 0; j < boxes[1]; j++) {
                for (int k = 0; k < boxes[2]; k++) {
                    double change = 0.0;

                    if (velocity.x > 0) {
                        double leftNeighbor;
                        if (i > 0) leftNeighbor = conc[i - 1][j][k];
                        else if (applyInletBC) leftNeighbor = cInlet;
                        else leftNeighbor = conc[i][j][k];
                        change -= velocity.x * dt / boxSize[0] * (conc[i][j][k] - leftNeighbor);
                    } else if (velocity.x < 0 && i < boxes[0] - 1) {
                        change -= (-velocity.x) * dt / boxSize[0]
                                  * (conc[i][j][k] - conc[i + 1][j][k]);
                    }

                    if (velocity.y > 0 && j > 0) {
                        change -= velocity.y * dt / boxSize[1]
                                  * (conc[i][j][k] - conc[i][j - 1][k]);
                    } else if (velocity.y < 0 && j < boxes[1] - 1) {
                        change -= (-velocity.y) * dt / boxSize[1]
                                  * (conc[i][j][k] - conc[i][j + 1][k]);
                    }

                    double newConc = conc[i][j][k] + change;
                    field.setConc(i, j, k, Math.max(0.0, newConc));
                }
            }
        }
    }
}
