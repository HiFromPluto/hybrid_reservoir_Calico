package bsim.a1;

import bsim.BSimRandom;
import bsim.BSimTicker;
import bsim.a0.A0IdealSource;
import bsim.c1c.DaninoSIC1cFilledPocket;
import bsim.circuit.DaninoSIPeriodCheck;
import bsim.draw.BSimP3DDrawer;
import bsim.p0.DaninoPocketGeometry;
import processing.core.PGraphics3D;

import javax.vecmath.Vector3d;
import java.awt.Color;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

/**
 * Live BSim Preview window for A1 on C1c_FILLED_POCKET.
 * Viewer only. Not a gate. Not A2. NOT_FIG4B.
 *
 * Occupancy is voxel-weighted (mechanics OFF): one marker per voxel, not
 * 11000 drawn rods. Close the window to stop.
 */
public final class DaninoSIA1LivePreview {

    private static final int STEPS_PER_FRAME = 20;
    private static final double AC_VISUAL_RADIUS_UM = 4.0;
    private static final Color AC_COLOR = new Color(40, 220, 255);
    private static final Color POCKET_COLOR = new Color(80, 220, 255);

    private DaninoSIA1LivePreview() { }

    public static void main(String[] args) {
        Locale.setDefault(Locale.US);
        System.setOut(new java.io.PrintStream(System.out, true, StandardCharsets.UTF_8));
        System.setErr(new java.io.PrintStream(System.err, true, StandardCharsets.UTF_8));
        DaninoSIPeriodCheck.refuseNarma(args);
        A1BoundedTransducer.freezeCheck();

        A1BoundedTransducer.Command command = A1BoundedTransducer.Command.U_PULSE;
        if (args.length > 0 && !args[0].isBlank()) {
            command = A1BoundedTransducer.Command.valueOf(args[0].trim());
        }

        BSimRandom rng = new BSimRandom(0L);
        DaninoSIC1cFilledPocket device = new DaninoSIC1cFilledPocket(rng);
        A1BoundedTransducer ac = new A1BoundedTransducer(command, false);
        DaninoSIC1cFilledPocket.LiveView live = device.openLive(
                0.40, "AHL_KICK_005", A1BoundedTransducer.RK_DT, ac::deposit);

        System.out.println("A1 LIVE VIEW  A1_BOUNDED_TRANSDUCER  HYPOTHETICAL_DESIGN_ENVELOPE  NOT_FIG4B");
        System.out.println("C1c_FILLED_POCKET  mechanics=OFF  voxels not rods  close the window to stop.");
        System.out.printf(Locale.US, "command=%s  AC=(%.3g, %.3g, %.3g) um  voxel=(%d,%d,%d)%n",
                command.name(), A0IdealSource.X_UM, A0IdealSource.Y_UM, A0IdealSource.Z_UM,
                A0IdealSource.I_AC, A0IdealSource.J_AC, A0IdealSource.K_AC);

        live.sim.setTicker(new BSimTicker() {
            private int frames;

            @Override
            public void tick() {
                for (int s = 0; s < STEPS_PER_FRAME; s++) {
                    live.step();
                }
                frames++;
                if (frames % 25 == 0) {
                    System.err.printf(Locale.US,
                            "  live t=%.2f  I_ac=%.4g  He_ac=%.4g  J_S=%.4g  M=%.4g  NOT_FIG4B%n",
                            live.time(),
                            live.luxI[A0IdealSource.I_AC + A0IdealSource.J_AC * DaninoSIC1cFilledPocket.NX],
                            live.field.getConc(A0IdealSource.I_AC, A0IdealSource.J_AC, A0IdealSource.K_AC),
                            ac.lastJ(),
                            ac.payload());
                }
            }
        });

        live.sim.setDrawer(new BSimP3DDrawer(live.sim, 900, 700) {
            private final Vector3d acPos = new Vector3d(
                    A0IdealSource.X_UM, A0IdealSource.Y_UM, A0IdealSource.Z_UM);
            private final Vector3d voxelPos = new Vector3d();

            @Override
            public void scene(PGraphics3D p3d) {
                float lx = (float) DaninoSIC1cFilledPocket.BOUND_X;
                float ly = (float) DaninoSIC1cFilledPocket.BOUND_Y;
                p3d.ortho(0, lx, ly, 0, -1000, 10000);
                p3d.camera(lx / 2f, ly / 2f, ly,
                        lx / 2f, ly / 2f, 0f,
                        0f, 1f, 0f);
                p3d.perspective((float) Math.PI / 2f, lx / ly, 0.1f, 10000f);

                draw(live.field, Color.ORANGE, (float) (255.0 / 0.8));
                drawPocket(p3d);

                int n = live.nVoxels();
                for (int c = 0; c < n; c++) {
                    voxelPos.set(
                            (live.voxelI[c] + 0.5) * DaninoPocketGeometry.DX_UM,
                            (live.voxelJ[c] + 0.5) * DaninoPocketGeometry.DY_UM,
                            0.825);
                    double u = live.luxI[c] / 800.0;
                    if (u < 0.0) {
                        u = 0.0;
                    }
                    if (u > 1.0) {
                        u = 1.0;
                    }
                    Color cell = new Color(
                            (int) (20 + 40 * u),
                            (int) (80 + 175 * u),
                            (int) (40 + 20 * u));
                    sphere(voxelPos, 0.85, cell, 220);
                }
                sphere(acPos, AC_VISUAL_RADIUS_UM, AC_COLOR, 255);
            }

            void drawPocket(PGraphics3D p3d) {
                p3d.noFill();
                p3d.stroke(POCKET_COLOR.getRed(), POCKET_COLOR.getGreen(), POCKET_COLOR.getBlue());
                p3d.pushMatrix();
                p3d.translate(
                        (float) (DaninoSIC1cFilledPocket.BOUND_X / 2.0),
                        (float) (DaninoSIC1cFilledPocket.BOUND_Y / 2.0),
                        (float) (DaninoSIC1cFilledPocket.BOUND_Z / 2.0));
                p3d.box(
                        (float) DaninoSIC1cFilledPocket.BOUND_X,
                        (float) DaninoSIC1cFilledPocket.BOUND_Y,
                        (float) DaninoSIC1cFilledPocket.BOUND_Z);
                p3d.popMatrix();
                p3d.noStroke();
            }
        });

        live.sim.preview();
    }
}
