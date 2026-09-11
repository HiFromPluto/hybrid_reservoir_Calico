package bsim.transport;

/**
 * Face-normal velocities for a regular finite-volume grid (micrometres/s).
 * Positive values point in the positive coordinate direction.
 */
public final class BSimFaceVelocity {
    private final int nx;
    private final int ny;
    private final int nz;
    private final double[][][] x;
    private final double[][][] y;
    private final double[][][] z;

    public BSimFaceVelocity(int[] boxes) {
        if (boxes == null || boxes.length != 3
                || boxes[0] <= 0 || boxes[1] <= 0 || boxes[2] <= 0) {
            throw new IllegalArgumentException("boxes must contain three positive dimensions");
        }
        nx = boxes[0];
        ny = boxes[1];
        nz = boxes[2];
        x = new double[nx + 1][ny][nz];
        y = new double[nx][ny + 1][nz];
        z = new double[nx][ny][nz + 1];
    }

    public static BSimFaceVelocity constant(int[] boxes, double vx, double vy, double vz) {
        requireFinite(vx);
        requireFinite(vy);
        requireFinite(vz);
        BSimFaceVelocity velocity = new BSimFaceVelocity(boxes);
        for (int i = 0; i <= velocity.nx; i++)
            for (int j = 0; j < velocity.ny; j++)
                for (int k = 0; k < velocity.nz; k++)
                    velocity.x[i][j][k] = vx;
        for (int i = 0; i < velocity.nx; i++)
            for (int j = 0; j <= velocity.ny; j++)
                for (int k = 0; k < velocity.nz; k++)
                    velocity.y[i][j][k] = vy;
        for (int i = 0; i < velocity.nx; i++)
            for (int j = 0; j < velocity.ny; j++)
                for (int k = 0; k <= velocity.nz; k++)
                    velocity.z[i][j][k] = vz;
        return velocity;
    }

    public void setX(int faceI, int j, int k, double value) {
        requireFinite(value);
        x[faceI][j][k] = value;
    }

    public void setY(int i, int faceJ, int k, double value) {
        requireFinite(value);
        y[i][faceJ][k] = value;
    }

    public void setZ(int i, int j, int faceK, double value) {
        requireFinite(value);
        z[i][j][faceK] = value;
    }

    public double getX(int faceI, int j, int k) {
        return x[faceI][j][k];
    }

    public double getY(int i, int faceJ, int k) {
        return y[i][faceJ][k];
    }

    public double getZ(int i, int j, int faceK) {
        return z[i][j][faceK];
    }

    public int getNx() { return nx; }
    public int getNy() { return ny; }
    public int getNz() { return nz; }

    private static void requireFinite(double value) {
        if (!Double.isFinite(value)) {
            throw new IllegalArgumentException("velocity must be finite: " + value);
        }
    }
}
