package bsim.laneb;

import bsim.transport.BSimTransportField;

/**
 * Coarse 20×10 maps, left–right contrast, Pearson autocorrelation.
 * λ = 50 µm. B0_MOTILITY_MEMORY.
 */
public final class LaneBMaps {

    private LaneBMaps() { }

    public static double[][] zeros() {
        return new double[LaneBIdentity.RX][LaneBIdentity.RY];
    }

    public static void addCell(double[][] map, double x, double y) {
        map[LaneBIdentity.readI(x)][LaneBIdentity.readJ(y)] += 1.0;
    }

    public static double[][] copy(double[][] map) {
        double[][] out = zeros();
        for (int i = 0; i < LaneBIdentity.RX; i++) {
            System.arraycopy(map[i], 0, out[i], 0, LaneBIdentity.RY);
        }
        return out;
    }

    public static double cellContrast(double[][] density) {
        double left = 0.0;
        double right = 0.0;
        int mid = LaneBIdentity.RX / 2;
        for (int i = 0; i < LaneBIdentity.RX; i++) {
            double col = 0.0;
            for (int j = 0; j < LaneBIdentity.RY; j++) {
                col += density[i][j];
            }
            if (i < mid) left += col;
            else right += col;
        }
        return (left - right) / (left + right + 1e-12);
    }

    public static double fieldSceneContrast(BSimTransportField att, BSimTransportField rep) {
        double voxelVol = LaneBIdentity.DX * LaneBIdentity.DY * LaneBIdentity.DZ;
        double sL = 0.0;
        double sR = 0.0;
        int mid = LaneBIdentity.NX / 2;
        for (int i = 0; i < LaneBIdentity.NX; i++) {
            for (int j = 0; j < LaneBIdentity.NY; j++) {
                double a = att.getConc(i, j, 0) * voxelVol;
                double r = rep.getConc(i, j, 0) * voxelVol;
                if (i < mid) {
                    sL += a;
                    sR += r;
                } else {
                    sL += r;
                    sR += a;
                }
            }
        }
        return (sL - sR) / (sL + sR + 1e-12);
    }

    public static double cellStripeContrast(double[] xs, double slabUm) {
        double sAtt = 0.0;
        double sRep = 0.0;
        for (int n = 0; n < xs.length; n++) {
            if (Math.floor(xs[n] / slabUm) % 2.0 == 0.0) sAtt += 1.0;
            else sRep += 1.0;
        }
        return (sAtt - sRep) / (sAtt + sRep + 1e-12);
    }

    public static double fieldStripeContrast(BSimTransportField att, BSimTransportField rep, double slabUm) {
        double voxelVol = LaneBIdentity.DX * LaneBIdentity.DY * LaneBIdentity.DZ;
        double sAtt = 0.0;
        double sRep = 0.0;
        for (int i = 0; i < LaneBIdentity.NX; i++) {
            double xc = (i + 0.5) * LaneBIdentity.DX;
            boolean even = Math.floor(xc / slabUm) % 2.0 == 0.0;
            for (int j = 0; j < LaneBIdentity.NY; j++) {
                double a = att.getConc(i, j, 0) * voxelVol;
                double r = rep.getConc(i, j, 0) * voxelVol;
                if (even) {
                    sAtt += a;
                    sRep += r;
                } else {
                    sAtt += r;
                    sRep += a;
                }
            }
        }
        return (sAtt - sRep) / (sAtt + sRep + 1e-12);
    }

    public static double medianAbsDeltaL(BSimTransportField att, BSimTransportField rep, double ell) {
        double[] lCent = new double[LaneBIdentity.NX];
        int jMid = LaneBIdentity.NY / 2;
        for (int i = 0; i < LaneBIdentity.NX; i++) {
            lCent[i] = att.getConc(i, jMid, 0) - rep.getConc(i, jMid, 0);
        }
        int n = (int) Math.floor((LaneBIdentity.LX - ell) / 1.0) + 1;
        double[] d = new double[n];
        for (int s = 0; s < n; s++) {
            double x = s * 1.0;
            d[s] = Math.abs(interpCentreL(lCent, x + ell) - interpCentreL(lCent, x));
        }
        java.util.Arrays.sort(d);
        if (n % 2 == 1) return d[n / 2];
        return 0.5 * (d[n / 2 - 1] + d[n / 2]);
    }

    static double interpCentreL(double[] lCent, double x) {
        double pos = (x - 0.5 * LaneBIdentity.DX) / LaneBIdentity.DX;
        if (pos <= 0.0) return lCent[0];
        if (pos >= LaneBIdentity.NX - 1) return lCent[LaneBIdentity.NX - 1];
        int i = (int) Math.floor(pos);
        double f = pos - i;
        return lCent[i] + f * (lCent[i + 1] - lCent[i]);
    }

    public static double meanAbsLigand(BSimTransportField att, BSimTransportField rep) {
        double s = 0.0;
        int n = LaneBIdentity.NX * LaneBIdentity.NY * LaneBIdentity.NZ;
        for (int i = 0; i < LaneBIdentity.NX; i++) {
            for (int j = 0; j < LaneBIdentity.NY; j++) {
                s += Math.abs(att.getConc(i, j, 0) - rep.getConc(i, j, 0));
            }
        }
        return s / n;
    }

    public static double[][] coarsenField(BSimTransportField att, BSimTransportField rep) {
        double[][] map = zeros();
        int sx = LaneBIdentity.NX / LaneBIdentity.RX;
        int sy = LaneBIdentity.NY / LaneBIdentity.RY;
        for (int i = 0; i < LaneBIdentity.NX; i++) {
            for (int j = 0; j < LaneBIdentity.NY; j++) {
                int ii = i / sx;
                int jj = j / sy;
                map[ii][jj] += att.getConc(i, j, 0) - rep.getConc(i, j, 0);
            }
        }
        return map;
    }

    public static double pearson(double[][] a, double[][] b) {
        int n = LaneBIdentity.RX * LaneBIdentity.RY;
        double ma = 0.0;
        double mb = 0.0;
        for (int i = 0; i < LaneBIdentity.RX; i++) {
            for (int j = 0; j < LaneBIdentity.RY; j++) {
                ma += a[i][j];
                mb += b[i][j];
            }
        }
        ma /= n;
        mb /= n;
        double num = 0.0;
        double va = 0.0;
        double vb = 0.0;
        for (int i = 0; i < LaneBIdentity.RX; i++) {
            for (int j = 0; j < LaneBIdentity.RY; j++) {
                double da = a[i][j] - ma;
                double db = b[i][j] - mb;
                num += da * db;
                va += da * da;
                vb += db * db;
            }
        }
        double den = Math.sqrt(va * vb);
        if (den < 1e-18) {
            return Double.NaN;
        }
        return num / den;
    }
}
