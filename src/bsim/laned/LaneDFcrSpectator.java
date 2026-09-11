package bsim.laned;

/**
 * FCR spectator on the D0 bath. Eq. 5 + ribosome dilution.
 * Does not set elongation. lambda_C is not a Monod knob.
 */
public final class LaneDFcrSpectator {

    private double sigma;
    private double phiRb;
    private boolean frozen;

    public LaneDFcrSpectator(double lambdaIPerHour) {
        this.sigma = LaneDD4Identity.sigmaStar(lambdaIPerHour);
        this.phiRb = LaneDD4Identity.phiRbStar(lambdaIPerHour);
        this.frozen = false;
    }

    public double sigma() {
        return sigma;
    }

    public double phiRb() {
        return phiRb;
    }

    public boolean frozen() {
        return frozen;
    }

    public void freeze() {
        frozen = true;
    }

    public double muEff() {
        return LaneDD4Identity.muEffFromPhiRb(phiRb);
    }

    /**
     * One RK4 step in hours. Dilution and mu_f use Warren lambda_monod.
     */
    public void stepHours(double dtH, double lambdaMonodPerHour) {
        if (frozen || dtH <= 0.0) {
            return;
        }
        double muF = LaneDD4Identity.muFcr(lambdaMonodPerHour);
        double k1s = dSigma(sigma, muF);
        double k1r = dPhiRb(sigma, phiRb, lambdaMonodPerHour);
        double k2s = dSigma(sigma + 0.5 * dtH * k1s, muF);
        double k2r = dPhiRb(sigma + 0.5 * dtH * k1s,
                phiRb + 0.5 * dtH * k1r, lambdaMonodPerHour);
        double k3s = dSigma(sigma + 0.5 * dtH * k2s, muF);
        double k3r = dPhiRb(sigma + 0.5 * dtH * k2s,
                phiRb + 0.5 * dtH * k2r, lambdaMonodPerHour);
        double k4s = dSigma(sigma + dtH * k3s, muF);
        double k4r = dPhiRb(sigma + dtH * k3s,
                phiRb + dtH * k3r, lambdaMonodPerHour);
        sigma += (dtH / 6.0) * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
        phiRb += (dtH / 6.0) * (k1r + 2.0 * k2r + 2.0 * k3r + k4r);
    }

    private static double chiRb(double s) {
        return LaneDD4Identity.PHI_RB0
                / (1.0 - s / LaneDD4Identity.GAMMA_TR_PER_H);
    }

    private static double chiCat(double s) {
        return 1.0 - (s / LaneDD4Identity.LAMBDA_C_PER_H) * chiRb(s);
    }

    private static double dSigma(double s, double muF) {
        return s * (muF * chiCat(s) - s * chiRb(s));
    }

    private static double dPhiRb(double s, double rb, double lambdaMonod) {
        return lambdaMonod * (chiRb(s) - rb);
    }
}
