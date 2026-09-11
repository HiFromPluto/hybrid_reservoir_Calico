package BacteriumFromScratch;

import java.util.List;

import javax.vecmath.Vector3d;

/**
 * Valdez 2025 eqs. (6)–(10) on sphero-cylinders.
 * DOI 10.1038/s42005-025-02078-1.
 * <p>
 * Not {@code RelaxationMover}: that integrates
 * {@code BSimCapsuleBacterium.computeNeighbourForce}
 * ({@code 0.4 k_cell (2R−d)^{2.5}}), not Hertzian δ^{3/2} with paper k_cc.
 */
public final class ValdezHertzian {

    private static final double EPS = 1e-12;

    private ValdezHertzian() {}

    /** Packing diagnostics. Public so ChassisPocket (and later dishes) can import the chassis. */
    public static final class PackingStats {
        public double maxDeltaCc;
        public double minCentreDist;
        public double maxForceMag;
        public int nPairs;
    }

    /**
     * Particle walls. Default is the closed box Jobs 2/3/3b/3c/I0 use.
     * ChassisPocket I0b passes {@link #openNeckYPlus} so the T4/T5 mouth
     * can omit the +y Hertzian in a frozen gap. Do not change the
     * no-argument {@link #relaxContacts} path.
     */
    public static final class WallSpec {
        public final boolean neckYPlus;
        public final double neckHalfWidthUm;
        public final double neckCenterXUm;

        private WallSpec(boolean neckYPlus, double neckHalfWidthUm, double neckCenterXUm) {
            this.neckYPlus = neckYPlus;
            this.neckHalfWidthUm = neckHalfWidthUm;
            this.neckCenterXUm = neckCenterXUm;
        }

        public static WallSpec closedBox() {
            return new WallSpec(false, 0.0, 0.0);
        }

        /** Gap in the +y face, width {@code widthUm}, centred at {@code centerXUm}. */
        public static WallSpec openNeckYPlus(double widthUm, double centerXUm) {
            if (!(widthUm > 0.0) || !Double.isFinite(widthUm) || !Double.isFinite(centerXUm)) {
                throw new IllegalArgumentException("neck width and centre must be finite and width > 0");
            }
            return new WallSpec(true, 0.5 * widthUm, centerXUm);
        }

        /**
         * Entire +y face open (Danino-class pocket mouth). Additive; does not
         * change {@link #closedBox()} or {@link #openNeckYPlus}.
         */
        public static WallSpec openYPlus() {
            return new WallSpec(true, Double.POSITIVE_INFINITY, 0.0);
        }

        boolean yPlusWallApplies(double x) {
            if (!neckYPlus) {
                return true;
            }
            if (Double.isInfinite(neckHalfWidthUm)) {
                return false;
            }
            return x < neckCenterXUm - neckHalfWidthUm
                    || x > neckCenterXUm + neckHalfWidthUm;
        }
    }

    /**
     * Quasi-static contact solve for one growth tick.
     * <p>
     * Job 2 Δt = 1 s is far above Hertzian contact time, so each tick
     * relaxes elastic Valdez (8)/(10) with overdamped (6)–(7) under a CFL
     * until overlap ≤ named residual. Dissipation uses v = 0 (static
     * contact solve). k_cc is not raised.
     */
    public static void relaxContacts(List<EcoliRodCell> cells, Vector3d bound) {
        relaxContacts(cells, bound, WallSpec.closedBox());
    }

    public static void relaxContacts(List<EcoliRodCell> cells, Vector3d bound,
            WallSpec walls) {
        if (walls == null) {
            throw new IllegalArgumentException("WallSpec must not be null");
        }
        for (EcoliRodCell c : cells) {
            c.velocity.set(0.0, 0.0, 0.0);
            c.angularVelocity.set(0.0, 0.0, 0.0);
        }
        for (int guard = 0; guard < ChassisParameters.MECH_SUBSTEP_CAP; guard++) {
            accumulateForces(cells, bound, walls);
            double maxDelta = maxPairDelta(cells);
            if (maxDelta <= ChassisParameters.TWO_BODY_GAP_RESIDUAL_UM) {
                return;
            }
            double stepUm = Math.min(
                    ChassisParameters.MECH_CFL_UM,
                    0.5 * Math.max(maxDelta, ChassisParameters.MECH_CFL_UM));
            double vmaxUmPerH = 0.0;
            for (EcoliRodCell cell : cells) {
                double[] drag = tiradoDrag(cell);
                double vH = cell.netForce.length() / drag[0];
                double w = cell.netTorque.length() / drag[0];
                double rPole = 0.5 * cell.L + ChassisParameters.RADIUS_UM;
                vmaxUmPerH = Math.max(vmaxUmPerH, vH + w * rPole);
            }
            if (vmaxUmPerH < 1e-18) {
                return;
            }
            double dtH = stepUm / vmaxUmPerH;
            for (EcoliRodCell cell : cells) {
                applyOverdamped(cell, cell.netForce, cell.netTorque, dtH);
            }
            if (!allFinite(cells)) {
                return;
            }
            for (EcoliRodCell c : cells) {
                c.velocity.set(0.0, 0.0, 0.0);
                c.angularVelocity.set(0.0, 0.0, 0.0);
            }
        }
    }

    public static PackingStats stats(List<EcoliRodCell> cells, Vector3d bound) {
        return stats(cells, bound, WallSpec.closedBox());
    }

    public static PackingStats stats(List<EcoliRodCell> cells, Vector3d bound,
            WallSpec walls) {
        PackingStats s = new PackingStats();
        s.maxDeltaCc = 0.0;
        s.minCentreDist = Double.POSITIVE_INFINITY;
        accumulateForces(cells, bound, walls);
        s.maxForceMag = 0.0;
        for (EcoliRodCell c : cells) {
            s.maxForceMag = Math.max(s.maxForceMag, c.netForce.length());
        }
        int n = cells.size();
        for (int i = 0; i < n; i++) {
            EcoliRodCell a = cells.get(i);
            Vector3d ca = a.centre();
            for (int j = i + 1; j < n; j++) {
                EcoliRodCell b = cells.get(j);
                s.maxDeltaCc = Math.max(s.maxDeltaCc, pairDelta(a, b));
                s.minCentreDist = Math.min(s.minCentreDist, dist(ca, b.centre()));
                s.nPairs++;
            }
        }
        if (n < 2) {
            s.minCentreDist = Double.POSITIVE_INFINITY;
        }
        return s;
    }

    static double pairDelta(EcoliRodCell a, EcoliRodCell b) {
        Vector3d c1 = new Vector3d();
        Vector3d c2 = new Vector3d();
        return overlapDelta(segmentDistance(a.x1, a.x2, b.x1, b.x2, c1, c2));
    }

    /** Surface gap = segment distance − w0. Negative ⇒ overlap. */
    static double surfaceGap(EcoliRodCell a, EcoliRodCell b) {
        Vector3d c1 = new Vector3d();
        Vector3d c2 = new Vector3d();
        return segmentDistance(a.x1, a.x2, b.x1, b.x2, c1, c2)
                - ChassisParameters.W0_UM;
    }

    static double dist(Vector3d a, Vector3d b) {
        Vector3d d = new Vector3d();
        d.sub(a, b);
        return d.length();
    }

    private static double maxPairDelta(List<EcoliRodCell> cells) {
        double maxDelta = 0.0;
        int n = cells.size();
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                maxDelta = Math.max(maxDelta, pairDelta(cells.get(i), cells.get(j)));
            }
        }
        return maxDelta;
    }

    private static boolean allFinite(List<EcoliRodCell> cells) {
        for (EcoliRodCell c : cells) {
            if (!finite(c.x1) || !finite(c.x2)) {
                return false;
            }
        }
        return true;
    }

    private static boolean finite(Vector3d v) {
        return Double.isFinite(v.x) && Double.isFinite(v.y) && Double.isFinite(v.z);
    }

    static void accumulateForces(List<EcoliRodCell> cells, Vector3d bound) {
        accumulateForces(cells, bound, WallSpec.closedBox());
    }

    static void accumulateForces(List<EcoliRodCell> cells, Vector3d bound,
            WallSpec walls) {
        for (EcoliRodCell c : cells) {
            c.netForce.set(0.0, 0.0, 0.0);
            c.netTorque.set(0.0, 0.0, 0.0);
        }
        Vector3d contactA = new Vector3d();
        Vector3d contactB = new Vector3d();
        int n = cells.size();
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                addCellCell(cells.get(i), cells.get(j), contactA, contactB);
            }
        }
        for (EcoliRodCell c : cells) {
            addCellWalls(c, bound, walls);
        }
    }

    /**
     * Valdez (8). Elastic Hertzian only here (v = 0 in the contact solve).
     * δ = d0 − d if d0 > d else 0, d0 = w0.
     */
    private static void addCellCell(EcoliRodCell a, EcoliRodCell b,
            Vector3d contactA, Vector3d contactB) {
        double d = segmentDistance(a.x1, a.x2, b.x1, b.x2, contactA, contactB);
        double delta = overlapDelta(d);
        if (delta <= 0.0) {
            return;
        }
        Vector3d n = new Vector3d();
        n.sub(contactA, contactB);
        if (n.lengthSquared() < EPS) {
            n.set(0.0, 1.0, 0.0);
        } else {
            n.normalize();
        }
        double d0 = ChassisParameters.W0_UM;
        double fn = (2.0 / 3.0) * ChassisParameters.K_CC * Math.sqrt(d0)
                * Math.pow(delta, 1.5);
        Vector3d f = new Vector3d();
        f.scale(fn, n);
        addForceAt(a, contactA, f);
        Vector3d fOpp = new Vector3d();
        fOpp.scale(-1.0, f);
        addForceAt(b, contactB, fOpp);
    }

    /** Valdez (10) on each face. F_s OFF. +y may be gapped by {@link WallSpec}. */
    private static void addCellWalls(EcoliRodCell cell, Vector3d bound, WallSpec walls) {
        addWallPole(cell, cell.x1, new Vector3d(1, 0, 0), cell.x1.x);
        addWallPole(cell, cell.x2, new Vector3d(1, 0, 0), cell.x2.x);
        addWallPole(cell, cell.x1, new Vector3d(-1, 0, 0), bound.x - cell.x1.x);
        addWallPole(cell, cell.x2, new Vector3d(-1, 0, 0), bound.x - cell.x2.x);

        addWallPole(cell, cell.x1, new Vector3d(0, 1, 0), cell.x1.y);
        addWallPole(cell, cell.x2, new Vector3d(0, 1, 0), cell.x2.y);
        if (walls.yPlusWallApplies(cell.x1.x)) {
            addWallPole(cell, cell.x1, new Vector3d(0, -1, 0), bound.y - cell.x1.y);
        }
        if (walls.yPlusWallApplies(cell.x2.x)) {
            addWallPole(cell, cell.x2, new Vector3d(0, -1, 0), bound.y - cell.x2.y);
        }

        addWallPole(cell, cell.x1, new Vector3d(0, 0, 1), cell.x1.z);
        addWallPole(cell, cell.x2, new Vector3d(0, 0, 1), cell.x2.z);
        addWallPole(cell, cell.x1, new Vector3d(0, 0, -1), bound.z - cell.x1.z);
        addWallPole(cell, cell.x2, new Vector3d(0, 0, -1), bound.z - cell.x2.z);
    }

    private static void addWallPole(EcoliRodCell cell, Vector3d pole,
            Vector3d nOut, double distInside) {
        double delta = ChassisParameters.RADIUS_UM - distInside;
        if (delta <= 0.0) {
            return;
        }
        // Cap penetration so a single bad step cannot Inf the Hertzian.
        if (delta > ChassisParameters.W0_UM) {
            delta = ChassisParameters.W0_UM;
        }
        Vector3d n = new Vector3d(nOut);
        n.normalize();
        double d0 = ChassisParameters.W0_UM;
        double fn = (2.0 * Math.sqrt(2.0) / 3.0) * ChassisParameters.K_AC
                * Math.sqrt(d0) * Math.pow(delta, 1.5);
        Vector3d f = new Vector3d();
        f.scale(fn, n);
        addForceAt(cell, pole, f);
    }

    /**
     * Overdamped (6)–(7): F = η_t ú, T = η_t θ̇ (RodCellVSPhage Integrate.cpp
     * uses η_t for both). Displacement CFL-capped.
     */
    private static void applyOverdamped(EcoliRodCell cell, Vector3d f, Vector3d torque,
            double dtH) {
        double[] drag = tiradoDrag(cell);
        double etaT = drag[0];
        cell.velocity.scale(1.0 / etaT, f);
        cell.angularVelocity.scale(1.0 / etaT, torque);

        Vector3d cm = cell.centre();
        Vector3d r1 = new Vector3d();
        r1.sub(cell.x1, cm);
        Vector3d r2 = new Vector3d();
        r2.sub(cell.x2, cm);

        Vector3d wXr1 = new Vector3d();
        wXr1.cross(cell.angularVelocity, r1);
        Vector3d wXr2 = new Vector3d();
        wXr2.cross(cell.angularVelocity, r2);

        Vector3d v1 = new Vector3d(cell.velocity);
        v1.add(wXr1);
        Vector3d v2 = new Vector3d(cell.velocity);
        v2.add(wXr2);
        double vPole = Math.max(v1.length(), v2.length());
        double cfl = ChassisParameters.MECH_CFL_UM;
        if (vPole * Math.abs(dtH) > cfl && vPole > 0.0) {
            dtH = Math.copySign(cfl / vPole, dtH);
        }

        cell.x1.scaleAdd(dtH, cell.velocity, cell.x1);
        cell.x1.scaleAdd(dtH, wXr1, cell.x1);
        cell.x2.scaleAdd(dtH, cell.velocity, cell.x2);
        cell.x2.scaleAdd(dtH, wXr2, cell.x2);
        cell.stretchPolesToLength();
        cell.centre();
    }

    /**
     * Tirado drag. eta = 1e6 * μ_liq as in RodCellVSPhage Integrate.cpp
     * (TAKEN BUILD integrator units; lengths μm, time hours).
     */
    static double[] tiradoDrag(EcoliRodCell cell) {
        double eta = ChassisParameters.etaScale();
        double l = Math.max(cell.L, 1e-6);
        double r = ChassisParameters.RADIUS_UM;
        double pc = (l + r) / r;
        double ct = 0.312 + (0.565 / pc) - (0.1 / (pc * pc));
        double cr = -0.662 + (0.917 / pc) - (0.05 / (pc * pc));
        double logp = Math.log(pc);
        double etaT = (3.0 * Math.PI * eta * (l + r)) / (logp + ct);
        double etaR = (Math.PI * eta * l * l * l) / (3.0 * (logp + cr));
        if (etaR < EPS) {
            etaR = etaT;
        }
        return new double[] {etaT, etaR};
    }

    private static void addForceAt(EcoliRodCell cell, Vector3d point, Vector3d f) {
        cell.netForce.add(f);
        Vector3d r = new Vector3d();
        r.sub(point, cell.centre());
        Vector3d torque = new Vector3d();
        torque.cross(r, f);
        cell.netTorque.add(torque);
    }

    static double overlapDelta(double segmentDistance) {
        double d0 = ChassisParameters.W0_UM;
        return d0 > segmentDistance ? d0 - segmentDistance : 0.0;
    }

    /**
     * Closest points on two segments (geomalgorithms.com), for Valdez
     * overlap only — not the BSim spring.
     */
    static double segmentDistance(Vector3d a1, Vector3d a2, Vector3d b1, Vector3d b2,
            Vector3d closestA, Vector3d closestB) {
        Vector3d u = new Vector3d();
        u.sub(a2, a1);
        Vector3d v = new Vector3d();
        v.sub(b2, b1);
        Vector3d w = new Vector3d();
        w.sub(a1, b1);
        double aa = u.dot(u);
        double bb = u.dot(v);
        double cc = v.dot(v);
        double dd = u.dot(w);
        double ee = v.dot(w);
        double D = aa * cc - bb * bb;
        double sN;
        double sD = D;
        double tN;
        double tD = D;
        if (D < EPS) {
            sN = 0.0;
            sD = 1.0;
            tN = ee;
            tD = cc;
        } else {
            sN = bb * ee - cc * dd;
            tN = aa * ee - bb * dd;
            if (sN < 0.0) {
                sN = 0.0;
                tN = ee;
                tD = cc;
            } else if (sN > sD) {
                sN = sD;
                tN = ee + bb;
                tD = cc;
            }
        }
        if (tN < 0.0) {
            tN = 0.0;
            if (-dd < 0.0) {
                sN = 0.0;
            } else if (-dd > aa) {
                sN = sD;
            } else {
                sN = -dd;
                sD = aa;
            }
        } else if (tN > tD) {
            tN = tD;
            if ((-dd + bb) < 0.0) {
                sN = 0.0;
            } else if ((-dd + bb) > aa) {
                sN = sD;
            } else {
                sN = -dd + bb;
                sD = aa;
            }
        }
        double sc = Math.abs(sN) < EPS ? 0.0 : sN / sD;
        double tc = Math.abs(tN) < EPS ? 0.0 : tN / tD;
        closestA.scaleAdd(sc, u, a1);
        closestB.scaleAdd(tc, v, b1);
        Vector3d dP = new Vector3d();
        dP.sub(closestA, closestB);
        return dP.length();
    }
}
