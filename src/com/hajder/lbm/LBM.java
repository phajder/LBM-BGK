package com.hajder.lbm;

import java.io.IOException;

import static com.hajder.lbm.PropertyReader.*;

/**
 * Lattice Boltzmann Method implementation with D2Q9 velocity model used to solve lid driven cavity problem.
 *
 *         y ^
 *           |       TOP
 *           |    u=u0, v=0
 *           -------------------
 *         L |                 | R
 *         E |                 | I
 *  u=v=0  F |                 | G  u=v=0
 *         T |                 | H
 *           |                 | T
 *           ------------------- --->
 *                  u=v=0           x
 *                 BOTTOM
 *
 * Along the y-axis i-pointer increases and along the x-axis increases j-pointer.
 * This determines described below base velocity vectors.
 * @author Piotr Hajder
 */
class LBM {
    private static final int Q=9;

    private static final double[] w = new double[] {
            4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.
    };
    private static final int[][] e = new int[][] {
            { 0, 0}, // (1,-1)        (1,0)         (1,1)
            { 0, 1}, //   6             2             5
            { 1, 0}, //                 ^
            { 0,-1}, //                 |
            {-1, 0}, // (0,-1)          |           (0,1)
            { 1, 1}, //   3 <---------- 0 ----------> 1
            { 1,-1}, //                 |
            {-1,-1}, //                 |
            {-1, 1}  // (-1,-1)         v           (-1,1)
    };               //   7             4             8

    private final BoundaryCondition bc;

// =============== PARAMETERS =============== //
    private boolean debugMode = false;
    private int maxStep = 100000;
    private double ep = 1e-4;
    private int Re = 1000;
    private double tau = 2./3.;
// ========================================== //

    private double u_lid;
    private double[][][] f, f_eq, f_tmp;
    private double[][] rho, u, v;

    LBM() throws IOException {
        int sizeX = Integer.valueOf(readProperty(SIZE_X)),
                sizeY = Integer.valueOf(readProperty(SIZE_Y));
        f = new double[sizeX][sizeY][Q];
        f_eq = new double[sizeX][sizeY][Q];
        f_tmp = new double[sizeX][sizeY][Q];

        rho = new double[sizeX][sizeY];
        u = new double[sizeX][sizeY];
        v = new double[sizeX][sizeY];

        Re = Integer.valueOf(readProperty(REYNOLDS));
        tau = Double.valueOf(readProperty(TAU));

        u_lid = Re * (tau - 0.5) / (3.0 * f[0].length);

        bc = new BoundaryCondition(f, u_lid);

        ep = Double.valueOf(readProperty(EP));
        maxStep = Integer.valueOf(readProperty(MAX_STEP));
        debugMode = Boolean.valueOf(readProperty(DEBUG));
    }

    /**
     * Calculates matrix norm.
     * Matrices have to be same size to perform this action.
     * @param m1 First matrix
     * @param m2 Second matrix
     * @return Norm value
     */
    private static double norm(double[][] m1, double[][] m2) {
        assert(m1.length == m2.length);
        assert(m1[0].length == m2[0].length);

        double sum = .0, tmp;
        for(int i=0; i<m1.length; i++) {
            for(int j=0; j<m1[i].length; j++) {
                tmp = m1[i][j] - m2[i][j];
                sum += tmp*tmp;
            }
        }
        return Math.sqrt(sum);
    }

    /**
     * LBM variables initialization.
     * Initial values of distribution function f are corresponding base weights.
     * This determines, that initial density rho is 1.0 (in lattice units).
     * Velocity u_x at the top of cavity is set to calculated before u_lid value.
     */
    private void initialize() {
        for(int i=0; i<f.length; i++) {
            for(int j=0; j<f[i].length; j++) {
                System.arraycopy(w, 0, f[i][j], 0, Q);
                rho[i][j] = 1.0;
                if(i == f.length-1)
                    u[i][j] = u_lid;
            }
        }
        bc.apply();
    }

    /**
     * Performs streaming LBM step.
     * Moves direction-specific densities f_a to the nearest neighbours.
     * At the end of this step boundary conditions are applied.
     */
    private void streamingStep() {
        final int LX = f[0].length, LY = f.length;
        int iS, jS;
        for(int k=0; k<Q; k++) {
            for (int i = 0; i < LY; i++) {
                iS = i + e[k][0];
                if(iS >= 0 && iS < LY) {
                    for (int j = 0; j < LX; j++) {
                        jS = j + e[k][1];
                        if (jS >= 0 && jS < LX)
                            f_tmp[iS][jS][k] = f[i][j][k];
                    }
                }
            }
        }

        //Copying f_tmp to f. Totally nonsense, but actually there is no way to apply BCs.
        for(int i=0; i<LY; i++) {
            for(int j=0; j<LX; j++) {
                System.arraycopy(f_tmp[i][j], 0, f[i][j], 0, Q);
            }
        }

        bc.apply();
    }

    /**
     * Calculates macroscopic parameters in lattice.
     * These parameters are:
     *  - density rho,
     *  - velocity u (directed towards x-axis),
     *  - velocity v (directed towards y-axis).
     */
    private void calculateMacroscopicParameters() {
        for(int i=0; i<f.length; i++) {
            for(int j=0; j<f[i].length; j++) {
                rho[i][j] = 0.0;
                u[i][j] = 0.0;
                v[i][j] = 0.0;

                for(int k=0; k<Q; k++) {
                    rho[i][j] += f[i][j][k];
                    u[i][j] += f[i][j][k] * e[k][0];
                    v[i][j] += f[i][j][k] * e[k][1];
                }
                u[i][j] /= rho[i][j];
                v[i][j] /= rho[i][j];
            }
        }
    }

    /**
     * Calculates equilibrium distribution function.
     * Performed for new velocity u_eq.
     * In this model there are no external forces, thus assumed:
     *  - u_eq = u_ij,
     *  - v_eq = v_ij.
     */
    private void calculateEquilibrium() {
        final double c = 1.0;
        double u_eq, v_eq, u_dot_u, c_square = c * c, e_dot_u;
        for(int i=0; i<f_eq.length; i++) {
            for(int j=0; j<f_eq[i].length; j++) {
                u_eq = u[i][j];
                v_eq = v[i][j];
                u_dot_u = u_eq * u_eq + v_eq * v_eq;
                f_eq[i][j][0] = w[0] * rho[i][j] * (1.0 - 1.5 * u_dot_u / c_square);
                for(int k=1; k<Q; k++) {
                    e_dot_u = e[k][0] * u_eq + e[k][1] * v_eq;
                    f_eq[i][j][k] = 1.0
                            + 3.0 * e_dot_u / c_square
                            + 4.5 * e_dot_u * e_dot_u / (c_square * c_square)
                            - 1.5 * u_dot_u / c_square;
                    f_eq[i][j][k] *= w[k] * rho[i][j];
                }
            }
        }
    }

    /**
     * Performs collision in entire domain.
     */
    private void collisionStep() {
        for(int i=0; i<f.length; i++) {
            for(int j=0; j<f[i].length; j++) {
                for(int k=0; k<Q; k++) {
                    f[i][j][k] -= (f[i][j][k] - f_eq[i][j][k]) / tau;
                }
            }
        }
    }

    private void printInfo() {
        System.out.println("********************************************************************************");
        System.out.println("Initial conditions:");
        System.out.println("   - Grid diameter = " + f.length + "x" + f[0].length);
        System.out.println("   - Reynolds number Re = " + Re);
        System.out.println("   - Single relaxation time tau = " + tau);
        System.out.println("   - Lid velocity u_lid = " + u_lid);
        System.out.println();
        System.out.println("Stop conditions:");
        System.out.println("   - Maximum number of iterations = " + maxStep);
        System.out.println("   - Root mean square change ep = " + ep);
        System.out.println("********************************************************************************");
        System.out.println("Starting computation...");
        System.out.println();
    }

    /**
     * Main loop.
     * Each 1000 step norm is calculated.
     */
    void compute() {
        if(debugMode) printInfo();
        double[][] previousU = new double[u.length][u[0].length];
        double norm;
        initialize();
        for(int t = 0; t < maxStep; t++) {
            streamingStep();
            calculateMacroscopicParameters();
            calculateEquilibrium();
            collisionStep();
            if((t+1) % 1000 == 0) {
                for(int i=0; i<u.length; i++) {
                    System.arraycopy(u[i], 0, previousU[i], 0, u[i].length);
                }
            }
            if(t%1000 == 0) {
                norm = norm(u, previousU);
                if(debugMode) {
                    System.out.println("Iteration: " + t);
                    System.out.println("Norm: " + norm);
                }
                if(t>0 && norm < ep) break;
            }
        }
    }

    double[][] getRho() {
        return rho;
    }

    double[][] getU() {
        return u;
    }

    double[][] getV() {
        return v;
    }
}
