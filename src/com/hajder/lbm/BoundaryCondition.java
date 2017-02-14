package com.hajder.lbm;

/**
 * Applies boundary conditions in lattice domain.
 * @author Piotr Hajder
 */
class BoundaryCondition {
    private double[][][] f;
    private double u_lid;

    BoundaryCondition(double[][][] f, double u_lid) {
        this.f = f;
        this.u_lid = u_lid;
    }

    /**
     * Applies on the top side Von Neumann boundary condition.
     * @param node Lattice node
     */
    private void topBC(double[] node) {
        double rho = node[0] + node[1] + node[3] + 2 * (node[2] + node[5] + node[6]);
        node[4] = node[2];
        node[7] = node[5] - 0.5*(node[3] - node[1] + rho*u_lid);
        node[8] = node[6] + 0.5*(node[3] - node[1] + rho*u_lid);
    }

    /**
     * Applies on the bottom side bounce-back boundary condition.
     * @param node Lattice node
     */
    private void bottomBC(double[] node) {
        node[2] = node[4];
        node[5] = node[7];
        node[6] = node[8];
    }

    /**
     * Applies on the left side bounce-back boundary condition.
     * @param node Lattice node
     */
    private void leftBC(double[] node) {
        node[1] = node[3];
        node[5] = node[7];
        node[8] = node[6];
    }

    /**
     * Applies on the right side bounce-back boundary condition.
     * @param node Lattice node
     */
    private void rightBC(double[] node) {
        node[3] = node[1];
        node[6] = node[8];
        node[7] = node[5];
    }

    /**
     * Performs all boundary conditions at once:
     *  - y=f.length-1 - top BC,
     *  - y=0 - bottom BC,
     *  - x=0 - left BC,
     *  - y=f[0].length-1 - right BC.
     */
    void apply() {
        for(int i=0; i<f[0].length; i++) {
            bottomBC(f[0][i]);
            topBC(f[f.length-1][i]);
        }
        for(int i=1; i<f.length-1; i++) {
            leftBC(f[i][0]);
            rightBC(f[i][f[0].length-1]);
        }
    }
}
