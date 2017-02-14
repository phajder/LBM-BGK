package com.hajder.lbm;

import java.io.*;

/**
 * Main function.
 * @author Piotr Hajder
 */
public class Main {
    private static final String OUTPUT = "output/";
    private static final String NORM = OUTPUT + "norm.data";
    private static final String U = OUTPUT + "u.data";
    private static final String V = OUTPUT + "v.data";
    private static final String RHO = OUTPUT + "rho.data";

    public static void main(String[] args) {
        try {
            LBM lbm = new LBM();
            lbm.compute();
            double[][] u = lbm.getU();
            double[][] v = lbm.getV();
            double[][] matrix = new double[u.length][u[0].length];
            for (int i = 0; i < u.length; i++) {
                for (int j = 0; j < u[i].length; j++) {
                    matrix[i][j] = Math.sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]);
                }
            }
            printData(matrix, NORM);
            printData(u, U);
            printData(v, V);
            printData(lbm.getRho(), RHO);
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.err.println("Error while reading property file. " +
                    "Make sure property file is in the same location as application.");
        } catch (NumberFormatException e1) {
            System.err.println("Invalid parameter in property file. " +
                    "Invalid number conversion " + e1.getMessage().toLowerCase());
        }
    }

    private static void printData(double[][] matrix, String filename) {
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new FileWriter(filename));
            for (double[] tab : matrix) {
                for(int i=0; i<tab.length; i++) {
                    pw.print(tab[i]);
                    if(i < tab.length - 1) pw.print(",");
                }
                pw.println();
            }
            pw.flush();
            pw.close();
        } catch (IOException e) {
            System.err.println("Error while saving output data to file.");
        } finally {
            if(pw != null) {
                pw.flush();
                pw.close();
            }
        }
    }
}
