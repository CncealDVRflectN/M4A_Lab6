public class Solution {
    private static double calcF(double x) {
        return 1.5 * Math.pow(Math.E, x) - 0.5 * Math.cos(x);
    }

    private static void printArr(double[] arr) {
        for (int i = 0; i < arr.length; i++) {
            System.out.print(arr[i] + " ");
        }
        System.out.println();
    }

    private static void printApprox(double[] coefs) {
        for (int i = 0; i < coefs.length; i++) {
            System.out.print(coefs[i] + " * x^" + i);
            if (i != coefs.length - 1) {
                System.out.print(" + ");
            }
        }
        System.out.println();
    }

    private static double calcApprox(double[] coefs, double x) {
        double result = 0;
        for (int i = 0; i < coefs.length; i++) {
            result += coefs[i] * Math.pow(x, i);
        }
        return result;
    }

    private static double[][] calcTransposeMatrix(double[][] mtr) {
        double[][] result = new double[mtr[0].length][mtr.length];
        for (int i = 0; i < mtr[0].length; i++) {
            for (int j = 0; j < mtr.length; j++) {
                result[i][j] = mtr[j][i];
            }
        }
        return result;
    }

    private static double[][] calcMatrixMultiplyMatrix(double[][] left, double[][] right) {
        double[][] result = new double[left.length][right[0].length];
        for (int i = 0; i < result.length; i++) {
            for(int j = 0; j < result[0].length; j++) {
                result[i][j] = 0;
                for (int k = 0; k < left[0].length; k++) {
                    result[i][j] += left[i][k] * right[k][j];
                }
            }
        }
        return result;
    }

    private static double[] calcMatrixMultiplyVector(double[][] mtr, double[] vector) {
        double[] result = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            result[i] = 0;
            for (int j = 0; j < mtr[i].length; j++) {
                result[i] += mtr[i][j] * vector[j];
            }
        }
        return result;
    }

    private static double[][] fillSForSquareRoot(double[][] mtr) {
        double[][] s = new double[mtr.length][mtr[0].length];
        double sum;
        for (int i = 0; i < s.length; i++) {
            sum = 0;
            for (int k = 0; k < i; k++) {
                sum += Math.pow(s[k][i], 2);
            }
            s[i][i] = Math.sqrt(mtr[i][i] - sum);
            for (int j = i + 1; j < s.length; j++) {
                sum = 0;
                for (int k = 0; k < i; k++) {
                    sum += s[k][i] * s[k][j];
                }
                s[i][j] = (mtr[i][j] - sum) / s[i][i];
            }
        }
        return s;
    }

    private static double[] squareRootMethod(double[][] s, double[] b) {
        double[] result = new double[s.length];
        double[] y = new double[s.length];
        double sum;
        for (int i = 0; i < s.length; i++) {
            sum = 0;
            for (int k = 0; k < i; k++) {
                sum += s[k][i] * y[k];
            }
            y[i] = (b[i] - sum) / s[i][i];
        }
        for (int i = s.length - 1; i >= 0; i--) {
            sum = 0;
            for (int k = i + 1; k < s.length; k++) {
                sum += s[i][k] * result[k];
            }
            result[i] = (y[i] - sum) / s[i][i];
        }
        return result;
    }

    private static void lsm(double[][] table, int N, int n) {
        double[][] left = new double[n + 1][n + 1];
        double[] right = new double[n + 1];
        double[] sol;
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n; j++) {
                for (int k = 0; k <= N; k++) {
                    left[i][j] += Math.pow(table[0][k], i + j);
                }
            }
            for (int j = 0; j <= N; j++) {
                right[i] += table[1][j] * Math.pow(table[0][j], i);
            }
        }
        sol = squareRootMethod(fillSForSquareRoot(calcMatrixMultiplyMatrix(calcTransposeMatrix(left), left)),
                calcMatrixMultiplyVector(calcTransposeMatrix(left), right));
        printApprox(sol);
        System.out.println("r*: " + Math.abs(calcApprox(sol, 0.1 / 3) - calcF(0.1 / 3)));
        System.out.println("r**: " + Math.abs(calcApprox(sol, 0.5 + 0.1 / 3) - calcF(0.5 +0.1 / 3)));
        System.out.println("r***: " + Math.abs(calcApprox(sol, 1 - 0.1 / 3) - calcF(1 - 0.1 / 3)));
    }

    public static void main(String[] args) {
        double[][] table = new double[2][11];
        int N = 10;
        for (int i = 0; i <= N; i++) {
            table[0][i] = (double)i / N;
            table[1][i] = calcF(table[0][i]);
        }
        System.out.print("X: ");
        printArr(table[0]);
        System.out.print("Y: ");
        printArr(table[1]);
        System.out.println("n = 3:");
        lsm(table, 10, 3);
        System.out.println("n = 5:");
        lsm(table, 10, 5);
    }
}
