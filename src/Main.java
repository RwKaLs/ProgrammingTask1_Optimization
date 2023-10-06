import org.apache.commons.math3.linear.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class Main {
    private static RealMatrix c, a, b;
    private static double eps;

    public static void main(String[] args) {
        initializeTask();
        //System.out.println("c: " + c);
        //System.out.println("a: " + a);
        //System.out.println("b: " + b);
    }

    public static int[] im(double[][] A) {
        int[] b0 = new int[A.length];
        for (int j = 0; j < A[0].length; j++) {
            int count = 0;
            int col = -1;
            int row = -1;
            for (int i = 0; i < A.length; i++) {
                if (A[i][j] != 0 && A[i][j] != 1) {
                    count = -1;
                    break;
                }
                if (A[i][j] == 1) {
                    count++;
                    col = j;
                    row = i;
                }
            }
            if (count == 1) {
                b0[row] = col;
            }
        }
        return b0;
    }

    public static double[][] transpose(double[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;

        double[][] result = new double[cols][rows];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[j][i] = matrix[i][j];
            }
        }

        return result;
    }

    public static void solve(RealMatrix C, RealMatrix A, RealMatrix b, double[] Cm, double[][] Am, double[] bm) {
        int[] bi = im(Am);
        System.out.println("basic coefs: " + Arrays.toString(bi));
        int[] nonbasic = new int[Cm.length - bi.length];
        int c = 0;
        for (int i = 0; i < Cm.length; i++) {
            int ans = -1;
            for (int el : bi) {
                if (i == el) {
                    ans = -1;
                    break;
                }
                ans = i;
            }
            if (ans != -1) {
                nonbasic[c] = ans;
                c++;
            }
        }
        System.out.println("nonbasic coefs: " + Arrays.toString(nonbasic));
        int[] nonbasic_c = nonbasic;
        int counter = 0;

        double[][] cbi_ = new double[1][bi.length];
        double[][] cnbi_ = new double[1][nonbasic.length];
        for (int i = 0; i < bi.length; i++) {
            cbi_[0][i] = Cm[bi[i]];
        }
        for (int i = 0; i < nonbasic.length; i++) {
            cnbi_[0][i] = Cm[nonbasic[i]];
        }
        double[][] Bi_ = new double[bi.length][bi.length];
        double[][] nBi_ = new double[nonbasic.length][nonbasic.length];
        for (int i = 0; i < bi.length; i++) {
            Bi_[i] = transpose(Am)[bi[i]];
        }
        for (int i = 0; i < nonbasic.length; i++) {
            nBi_[i] = transpose(Am)[nonbasic[i]];
        }


        RealMatrix cbi = MatrixUtils.createRealMatrix(cbi_);
        RealMatrix cnbi = MatrixUtils.createRealMatrix(cnbi_);

        RealMatrix Bi = MatrixUtils.createRealMatrix(Bi_);
        RealMatrix nBi = MatrixUtils.createRealMatrix(nBi_);

        while (true) {
            System.out.println("Iteration " + counter);
            counter++;

            System.out.println("    cb" + counter + ": " + cbi);
            System.out.println("    B" + counter + ": " + Bi);
            System.out.println("    nB" + counter + ": " + nBi);
            RealMatrix Bi_inverse = MatrixUtils.inverse(Bi);
            System.out.println("    B^-1: " + Bi_inverse);
            RealMatrix XBi = Bi_inverse.transpose().multiply(b.transpose());
            System.out.println("    Xb" + counter + ": " + XBi);
            double z = cbi.multiply(XBi).getData()[0][0];
            System.out.println("    z: " + z);
            System.out.println("    Optimality computations:");
            RealMatrix CbiBi_1 = cbi.multiply(Bi_inverse.transpose());
            RealMatrix zj_cj = CbiBi_1.multiply(nBi.transpose()).subtract(cnbi);
            System.out.println("    zj-cj: " + zj_cj);
            double[] zj_cj_ = zj_cj.getData()[0];

            double min = 1000000000;
            int index_entering = -1;
            int index_leaving = -1;
            //maximize
            int stop = 1;
            int cnt = 0;
            for (double el :
                    zj_cj_) {
                if (el < 0) {
                    stop = 0;
                }
                if (el < min && el < 0) {
                    index_entering = nonbasic[cnt];
                    min = el;
                }
                cnt++;
            }
            if (stop == 1) {
                System.out.print("ANSWER:\nF(x) = " + z + "\nx = ");
                int cn = 0;
                ArrayList<Double> ans = new ArrayList<>();
                for (int k = 0; k < Am[0].length; k++) {
                    int flag = -1;
                    for (int el : bi) {
                        if (k == el) {
                            ans.add(XBi.getEntry(cn, 0));
                            cn++;
                            flag = 1;
                            break;
                        }
                    }
                    if (flag == -1) {
                        ans.add(0.0);
                    }
                }
                System.out.println(ans);
                return;
            }
            System.out.println("    index of entering: " + index_entering);
            double[][] enterMx_ = new double[1][bi.length];
            enterMx_[0] = transpose(Am)[index_entering];
            RealMatrix Pe = MatrixUtils.createRealMatrix(enterMx_);
            System.out.println("    entering vec: " + Pe);
            System.out.println("    Feasibility computations:");
            System.out.println("    XB" + counter + ": " + XBi);
            RealMatrix Bi_1P1 = Bi_inverse.preMultiply(Pe);
            System.out.println("    B0^-1Pe: " + Bi_1P1);

            double x1 = 1000000000;
            double[] ratios = new double[bi.length];
            double[] xbi = new double[bi.length];
            xbi = XBi.transpose().getData()[0];
            double[] bipj = new double[bi.length];
            bipj = Bi_1P1.getData()[0];
            for (int j = 0; j < bi.length; j++) {
                if (xbi[j] / bipj[j] < x1 && bi[j] / bipj[j] > 0) {
                    x1 = xbi[j] / bipj[j];
                    index_leaving = bi[j];
                }
                ratios[j] = xbi[j] / bipj[j];

                int contains_pos = 0;
                for (double el : ratios) {
                    if (el > 0) {
                        contains_pos = 1;
                        break;
                    }
                }
                if (contains_pos == 0) {
                    System.out.println("solution is unbounded");
                    return;
                }
            }
            System.out.println("    ratios: " + Arrays.toString(ratios));
            System.out.println("    index of leaving: " + index_leaving);
            System.out.println("    leaving vec: " + Arrays.toString(A.getColumn(index_leaving)));
            bi[index_leaving - nonbasic.length] = index_entering;
            cbi.setEntry(0, index_leaving - nonbasic.length, Cm[index_entering]);
            Bi.setRow(index_leaving - nonbasic.length, transpose(Am)[index_entering]);

            c = 0;
            for (int i = 0; i < Cm.length; i++) {
                int ans = -1;
                for (int el : bi) {
                    if (i == el) {
                        ans = -1;
                        break;
                    }
                    ans = i;
                }
                if (ans != -1) {
                    nonbasic[c] = ans;
                    c++;
                }
            }
            for (int i = 0; i < nonbasic.length; i++) {
                cnbi_[0][i] = Cm[nonbasic[i]];
            }
            for (int i = 0; i < nonbasic.length; i++) {
                nBi_[i] = transpose(Am)[nonbasic[i]];
            }
            cnbi = MatrixUtils.createRealMatrix(cnbi_);
            nBi = MatrixUtils.createRealMatrix(nBi_);
        }
    }

    public static void maximize(RealMatrix C, RealMatrix A, RealMatrix b, double[] Cm, double[][] Am, double[] bm) {
        solve(C, A, b, Cm, Am, bm);
    }

    public static void minimize(RealMatrix C, RealMatrix A, RealMatrix b, double[] Cm, double[][] Am, double[] bm) {
        for (int i = 0; i < Cm.length; i++) {
            Cm[i] = -Cm[i];
        }
        C.setRow(0, Cm);
        solve(C, A, b, Cm, Am, bm);
    }

    static void initializeTask() {
        System.out.print("Input number of variables: ");
        Scanner in = new Scanner(System.in);
        int n = Integer.parseInt(in.nextLine());
        System.out.print("Input number of constraints: ");
        int m = Integer.parseInt(in.nextLine());
        System.out.print("Input vector c, separated by spaces: ");
        double[][] inC = new double[1][n];
        inC[0] = Arrays.stream(in.nextLine().split(" ")).mapToDouble(Double::valueOf).toArray();
        double[][] inA = new double[m][n];
        System.out.println("Input matrix line by line separated by spaces: ");
        for (int i = 0; i < m; ++i) {
            inA[i] = Arrays.stream(in.nextLine().split(" ")).mapToDouble(Double::valueOf).toArray();
        }
        System.out.print("Input vector b, separated by spaces: ");
        double[][] inB = new double[1][n];
        inB[0] = Arrays.stream(in.nextLine().split(" ")).mapToDouble(Double::valueOf).toArray();
        a = MatrixUtils.createRealMatrix(inA);
        b = MatrixUtils.createRealMatrix(inB);
        c = MatrixUtils.createRealMatrix(inC);
        maximize(c, a, b, inC[0], inA, inB[0]);
        System.out.println("\n\n\n");
        minimize(c, a, b, inC[0], inA, inB[0]);
    }
}
