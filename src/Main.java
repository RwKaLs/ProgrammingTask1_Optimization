import org.apache.commons.math3.linear.*;

import java.util.Arrays;
import java.util.Scanner;

public class Main {
    private static RealMatrix c, a, b;
    private static double eps;
    public static void main(String[] args) {
        initializeTask();
        System.out.println("c: " + c);
        System.out.println("a: " + a);
        System.out.println("b: " + b);
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
        System.out.print("Input vector c, separated by spaces: ");
        double[][] inB = new double[1][n];
        inB[0] = Arrays.stream(in.nextLine().split(" ")).mapToDouble(Double::valueOf).toArray();
        a = MatrixUtils.createRealMatrix(inA);
        b = MatrixUtils.createRealMatrix(inB);
        c = MatrixUtils.createRealMatrix(inC);
        System.out.print("Input approximation accuracy: ");
        eps = Double.parseDouble(in.nextLine());
    }
}
