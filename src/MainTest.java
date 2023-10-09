import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Test;

public class MainTest {
    @Test
    public void testSolving0() {
        double[][] inC = new double[][]{
                {5, 6, 4, 7, 0, 0, 0, 0}
        };
        double[][] inB = new double[][]{
                {120, 100, 90, 80}
        };
        double[][] inA = new double[][]{
                {3, 0, 2, 0 ,1, 0, 0, 0},
                {0, 2, 0, 4, 0, 1, 0, 0},
                {2, 0, 3, 0, 0, 0, 1, 0},
                {0, 4, 0, 2, 0, 0, 0, 1}
        };
        RealMatrix c = MatrixUtils.createRealMatrix(inC);
        RealMatrix b = MatrixUtils.createRealMatrix(inB);
        RealMatrix a = MatrixUtils.createRealMatrix(inA);
        double eps = 0.01;
        Assert.assertEquals(404.00, Main.maximize(c, a, b, inC[0], inA, inB[0], eps, true), eps);
    }

    @Test
    public void testSolving1() {
        double[][] inC = new double[][]{
                {1, -3, 2, 0, 0, 0}
        };
        double[][] inB = new double[][]{
                {7, 12, 10}
        };
        double[][] inA = new double[][]{
                {3, -1, 2, 1, 0, 0},
                {-2, 4, 0, 0, 1, 0},
                {-4, 3, 8, 0, 0, 1}
        };
        RealMatrix c = MatrixUtils.createRealMatrix(inC);
        RealMatrix b = MatrixUtils.createRealMatrix(inB);
        RealMatrix a = MatrixUtils.createRealMatrix(inA);
        double eps = 0.01;
        Assert.assertEquals(-11.00, Main.minimize(c, a, b, inC[0], inA, inB[0], eps, false), eps);
    }
}
