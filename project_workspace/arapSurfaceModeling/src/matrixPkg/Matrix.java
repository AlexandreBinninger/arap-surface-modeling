package matrixPkg;
import Utils.Pair;

public interface Matrix {
	
	public Pair<Matrix, Matrix> getSVD();
	
	public double get(int i, int j);
	
	public Matrix getMatrix(double[][] array);

}
