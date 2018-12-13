package matrixPkg;
import Utils.Pair;

public interface Matrix {
	
	public Pair<Matrix, Matrix> getSVD();
	
	public int getRowDimension();
	
	public int getColumnDimension();
	
	public double get(int i, int j);
	
	public void set(int i, int j, double s);
	
	public Matrix getMatrix(double[][] array);

	public Matrix times(double s);
	
	public Matrix times(Matrix B);
	
	public Matrix plus(Matrix B);
	
	public Matrix getTranspose();
	
	public Matrix solve(Matrix B);
}
