package matrixPkg;

import Jama.*;
import Utils.Pair;

public class Jama_Matrix implements Matrix {

	Jama.Matrix M;

	public Jama_Matrix(Jama.Matrix MM) {
		this.M = MM;
	}

	public Jama_Matrix(double[][] array) {
		this.M = new Jama.Matrix(array);
	}

	public Matrix clone() {
		return new Jama_Matrix(M);
	}

	public int getRowDimension() {
		return M.getRowDimension();
	}
	
	public int getColumnDimension() {
		return M.getColumnDimension();
	}

	public Jama_Matrix getMatrix(double[][] array) {
		Jama_Matrix result = new Jama_Matrix(array);
		return result;
	}

	public Pair<Matrix, Matrix> getSVD() {
		SingularValueDecomposition svd = this.M.svd();
		Matrix u = (Matrix) svd.getU();
		Matrix v = (Matrix) svd.getV();
		Pair<Matrix, Matrix> uv = new Pair<Matrix, Matrix>(u, v);
		return uv;
	}
	
	@Override
	public double get(int i, int j){
		return M.get(i,  j);
	}
	
	public void set(int i, int j, double s) {
		M.set(i, j, s);
	}
	
	public Jama_Matrix times(double s) {
		return (new Jama_Matrix(this.M.times(s)));
	}

	public Matrix times(Matrix B) {
		if (B instanceof Jama_Matrix) {
			return (new Jama_Matrix(this.M.times(((Jama_Matrix)B).M)));
		}
		System.err.println("Jama_Matrix times (not Jama_Matrix)");
		return null;
	}

	public Matrix plus(Matrix B) {
		if (B instanceof Jama_Matrix) {
			return (new Jama_Matrix(this.M.plus(((Jama_Matrix)B).M)));
		}
		System.err.println("Jama_Matrix plus (not Jama_Matrix)");
		return null;
	}
	
	public Jama_Matrix getTranspose() {
		return (new Jama_Matrix(this.M.transpose()));
	}
	
	public Matrix solve(Matrix B) {
		if (B instanceof Jama_Matrix) {
			Jama.CholeskyDecomposition A = new Jama.CholeskyDecomposition(M);
			return (new Jama_Matrix(A.solve(((Jama_Matrix)B).M)));
		}
		System.err.println("solve not Jama_Matrix");
		return null;
	}
}
