package matrixPkg;

import java.util.Arrays;

import Jama.*;
import Jcg.geometry.Vector_3;
import Utils.Pair;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.SparseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleCholeskyDecomposition;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleSingularValueDecomposition;
import cern.colt.matrix.tdouble.algo.decomposition.SparseDoubleCholeskyDecomposition;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.io.MatrixVectorReader;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;

public class Parallel_Colt implements Matrix {

	DoubleMatrix2D M;

	public Parallel_Colt(DoubleMatrix2D MM) {
		this.M = MM;
	}

	public Parallel_Colt(double[][] array) {
		this.M = new SparseDoubleMatrix2D(array);
	}
	
	
	public DoubleMatrix2D getM() {
		return M;
	}

	public Matrix clone() {
		return (Matrix) (new Parallel_Colt(M));
	}

	public int getRowDimension() {
		return M.rows();
	}
	
	public int getColumnDimension() {
		return M.columns();
	}

	public Matrix getMatrix(double[][] array) {
		SparseDoubleMatrix2D result = new SparseDoubleMatrix2D(array);
		return (Matrix)result;
	}

	public Pair<Matrix, Matrix> getSVD() {
//		System.out.println("before svd1");
//		System.out.println(Arrays.deepToString(this.M.getArray()));
		DenseDoubleSingularValueDecomposition svd = new DenseDoubleSingularValueDecomposition(M, true, true);
//		System.out.println("after svd1");
		Matrix u = (Matrix) new Parallel_Colt(new SparseDoubleMatrix2D(svd.getU().toArray()));
		Matrix v = (Matrix) new Parallel_Colt(new SparseDoubleMatrix2D(svd.getV().toArray()));
		Pair<Matrix, Matrix> uv = new Pair<Matrix, Matrix>(u, v);
		return uv;
	}
	
	public Matrix getS() {
		DenseDoubleSingularValueDecomposition svd = new DenseDoubleSingularValueDecomposition(M, false, false);
		Matrix s = (Matrix) new Parallel_Colt(new SparseDoubleMatrix2D(svd.getS().toArray()));
		return s;
	}
	
	@Override
	public double get(int i, int j){
		return M.get(i,  j);
	}
	
	public void set(int i, int j, double s) {
		M.set(i, j, s);
	}
	
	public Matrix times(double s) {
		SparseDoubleMatrix2D A = (SparseDoubleMatrix2D) (M.copy());
		for(int i=0; i<this.getRowDimension(); i++) {
			for(int j=0; j<this.getColumnDimension(); j++) {
				A.set(i, j, M.get(i, j) * s);
			}
		}
		return (Matrix) (new Parallel_Colt(A));
	}
	
	public Matrix times(Matrix B) {
		if (B instanceof Parallel_Colt) {
			DenseDoubleAlgebra alg = new DenseDoubleAlgebra();
			return (Matrix) (new Parallel_Colt(new SparseDoubleMatrix2D(alg.mult(M, ((Parallel_Colt)B).M).toArray())));
		}
		System.err.println("Parallel_Colt times (not Parallel_Colt)");
		return null;
	}

	public Matrix plus(Matrix B) {
		if (B instanceof Parallel_Colt) {
			Matrix C = ((Parallel_Colt)B).clone();
			for (int i = 0; i < B.getRowDimension(); i++) {
				for (int j = 0; j<B.getColumnDimension(); j++) {
					C.set(i, j, C.get(i, j) + M.get(i, j));
				}
			}
			return (Matrix) C;
		}
		System.err.println("Parallel_Colt times (not Parallel_Colt)");
		return null;
	}
	
	public Matrix getTranspose() {
		DenseDoubleAlgebra alg = new DenseDoubleAlgebra();
		return (Matrix) (new Parallel_Colt(new SparseDoubleMatrix2D(alg.transpose(this.M).toArray())));
	}
	
	public Matrix solve(Matrix B) {
		if (B instanceof Parallel_Colt) {
			SparseDoubleMatrix2D X = new SparseDoubleMatrix2D(B.getRowDimension(), B.getColumnDimension());
			DenseDoubleCholeskyDecomposition chol = new DenseDoubleCholeskyDecomposition(new DenseDoubleMatrix2D(M.toArray()));
//			for(int j = 0; j < B.getColumnDimension(); j++) {
//				SparseDoubleMatrix1D bColumn = new SparseDoubleMatrix1D(B.getRowDimension());
//				for(int i = 0; i < B.getRowDimension(); i++) {
//					bColumn.set(i, B.get(i, j));
//				}
//				chol.solve(bColumn);
//				for(int i = 0; i < B.getRowDimension(); i++) {
//					X.set(i, j, bColumn.get(i));
//				}
//			}
			DenseDoubleMatrix2D C = new DenseDoubleMatrix2D(((Parallel_Colt)B).M.toArray());
			chol.solve(C);
			return (new Parallel_Colt(new SparseDoubleMatrix2D(C.toArray())));
		}
		System.err.println("solve not Parallel_Colt");
		return null;
	}
	
	public double determinant(){
		if (this.getColumnDimension()!=this.getRowDimension()){
			System.err.println("computing determinant of non-squared matrix");
			return 0;
		} else {
			SparseDoubleAlgebra alg = new SparseDoubleAlgebra();
			return alg.det(M);
		}
	}
	
	public void toPrint(){
		int imax = this.getRowDimension();
		int jmax=this.getColumnDimension();
		System.out.print("[");
		for (int i=0; i<imax; i++){
			System.out.print("[");
			for (int j=0; j<jmax-1; j++){
				System.out.print((this.get(i, j))+", ");
			}
			if (i!=imax-1){
				System.out.println((this.get(i, jmax-1))+"]");
			} else{
				System.out.print((this.get(i, jmax-1))+"]");
			}
		}
		System.out.println("]");
	}

	public Matrix vectorToMatrix(Vector_3 vector){
		double[][] array = {{vector.x}, {vector.y}, {vector.z}};
		return new Parallel_Colt(array);
	}


}
