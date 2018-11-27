package matrixPkg;

import Utils.Pair;

public class Jama_Matrix implements Matrix {

	@Override
	public Pair<Matrix, Matrix> SVD() {
		// TODO Auto-generated method stub
		SingularValueDecomposition h_svd = this.svd();
		Matrix u = h_svd.getU();
		Matrix v = h_svd.getV();
		Matrix r = v.times(u.transpose());
		return null;
	}
	
}
