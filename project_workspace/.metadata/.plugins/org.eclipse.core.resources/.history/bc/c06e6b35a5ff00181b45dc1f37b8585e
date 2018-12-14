package core;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

import Jcg.geometry.*;
import Jcg.polyhedron.*;
import Utils.Rotation_3;
import matrixPkg.Jama_Matrix;
import matrixPkg.Matrix;
import Utils.Pair;

public class Computations {

	static double getOppositeAngle(Halfedge<Point_3> h) {
		Vertex<Point_3> va = h.getVertex();
		h = h.getNext();
		Vertex<Point_3> vb = h.getVertex();
		h = h.getNext();
		Vertex<Point_3> vc = h.getVertex();
		Point_3 a = va.getPoint();
		Point_3 b = vb.getPoint();
		Point_3 c = vc.getPoint();
		Vector_3 e1 = new Vector_3(b, a);
		Vector_3 e2 = new Vector_3(b, c);
		return Math.acos((double)e1.innerProduct(e2) / Math.sqrt(((double)e1.squaredLength() * (double)e2.squaredLength())));
	}

	static double getWeight(Halfedge<Point_3> h) {
		// TODO : if boundary edge, opposite halfedge exists but face is null, fix the
		// method if needed
		double alpha = getOppositeAngle(h);
		h = h.getOpposite();
		double beta = getOppositeAngle(h);
		double weight = 0.5 * ((1 / Math.atan(alpha)) + (1 / Math.atan(beta)));
		if (Double.isInfinite(weight) || Double.isNaN(weight)) {
			return 0;
		} else {
			return weight;
		}
	}

	static ArrayList<Halfedge<Point_3>> getNeighbors(Halfedge<Point_3> h) {
		Halfedge<Point_3> e = h.getNext();
		Halfedge<Point_3> o = h.getOpposite();
		ArrayList<Halfedge<Point_3>> neighbors = new ArrayList<Halfedge<Point_3>>();
		while (e != o) {
			neighbors.add(e);
			e = e.getOpposite().getNext();
		}
		neighbors.add(o);
		return (neighbors);
	}
	
//	static ArrayList<Vertex<Point_3>> getNeighbors(Vertex<Point_3> v){
//		Halfedge<Point_3> h = v.getHalfedge();
//		Halfedge<Point_3> e = h.getNext();
//		Halfedge<Point_3> o = h.getOpposite();
//		ArrayList<Vertex<Point_3>> neighbors = new ArrayList<Vertex<Point_3>>();
//		while (e != o) {
//			neighbors.add(e.getVertex());
//			e = e.getOpposite().getNext();
//		}
//		neighbors.add(o.getVertex());
//		return (neighbors);
//	}

	static HashMap<Halfedge<Point_3>, Double> getWeightsMap(Halfedge<Point_3> h, ArrayList<Halfedge<Point_3>> neighbors) {
		HashMap<Halfedge<Point_3>, Double> weights = new HashMap<Halfedge<Point_3>, Double>();
		Iterator<Halfedge<Point_3>> iter = neighbors.iterator();
		while (iter.hasNext()) {
			Halfedge<Point_3> e = iter.next();
			weights.put(e, (Double)getWeight(e));
		}
		return weights;
	}
	
//	HashMap<Vertex<Point_3>, Double> getWeightsArray(Vertex v, ArrayList<Vertex<Point_3>> neighbors){
//		HashMap<Vertex<Point_3>, Double> weights = new HashMap<Vertex<Point_3>, Double>();
//		for (Vertex<Point_3> n : neighbors){
//			
//			weights.put(v, )
//		}
//	}
	
//	Matrix getColumnMatrix(Halfedge<Point_3> h, ArrayList<Halfedge<Point_3>> neighbors) {
//		int n = neighbors.size();
//		Matrix P = new Jama_Matrix(new double[3][n]);
//		Vertex<Point_3> hVertex = h.getVertex();
//		Point_3 i = hVertex.getPoint();
//		for (int k = 0; k < n; k++) {
//			Halfedge<Point_3> e = neighbors.get(k);
//			Vertex<Point_3> eVertex = e.getVertex();
//			Point_3 j = eVertex.getPoint();
//			Vector_3 eij = new Vector_3(i, j);
//			P.set(k, 0, (double) eij.getX());
//			P.set(k, 1, (double) eij.getY());
//			P.set(k, 2, (double) eij.getZ());
//		}
//		return P;	
//	}
	
	static Matrix getColumnMatrix(Matrix points, int i, ArrayList<Halfedge<Point_3>> neighbors) {
		// points contains the vertices, P is the Pi of one vertex
		int n = neighbors.size();
		Matrix P = new Jama_Matrix(new double[n][3]);
		for (int k=0; k < n; k++) {
			int j = neighbors.get(k).getVertex().index;
			P.set(k, 0, (double) (points.get(i, 0) - points.get(j, 0))); // eij = pi - pj
			P.set(k, 1, (double) (points.get(i, 1) - points.get(j, 1)));
			P.set(k, 2, (double) (points.get(i, 2) - points.get(j, 2)));
		}
		return P;	
	}
	
	static Matrix getWeightsMatrix(ArrayList<Halfedge<Point_3>> neighbors, HashMap<Vertex<Point_3>, Double> weights) {
		int n = neighbors.size();
		Matrix D = new Jama_Matrix(new double[n][n]);
		for (int k = 0; k < n; k++) {
			Halfedge<Point_3> e = neighbors.get(k);
			Vertex<Point_3> eVertex = e.getVertex();
			D.set(k, k, (double) weights.get(eVertex));
		}
		return D;
	}
	
//	Matrix getCovarianceMatrix(Matrix P, Matrix PPrime, Matrix D) {
	static Matrix getCovarianceMatrix(Matrix points, Matrix pointsPrime, int i, ArrayList<Halfedge<Point_3>> neighbors, HashMap<Vertex<Point_3>, Double> weights) {
		Matrix P = getColumnMatrix(points, i, neighbors);
//		System.out.println(Arrays.deepToString(((Jama_Matrix)P).getM().getArray()));
		Matrix PPrime = getColumnMatrix(pointsPrime, i, neighbors);
//		System.out.println(Arrays.deepToString(((Jama_Matrix)PPrime).getM().getArray()));
		Matrix D = getWeightsMatrix(neighbors, weights);
//		System.out.println(Arrays.deepToString(((Jama_Matrix)D).getM().getArray()));
//		System.out.println(Arrays.deepToString((((Jama_Matrix)(D.times(PPrime))).getM().getArray())));
		Matrix S = P.getTranspose().times(D.times(PPrime));
//		System.out.println(Arrays.deepToString(((Jama_Matrix)S).getM().getArray()));
		return S;
	}
	
//	Matrix getCovarianceMatrix(Matrix M, Halfedge<Point_3> h, Halfedge<Point_3> hPrime) {
//		ArrayList<Halfedge<Point_3>> neighbors = getNeighbors(h);
//		ArrayList<Halfedge<Point_3>> neighborsPrime = getNeighbors(hPrime); // hopefully the order is the same
//		int n = neighbors.size();
//		Matrix S = M.getMatrix(new double[n][n]);
//		HashMap<Halfedge<Point_3>, Double> weights = getWeightsMap(h, neighbors);
//		Vertex<Point_3> hVertex = h.getVertex();
//		Point_3 i = hVertex.getPoint();
//		Vertex<Point_3> hVertexPrime = hPrime.getVertex();
//		Point_3 iPrime = hVertexPrime.getPoint();
//		for (int k = 0; k < n; k++) {
//			Halfedge<Point_3> e = neighbors.get(k);
//			Halfedge<Point_3> ePrime = neighborsPrime.get(k);
//			Vertex<Point_3> eVertex = e.getVertex();
//			Point_3 j = eVertex.getPoint();
//			Vertex<Point_3> eVertexPrime = ePrime.getVertex();
//			Point_3 jPrime = eVertexPrime.getPoint();
//			Vector_3 eij = new Vector_3(i, j);
//			Vector_3 eijPrime = new Vector_3(iPrime, jPrime);
//			double[][] eij_column = {{(double) eij.getX()}, {(double) eij.getY()}, {(double) eij.getZ()}};
//			double[][] eijPrime_row = {{(double) eijPrime.getX(), (double) eijPrime.getY(), (double) eijPrime.getZ()}};
//			Matrix eij_column_M = M.getMatrix(eij_column);
//			Matrix eijPrime_row_M = M.getMatrix(eijPrime_row);
//			S = S.plus(eij_column_M.times(eijPrime_row_M).times(weights.get(e)));
//		}
//		return S;
//	}
	
//	Matrix getWeightsMatrix(Matrix M, Halfedge<Point_3> h) {
//		ArrayList<Double> weights = new ArrayList<Double>();
//		ArrayList<Halfedge<Point_3>> neighbors = getNeighbors(h);
//		Iterator<Halfedge<Point_3>> iter = neighbors.iterator();
//		while(iter.hasNext()) {
//			Halfedge<Point_3> e = iter.next();
//			weights.add(getWeight(e));
//		}
//		double[][] weights_array = new double[weights.size()][weights.size()];
//		for (int i = 0; i < weights_array.length; i++) {
//			weights_array[i][i] = weights.get(i);
//		}
//		return M.getMatrix(weights_array);
//	}
	
//	Rotation_3 getHalfedgeRotation(Matrix S, Halfedge<Point_3> h) {
	static Rotation_3 getVertexRotation(Matrix points, Matrix pointsPrime, int i, ArrayList<Halfedge<Point_3>> neighbors, HashMap<Vertex<Point_3>, Double> weights) {
		//System.out.println("coucou");
		Matrix S = getCovarianceMatrix(points, pointsPrime, i, neighbors, weights);
//		System.out.println(Arrays.deepToString(((Jama_Matrix)S).getM().getArray()));
//		System.out.println("before svd");
		Pair<Matrix, Matrix> UV = S.getSVD();
//		System.out.println("after svd");
		Matrix U = UV.getFirst();
		Matrix V = UV.getSecond();
		Rotation_3 R = new Rotation_3(V.times(U.getTranspose()));
		if (R.getMatrix().determinant() < 0) {
			Matrix Sigma = S.getS();
			int lambda_index = 0;
			double lambda = Math.abs(S.get(0, 0));
			for(int k = 1; k<S.getColumnDimension(); k++) {
				if(Math.abs(S.get(k, k)) < lambda) {
					lambda = Math.abs(S.get(k, k));
					lambda_index = k;
				}
			}
			for(int k = 0; k < U.getRowDimension(); k++) {
				U.set(k, lambda_index, -U.get(k, lambda_index));
			}
			R = new Rotation_3(V.times(U.getTranspose()));
		}
		return R;
	}
	
	static Point_3 getBi(Vertex<Point_3> v, HashMap<Vertex<Point_3>, Rotation_3> VertRotMap, ArrayList<Halfedge<Point_3>> neighbors, HashMap<Vertex<Point_3>, Double> weights){
		Point_3 b = new Point_3(0., 0., 0.);
		Rotation_3 R_i = VertRotMap.get(v);
		for (Halfedge<Point_3> h : neighbors){
			Rotation_3 R_j = VertRotMap.get(h.getVertex());
			Matrix R_ij = R_i.getMatrix().plus(R_j.getMatrix());
			Rotation_3 transformation = new Rotation_3(R_ij.times(0.5 * weights.get(h.getVertex())));
			Point_3 tmp = Point_3.linearCombination(new Point_3[] {v.getPoint(),  h.getVertex().getPoint()}, new Number[] {1, -1});
//			Point_3 tmp = new Point_3(0., 0., 0.);
//			tmp.setX((double) v.getPoint().getX() - (double) h.getVertex().getPoint().getX());
//			tmp.setY((double) v.getPoint().getY() - (double) h.getVertex().getPoint().getY());
//			tmp.setZ((double) v.getPoint().getZ() - (double) h.getVertex().getPoint().getZ());
			tmp = transformation.transform(tmp);
			b = Point_3.linearCombination(new Point_3[] {b, tmp}, new Number[] {1, 1});
//			b.setX((double)b.getX() + (double)tmp.getX());
//			b.setY((double)b.getY() + (double)tmp.getY());
//			b.setZ((double)b.getZ() + (double)tmp.getZ());
		}
		return b;
	}
	
}
