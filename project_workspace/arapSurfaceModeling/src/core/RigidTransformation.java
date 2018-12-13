package core;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
//import java.util.Iterator;

import Jcg.geometry.*;
import Jcg.polyhedron.*;
import Utils.Rotation_3;
import matrixPkg.Jama_Matrix;
import matrixPkg.Matrix;

public class RigidTransformation {
	
	public Polyhedron_3<Point_3> polyhedron3D;
	public HashMap<Vertex<Point_3>, Rotation_3> VertRotMap; //Map the vertex i to the rotation R_i
	public HashMap<Vertex<Point_3>, ArrayList<Halfedge<Point_3>>> globalNeighbors; // Hashmap<i, neighborsOfI>
	public HashMap<Vertex<Point_3>, HashMap<Vertex<Point_3>, Double>> weightij; // Hashmap<i, Hashmap<j, wij>>
	public Matrix L; //Laplace-Beltrami operator
	public ArrayList<Integer> mobilePoints; // the points the user is allowed to move
	public ArrayList<Integer> fixedPoints; // the points the user wants to stay at a given position
	public Matrix p; //Matrix of all points
	public Matrix pPrime; //Matrix of new points
	public Matrix b; //eq (9) of the paper

	/*
	 * 1) Precompute the weight coefficients w_ij
	 * We want to solve Equation 9 from the paper
	 * 2) Prefactorization of the System
	 * 3) Initial guess p0 : last position
	 * 4) Estimate local rotation R_i (sec 2.1)
	 * 5) Solve 9 to obtain new positions
	 * 6) Goto 4)
	 * 
	 */
	
	public RigidTransformation(Polyhedron_3<Point_3> polyhedron) {//some constraints in parameters
		// Step 1 & 2
		polyhedron3D = polyhedron;
		VertRotMap = new HashMap<Vertex<Point_3>, Rotation_3>();
		globalNeighbors = new HashMap<Vertex<Point_3>, ArrayList<Halfedge<Point_3>>>();
		weightij = new HashMap<Vertex<Point_3>, HashMap<Vertex<Point_3>, Double>>();
		L = new Jama_Matrix(new Jama.Matrix(polyhedron3D.vertices.size(), polyhedron3D.vertices.size()));
		mobilePoints = new ArrayList<Integer>();
		fixedPoints = new ArrayList<Integer>();
		p = new Jama_Matrix(new Jama.Matrix(polyhedron3D.vertices.size(), 3));
		pPrime = new Jama_Matrix(new Jama.Matrix(polyhedron3D.vertices.size(), 3));
		b = new Jama_Matrix(new Jama.Matrix(polyhedron3D.vertices.size(), 3));
		
		for (Vertex<Point_3> v : polyhedron3D.vertices){
			weightij.put(v, new HashMap<Vertex<Point_3>, Double>());
			Point_3 vPoint = v.getPoint();
			p.set(v.index, 0, (Double) vPoint.getX());
			p.set(v.index, 1, (Double) vPoint.getY());
			p.set(v.index, 2, (Double) vPoint.getZ());
			pPrime.set(v.index, 0, (Double) vPoint.getX());
			pPrime.set(v.index, 1, (Double) vPoint.getY());
			pPrime.set(v.index, 2, (Double) vPoint.getZ());
		}
		for (Halfedge<Point_3> e : polyhedron3D.halfedges){
			globalNeighbors.put(e.getVertex(), Computations.getNeighbors(e));
			HashMap<Halfedge<Point_3>, Double> tmp = Computations.getWeightsMap(e, globalNeighbors.get(e.getVertex()));
			int i = e.getVertex().index; // hopefully, it is polyhedron3D.vertices.indexOf(e.getVertex()). If it doesn't work, just initialize a global array
			L.set(i, i, 0);
			for (Halfedge<Point_3> f : tmp.keySet()){ // for all incident halfedge f of e
				int j = f.getVertex().index;
				weightij.get(e.getVertex()).put(f.getVertex(), tmp.get(f));
				L.set(i, j, -tmp.get(f));
				L.set(j, i, -tmp.get(f));
				L.set(i, i, L.get(i, i)+tmp.get(f));
			}
		}
	}
	
	public void updateEverything() {
		p = pPrime.clone();
		for (Vertex<Point_3> v : polyhedron3D.vertices){
			weightij.put(v, new HashMap<Vertex<Point_3>, Double>());
			int i = v.index;
			Point_3 pi = new Point_3(p.get(i, 0), p.get(i, 1), p.get(i, 2));
			v.setPoint(pi);
		}
		for (Halfedge<Point_3> e : polyhedron3D.halfedges){
			HashMap<Halfedge<Point_3>, Double> tmp = Computations.getWeightsMap(e, globalNeighbors.get(e.getVertex()));
			int i = e.getVertex().index; // hopefully, it is polyhedron3D.vertices.indexOf(e.getVertex()). If it doesn't work, just initialize a global array
			if (mobilePoints.contains(i) || fixedPoints.contains(i)) {
				for (int k = 0; k<L.getColumnDimension(); k++) {
					L.set(i, k, 0);
					L.set(k, i, 0);
				}
				L.set(i, i, 1);
			} else {
				L.set(i, i, 0);
				for (Halfedge<Point_3> f : tmp.keySet()){
					int j = f.getVertex().index;
					weightij.get(e.getVertex()).put(f.getVertex(), tmp.get(f));
					L.set(i, j, -tmp.get(f));
					L.set(j, i, -tmp.get(f));
					L.set(i, i, L.get(i, i)+tmp.get(f));
				}
			}
		}	
	}
	
	public void arapIteration() {
		// Step 3
		//TODO on devrait pas faire pareil pour fixedPoints ?
		for (Integer index : mobilePoints){
			Vertex<Point_3> v = polyhedron3D.vertices.get(index);
			Point_3 vPoint = v.getPoint();
			pPrime.set(v.index, 0, (Double) vPoint.getX());
			pPrime.set(v.index, 1, (Double) vPoint.getY());
			pPrime.set(v.index, 2, (Double) vPoint.getZ());
			//TODO : bizarre non ? i==index non ?
			int i = v.index;
			if (mobilePoints.contains(i) || fixedPoints.contains(i)) {
				for (int k = 0; k<L.getColumnDimension(); k++) {
					L.set(i, k, 0);
					L.set(k, i, 0);
				}
				L.set(i, i, 1);
			}
			// the few next lines should be uncommented if we need to update the weights too
//			Halfedge<Point_3> e = v.getHalfedge();
//			HashMap<Halfedge<Point_3>, Double> tmp = Computations.getWeightsMap(e, globalNeighbors.get(e.getVertex()));
//			int i = v.index; // hopefully, it is polyhedron3D.vertices.indexOf(e.getVertex())
//			L.set(i, i, 0);
//			for (Halfedge<Point_3> f : tmp.keySet()){
//				int j = f.getVertex().index;
//				weightij.get(e.getVertex()).put(f.getVertex(), tmp.get(f));
//				L.set(i, j, -tmp.get(f));
//				L.set(i, i, L.get(i, i)+tmp.get(f));
//			}
		}
		
		// Step 4
		System.out.println("step 4");
		int count = 0;
		for (Vertex<Point_3> v : polyhedron3D.vertices) {
//			System.out.println("step 4");
			int i = v.index;

//			if (mobilePoints.contains(i) || fixedPoints.contains(i)) {
//				VertRotMap.put(v, Rotation_3.rotationAxisX(0.));
//			} else {
				VertRotMap.put(v, Computations.getVertexRotation(p, pPrime, i, globalNeighbors.get(v), weightij.get(v)));
//			}
				System.out.println(count);
				count++;
		}
		
		// Step 5
		for (Vertex<Point_3> v : polyhedron3D.vertices) {
//			System.out.println("step 5");
			int i = v.index;
			Point_3 bi = Computations.getBi(v, VertRotMap, globalNeighbors.get(v), weightij.get(v));
			b.set(i, 0, (Double) bi.getX());
			b.set(i, 1, (Double) bi.getY());
			b.set(i, 2, (Double) bi.getZ());
		}
//		System.out.println(Arrays.deepToString(((Jama_Matrix)L).getM().getArray()));
		Matrix pSecond = L.solve(b);
		
		for (Integer index : fixedPoints){
			pSecond.set(index, 0, pPrime.get(index, 0));
			pSecond.set(index, 1, pPrime.get(index, 1));
			pSecond.set(index, 2, pPrime.get(index, 2));
		}

		for (Integer index1 : mobilePoints){
			pSecond.set(index1, 0, pPrime.get(index1, 0));
			pSecond.set(index1, 1, pPrime.get(index1, 1));
			pSecond.set(index1, 2, pPrime.get(index1, 2));
		}

		pPrime = pSecond;
		
	}
	
	public boolean isSymmetric(Matrix M){
		int imax = M.getColumnDimension();
		int jmax = M.getRowDimension();
		for (int i =0; i<imax; i++){
			for (int j=0; j<jmax; j++){
				if (M.get(i, j) != M.get(j, i)){
					System.out.println("Matrix is not symmetric at indexes (i, j) = ("+(i)+", "+(j)+" --> M(i, j) = "+(M.get(i, j))+" and M(j, i) = "+M.get(j, i));
					return false;
				}
			}
		}
		return true;
	}
	
}
