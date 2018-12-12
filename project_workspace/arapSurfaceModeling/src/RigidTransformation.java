import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import Jcg.geometry.*;
import Jcg.polyhedron.*;
import Utils.Rotation_3;
import matrixPkg.Jama_Matrix;
import matrixPkg.Matrix;

public class RigidTransformation {
	
	public Polyhedron_3<Point_3> polyhedron3D;
	public HashMap<Vertex, Rotation_3> VertRotMap;
	public HashMap<Vertex<Point_3>, ArrayList<Halfedge<Point_3>>> globalNeighbors; // Hashmap<i, neighborsOfI>
	public HashMap<Vertex<Point_3>, HashMap<Vertex<Point_3>, Double>> weightij; // Hashmap<i, Hashmap<j, wij>>
	public Matrix L;
	public ArrayList<Integer> mobilePoints; // the points the user is allowed to move
	public ArrayList<Integer> fixedPoints; // the points the user wants to stay at a given position
	public Matrix p;
	public Matrix pPrime;

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
	
	public RigidTransformation() {//some constraints in parameters
		// Step 1 & 2
		globalNeighbors = new HashMap<Vertex<Point_3>, ArrayList<Halfedge<Point_3>>>();
		L = new Jama_Matrix(new Jama.Matrix(polyhedron3D.vertices.size(), polyhedron3D.vertices.size()));
		p = new Jama_Matrix(new Jama.Matrix(polyhedron3D.vertices.size(), 3));
		pPrime = new Jama_Matrix(new Jama.Matrix(polyhedron3D.vertices.size(), 3));
		
		for (Vertex<Point_3> v : polyhedron3D.vertices){
			weightij.put(v, new HashMap<Vertex<Point_3>, Double>());
			Point_3 vPoint = v.getPoint();
			p.set(v.index, 0, (double) vPoint.getX());
			p.set(v.index, 1, (double) vPoint.getY());
			p.set(v.index, 2, (double) vPoint.getZ());
			pPrime.set(v.index, 0, (double) vPoint.getX());
			pPrime.set(v.index, 1, (double) vPoint.getY());
			pPrime.set(v.index, 2, (double) vPoint.getZ());
		}
		for (Halfedge<Point_3> e : polyhedron3D.halfedges){
			globalNeighbors.put(e.getVertex(), Computations.getNeighbors(e));
			HashMap<Halfedge<Point_3>, Double> tmp = Computations.getWeightsArray(e, globalNeighbors.get(e.getVertex()));
			int i = e.getVertex().index; // hopefully, it is polyhedron3D.vertices.indexOf(e.getVertex()). If it doesn't work, just initialize a global array
			for (Halfedge<Point_3> f : tmp.keySet()){
				int j = f.getVertex().index;
				weightij.get(e.getVertex()).put(f.getVertex(), tmp.get(f));
				L.set(i, j, -tmp.get(f));
				L.set(i, i, L.get(i, i)+tmp.get(f));
			}
		}
	}
	
	public void arapIteration() {
		// Step 1 & 2
		for (Integer index : mobilePoints){
			Vertex<Point_3> v = polyhedron3D.vertices.get(index);
			Point_3 vPoint = v.getPoint();
			pPrime.set(v.index, 0, (double) vPoint.getX());
			pPrime.set(v.index, 1, (double) vPoint.getY());
			pPrime.set(v.index, 2, (double) vPoint.getZ());
			Halfedge<Point_3> e = v.getHalfedge();
			HashMap<Halfedge<Point_3>, Double> tmp = Computations.getWeightsArray(e, globalNeighbors.get(e.getVertex()));
			int i = v.index; // hopefully, it is polyhedron3D.vertices.indexOf(e.getVertex())
			L.set(i, i, 0);
			for (Halfedge<Point_3> f : tmp.keySet()){
				int j = f.getVertex().index;
				weightij.get(e.getVertex()).put(f.getVertex(), tmp.get(f));
				L.set(i, j, -tmp.get(f));
				L.set(i, i, L.get(i, i)+tmp.get(f));
			}
		}
		
		
	}
	
	// TODO
	public Point_3 getBi(Vertex i){
		Point_3 b = new Point_3(0., 0., 0.);
		Rotation_3 R_i = VertRotMap.get(i);
		for (Halfedge<Point_3> h : Computations.getNeighbors(i.getHalfedge())){
			Rotation_3 R_j = VertRotMap.get(h.vertex);
			Matrix R_ij = R_i.getMatrix().plus(R_j.getMatrix());
//			R_ij.times(weightij.get(i).);
		}
		return b;
	}
	
	
}
