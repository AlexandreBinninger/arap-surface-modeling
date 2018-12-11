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
	public HashMap<Vertex<Point_3>, HashMap<Vertex<Point_3>, Double>> weightij; // Hashmap<i, Hashmap<j, wij>>
	public Jama_Matrix L;
	
	
	public void RigidTransformation(){//some constraints in parameters
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
		
		L = new Jama_Matrix(new Jama.Matrix(polyhedron3D.vertices.size(), polyhedron3D.vertices.size()));
		
		// Step 1 & 2
		
		for (Vertex<Point_3> v : polyhedron3D.vertices){
			weightij.put(v, new HashMap<Vertex<Point_3>, Double>());
		}
		for (Halfedge<Point_3> e : polyhedron3D.halfedges){
			HashMap<Halfedge<Point_3>, Double> tmp = Computations.getWeightsArray(e, Computations.getNeighbors(e));
			int i = e.getVertex().index;
			for (Halfedge<Point_3> f : tmp.keySet()){
				int j = f.getVertex().index;
				weightij.get(e.getVertex()).put(f.getVertex(), tmp.get(f));
				L.set(i, j, -tmp.get(f));
				L.set(i, i, L.get(i, i)+tmp.get(f));
			}
		}
		
		//Step 3
		
		//Step 4
		
		//Step 5
		
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
