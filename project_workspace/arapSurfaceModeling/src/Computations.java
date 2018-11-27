import java.util.ArrayList;
import Jcg.geometry.*;
import Jcg.polyhedron.*;
import Utils.Rotation_3;

public class Computations {
	
	double getOppositeAngle(Halfedge<Point_3> h) {
		Vertex<Point_3> va = h.getVertex();
		h = h.getNext();
		Vertex<Point_3> vb = h.getVertex();
		h = h.getNext();
		Vertex<Point_3> vc = h.getVertex();
		Point_3 a = va.getPoint();
		Point_3 b = vb.getPoint();
		Point_3 c = vc.getPoint();
		Vector_3 e1 = new Vector_3(b,a);
		Vector_3 e2 = new Vector_3(b,c);
		return (double) e1.innerProduct(e2);
	}
	
	double getWeight(Halfedge<Point_3> h) {
		// if boundary edge, opposite halfedge exists but face is null, fix the method if needed
		double alpha = getOppositeAngle(h);
		h = h.getOpposite();
		double beta = getOppositeAngle(h);
		return(0.5 * (1/Math.atan(alpha))* (1/Math.atan(beta)));
	}
	
	ArrayList<Vertex<Point_3>> getNeighbors(Halfedge<Point_3> h) {
		Halfedge<Point_3> e = h.getNext();
		Halfedge<Point_3> o = h.getOpposite();
		ArrayList<Vertex<Point_3>> neighbors = new ArrayList<Vertex<Point_3>>();
		while (e != o) {
			neighbors.add(e.getVertex());
			e = e.getOpposite().getNext();
		}
		return(neighbors);
	}

	Rotation_3 getHalfedgeRotation(Halfedge<Point_3> h) {
		return null;
	}
	
}