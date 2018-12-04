import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import Jcg.geometry.*;
import Jcg.polyhedron.*;
import Utils.Rotation_3;
import matrixPkg.Matrix;

public class RigidTransformation {
	
	public Polyhedron_3<Point_3> polyhedron3D;
	public HashMap<Vertex, Rotation_3> VertRotMap;
	
	
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
		
		
	}
}
