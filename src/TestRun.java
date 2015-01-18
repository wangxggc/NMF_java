

import java.util.Date;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SmpBlas;

public class TestRun {
	public static void main(String[] args) {
		NMF decompose = new NMF(10, 8, 4);
		decompose.TestDecompose();
	}
}
