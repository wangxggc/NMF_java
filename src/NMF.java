import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import model.SparseVector;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.SeqBlas;
import cern.colt.matrix.linalg.SmpBlas;

/**
 * Non-negative Matrix Factorization(NMF)</br>
 * Decompose a matrix D to a product of two matrices U and V</br>
 * format this in mathematics is D=UV</br>
 * which U and V are all non-negative</br>
 * in other words, all elements in U and V are positive</br>
 * 2015-01-18</br>
 * 
 * @author wxg
 * @address BUAA BeiJing HaiDian China
 * @thanks clot.jar for matrix multiplication
 */
public class NMF {
	public DenseDoubleMatrix2D D;
	public DenseDoubleMatrix2D U;
	public DenseDoubleMatrix2D V;
	
	//precision for iterations of UpdateV
	public Double vEps;
	//precision for iterations of UpdateU
	public Double uEps;
	//define value that to be taken as zero
	public Double zero;
	//iteration times for UpdateU and UpdateV
	public int outerLoopMax;
	//inner iterations times of update of U and V
	public int innerLoopMax;
	//define number of rows in D and U
	public int mMax;
	//define number of columns in D and V
	public int nMax;
	//define number of columns in U and number of rows in V
	public int kMax;
	//define cores of cpu to be used for matrix multiplication
	public int cpu_s;
	//define step length for iterations that used gradient descent
	public double c;
	
	/**
	 * Constructor  
	 * @param m number of rows in D and U
	 * @param n number of columns in D and V
	 * @param k number of columns in U and number of rows in V
	 */
	public NMF(int m,int n,int k){
		this.D = new DenseDoubleMatrix2D(m, n);
		this.U = new DenseDoubleMatrix2D(m, k);
		this.V = new DenseDoubleMatrix2D(k, n);
		this.vEps = new Double(0.000001);
		this.uEps = this.vEps * Math.sqrt(this.mMax+this.nMax);
		this.zero = new Double(0.00000001);
		this.innerLoopMax = 100;
		this.outerLoopMax = 100;
		this.mMax = m;
		this.nMax = n;
		this.kMax = k;
		this.c = 0.0;
		this.cpu_s = 4;
	}
	/**
	 * Constructor with data file </br>
	 * File format:</br>
	 * 1  2  3</br>
	 * 2  3  4</br>
	 * 2 spaces between each two number</br>
	 * @param datafile data file path
	 * @param k number of columns in U and number of rows in V
	 */
	public NMF(String datafile,int k){
		this.vEps = new Double(0.000001);
		this.uEps = this.vEps;//this.vEps * Math.sqrt(this.mMax+this.nMax);
		this.zero = new Double(0.00000001);
		this.innerLoopMax = 100;
		this.outerLoopMax = 100;
		this.c = 0.0;
		List<double[]> doubleList = new ArrayList<double[]>();
		this.nMax = 0;
		try {
			BufferedReader reader = new BufferedReader(new FileReader(new File(datafile)));
			String line = null;
			while((line = reader.readLine())!=null){
				String[] values = line.split("  ");
				if(this.nMax<values.length-1){
					this.nMax = values.length-1;
				}
				double[] doubleLine = new double[values.length-1];
				for(int y=1;y<values.length;y++){
					double value = Double.parseDouble(values[y].trim());
					if(this.validate(value)){
						doubleLine[y-1]=value;
					}else{
						doubleLine[y-1]=0.0;
					}
				}
				doubleList.add(doubleLine);
			} 
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		double[][] doubleMatric = new double[doubleList.size()][];
		for(int i=0;i<doubleList.size();i++){
			doubleMatric[i] = doubleList.get(i);
		}
		this.D = new DenseDoubleMatrix2D(doubleMatric);
		this.mMax = doubleList.size();
		this.kMax = k;
		this.U = new DenseDoubleMatrix2D(mMax, kMax);
		this.V = new DenseDoubleMatrix2D(kMax, nMax);
	}
	/**
	 * Random initialize matrices A with columns normalized
	 * @param A
	 */
	public void RandomInit(DoubleMatrix2D A){
		for(int x=0;x<A.rows();x++){
			for(int y=0;y<A.columns();y++){
				int value = new Double(Math.random()*10).intValue();
				A.set(x, y, new Double(value));
			}
		}
		NormalizationByColumn(A);
	}
	/**
	 * Random initialize matrices U and V with columns normalized
	 */
	public void RandomInit(){
		RandomInit(U);
		RandomInit(V);
	}
	/**
	 * check element's validation in matrix </br>
	 * bigger than this.zero return true, else return false
	 * @param value value of element
	 * @return true for bigger than this.zero
	 */
	public boolean validate(double value){
		if(Math.abs(value)>this.zero){
			return true;
		}else{
			return false;
		}
	}
	/**
	 * multiply a matrix A with a number times
	 * @param A matrix
	 * @param times number
	 * @return A*times
	 */
	public DoubleMatrix2D Multi(DoubleMatrix2D A,double times){
		DoubleMatrix2D matric = new DenseDoubleMatrix2D(A.rows(), A.columns());
		for(int i=0;i<A.rows();i++){
			for(int j=0;j<A.columns();j++){
				double value = A.get(i, j)*times;
				if(this.validate(value))
					matric.set(i, j, value);
				else
					matric.set(i, j, 0.0);
			}
		}
		return matric;
	}
	/**
	 * matrices minus A-B
	 * @param A
	 * @param B
	 * @return
	 */
	public DoubleMatrix2D Minus(DoubleMatrix2D A,DoubleMatrix2D B){
		return this.LineCal(A, B, true);
	}
	/**
	 * matrices add A+B
	 * @param A
	 * @param B
	 * @return
	 */
	public DoubleMatrix2D Add(DoubleMatrix2D A,DoubleMatrix2D B){
		return this.LineCal(A, B, false);
	}
	/**
	 * used in Add and Minus
	 * @param A
	 * @param B
	 * @param sub
	 * @return
	 */
	private DoubleMatrix2D LineCal(DoubleMatrix2D A,DoubleMatrix2D B,boolean sub){
		DoubleMatrix2D matric = new DenseDoubleMatrix2D(A.rows(), A.columns());
		for(int i=0;i<A.rows();i++){
			for(int j=0;j<A.columns();j++){
				double value = 0.0;
				if(sub){
					value = A.get(i, j)-B.get(i, j);
				}else{
					value = A.get(i, j)+B.get(i, j);
				}
				if(this.validate(value))
					matric.set(i, j, value);
				else
					matric.set(i, j, 0.0);
			}
		}
		return matric;
	}
	/**
	 * perform non-negative mapping for matrix A
	 * @param A 
	 * @return A this all negative values been set to be 0 
	 */
	public DoubleMatrix2D NonNegativeMapping(DoubleMatrix2D A){
		for(int i=0;i<A.rows();i++){
			for(int j=0;j<A.columns();j++){
				if(A.get(i, j)<0){
					A.set(i, j, 0.0);
				}
			}
		}
		return A;
	}
	/**
	 * normalize columns in matrix A
	 * @param A
	 */
	public void NormalizationByColumn(DoubleMatrix2D A){
		for(int y=0;y<A.columns();y++){
			double value = 0.0;
			for(int x=0;x<A.rows();x++){
				value += A.get(x, y)*A.get(x, y);
			}
			if(!this.validate(value)){
				continue;
			}
			value = Math.sqrt(value);
			for(int x=0;x<A.rows();x++){
				if(this.validate(A.get(x, y)/value))
					A.set(x, y, A.get(x, y)/value);
				else
					A.set(x, y, 0);
			}
		}
	}
	/**
	 * Perform UpdateU
	 */
	public void UpdateU(){
		SmpBlas.allocateBlas(cpu_s, SeqBlas.seqBlas);

		DoubleMatrix2D VVT2 = init(V.rows(), V.rows());
		DoubleMatrix2D DVT2 = init(D.rows(), V.rows());
		SmpBlas.smpBlas.dgemm(false, true, 2, V, V, 0, VVT2);
		SmpBlas.smpBlas.dgemm(false, true, 2, D, V, 0, DVT2);
		
		DoubleMatrix2D UU;
		
		int t = 1;
		do{
			DoubleMatrix2D gradientU = init(D.rows(), V.rows());
			gradientU.assign(DVT2);
			SmpBlas.smpBlas.dgemm(false, false, 1, U, VVT2, -1,gradientU);
			
			double tt = this.c/Math.sqrt(t);
			DoubleMatrix2D SU = this.Minus(U, this.Multi(gradientU, tt));
			UU = this.Minus(SU, U);
			U.assign(NonNegativeMapping(SU));
			t++;
		}while(t<this.innerLoopMax && Algebra.DEFAULT.normF(UU)>this.uEps);
		System.out.print("UpdateU() iter times:"+t+"\t");
	}
	/**
	 * Perform UpdateV in Multi Thread
	 */
	public void UpdateVMultiThread(){
		DoubleMatrix2D R = new DenseDoubleMatrix2D(U.columns(), D.columns());
		SmpBlas.smpBlas.dgemm(true, false, 1, U, D, 0, R);
		DoubleMatrix2D S = new DenseDoubleMatrix2D(U.columns(), U.columns());
		SmpBlas.smpBlas.dgemm(true, false, 1, U, U, 0, S);
		
		int y;
		for(y=0;y<this.nMax;y++){
			//初始化
			UpdateVThread thread = new UpdateVThread();
			thread.y = y;
			thread.S = S;
			thread.R = R;
			thread.run();
		}
		System.out.print("UpdateV() number of Threads:"+(y)+" \t");
	}
	/**
	 * Perform decompose
	 */
	public void Decompose(){
		//初始化
		RandomInit(U);
		RandomInit(V);
		//迭代
		for(int i=0;i<this.outerLoopMax;i++){
			System.out.print("Iteration:"+i+"\t");
			UpdateU();
			UpdateVMultiThread();
			//System.gc();
			System.out.println();
		}
		save("result");
	}
	/**
	 * test validation of decompose
	 */
	public void TestDecompose(){
		//初始化
		DoubleMatrix2D UU;
		RandomInit(U);
		UU = U.copy();
		RandomInit(V);
		D = (DenseDoubleMatrix2D) Algebra.DEFAULT.mult(U,V);
		this.c = 0.1;
		//迭代
		RandomInit(U);
		RandomInit(V);
		for(int i=0;i<this.outerLoopMax;i++){
			System.out.print("Iteration:"+i+"\t");
			UpdateU();
			UpdateVMultiThread();
			//System.gc();
			System.out.println();
		}
		print(D);
		print(Algebra.DEFAULT.mult(U, V));
	}
	/**
	 * Save a matrix with filename
	 * @param A
	 * @param fileName
	 */
	public void Save(DoubleMatrix2D A,String fileName) {
		// TODO Auto-generated method stub
		FileWriter writer = null;
		try {
			writer = new FileWriter(new File(fileName),true);
			String output = "";
			writer.write(output);
			for(int x=0;x<A.rows();x++){
				output = "";
				for(int y=0;y<A.columns();y++){
					output += A.get(x, y)+"\t";
				}
				output += "\n";
				writer.write(output);
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally{
			if(writer!=null){
				try {
					writer.flush();
					writer.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	/**
	 * Save matrices D, U, V with filenames fileName+"D or U or V"
	 * @param fileName
	 */
	public void save(String fileName){
		Save(D,fileName+"D");
		Save(U,fileName+"U");
		Save(V,fileName+"V");
	}
	/**
	 * print matirx 
	 * @param A
	 */
	public void print(DoubleMatrix2D A){
		String output = "\n------------------------------------------Matric (xMax="+A.rows()+";yMax="+A.columns()+")------------------------------------------\n\n";
		for(int x=0;x<A.rows();x++){
			for(int y=0;y<A.columns();y++){
				output += A.get(x, y)+"\t";
			}
			output += "\n";
		}
		output += "\n";
		System.out.println(output);
	}
	/**
	 * thread for UpdateV
	 * @author wxg
	 *
	 */
	class UpdateVThread implements Runnable{
		int y;
		DoubleMatrix2D R;
		DoubleMatrix2D S;
		@Override
		public void run() {
			// TODO Auto-generated method stub
			SparseVector v;
			double tt = 1;
			do{
				v = getColumnsVector(V,y);
				for(int x=0;x<kMax;x++){
					double value = R.get(x, y);
					for(int k=0;k<kMax;k++){
						if(k!=x){
							value-=S.get(x, k)*V.get(k,y);
						}
					}
					value /= S.get(x, x);
					if(value>zero){
						V.set(x, y, value);
					}else{
						V.set(x, y,0.0);
					}
				}
				tt++;
				//V.print();
			}while(tt<=innerLoopMax&&v.SubBy(getColumnsVector(V,y)).getNorm_Max()>vEps);
		}		
	}
	public SparseVector getColumnsVector(DoubleMatrix2D A,int collumn){
		SparseVector sv = new SparseVector(A.rows());
		for(int x=0;x<A.rows();x++){
			if(validate(A.get(x, collumn))){
				sv.set(x+1, A.get(x,collumn));
			}
		}
		return sv;
	}
	public DoubleMatrix2D init(int m,int n){
		DoubleMatrix2D M = new DenseDoubleMatrix2D(m,n);
		for(int i=0;i<m;i++){
			for(int j =0;j<n;j++){
				M.set(i,j,0);
			}
		}
		return M;
	}
}
