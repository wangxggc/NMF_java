package model;

import java.util.HashMap;
import java.util.Map;

public class SparseVector {
	private Map<Index, Double> vector;
	private long nMax;
	private double zero;
	
	public SparseVector(long nMax){
		this.nMax = nMax;
		this.zero = 0.00000000000001;
		this.vector = new HashMap<Index, Double>();
	}
	
	public Double get(long x){
		if(x>this.nMax||x<=0){
			System.err.println("数组下标越界，请校对数组下标！");
			return null;
		}
		Double d = vector.get(new Index(x,x));
		if(d!=null){
			return d;
		}else{
			return 0.0;
		}
	}
	public void set(long x,Double value){
		if(x>this.nMax||x<=0){
			System.err.println("数组下标越界，请校对数组下标！");
			return;
		}
//		if(Math.abs(value)>this.eps){
			//保持稀疏性
			this.vector.put(new Index(x,x), value);
//		}
	}
	public Double remove(long x){
		return this.vector.remove(new Index(x,x));
	}
	//满足的才放入矩阵中存储，维持矩阵的稀疏性
	private boolean validate(double value){
		if(Math.abs(value)>this.zero){
			return true;
		}else{
			return false;
		}
	}
	public void clear(){
		this.vector.clear();
	}
	/**
	 * 随即初始化，0~1
	 */
	public void randomInit(){
		for(long x=1;x<=this.nMax;x++){
			int value = new Double(Math.random()*10).intValue();
			this.set(x, new Double(value));
		}
		this.Normalization();
	}
	/**
	 * 归一化
	 */
	public void Normalization(){
		double value = 0.0;
		for(long x=1;x<=this.nMax;x++){
			value += get(x)*get(x);
		}
		for(long x=1;x<=this.nMax;x++){
			if(this.validate(this.get(x)/value)){
				this.set(x, this.get(x)/value);
			}
		}
	}
	/**
	 * 获取拷贝
	 * @return
	 */
	public SparseVector getCopy(){
		return this.MultiplyBy(1.0);
	}
	/**
	 * 相加
	 * @param v
	 * @return
	 */
	public SparseVector AddBy(SparseVector v){
		return this.AddOrSub(v, false);
	}
	/**
	 * 相减
	 * @param v
	 * @return
	 */
	public SparseVector SubBy(SparseVector v){
		return this.AddOrSub(v, true);
	}
	/**
	 * 乘以一个数
	 * @param T
	 * @return
	 */
	public SparseVector MultiplyBy(double T){
		SparseVector sv = new SparseVector(this.nMax);
		for(long x=1;x<=this.nMax;x++){
			if(this.validate(this.get(x)*T)){
				sv.set(x, this.get(x)*T);
			}
		}
		return sv;
	}
	/**
	 * 乘以一个向量
	 * @param v
	 * @return
	 */
	public SparseVector MultiplyBy(SparseVector v){
		SparseVector sv = new SparseVector(this.nMax);
		if(this.nMax!=v.nMax){
			System.err.println("向量维数不同无法相乘！");
		}
		for(long x=1;x<=this.nMax;x++){
			if(this.validate(this.get(x)*v.get(x))){
				sv.set(x, this.get(x)*v.get(x));
			}
		}
		return sv;
	}
	private SparseVector AddOrSub(SparseVector v,boolean sub){
		if(v.nMax!=this.nMax){
			System.err.println("向量维数不同无法加减！");
		}
		SparseVector sv = new SparseVector(nMax);
		for(long x=1;x<=this.nMax;x++){
			double value = 0.0;
			if(sub){
				value = this.get(x)-v.get(x);
			}else{
				value = this.get(x)+v.get(x);
			}
			if(this.validate(value)){
				sv.set(x, value);
			}else{
				sv.remove(x);
			}
		}
		return sv;
	}
	/**
	 * 二范式
	 * @return
	 */
	public double getNorm_2(){
		double value = 0.0;
		for(long x=1;x<=this.nMax;x++){
			value += this.get(x)*this.get(x);
		}
		return Math.sqrt(value);
	}
	/**
	 * 无穷范式
	 * @return
	 */
	public double getNorm_Max(){
		double value = 0.0;
		for(long x=1;x<=this.nMax;x++){
			if(value<Math.abs(this.get(x))){
				value = Math.abs(this.get(x));
			}
		}
		return value;
	}
	
	public Map<Index, Double> getVector() {
		return vector;
	}
	public void setVector(Map<Index, Double> vector) {
		this.vector = vector;
	}
	public long getnMax() {
		return nMax;
	}
	public void setnMax(long nMax) {
		this.nMax = nMax;
	}
	public double getZero() {
		return zero;
	}
	public void setZero(double zero) {
		this.zero = zero;
	}
}
