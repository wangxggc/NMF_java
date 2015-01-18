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
			System.err.println("�����±�Խ�磬��У�������±꣡");
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
			System.err.println("�����±�Խ�磬��У�������±꣡");
			return;
		}
//		if(Math.abs(value)>this.eps){
			//����ϡ����
			this.vector.put(new Index(x,x), value);
//		}
	}
	public Double remove(long x){
		return this.vector.remove(new Index(x,x));
	}
	//����Ĳŷ�������д洢��ά�־����ϡ����
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
	 * �漴��ʼ����0~1
	 */
	public void randomInit(){
		for(long x=1;x<=this.nMax;x++){
			int value = new Double(Math.random()*10).intValue();
			this.set(x, new Double(value));
		}
		this.Normalization();
	}
	/**
	 * ��һ��
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
	 * ��ȡ����
	 * @return
	 */
	public SparseVector getCopy(){
		return this.MultiplyBy(1.0);
	}
	/**
	 * ���
	 * @param v
	 * @return
	 */
	public SparseVector AddBy(SparseVector v){
		return this.AddOrSub(v, false);
	}
	/**
	 * ���
	 * @param v
	 * @return
	 */
	public SparseVector SubBy(SparseVector v){
		return this.AddOrSub(v, true);
	}
	/**
	 * ����һ����
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
	 * ����һ������
	 * @param v
	 * @return
	 */
	public SparseVector MultiplyBy(SparseVector v){
		SparseVector sv = new SparseVector(this.nMax);
		if(this.nMax!=v.nMax){
			System.err.println("����ά����ͬ�޷���ˣ�");
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
			System.err.println("����ά����ͬ�޷��Ӽ���");
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
	 * ����ʽ
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
	 * ���ʽ
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
