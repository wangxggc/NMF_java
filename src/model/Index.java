package model;

/**
 * Ë÷Òý£¬·¶Î§10^9*10^9
 * @author wxg
 */
public class Index {
	public long x;
	public long y;
	
	public Index(long x,long y){
		this.x = x;
		this.y = y;
	}
	@Override
	public String toString() {
		// TODO Auto-generated method stub
		String sx = new Long(x+1000000).toString().substring(1);
		String sy = new Long(y+1000000).toString().substring(1);
		return sx+sy;
	}
	@Override
	public boolean equals(Object obj) {
		// TODO Auto-generated method stub
		return toString().equals(obj.toString());
	}
	@Override
	public int hashCode() {
		// TODO Auto-generated method stub
		return toString().hashCode();
	}
}
