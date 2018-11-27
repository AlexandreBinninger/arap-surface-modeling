package Utils;

public class Pair<K, V>{
	
	private K first;
	private V second;
	
	public Pair(K fst, V snd){
		this.first=fst;
		this.second=snd;
	}
	
	public void setFirst(K fst){
		this.first=fst;
	}
	
	public void setSecond(V snd){
		this.second=snd;
	}
	
	public K getFirst(){
		return this.first;
	}
	
	public V getSecond(){
		return this.second;
	}
	
}
