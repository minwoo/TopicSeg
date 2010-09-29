/*
 * Copyright (C) 2010 Cluster of Excellence, Univ. of Saarland
 * Minwoo Jeong (minwoo.j@gmail.com) is a main developer.
 * This file is part of "TopicSeg" package.
 * This software is provided under the terms of LGPL.
 */

package topicseg.document;

import gnu.trove.map.hash.TIntIntHashMap;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.log4j.Logger;

/**
 * unisaar.topicseg.document::SparseWordVector.java
 *
 * @author minwoo
 */
public class SparseWordVector implements Serializable, Cloneable {
	
	private static final long serialVersionUID = 1L;
	private transient Logger logger = Logger.getLogger(SparseWordVector.class);

	protected TIntIntHashMap map = new TIntIntHashMap();
	protected ArrayList<WordId> wordVec;
	protected double totalWeight;
	
	public SparseWordVector() {
		wordVec = new ArrayList<WordId>();
		totalWeight = 0;
	}
	
	public SparseWordVector(int[] ids) {
		this();
		this.add(ids);
	}
	
	/**
	 * Add a term into vector with its weight
	 * @param id
	 * @param val
	 */
	public void add(int id, double val) {
		this.update(id, val, 1);
	}
	
	public void trim() {
		TIntIntHashMap newMap = new TIntIntHashMap();
		ArrayList<WordId> newWordVec = new ArrayList<WordId>();
		
		for (int id : map.keys()) {
			int index = map.get(id);
			double val = wordVec.get(index).getVal();
			if (val != 0) {
				WordId elem = new WordId(id, val);
				newWordVec.add(elem);
				newMap.put(id, newWordVec.size()-1);
			}
		}
		
		map = newMap;
		wordVec = newWordVec;
	}
	
	public void update(int id, double val, int sign) {
		if (map.containsKey(id)) {
			int index = map.get(id);
			double v = wordVec.get(index).update(sign * val);
		}
		else {
			WordId elem = new WordId(id, sign * val);
			wordVec.add(elem);
			map.put(id, wordVec.size()-1);
		}
		totalWeight += sign * val;
	}

	/**
	 * Add term array into vector with its weights
	 * @param ids
	 * @param vals
	 */
	public void add(int[] ids, double[] vals) {
		assert(ids.length == vals.length);
		for (int i = 0; i < ids.length; i++) 
			add(ids[i], vals[i]);
	}
	
	/**
	 * Add term count into vector
	 * @param ids
	 */
	public void add(int[] ids) {
		for (int id : ids) 
			add(id, 1);
	}
	
	public void add(SparseWordVector wordVector) {
		this.update(wordVector, 1);
	}
	
	public void update(SparseWordVector wordVector, int sign) {
		ArrayList<WordId> wordVec = wordVector.getList();
		for (WordId wordId : wordVec) {
			this.update(wordId.getId(), wordId.getVal(), sign);
		}
	}

	/**
	 * Check whether the current vector contains given term id 
	 * @param id
	 * @return
	 */
	public boolean contains(int id) {
		if (map.containsKey(id))
			return true;
		else
			return false;
	}
	
	/**
	 * Get term weight (freq) given term id 
	 * @param id
	 * @return term weight value
	 */
	public double getVal(int id) {
		if (map.containsKey(id))
			return wordVec.get(map.get(id)).getVal();
		else
			return 0;
	}
	
	/**
	 * Get term weight vector given term ids
	 * No index holds.
	 * @param ids
	 * @return term weight vector
	 */
	public double[] getVals(int[] ids) {
		double[] ret = new double[ids.length];
		
		for (int i = 0 ; i < ids.length; i++) {
			if (map.containsKey(ids[i]))
				ret[i] = wordVec.get(map.get(ids[i])).getVal();
			else
				ret[i] = 0;
		}
		
		return ret;
	}
	
	public int[] getIds() {
		int[] ret = new int[wordVec.size()];
		
		for (int i = 0; i < wordVec.size(); i++)
			ret[i] = wordVec.get(i).getId();
		
		return ret;
	}
	
	public double[] getVals() {
		double[] ret = new double[wordVec.size()];
		
		for (int i = 0; i < wordVec.size(); i++) 
			ret[i] = wordVec.get(i).getVal();
		
		return ret;
	}
	
	public ArrayList<Double> getVals2() {
		ArrayList<Double> vals = new ArrayList<Double>();
		
		for (int i = 0; i < wordVec.size(); i++) { 
			double v = wordVec.get(i).getVal();
			if (v > 0)
				vals.add(v);
		}
		return vals;
	}	
	
	/**
	 * Get term weight vector that is a dense vector and has T elements 
	 * @param total
	 * @return
	 */
	public double[] toArray(int T) {
		double[] ret = new double[T];
		Arrays.fill(ret, 0);
		
		for (int i = 0 ; i < wordVec.size(); i++) {
			int index = wordVec.get(i).getId();
			ret[index] = wordVec.get(i).getVal();
		}
		
		return ret;
	}
	
    public String toString() {
        StringBuffer sb = new StringBuffer();
        
        sb.append((int)sum()+"\n");
		for (int i = 0; i < wordVec.size(); i++) { 
			WordId w = wordVec.get(i);
			for (int j = 0; j < w.getVal(); j++) {
				sb.append((w.getId()+1)+"\n");
			}
		}
        
        return sb.toString();
    }
	
	public ArrayList<WordId> getList() { return wordVec; }
	public WordId getWordId(int index) { return wordVec.get(index);	}
	public double totalCount() { return totalWeight; }
	public double sum() {
		double ret = 0;
		for (WordId word : wordVec)
			ret += word.getVal();
		return ret;
	}

	public int size() { return wordVec.size(); }
	public int nonzero() {
		int ret = 0;
		for (WordId word : wordVec)
			if (word.getVal() > 0)
				ret++;
		return ret;
	}
	
	/**
	 * unisaar.topicseg.document::WordId
	 */
	public class WordId {
		public int id;
		public double val;	// val >= 0 (non-negative weight)
		public WordId() { id = -1; val = 0; }
		public WordId(int id, double val) { 
			this.id = id; 
			this.val = val; //Math.max(0, val);
		}
		public int getId() { return id; }
		public double getVal() { return val; }
		public double update(double diff) { 
			this.val += diff; 
			//this.val = Math.max(0, this.val);
			return this.val;
		}
	}
}
