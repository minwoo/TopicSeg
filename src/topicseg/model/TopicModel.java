/*
 * Copyright (C) 2010 Cluster of Excellence, Univ. of Saarland
 * Minwoo Jeong (minwoo.j@gmail.com) is a main developer.
 * This file is part of "TopicSeg" package.
 * This software is provided under the terms of LGPL.
 */

package topicseg.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Set;

import org.apache.log4j.Logger;

import topicseg.document.Alphabet;
import topicseg.document.SparseWordVector;

/**
 * unisaar.topicseg.model::TopicModel.java
 *
 * @author minwoo
 */
public class TopicModel {

	private transient Logger logger = Logger.getLogger(TopicModel.class);
	protected DirchletMultinomial dcm;
	protected double lmPrior;
	protected double tpPrior;

	protected class Table {
		public SparseWordVector wordVec;
		public int nDish;
		public double logProb;
		public int nRestaurant;
		
		public Table() {
			wordVec = new SparseWordVector();
			nDish = 0; logProb = 0; nRestaurant = 0;
		}
	}

	protected HashMap<Integer, Table> tableMap;
	protected PriorityQueue<Integer> emptyTableQueue;
	protected int nTable;
	protected int nTotalDish;
	protected int nVocab;
	
	public TopicModel() {
		dcm = new DirchletMultinomial();
		tableMap = new HashMap<Integer, Table>();
		emptyTableQueue =  new PriorityQueue<Integer>();
		nTable = 0;
	}
	
	public TopicModel(int W, double lmPrior, double tpPrior) {
		this();
		initialize(W, lmPrior, tpPrior);
	}
	
	public void initialize(int W, double lmPrior, double tpPrior) {
		tableMap.clear();
		emptyTableQueue.clear();
		nTable = 0;
		nVocab = W;
		this.lmPrior = lmPrior;
		this.tpPrior = tpPrior;
	}
	
	/*
	public void update(SparseWordVector start_cum_vec, SparseWordVector end_cum_vec, int id, int sign) {
		Table table = table_map.get(id);
		table.word_vec.update(end_cum_vec, sign);
		table.word_vec.update(start_cum_vec, -1 * sign);
		table.n_dish += sign;
		n_dish += sign;
	}
	*/
	
	public void update(SparseWordVector vec, int key, int val) {
		if (!tableMap.containsKey(key)) 
			return;
		
		Table table = tableMap.get(key);
		table.wordVec.update(vec, val);
		table.nDish += val;
		nTotalDish += val;
	}
	
	public void update(int val) {
		nTotalDish += val;
	}
	
	public double recomputeLogProb(int id, boolean restore) {
		double v = dcm.logDCM(this.getWordVector(id), nVocab, lmPrior);
		if (restore)
			tableMap.get(id).logProb = v;
		return v;
	}
	
	public double recomputeLogProb(int id) {
		return recomputeLogProb(id, false);
	}

	public void storeLogProb(int key, double val) {
		tableMap.get(key).logProb = val;
		//stampedLogProb.set(id, val); 
	}
	
	public double getLogProb(int key) {
		if (tableMap.containsKey(key))
			return tableMap.get(key).logProb; 
		return 0; 
	}
	
	public double getLogProb() {
		double v = 0;
		for (Table table : tableMap.values())
			v += table.logProb;
		return v;
	}
	
	public class WordId {
		public int id;
		public double val;	// val >= 0 (non-negative weight)
		public WordId() { id = -1; val = 0; }
		public WordId(int id, double val) { 
			this.id = id; 
			this.val = val; //Math.max(0, val);
		}
	}
	
	class WordIdCompare implements Comparator<WordId> {
		public int compare(WordId o1, WordId o2) {
			if (o1.val > o2.val) return -1;
			else if (o1.val < o2.val) return 1;
			else return 0;
		}
	}

	public String getTopicWord(int N, double prior, Alphabet dict) {
        StringBuffer sb = new StringBuffer();
        
        for (int tId : tableMap.keySet()) {
        	Table table = tableMap.get(tId);
        	int[] ids = table.wordVec.getIds();
        	double[] vals = table.wordVec.getVals();
        	double sum = table.wordVec.sum() + prior * table.wordVec.getVals2().size();
        	ArrayList<WordId> queue = new ArrayList<WordId>();
        	for (int i = 0; i < ids.length; i++) {
        		if (vals[i] <= 0.0) continue;
        		double prob = (vals[i] + prior) / sum;
        		
        		double norm = 0.0;
        		int K = 0;
        		for (int tId2 : tableMap.keySet()) {
        			Table table2 = tableMap.get(tId2);
        			double prob2 = Math.log((table2.wordVec.getVal(ids[i]) + prior) / (table2.wordVec.sum() + table2.wordVec.getVals2().size() * prior));
        			norm += prob2;
        			K++;
        		}
        		norm *= 1.0 / (double)K;
        		
        		double score = prob * (Math.log(prob) - norm); 
        		queue.add(new WordId(ids[i], score));
        	}
        	WordIdCompare wordIdCompare = new WordIdCompare();
        	Collections.sort(queue, wordIdCompare);
        	int j = 0;
        	WordId w = null;
        	sb.append("TOPIC " + tId +"\n");
        	while ( queue.size() > j && j < N) {
        		w = queue.get(j);
        		sb.append(dict.getObject(w.id) + " " + w.val + "\n");
        		j++;
        	}
        	sb.append("\n");
        }
        
        return sb.toString();
	}
	
	public int size() { 
		return nTable;
	}
	
	public SparseWordVector getWordVector(int key) {
		return tableMap.get(key).wordVec;
	}
	
	public void setWordVector(int key, SparseWordVector vec) {
		tableMap.get(key).wordVec = vec;
	}
	
	public int nVocab() { 
		return nVocab; 
	}
	
	public double nVocab(int key) {
		return tableMap.get(key).wordVec.size();
	}
	
	public int nDish() { 
		return nTotalDish; 
	}
	
	public int nDish(int key) {
		if (tableMap.containsKey(key))
			return tableMap.get(key).nDish;
		else
			return -1;
	}
	
	public double nWord(int key) {
		if (tableMap.containsKey(key))
			return tableMap.get(key).wordVec.sum();
		else
			return 0;
	}
	
	public double nWord() {
		int v = 0;
		
		for (Table table : tableMap.values()) 
			v += table.wordVec.sum();
		
		return v;
	}
	
	public int addTopic() {
		Table aNewTable = new Table();
		int tableId;
		
		if (emptyTableQueue.isEmpty())
			tableId = nTable;
		else
			tableId = emptyTableQueue.poll();
		
		tableMap.put(tableId, aNewTable);
		nTable++;
		
		return tableId;
	}
	
	public void removeTopic(int key) {
		if (!tableMap.containsKey(key))
			return;
		
		nTotalDish -= tableMap.get(key).nDish;
		tableMap.remove(key);
		emptyTableQueue.add(key);
		nTable--;
	}
	
	public void increaseRestaurant(int idx) {
		tableMap.get(idx).nRestaurant++;
	}
	
	public void decreseRestaurant(int idx) {
		tableMap.get(idx).nRestaurant--;
	}
	
	public int nRestaurant(int key) {
		if (tableMap.containsKey(key))
			return tableMap.get(key).nRestaurant;
		else
			return 0;
	}
	
	public Set<Integer> getIds() {
		return tableMap.keySet();
	}
	
	public double topicPropotional() {
		double[] counts = new double[tableMap.size()];
		if (counts.length <= 0) 
			return 0;
		
		int idx = 0;
		for (Table table : tableMap.values()) {
			counts[idx] = (double)table.nDish;
			idx++;
		}
			
		return dcm.logDCM(counts, tpPrior);
	}
	
	//-------------------------------------------------
	/*
	public void enterCustomer(SparseWordVector wordVec, int topicId) {
		if (topicId < wordCount.size())
			update(wordVec, topicId, 1);
		else {
			wordCount.add(new SparseWordVector());
			dishCount.add(1); nDish += 1;
			stampedLogProb.add(0.0);
		}
	}
	
	public void leaveCustomer(SparseWordVector wordVec, int topicId) {
		if (topicId < wordCount.size())
			update(wordVec, topicId, -1);
	}
	
	public void update(int topicId, int sign) {
		int nCustomer = dishCount.get(topicId);
		dishCount.set(topicId, nCustomer + sign);
		nDish += sign;
		assert(nDish >= 0);
	}*/

}
