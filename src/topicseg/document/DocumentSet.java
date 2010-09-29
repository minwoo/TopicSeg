/*
 * Copyright (C) 2010 Cluster of Excellence, Univ. of Saarland
 * Minwoo Jeong (minwoo.j@gmail.com) is a main developer.
 * This file is part of "TopicSeg" package.
 * This software is provided under the terms of LGPL.
 */

package topicseg.document;

import java.util.ArrayList;
import java.util.Iterator;

import topicseg.model.TopicModel;

/**
 * unisaar.topicseg.document::DocumentSet.java
 *
 * @author minwoo
 */
public class DocumentSet extends ArrayList<Document> {
	
	protected Alphabet dict;
	//protected ArrayList<Document> doc_set;
	protected Annotation annotation;
	protected TopicModel globalTopic;
	//protected ArrayList<TopicSegment> localTopic;
	
	protected int docset_id; 
 	protected String common_filename; 
	public int maxGlobalCluster = 0;
	
	public DocumentSet() {
		this(new Alphabet());
	}
	
	public DocumentSet(Alphabet dict) {
		super();
		this.dict = dict;
		//doc_set = new ArrayList<Document>();
		annotation = new Annotation();
		globalTopic = new TopicModel();
		//localTopic = new ArrayList<TopicSegment>();
	}
	
//	public boolean add(Document doc) { 
//		super.add(doc); 
//		//localTopic.add(new TopicSegment());
//		return true;
//	}
	
	//public Document get(int index) { return doc_set.get(index); }
	
	//public int size() { return doc_set.size(); }
	
	public Annotation getAnnotation() {	return annotation; }
	
	public TopicModel getGlobalTopic() { return globalTopic; }
	
//	public TopicSegment getLocalTopic(int idx) { return localTopic.get(idx);	}
	
	public void setInfo(int id, String filename) { this.docset_id = id; this.common_filename = filename; }
	
	public String getFilename() { return common_filename; }
	
	public int getId() { return docset_id; }

	public Alphabet getAlphabet() { return dict; }
	
	public double[][][] getTFIDF() {
		int W = dict.size();
		int N = 0;
		double[][][] tfidf = new double[this.size()][][];
		
		ArrayList<SparseWordVector> wordMatrix = new ArrayList<SparseWordVector>();
		for (int d = 0; d < size(); d++) {
			Document doc = get(d);
			for (int t = 0; t < doc.size(); t++) {
				SparseWordVector termVec = doc.getSentence(t);
				wordMatrix.add(termVec);
				N++;
			}
		}
		
		for (int d = 0; d < size(); d++) {
			Document doc = get(d);
			tfidf[d] = new double[doc.size()][W];
			
			for (int t = 0; t < doc.size(); t++) {
				SparseWordVector termVec = doc.getSentence(t);
				int[] ids = termVec.getIds();
				double[] vals = termVec.getVals();
				double total = termVec.sum();
				
				for (int i = 0; i < vals.length; i++) {
					tfidf[d][t][ids[i]] = vals[i] / total;
					
					int count = 0;
					for (SparseWordVector vec : wordMatrix)
						if (vec.contains(ids[i])) count++;
					double idf = Math.log(N / (1 + (double)count));
					tfidf[d][t][ids[i]] *= idf;
				}
				
			}
		}
		return tfidf;		
	}
	
	public void tfidf() {
		int W = dict.size();
		int N = 0;
		
		ArrayList<SparseWordVector> wordMatrix = new ArrayList<SparseWordVector>();
		for (int d = 0; d < size(); d++) {
			Document doc = get(d);
			for (int t = 0; t < doc.size(); t++) {
				SparseWordVector termVec = doc.getSentence(t);
				wordMatrix.add(termVec);
				N++;
			}
		}
		
		for (int d = 0; d < size(); d++) {
			Document doc = get(d);
			double tfidf = 0;
			
			for (int t = 0; t < doc.size(); t++) {
				SparseWordVector termVec = doc.getSentence(t);
				int[] ids = termVec.getIds();
				double[] vals = termVec.getVals();
				double total = termVec.sum();
				
				for (int i = 0; i < vals.length; i++) {
					tfidf = vals[i] / total;
					
					int count = 0;
					for (SparseWordVector vec : wordMatrix)
						if (vec.contains(ids[i])) count++;
					double idf = Math.log(N / (1 + (double)count));
					tfidf *= idf;
					termVec.add(ids[i], tfidf);
				}
			}
		}
	}

	//public Iterator<Document> iterator() { return doc_set.iterator(); }
	
}
