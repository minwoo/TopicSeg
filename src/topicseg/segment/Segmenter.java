/*
 * Copyright (C) 2010 Cluster of Excellence, Univ. of Saarland
 * Minwoo Jeong (minwoo.j@gmail.com) is a main developer.
 * This file is part of "TopicSeg" package.
 * This software is provided under the terms of LGPL.
 */

package topicseg.segment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import topicseg.document.*;
import topicseg.document.Annotation.Label;
import topicseg.document.Document.SegmentBoundary;
import topicseg.utils.Eval;
import topicseg.utils.Option;

/**
 * unisaar.topicseg.algorithm::TopicSeg.java
 *
 * @author minwoo
 */
public abstract class Segmenter {

	public Segmenter() {
	}
	
	public abstract void initialize(Option config);
	public abstract void initializeRandom(int initSeed);
	public abstract void run(DocumentSet docset, double alpha, double beta);
	public abstract void run(Corpus corpus, double alpha, double beta);
	
    public double[] evaluate(Corpus corpus) {
    	double[] eval_metric = new double[10];
    	Arrays.fill(eval_metric, 0.0);
    	
    	int n_docs = 0;
        for (DocumentSet docSet : corpus) {
        	Annotation reference = docSet.getAnnotation();
        	int[][] refs = reference.toArray();
        	int K = reference.size();
        	Annotation hypothesis = new Annotation();
        	int nData = 0;
        	for (Document doc : docSet)
        		nData += doc.size();
        	
        	for (Document doc : docSet) {
        		int d = doc.getId();
        		int[] hyp = doc.getTopicLabel();
        		refs[d] = Eval.transform(refs[d], K);
        		hyp = Eval.transform(hyp, K);
        		eval_metric[0] += Eval.accuracy(refs[d], hyp);
        		eval_metric[1] += Eval.Pk(refs[d], hyp, K);
        		eval_metric[2] += Eval.WinDiff(refs[d], hyp, K);
        		eval_metric[4] += Eval.lengthRatio(refs[d], hyp);
        		n_docs++;
        		
        		for (int seg = 0; seg < doc.nSegment(); seg++) {
        			int start = doc.getSegmentStartPos(seg);
        			int end = doc.getSegment(seg) - 1;
        			int label = doc.getLabel(seg);
        			if (label >= 0)
        				hypothesis.add(label, d, new int[] {start, end});
        		}
        		
        	}
        	
        	ArrayList<ArrayList<Label>> refClusters = reference.getClusters();
        	ArrayList<ArrayList<Label>> hypClusters = hypothesis.getClusters();
        	double scores[] = Eval.fscore(refClusters, hypClusters);
        	double miscore = Eval.varInfo(refClusters, hypClusters, nData);
        	double rindex = Eval.randIndex(refClusters, hypClusters, nData);
        	eval_metric[5] += scores[0]; eval_metric[6] += scores[1]; eval_metric[7] += scores[2];
        	eval_metric[8] += miscore;
        	eval_metric[9] += rindex;
        }
    	
        for (int i = 0; i < 5; i++)
        	eval_metric[i] /= (double) n_docs;
        for (int i = 5; i < eval_metric.length; i++)
        	eval_metric[i] /= (double) corpus.size();
        
        double f1score = 2 * (eval_metric[5] * eval_metric[6]) / (eval_metric[5] + eval_metric[6]);
        eval_metric[7] = f1score;

    	return eval_metric;
    }
    
    // This method aims to produce the result in separate for each document in order to use for significance test.
    public void writeResult(Corpus corpus, String filename, boolean usePerDocEval) throws IOException {
    	
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));
    	
    	double[] eval_metric = new double[10];
    	Arrays.fill(eval_metric, 0.0);
    	
    	int n_docs = 0; int docId = 0;
        for (DocumentSet docSet : corpus) {
        	Annotation reference = docSet.getAnnotation();
        	int[][] refs = reference.toArray();
        	int K = reference.size();
        	Annotation hypothesis = new Annotation();
        	int nData = 0;
        	for (Document doc : docSet)
        		nData += doc.size();
        	
        	double[] docset_eval = new double[10];
        	Arrays.fill(docset_eval, 0.0);
        	for (Document doc : docSet) {
        		int d = doc.getId();
        		int[] hyp = doc.getTopicLabel();
        		refs[d] = Eval.transform(refs[d], K);
        		hyp = Eval.transform(hyp, K);
        		
        		double acc = Eval.accuracy(refs[d], hyp);
        		double pk = Eval.Pk(refs[d], hyp, K);
        		double wd = Eval.WinDiff(refs[d], hyp, K);
        		double len = Eval.lengthRatio(refs[d], hyp);
        		// per doc. evaluation
        		Annotation perDocRef = new Annotation();
        		Annotation perDocHyp = new Annotation();

        		eval_metric[0] += acc;
        		eval_metric[1] += pk;
        		eval_metric[2] += wd;
        		eval_metric[4] += len;
        		n_docs++;

        		for (int seg = 0; seg < doc.nSegment(); seg++) {
        			int start = doc.getSegmentStartPos(seg);
        			int end = doc.getSegment(seg) - 1;
        			int label = doc.getLabel(seg);
        			if (label >= 0) {
        				hypothesis.add(label, d, new int[] {start, end});
        				perDocHyp.add(label, d, new int[] {start, end});
        			}
        		}
        		
        		if (usePerDocEval) {
	        		for (ArrayList<Label> refLabels : reference.getClusters()) {
	        			for (Label refLabel : refLabels) {
	        				if (refLabel.docId == d) {
	        					perDocRef.add(refs[d][refLabel.start], d, new int[] {refLabel.start, refLabel.end});
	        				}
	        			}
	        		}
	            	double perDocScores[] = Eval.fscore(perDocRef.getClusters(), perDocHyp.getClusters());
	            	double perDocNvi = Eval.varInfo(perDocRef.getClusters(), perDocHyp.getClusters(), doc.size());
	            	double perDocRi = Eval.randIndex(perDocRef.getClusters(), perDocHyp.getClusters(), doc.size());
	
	    			writer.write(String.format("docName=%s docid=%d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.2f %.4f\n",
					docSet.getFilename(), doc.getId(),
					perDocScores[0], perDocScores[1], perDocScores[2], perDocNvi, perDocRi,  
					pk, wd,	len, acc));
        		}
        		else {
            		docset_eval[0] += Eval.accuracy(refs[d], hyp);
            		docset_eval[1] += Eval.Pk(refs[d], hyp, K);
            		docset_eval[2] += Eval.WinDiff(refs[d], hyp, K);
            		docset_eval[4] += Eval.lengthRatio(refs[d], hyp);
        		}
    			
        	}
        	
        	ArrayList<ArrayList<Label>> refClusters = reference.getClusters();
        	ArrayList<ArrayList<Label>> hypClusters = hypothesis.getClusters();
        	double scores[] = Eval.fscore(refClusters, hypClusters);
        	double nvi = Eval.varInfo(refClusters, hypClusters, nData);
        	double ri = Eval.randIndex(refClusters, hypClusters, nData);
        	eval_metric[5] += scores[0]; eval_metric[6] += scores[1]; eval_metric[7] += scores[2];
        	eval_metric[8] += nvi;
        	eval_metric[9] += ri;
        	
        	if (!usePerDocEval) {
	            for (int i = 0; i < docset_eval.length; i++)
	            	docset_eval[i] /= (double) docSet.size();
	            
				writer.write(String.format("docName=%s docid=%d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.2f %.4f\n",
						docSet.getFilename(), docId,
						scores[0], scores[1], scores[2], nvi, ri,  
						docset_eval[1], docset_eval[2],	docset_eval[4], docset_eval[0]));
        	}
			docId++;
        }
    	
        for (int i = 0; i < 5; i++)
        	eval_metric[i] /= (double) n_docs;
        for (int i = 5; i < eval_metric.length; i++)
        	eval_metric[i] /= (double) corpus.size();
        
        double f1score = 2 * (eval_metric[5] * eval_metric[6]) / (eval_metric[5] + eval_metric[6]);
		writer.write(String.format("!! average %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.2f %.4f\n",
				eval_metric[5], eval_metric[6], f1score, eval_metric[8], eval_metric[9], 
				eval_metric[1], eval_metric[2],	eval_metric[4], eval_metric[0]));
		
		writer.close();
	}

}
