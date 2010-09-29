/*
 * Copyright (C) 2010 Cluster of Excellence, Univ. of Saarland
 * Minwoo Jeong (minwoo.j@gmail.com) is a main developer.
 * This file is part of "TopicSeg" package.
 * This software is provided under the terms of LGPL.
 */

package topicseg.segment;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import org.apache.commons.math.stat.clustering.*;
import org.apache.log4j.Logger;

import cern.jet.random.Uniform;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

import topicseg.document.*;
import topicseg.model.DirchletMultinomial;
import topicseg.utils.Eval;
import topicseg.utils.Option;

/**
 * unisaar.topicseg.algorithm::SingleSeg.java
 * Many stuffs in this code is based on BayesSeg package (http://groups.csail.mit.edu/rbg/code/bayesseg/).
 * 
 * A pipeline system that first performs BayesSeg for each document in isolation and then aligns the predicted segments. 
 *
 * @author minwoo
 */
public class SingleSeg extends Segmenter {

	private transient Logger logger = Logger.getLogger(SingleSeg.class);
	
	protected DirchletMultinomial dcm;
	protected Corpus mCorpus;
	
	int maxKMeansIter;
	boolean useVerbose;
	int initRandSeed;
	RandomEngine randomEngine;
	
	public SingleSeg() {
		dcm = new DirchletMultinomial();
	}

	public void initialize(Option config) {
		this.initRandSeed = config.contains("rand.seed") ? config.getInteger("rand.seed") : 0;
		this.maxKMeansIter = config.contains("kmeans.maxIter") ? config.getInteger("kmeans.maxIter") : 1000;
		this.useVerbose = config.contains("print.verbose") ? config.getBoolean("print.verbose") : true;
		
        // initializing random engine
        if (initRandSeed > 0) 
	        randomEngine = new MersenneTwister(initRandSeed);
        else {
	        Date thedate = new Date();
	        randomEngine = new MersenneTwister(thedate);
        }
	    Uniform.staticSetRandomEngine(randomEngine);
	    
        if (config.contains("exp.log")) {
			try {
				Option.addFileLogger(logger, config.getString("exp.log"));
			}
			catch (IOException e) {
				e.printStackTrace();
			}
        }
	}

	@Override
	public void run(DocumentSet docset, double alpha, double beta) {
		Annotation annotation = docset.getAnnotation();
		int[][] refs = annotation.toArray();
		int nRefCluster = annotation.size();
		
		for (int d = 0; d < docset.size(); d++) {
			Document doc = docset.get(d);
			int K = Eval.numSeg(refs[doc.getId()], 0, refs[doc.getId()].length-1) + 1;
			HashSet<Integer> localSegIds = new HashSet<Integer>();
			for (int i = 0; i < refs[doc.getId()].length; i++)
				if (refs[doc.getId()][i] < 0) localSegIds.add(refs[doc.getId()][i]);
			nRefCluster += localSegIds.size();
			
			// segment
			dpseg(doc, alpha, K);
		}
		
		nRefCluster = annotation.size();
		// alignment
		kmeans(docset, alpha, nRefCluster);
	}

	@Override
	public void run(Corpus corpus, double alpha, double beta) {
    	long timer = System.currentTimeMillis(); // timer

		for (int s = 0; s < corpus.size(); s++) {
			DocumentSet docset = corpus.get(s);
			run(docset, alpha, beta);
		}
		
	    // print the reference and hypothesis
	    if (useVerbose) {
	        for (DocumentSet docSet : corpus) {
	        	Annotation annotation = docSet.getAnnotation();
	        	if (annotation.size() > 0) {
		        	int[][] refs = annotation.toArray();
		        	for (Document doc : docSet) {
		        		int[] hyp = doc.getTopicLabel();
		                StringBuffer sb = new StringBuffer();
		        		sb.append("[" + docSet.getId() + "," + doc.getId() + "] ");
		        		for (int x = 0; x < hyp.length; x++)
		        			sb.append(refs[doc.getId()][x] + ":" + hyp[x] + " ");
		        		logger.info(sb);
		        	}
	        	}
	        }
	    }
		
		print(corpus, timer);
	}
	
	private void print(Corpus corpus, long timer) {
	    // print out
		double[] evalMetrics = evaluate(corpus);
		double tookTime = ((long)System.currentTimeMillis() - timer);
		logger.info(String.format("SingleSeg || loglike=%e p=%.4f/r=%.4f/f1=%.4f vi=%.4f ri=%.4f pk=%.4f wd=%.4f time=%.2f nSeg=%.2f",
				0.,
				evalMetrics[5], evalMetrics[6], evalMetrics[7],
				evalMetrics[8], evalMetrics[9],
				evalMetrics[1], evalMetrics[2],
				tookTime / 1000,
				evalMetrics[4]
				));
	}

	/**
	 * Segmentation with dynamic programming for isolate documents
	 * @param docset
	 * @param prior
	 */
	public void dpseg(Document doc, double prior) {
		
		int nSent = doc.size();
		
        // computing score matrix
        double[][] segLLs = new double[nSent+1][nSent+1];
        for (int t = 0; t <= nSent; t++) {
            for (int t2 = 0; t2 < t; t2++) {
            	//seglls[t][t2] += logprob(cumSums, t2+1, t, prior);
            	segLLs[t][t2] += computeLogprob(doc, t2, t, prior);
            }
        }
        
        // dynamic programming
        double C[] = new double[nSent+1]; 
        int B[] = new int[nSent+1]; 
        C[0] = 0; B[0] = 0;
        for (int t = 1; t <= nSent; t++){
            double bestValue = -Double.MAX_VALUE; 
            int bestIndex = -1;
            for (int t2 = 0; t2 < t; t2++) {
                double score = C[t2] + segLLs[t][t2];
                if (score > bestValue) { 
                	bestValue = score; 
                	bestIndex = t2; 
                }
            }
            C[t] = bestValue; 
            B[t] = bestIndex;
        } 
        
        // back tracking
        int nextSegPt = B[B.length-1];
        ArrayList<Integer> segpts = new ArrayList<Integer>();
        segpts.add(0);
        while (nextSegPt > 0) {
            segpts.add(nextSegPt);
            nextSegPt = B[nextSegPt];
        }
        Collections.sort(segpts);
        int label = 0;
        int[] hiddentopic = new int[doc.size()];
		for (int i = 0; i < segpts.size() - 1; i++) {
			int start = segpts.get(i);
			int end = segpts.get(i+1);
			for (int x = start; x < end; x++)
				hiddentopic[x] = label;
			label++;
		}
		
		int start = segpts.get(segpts.size()-1) ;
		int end = nSent;
		for (int x = start; x < end; x++)
			hiddentopic[x] = label;
		doc.setTopicLabel(hiddentopic);
		//System.out.println(segpts);
	}
	
	public void kmeans(DocumentSet docSet, double prior, int nRefCluster) {
//		double[][][] phi = new double[docSet.size()][][];
		int W = docSet.getAlphabet().size();
		double[][][] tfidf = docSet.getTFIDF();
		Collection<SegmentPoint> data = new ArrayList<SegmentPoint>();

		int maxK = 0;
		
		for (int i = 0; i < docSet.size(); i++) {
			Document doc = docSet.get(i);
			int nSeg = doc.nSegment();
//			phi[i] = new double[nSeg][];

			for (int j = 0; j < nSeg; j++) {
				int start = doc.getSegmentStartPos(j);
				int end = doc.getSegment(j);
//				phi[i][j] = probLM(doc, start, end, prior);
				double[] phi = new double[W];
				for (int k = start; k < end; k++)
					for (int x = 0; x < W; x++)
						phi[x] += tfidf[i][k][x];
				SegmentPoint point = new SegmentPoint(i, j, start, end, phi);
//				SegmentPoint point = new SegmentPoint(i, j, start, end, phi[i][j]);
				data.add(point);
			}
			
			//if (nSeg > maxK)
			maxK += nSeg;
		}
		
		int k = Math.min(nRefCluster, maxK);
		KMeansPlusPlusClusterer<SegmentPoint> clusterer = new KMeansPlusPlusClusterer<SegmentPoint>(new Random(initRandSeed));
		List<Cluster<SegmentPoint>> clusters = clusterer.cluster(data, k, 1000);
		int clusterId = 0;
		for (Cluster<SegmentPoint> clu : clusters) {
			List<SegmentPoint> points = clu.getPoints();
			for (SegmentPoint p : points) {
				Document doc = docSet.get(p.docId);
				for (int pos = p.start; pos < p.end; pos++)
					doc.setTopicLabel(pos, points.size() > 1 ? clusterId : -1 * clusterId);
				//System.out.println(p.docId + " " + p.segId + " " + p.start + " " + p.end);
			}
			//System.out.println();
			clusterId++;
		}
		
		for (Document doc : docSet)
			doc.resetTopicLabel();
	}
	
	public class SegmentPoint implements Clusterable<SegmentPoint> {
		
		public int docId;
		public int segId;
		public int start;
		public int end;
		private double[] prob;
		
		public SegmentPoint(int docId, int segId, int start, int end, double[] prob) {
			this.docId = docId;
			this.segId = segId;
			this.start = start;
			this.end = end;
			this.prob = prob;
		}
		
		public double[] getProb() {
			return prob;
		}
		
		public SegmentPoint centroidOf(Collection<SegmentPoint> points) {
			double[] centroid = new double[getProb().length];
			for (SegmentPoint p : points) {
				for (int i = 0; i < centroid.length; i++) {
					centroid[i] += p.getProb()[i];
				}
			}
			for (int i = 0; i < centroid.length; i++)
				centroid[i] /= points.size();
			
			return new SegmentPoint(-1, -1, -1, -1, centroid);
		}

		public double distanceFrom(SegmentPoint p) {
//			double v = (klDiv(p.getProb(), getProb()) + klDiv(getProb(), p.getProb())) / 2.0;
//			// crude heuristics
//			if (p.docId == docId) {
//				if (Math.abs(p.segId - segId) < 2) 
//					v *= 1.5;
//				else
//					v *= 5.0;
//			}
//			return v;
			double val1 = 0., val2 = 0., val3 = 0.;
			double[] v1 = getProb();
			double[] v2 = p.getProb();
			for (int i = 0; i < v1.length; i++) {
				val1 += v1[i] * v2[i];
				val2 += v1[i]; val3 += v2[i];
			}
			
			double dist = val1 > 0 ? 1.0 - val1 / (Math.sqrt(val2) * Math.sqrt(val3)) : 1.0;
			return dist;
			
		}
		
	}
	
	/*----------------------------------------------------------------------
	 * Automatically learning of prior hyperparameters
	 * TODO: Remove these methods  
	 *---------------------------------------------------------------------- */ 

//	public double[] runEM(Corpus corpus, double[] priors) {
//		mCorpus = corpus; // keep the pointer of data set for EM iteration
//		
//		PriorOptimizer opt = new PriorOptimizer();
//		opt.m_eps = 1e-5; opt.m_num_corrections = 6;
//		double improvement = Double.MAX_VALUE;
//		double old_ll = Double.MAX_VALUE;
//		int ctr = 0;
//		double[] params = new double[priors.length];
//		for (int i = 0; i < params.length; i++) params[i] = Math.log(priors[i]);
//		segment(corpus, Math.exp(params[0]));
//		
//		do {
//			opt.setEstimate(params);
//			opt.setMaxIteration(200);
//			opt.setDebug(true);
//			try {
//				opt.findArgmin();
//			} catch (Exception e) {
//				e.printStackTrace();
//			}
//			params = opt.getVarbValues();
//			segment(corpus, Math.exp(params[0]));
//			
//			improvement = old_ll - opt.getMinFunction();
//			old_ll = opt.getMinFunction();
//			
//		} while (improvement > 0 && ++ctr < 20);
//		
//		double[] exp_params = new double[params.length];
//		for (int i = 0; i < params.length; i++) exp_params[i] = Math.exp(params[i]);
//		return exp_params;
//	}
	
	
	// Dense matrix version
	// TODO: Remove this method
    protected double logprob(double[][] cumSums, int start, int end, double prior) {
        double ret = 0;
        int W = cumSums[0].length;
        
        if (prior == 0) 
        	return -Double.MAX_VALUE;
        
        double[] counts = new double[W];
        double N = 0;
        for (int i = 0; i < W; i++) {
        	counts[i] = cumSums[end][i] - cumSums[start-1][i];
        	N += counts[i];
        }
        
        ret = dcm.logDCM(counts, prior);
        return ret;
    }

    protected double computeLogprob(Document doc, int start, int end, double prior) {
        double ret = 0;
        int W = doc.getAlphabet().size(); 
        
        if (prior == 0) 
        	return -Double.MAX_VALUE;
        
        SparseWordVector wordVector = doc.getMergedWordVector(start, end);
        ret = dcm.logDCM(wordVector, W, prior);
        
        return ret;
    }
    
    protected double[] probLM(Document doc, int start, int end, double prior) {
        int W = doc.getAlphabet().size();
        SparseWordVector wordVector = doc.getMergedWordVector(start, end);
        
        
    	double S = prior * W;
    	double[] prob = new double[W];
    	Arrays.fill(prob, 0);
    	
        double N = wordVector.sum();
        int[] ids = wordVector.getIds();
    	double[] vals = wordVector.getVals();
    	
        for (int i = 0; i < ids.length; i++) {
        	//double v = gammaCache.logGamma(vals[i] + prior) - gammaCache.logGamma(N+S);
        	prob[ids[i]] = (vals[i] + prior) / (S + N); 
        	//prob[ids[i]] = Math.exp(v);
        }

       	return prob;
    }
    
    public static double klDiv(double distrib1[], double distrib2[]) {
        double total = 0;
        double alpha = .99;
        for (int i = 0; i < distrib1.length; i++){
            double q = alpha * distrib1[i] + (1.0-alpha) * distrib2[i];
            double r = distrib2[i];
            if (r != 0)
                total += r * (Math.log(r) - Math.log(q));
        }
        return total;
    }

    protected double computeGradient(double[][] cumSums, int start, int end, double prior) {
        double ret = 0;
        int W = cumSums[0].length;
        
        if (prior == 0) 
        	return Double.MAX_VALUE;
        
        double[] counts = new double[W];
        double N = 0;
        for (int i = 0; i < W; i++) {
        	counts[i] = cumSums[end][i] - cumSums[start-1][i];
        	N += counts[i];
        }
        
        ret = dcm.gradientDCM(counts, N, prior);
        return ret;
    }

    
    protected double computeGradient(Document doc, int start, int end, double prior) {
    	double ret = 0;
    	int W = doc.getAlphabet().size();
    	
    	if (prior == 0)
    		return Double.MAX_VALUE;
    	
        SparseWordVector wordVector = doc.getMergedWordVector(start, end);
        ret = dcm.gradientDCM(wordVector, W, prior);
        
        return ret;
    }
    
    /**
     * Compute log-likelihood of all documents in corpus
     * @param logpriors
     * @return
     */
//    protected double computeDataLogLL(double[] logpriors) {
//    	if (logpriors.length <= 0 || logpriors[0] == 0)
//    		return -Double.MAX_VALUE;
//    	
//    	double ret = 0;
//    	double prior = Math.exp(logpriors[0]); // for now, we have only one parameter to be estimated
//    	
//		for (int s = 0; s < mCorpus.size(); s++) {
//			DocumentSet docset = mCorpus.get(s);
//			for (int d = 0; d < docset.size(); d++) {
//				Document doc = docset.get(d);
//				int nSent = doc.size();
//				ArrayList<Integer> segpts = doc.getSegment();
//
//				//ret += computeLogprob(doc, 0, segpts.get(0), prior);
//				for (int i = 0; i < segpts.size() - 1; i++) {
//					int start = segpts.get(i);
//					int end = segpts.get(i+1);
//					ret += computeLogprob(doc, start, end, prior);
//				}
//				int start = segpts.get(segpts.size()-1) ;
//				int end = nSent;
//				ret += computeLogprob(doc, start, end, prior);
//			}
//		}
//		
//    	
//    	return ret;
//    }
    
    /**
     * Compute the gradient values
     * It requests an array of return values because our LBFGS optimizer works with this format.  
     * @param logpriors
     * @return
     */
//    protected double[] computeDataGradient(double[] logpriors) {
//    	double[] ret = new double[logpriors.length];
//    	double prior = Math.exp(logpriors[0]);
//    	
//		for (int s = 0; s < mCorpus.size(); s++) {
//			DocumentSet docset = mCorpus.get(s);
//			for (int d = 0; d < docset.size(); d++) {
//				Document doc = docset.get(d);
//				int nSent = doc.size();
//				ArrayList<Integer> segpts = doc.getSegment();
//
//				//ret[0] += computeGradient(doc, 0, segpts.get(0), prior);
//				for (int i = 0; i < segpts.size() - 1; i++) {
//					int start = segpts.get(i);
//					int end = segpts.get(i+1);
//					ret[0] += computeGradient(doc, start, end, prior);
//				}
//				int start = segpts.get(segpts.size()-1);
//				int end = nSent;
//				ret[0] += computeGradient(doc, start, end, prior);
//			}
//		}    	
//    	
//    	return ret;
//    }
    
    /**
     * LBFGS Optimizer 
     */
//    protected class PriorOptimizer extends LBFGSWrapper {
//    	public PriorOptimizer() {
//    		super(1); // # parameters = 1
//    	}
//    		
//		public double objectiveFunction(double[] params) {
//			return -computeDataLogLL(params);
//		}
//		
//		public double[] evaluateGradient(double[] params) {
//			double[] ret = computeDataGradient(params);
//			for (int i = 0; i < ret.length; i++) ret[i] = -ret[i];
//			return ret;
//		}
//
//    }

	public void dpseg(Document doc, double prior, int N) {
		int nSent = doc.size();
		int nWord = doc.getAlphabet().size();
        
        // computing score matrix
        double[][] segLLs = new double[nSent+1][nSent+1];
        for (int t = 0; t <= nSent; t++) {
            for (int t2 = 0; t2 < t; t2++) {
            	segLLs[t][t2] += computeLogprob(doc, t2, t, prior);
            }
        }
        
        // dynamic programming - search
        double C[][] = new double[N+1][nSent+1]; 
        int B[][] = new int[N+1][nSent+1]; 
        C[0][0] = 0; B[0][0] = 0;
    	for (int t = 1; t <= nSent; t++) {
    		C[0][t] = -Double.MAX_VALUE;
    		B[0][t] = 1;
    	}
        for (int i = 1; i <= N; i++) {
        	for (int t = 0; t < i; t++) {
        		C[i][t] = -Double.MAX_VALUE;
        		B[i][t] = -1;
        	}
	        for (int t = 1; t <= nSent; t++){
	            double bestValue = -Double.MAX_VALUE; 
	            int bestIndex = -1;
	            for (int t2 = 0; t2 < t; t2++) {
	                double score = C[i-1][t2] + segLLs[t][t2];
	                if (score >= bestValue) { 
	                	bestValue = score; 
	                	bestIndex = t2; 
	                }
	            }
	            C[i][t] = bestValue; 
	            B[i][t] = bestIndex;
	        }
        }
        
        
        ArrayList<Integer> segpts = new ArrayList<Integer>();
        for (int k = 0; k < N; k++) segpts.add(-1);
        segpts.set(N-1, B[N][nSent]);
        for (int k = N-1; k > 0; k--) {
        	segpts.set(k-1, B[k][segpts.get(k)]);
        }

        int label = 0;
        int[] hiddentopic = new int[doc.size()];
		for (int i = 0; i < segpts.size() - 1; i++) {
			int start = segpts.get(i);
			int end = segpts.get(i+1);
			for (int x = start; x < end; x++)
				hiddentopic[x] = label;
			label++;
		}
		int start = segpts.get(segpts.size()-1) ;
		int end = nSent;
		for (int x = start; x < end; x++)
			hiddentopic[x] = label;
		doc.setTopicLabel(hiddentopic);
		//System.out.println(doc.getId() + " " + segpts);
	}

	public void initializeRandom(int initSeed) {
		initRandSeed = initSeed;		
	}


}
