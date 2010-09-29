/*
 * Copyright (C) 2010 Cluster of Excellence, Univ. of Saarland
 * Minwoo Jeong (minwoo.j@gmail.com) is a main developer.
 * This file is part of "TopicSeg" package.
 * This software is provided under the terms of LGPL.
 */

package topicseg.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.log4j.Logger;

import topicseg.document.SparseWordVector;

import cern.jet.stat.Gamma;

/**
 * unisaar.topicseg.model::LanguageModel.java
 * An implementation of multinomial language model
 * 
 * @author minwoo
 */
public class DirchletMultinomial {
	
	private transient Logger logger = Logger.getLogger(DirchletMultinomial.class);
	
	protected GammaCache gammaCache;
	
	public DirchletMultinomial() {
		gammaCache = new GammaCache();
	}

    /**
     * Compute log probability of document (in fact segment) using DCM
     * This method takes sparse word vector, and considers only non-zero terms.
     * @param wordVector
     * @param nWord
     * @param prior
     * @return
     */
    public double logDCM(SparseWordVector wordVector, int nWord, double prior) {
    	double S = prior * nWord;
        //int[] ids = wordVector.getIds();
        //double[] vals = wordVector.getVals();
    	ArrayList<Double> vals = wordVector.getVals2();
        
        double ret = gammaCache.logGamma(S) - vals.size() * gammaCache.logGamma(prior);

        double N = 0;
        for (int i = 0; i < vals.size(); i++) {
        	N += vals.get(i);
            ret += gammaCache.logGamma(vals.get(i) + prior);
        }
        ret -= gammaCache.logGamma(N + S);
        
        return ret;
    }
    
    /**
     * Compute gradient of log probability of document (in fact segment) using DCM
     * This method takes sparse word vector, and considers only non-zero terms.
     * @param wordVector
     * @param nWord
     * @param prior
     * @return
     */
    public double gradientDCM(SparseWordVector wordVector, int nWord, double prior) {
    	double S = prior * nWord;
        //int[] ids = wordVector.getIds();
        double[] vals = wordVector.getVals();
        
        double ret = nWord * gammaCache.digamma(S) - wordVector.nonzero() * gammaCache.digamma(prior);

        double N = 0;
        for (int i = 0; i < vals.length; i++) {
        	N += vals[i];
            ret += gammaCache.digamma(vals[i] + prior);
        }
        ret -= nWord * gammaCache.digamma(N + S);
        
        return ret;
    	
    }
    
    /**
     * Compute log probability of document (in fact segment) using DCM
     * This method takes two arrays for index and values (term counts) that represent non-zero terms.
     * @param wordVector
     * @param nWord
     * @param prior
     * @return
     */
    public double logDCM(int[] ids, double[] vals, int nWord, double prior) {
    	double S = prior * nWord;
        
        double ret = gammaCache.logGamma(S) - ids.length * gammaCache.logGamma(prior);

        double N = 0;
        for (int i = 0; i < ids.length; i++) {
        	N += vals[i];
            ret += gammaCache.logGamma(vals[i] + prior);
        }
        ret -= gammaCache.logGamma(N + S);
        
        return ret;
    }

    /**
     * Compute gradient of log probability of document (in fact segment) using DCM
     * This method takes two arrays for index and values (term counts) that represent non-zero terms.
     * @param wordVector
     * @param nWord
     * @param prior
     * @return
     */
    public double gradientDCM(int[] ids, double[] vals, int nWord, double prior) {
    	double S = prior * nWord;
        
        double ret = nWord * gammaCache.digamma(S) - ids.length * gammaCache.digamma(prior);

        double N = 0;
        for (int i = 0; i < ids.length; i++) {
        	N += vals[i];
            ret += gammaCache.digamma(vals[i] + prior);
        }
        ret -= nWord * gammaCache.digamma(N + S);
        
        return ret;
    	
    }
    
    
    /**
     * Compute log probability of document (in fact segment) using DCM
     * This method takes dense form of word vector, because it would be sometime efficient.
     * @param wordVector
     * @param nVocab
     * @param prior
     * @return
     */
    public double logDCM(double[] counts, double prior){
        int nWord = counts.length;
		double logGammaOfWPrior = gammaCache.logGamma(prior * nWord);
		double W_logGammaOfPrior = nWord * gammaCache.logGamma(prior);
        
        double ret = logGammaOfWPrior - W_logGammaOfPrior;
        if (Math.abs(W_logGammaOfPrior - W_logGammaOfPrior) > .0001){
            logger.error(String.format("believed: %.4e ; true %.4e", W_logGammaOfPrior, gammaCache.logGamma(prior)*nWord));
            logger.error(String.format("W = %d/%d",nWord,counts.length));
        }
        
        double N = 0;
        for (int i = 0; i < counts.length; i++) {
            N += counts[i] + prior;
            ret += gammaCache.logGamma(counts[i] + prior);
        }
        ret -= gammaCache.logGamma(N);
        
        return ret;
    }
    
    /**
     * Compute gradient of log probability of document (in fact segment) using DCM
     * This method takes dense form of word vector, because it would be sometime efficient.
     * @param wordVector
     * @param nVocab
     * @param prior
     * @return
     */
    public double gradientDCM(double[] counts, double N, double prior){
        int nWord = counts.length;
    	double S = prior * nWord;
        double ret = nWord * (gammaCache.digamma(S) - gammaCache.digamma(N + S) - gammaCache.digamma(prior));
        for (int i = 0; i < counts.length; i++) {
            ret += gammaCache.digamma(counts[i] + prior);
        }
        return ret;
    }
    
}
