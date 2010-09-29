/*
 * Copyright (C) 2010 Cluster of Excellence, Univ. of Saarland
 * Minwoo Jeong (minwoo.j@gmail.com) is a main developer.
 * This file is part of "TopicSeg" package.
 * This software is provided under the terms of LGPL.
 */

package topicseg.model;

import java.util.HashMap;

import cern.jet.stat.Gamma;

/**
 * unisaar.topicseg.model::GammaCache.java
 *
 * @author minwoo
 */
public class GammaCache {
	
	private HashMap<Double, Double> gammaCache;
	private HashMap<Double, Double> digammaCache;
    private long nHits;
    private long nMisses;
	
	public GammaCache() {
		nHits = 0; nMisses = 0;
		gammaCache = new HashMap<Double, Double>();
		digammaCache = new HashMap<Double, Double>();
	}
	
    public GammaCache(int capacity, float loadFactor) { 
		nHits = 0; nMisses = 0;
		gammaCache = new HashMap<Double, Double>(capacity, loadFactor);
		digammaCache = new HashMap<Double, Double>(capacity, loadFactor);
    }
    
    /**
     * Initialize the hash table
     */
    public void clear() {
		nHits = 0; nMisses = 0;
		gammaCache = new HashMap<Double, Double>();
		digammaCache = new HashMap<Double, Double>();
    }
    
    /**
     * Compute log gamma
     * Return cached value if stored
     * Gamma function is in Colt package
     * @param x
     * @return
     */
    public synchronized double logGamma(final double x) {
    	if (x == 0)
    		return 0;
    	
        Double ret = gammaCache.get(x);
        
        if (ret == null){
            double newVal = Gamma.logGamma(x);
            gammaCache.put(x, newVal);
            nMisses++;
            return newVal;
        } 
        else {
            nHits++;
            return ret;
        }
    }
    
    /**
     * Compute digamma
     * Return cached value if stored
     * Digamma function is in Alias' lingpipe package
     * @param x
     * @return
     */
    public synchronized double digamma(final double x) {
        Double ret = digammaCache.get(x);
        
        if (ret == null) {
            double newVal = com.aliasi.util.Math.digamma(x);
            digammaCache.put(x, newVal);
            return newVal;
        } 
        else {
            return ret;
        }
    }
}
