/*
 * Copyright (C) 2010 Cluster of Excellence, Univ. of Saarland
 * Minwoo Jeong (minwoo.j@gmail.com) is a main developer.
 * This file is part of "TopicSeg" package.
 * This software is provided under the terms of LGPL.
 */

package topicseg.utils;

/**
 * unisaar.topicseg.stats::Annealer.java
 * This code was from 'BayesSeg' package.
 */
public class Annealer {
	
    double burnin_duration;
    double cooling_duration;
    double max_burnin;
    int it_num;
    int num_its;
	
    public Annealer(double burnin_duration, double cooling_duration, double max_burnin, int num_its) {
        this.burnin_duration = burnin_duration;
        this.cooling_duration = cooling_duration;
        this.max_burnin = max_burnin;
        this.num_its = num_its;
        it_num = 0;
    }

    public double anneal(double prob) {
        double out = annealWithoutUpdate(prob);
        it_num++;
        return out;
    }
    
    public double annealWithoutUpdate(double prob) {
        double temperature = 1;
		int burnin_end = (int) (burnin_duration * num_its);
		int cooling_start = (int) ((1 - cooling_duration) * num_its );
		
		if (it_num < burnin_end) 
		    temperature = 1 + ((double) max_burnin - 1) * (1 - (double) it_num / (double) burnin_end);
		
		if (it_num > cooling_start)
		    temperature = ((double) num_its - it_num + 1) / ((double) num_its - cooling_start);
		
		return Math.pow(prob, 1 / temperature);
    }
    
    public void reset() { it_num = 0; }
    public void update(){ it_num++; }
    public double getHalfProbAnnealed() { return annealWithoutUpdate(.5); }
}
