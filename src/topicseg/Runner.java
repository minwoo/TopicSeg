/*
 * Copyright (C) 2010 Cluster of Excellence, Univ. of Saarland
 * Minwoo Jeong (minwoo.j@gmail.com) is a main developer.
 * This file is part of "TopicSeg" package.
 * This software is provided under the terms of LGPL.
 */

package topicseg;

import jargs.gnu.CmdLineParser;

import org.apache.log4j.*;

import topicseg.document.*;
import topicseg.segment.*;
import topicseg.utils.*;

/**
 * unisaar.topicseg::Runner.java
 * A command line running tool of topic segmentation 
 * 
 * @author minwoo
 */
public class Runner {
	
    public static void main(String[] args) {
    	
        Logger logger = Logger.getLogger(Runner.class);
        logger.info("===" + Runner.class.getName() + "===");

        // command line parsing
        CmdLineParser cmdParser = new CmdLineParser();
        CmdLineParser.Option debug = cmdParser.addBooleanOption('d', "debug");
        CmdLineParser.Option verbose = cmdParser.addBooleanOption('v', "verbose");
        CmdLineParser.Option configfile = cmdParser.addStringOption('c', "config");

        try {
           cmdParser.parse(args);
        }
        catch (CmdLineParser.OptionException e) {
            logger.error(e.getMessage());
            logger.error("Usage: java -cp ${CLASSPATH} unisaar.topicseg.Runner " +
                    "[-c,--config] config_file [{-v,--verbose}] [{-d,--debug}]");
            System.exit(2);
        }

        String configFileName = (String)cmdParser.getOptionValue(configfile);
        Boolean isDebug = (Boolean)cmdParser.getOptionValue(debug, Boolean.TRUE);
        Boolean isVerbose = (Boolean)cmdParser.getOptionValue(verbose, Boolean.TRUE);
        
        // running
        try {
            Option config = new Option(configFileName);
            
            if (config.contains("exp.log")) 
            	Option.addFileLogger(logger, config.getString("exp.log"));
            
            String classifierName = config.contains("model.class") ? config.getString("model.class") : "unisaar.topicseg.segment.MultiSeg";
            // hyperparameters for DCM
            double dcmGlobalPrior = config.contains("prior.dcm.global") ? config.getDouble("prior.dcm.global") : 0.1;
            double dcmLocalPrior = config.contains("prior.dcm.local") ? config.getDouble("prior.dcm.local") : 0.1;

            int nTrials = config.contains("exp.trial") ? config.getInteger("exp.trial") : 1;
            String outputFilename = config.contains("exp.output") ? config.getString("exp.output") : "";
            boolean usePerDocEval = config.contains("exp.perDocEval") ? config.getBoolean("exp.perDocEval") : false;
            
            for (int n = 0; n < nTrials; n++) {
	            // segmenter
	            Segmenter segmenter = (Segmenter) Class.forName(classifierName).getConstructor(new Class[]{}).newInstance(new Object[]{});
	            segmenter.initialize(config);
	            segmenter.initializeRandom(n+1);
	            
	        	// corpus processing
	            if (config.contains("corpus.datasetDir")) {
	            	// data loading
	            	Corpus corpus = new Corpus(config.getString("corpus.datasetDir"), config);
	            	
	            	// segmentation
	            	segmenter.run(corpus, dcmGlobalPrior, dcmLocalPrior);
	            	
	            	// write out
	            	if (!outputFilename.equals("")) {
	            		segmenter.writeResult(corpus, outputFilename + "." + n, usePerDocEval);
	            	}
	                if (config.contains("output.outputDir") && config.contains("output.rawDataDir")) {
	                	corpus.writeSegmentResult(config.getString("output.rawDataDir"), config.getString("output.outputDir"));
	                }
	            }
            }
            
        }
        catch (Exception e) {
        	logger.error("error " + e.getMessage());
        	e.printStackTrace();
            System.exit(2);
        }

    }
}
