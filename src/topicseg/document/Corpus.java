/*
 * Copyright (C) 2010 Cluster of Excellence, Univ. of Saarland
 * Minwoo Jeong (minwoo.j@gmail.com) is a main developer.
 * This file is part of "TopicSeg" package.
 * This software is provided under the terms of LGPL.
 */

package topicseg.document;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

import org.apache.log4j.Logger;

import opennlp.tools.lang.english.Tokenizer;
import topicseg.document.Document.SegmentBoundary;
import topicseg.model.TopicModel;
import topicseg.utils.Option;

import edu.mit.nlp.util.Utils;

/**
 * unisaar.topicseg.document::Corpus.java
 *
 * @author minwoo
 */
public class Corpus extends ArrayList<DocumentSet> {
	
	Alphabet sharedDict;
    Logger logger = Logger.getLogger(Corpus.class);
	private final String labelfileSuffix = "label";
	private boolean useSharedVocab = false;
    private boolean useStemmer;
    private boolean removeStopWord;
    private boolean useTokenizer;
    private static ArrayList<String> stopWords = new ArrayList<String>();
    private opennlp.tools.lang.english.Tokenizer tokenizer;
    
    private double dirPriorOfGlobalLM;
    private double globalTypePrior = 1;
    private double localTypePrior = 1;
    private double lengthPrior = 0;
    private HashMap<Integer, Double> globalTypePriors = new HashMap<Integer, Double>();
    private HashMap<Integer, Double> localTypePriors = new HashMap<Integer, Double>();
    private HashMap<Integer, Double> lengthPriors = new HashMap<Integer, Double>();
    private HashMap<String, Integer> docsetSpecificMaxCluster = new HashMap<String, Integer>();
    private int maxCluster;
    
    public Corpus(Option config) throws IOException {
    	super();
    	
    	this.useSharedVocab = config.contains("model.sharedDict") ? config.getBoolean("model.sharedDict") : false;
    	this.useStemmer = config.contains("corpus.useStemmer") ? config.getBoolean("corpus.useStemmer") : false;
    	this.removeStopWord = config.contains("corpus.removeStopword") ? config.getBoolean("corpus.removeStopword") : false;
    	this.useTokenizer = config.contains("corpus.useTokenizer") ? config.getBoolean("corpus.useTokenizer") : false;
    	String stopWordFilename = config.contains("corpus.stopwordList") ? config.getString("corpus.stopwordList") : "";
    	String tokenizerModelFilename = config.contains("corpus.tokenizerModel") ? config.getString("corpus.tokenizerModel") : "";
    	
		this.dirPriorOfGlobalLM = config.contains("prior.dcm.global") ? config.getDouble("prior.dcm.global") : 1; // just in case
    	this.globalTypePrior = config.contains("prior.docset.globalType") ? config.getDouble("prior.docset.globalType") : 1;
    	this.localTypePrior = config.contains("prior.docset.localType") ? config.getDouble("prior.docset.localType") : 1;
    	this.lengthPrior = config.contains("prior.docset.length") ? config.getDouble("prior.docset.length") : 1;
		this.maxCluster = config.contains("sampler.maxGlobalCluster") ? config.getInteger("sampler.maxGlobalCluster") : 0;

    	HashMap<String, Object> docSpecificPriors = config.getListStartWith("prior.doc.");
    	for (String key : docSpecificPriors.keySet()) {
    		String key2 = key.replace("prior.doc.", "").replace('.', ' ');
    		String[] tok = key2.split(" ");
    		int docNo = Integer.parseInt(tok[0]);
            double val = Double.parseDouble((String) docSpecificPriors.get(key));
    		if (tok[1].equals("globalType"))
    			globalTypePriors.put(docNo, val);
    		else if (tok[1].equals("localType"))
    			localTypePriors.put(docNo, val);
    		else if (tok[1].equals("length"))
    			lengthPriors.put(docNo, val);
    	}
    	
    	HashMap<String, Object> docsetSpecificParams = config.getListStartWith("sampler.docset.maxGlobalCluster");
    	for (String key : docsetSpecificParams.keySet()) {
    		String docsetName = key.replace("sampler.docset.maxGlobalCluster.", "").replace('.', ' ');
    		docsetSpecificMaxCluster.put(docsetName, Integer.parseInt((String)docsetSpecificParams.get(key)));
    	}
    	
    	
		if (useSharedVocab)
			sharedDict = new Alphabet();
		if (useTokenizer)
			tokenizer = new Tokenizer(tokenizerModelFilename);
		if (removeStopWord)
			this.loadStopWords(stopWordFilename);
		
        if (config.contains("exp.log")) 
        	Option.addFileLogger(logger, config.getString("exp.log"));
    }
    
	public Corpus(String dirpath, Option config) throws IOException {
		this(config);
		try {
			this.loadTextCorpus(dirpath);
		}
		catch (IOException e) {
        	logger.error("error " + e.getMessage());
			e.printStackTrace();
            System.exit(2);
		}
	}
	
	/**
	 * Loading the text in the directory
	 * The files that have same filename excluding extension will be stored in the same document set.
	 * TODO Implement built-in word tokenizer, stemmer, stop word remover, and tf-idf word selector    
	 * @param dirpath
	 * @param option
	 * @throws IOException 
	 */
	public void loadTextCorpus(String dirpath) throws IOException {
        File dirs = new File(dirpath);
        FilenameFilter firstFileFilter = new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(".0"); // every doc set has one doc at least 
            }
        };
        String[] filenames = dirs.list(firstFileFilter);
        
        for (int i = 0; i < filenames.length; i++) {
        	final String filename = filenames[i].replace(".0", "");
            
            // for each document set (filename.*)
            FilenameFilter myFilter = new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return name.startsWith(filename+"."); 
                }
            };
            String[] docsetFiles = dirs.list(myFilter);
            
            DocumentSet docset = useSharedVocab ? new DocumentSet(this.sharedDict) : new DocumentSet();
            docset.setInfo(i, filename);
            
            for (String docfile : docsetFiles) {
            	// annotation
            	if (docfile.endsWith(labelfileSuffix)) {
            		Annotation annotation = docset.getAnnotation();
            		
        			FileReader inFile = new FileReader(dirs+"/"+docfile);
        			BufferedReader inStream = new BufferedReader(inFile);
        			String line = null;
        			while ((line = inStream.readLine()) != null) {
        				StringTokenizer tokens = new StringTokenizer(line, " ");
        				String topicLabel = tokens.nextToken();
        				
        				while (tokens.hasMoreTokens()) {
        					String token = tokens.nextToken();
        					String[] segs = token.split("::", -1);
        					int docId = Integer.parseInt(segs[0]);
        					String[] segPoints = segs[1].split("-", -1);
        					int start = Integer.parseInt(segPoints[0]);
        					int end = Integer.parseInt(segPoints[1]);
        					annotation.add(topicLabel, docId, new int[] {start, end});
        				}
        			}
            	}
            	else { // document
            		// MJ: it would be Document() in which every separate documents have their own independent dictionary.
            		Document doc = useSharedVocab ? new Document(sharedDict) : new Document(docset.getAlphabet());
            		
        			FileReader inFile = new FileReader(dirs+"/"+docfile);
        			BufferedReader inStream = new BufferedReader(inFile);
        			String line = null;
        			while ((line = inStream.readLine()) != null) {
        				line = line.replaceAll("&[^;]*;","");
        				if (useTokenizer) {
        					String[] tokens = tokenizer.tokenize(line);
        					if (tokens.length > 0)
        						line = tokens[0];
        					for (int j = 1; j < tokens.length; j++)
        						line += " " + tokens[j];
        				}
        				
        				line = line.toLowerCase();
        				
        				if (useStemmer) {
        					String newLine = "";
            				StringTokenizer tokens = new StringTokenizer(line, " ");
            				while (tokens.hasMoreTokens()) {
            					String token = tokens.nextToken();
            					String stemmedWord = Utils.stemWord(token);
            					if (removeStopWord) {
            						if (!stopWords.contains(stemmedWord))
            							newLine += stemmedWord + " ";
            					}
            					else {
            						newLine += stemmedWord + " ";
            					}
            				}
            				doc.add(newLine);
        				}
        				else {
        					String newLine = "";
            				StringTokenizer tokens = new StringTokenizer(line, " ");
            				while (tokens.hasMoreTokens()) {
            					String token = tokens.nextToken();
            					if (removeStopWord) {
            						if (!stopWords.contains(token))
            							newLine += token + " ";
            					}
            					else {
            						newLine += token + " ";
            					}
            				}
            				doc.add(newLine);
        					//doc.add(line);
        				}
        				//System.out.println(line); 
        			}
            		String[] toks = docfile.replace('.', ' ').split(" ");
            		int docId = Integer.parseInt(toks[toks.length-1]);
        			doc.setId(docId);
        			double alpha = globalTypePrior, beta = localTypePrior, gamma = lengthPrior;
        			if (globalTypePriors.containsKey(docId))
        				alpha = globalTypePriors.get(docId);
        			if (localTypePriors.containsKey(docId))
        				beta = localTypePriors.get(docId);
        			if (lengthPriors.containsKey(docId))
        				gamma = lengthPriors.get(docId);
        			doc.initDocPrior(alpha, beta, gamma);
        			//doc.tfidf();
        			docset.add(doc);

        			inStream.close(); inFile.close();
            	}
            }
            
            //data.add(docset);
//            if (docset.size() < 2 || docset.getAlphabet().size() < 100) 
//            	continue;
			//docset.tfidf();
            this.add(docset);
            if (docsetSpecificMaxCluster.containsKey(docset.getFilename())) 
            	docset.maxGlobalCluster = docsetSpecificMaxCluster.get(docset.getFilename());
            else
            	docset.maxGlobalCluster = maxCluster;
            
            Annotation annotation = docset.getAnnotation();
            int[] T = new int[docset.size()];
            for (Document doc : docset) T[doc.getId()] = doc.size();
            //for (int d = 0; d < docset.size(); d++) T[d] = docset.get(d).size();
            annotation.set(docset.size(), T);
            
//            // temporary
//            int[][] refs = annotation.toArray();
//            for (Document doc : docset) {
//	    		ArrayList<SegmentBoundary> refBoundaries = doc.getSegmentBoundary(refs[doc.getId()]);
//	    		for (SegmentBoundary refBoundary : refBoundaries) {
//	    			if (refBoundary.id < 0) {
//	    				annotation.add(refBoundary.id + (-100 * doc.getId()), doc.getId(), new int[] {refBoundary.start, refBoundary.end-1});
//	    			}
//	    		}
//	    		
//            }
            //System.out.println(annotation.toString());
            
            logger.info("loading docset \"" + dirs+"/"+filename + ".*\", #doc=" + docset.size() + " #Vocab=" + docset.getAlphabet().size() + " annotation=" + (docset.getAnnotation().size() > 0 ? "T" : "F"));
            //System.out.println("[Set "+docset.getSetId() + "] " + docset.getFilename() + "\nAnnotation\n" + docset.getAnnotation().toString());
            //System.out.println("# of vocab = " + (useSharedVocab?sharedDict.size():docset.getAlphabet().size()) );
        }
	}

	/**
	 * Loading stop word list
	 * @param stopWordFilename
	 * @param useStemmer
	 */
    private void loadStopWords(String stopWordFilename) {
        try {
            stopWords.clear();
            BufferedReader br = new BufferedReader(new FileReader(stopWordFilename));
            String line = br.readLine();
            while ( (line = br.readLine()) != null) {
                line = line.trim().toLowerCase();
                if (useStemmer)
                	stopWords.add(Utils.stemWord(line));
                else
                	stopWords.add(line);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            throw new IllegalStateException();
        } catch (IOException e) {
            e.printStackTrace();
            throw new IllegalStateException();
        }
    }
	
    
    public void writeSegmentResult(String inputDirpath, String outputDirpath) throws NumberFormatException, IOException {
    	
    	for (DocumentSet docSet : this) {
    		String filename = docSet.getFilename();
    		//Annotation reference = docSet.getAnnotation();
    		//int[][] refs = reference.toArray();
    		
    		// output the top words
			BufferedWriter topicWordWriter = new BufferedWriter(new FileWriter(outputDirpath+"/"+filename+".topicWord"));
			TopicModel topicModel = docSet.getGlobalTopic();
			
			topicWordWriter.write(topicModel.getTopicWord(40, 0.1, docSet.getAlphabet()));
			
			topicWordWriter.close();

    		for (Document doc : docSet) {
    			BufferedReader inStream = new BufferedReader(new FileReader(inputDirpath+"/"+filename+"."+doc.getId()));
    			ArrayList<String> rawTexts = new ArrayList<String>();
    			String line = null;
    			while ((line = inStream.readLine()) != null) 
    				if (line.trim().length() > 0) rawTexts.add(line);
    			inStream.close();
    			
    			//int[] labels = refs[doc.getId()];
    			int[] labels = doc.getTopicLabel();
    			//assert(labels.length == rawTexts.size());
    			
    			BufferedWriter writer = new BufferedWriter(new FileWriter(outputDirpath+"/"+filename+"."+doc.getId()));
    			
    			int index = 0;
    	        int pre_label = labels[0];
    	        writer.write("====="); if (pre_label >= 0) writer.write(" " + labels[0]); writer.write("\n");
    	        if (rawTexts.get(index).startsWith("=====")) {
        	        writer.write(rawTexts.get(index).replace("=====", "-----") + "\n");
        	        index++;
    	        }
    	        writer.write(rawTexts.get(index) + "\n"); index++;
    	        int i = 1; boolean check = true;
    	        for (; index < rawTexts.size(); index++)
    	        {
    	    		if (i < labels.length && pre_label != labels[i] && check) { 
    	    	        writer.write("====="); if (labels[i] >= 0) writer.write(" " + labels[i]); writer.write("\n");
    	    	        check = false;
    	    		}
        	        if (rawTexts.get(index).startsWith("=====")) {
            	        writer.write(rawTexts.get(index).replace("=====", "-----") + "\n");
            	        //continue;
        	        }
        	        else {
        	        	writer.write(rawTexts.get(index) + "\n");
        	    		pre_label = labels[i];
        	    		i++;
        	    		check = true;
        	        }
    	        }
    	        writer.write("=====\n");
    			writer.close();
    		}
    	}
    }
}
