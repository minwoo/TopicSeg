/*
 * Copyright (C) 2010 Cluster of Excellence, Univ. of Saarland
 * Minwoo Jeong (minwoo.j@gmail.com) is a main developer.
 * This file is part of "TopicSeg" package.
 * This software is provided under the terms of LGPL.
 */

package topicseg.utils;

import java.io.*;
import java.util.HashMap;
import java.util.regex.*;

import org.apache.log4j.DailyRollingFileAppender;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;

/**
 * unisaar.topicseg.util::Configure
 *
 * @author minwoo
 */
public class Option implements Serializable {
	
	private static final long serialVersionUID = 1L;
	private HashMap<String, Object> map = new HashMap<String, Object>();

    public Option() {
    	map.clear();
    }

    public Option(String fileName) throws IOException {
    	map.clear();
        parse(fileName);
    }

    public void parse(String fileName) throws IOException {
        if (fileName == null)
            throw new IllegalArgumentException("Can't read \"null\" config file");

        FileReader fr = new FileReader(fileName);
        BufferedReader br = new BufferedReader(fr);
        Pattern p = Pattern.compile("[\t =]+");
        String line = null;
        while ( (line = br.readLine()) != null) {
            String[] tokens = p.split(line);
            if (tokens.length < 2)
            	continue;
                //throw new IllegalArgumentException("Invalid config at " + line);
            if (!map.containsKey(tokens[0]))
                map.put(tokens[0], tokens[1]);
        }
        
        br.close(); fr.close();
    }

    public Object get(String key) {
        if (!map.containsKey(key))
            throw new IllegalArgumentException("Invalid key");
        return map.get(key);
    }
    
    public int getInteger(String key) {
        if (!map.containsKey(key))
            throw new IllegalArgumentException("Invalid key");
        int r = Integer.parseInt((String) map.get(key));
        return r;
    }

    public double getDouble(String key) {
        if (!map.containsKey(key))
            throw new IllegalArgumentException("Invalid key");
        double r = Double.parseDouble((String) map.get(key));
        return r;
    }

    public String getString(String key) {
        if (!map.containsKey(key))
            throw new IllegalArgumentException("Invalid key");
        return (String) map.get(key);
    }

    public boolean getBoolean(String key) {
        if (!map.containsKey(key))
            throw new IllegalArgumentException("Invalid key");
        String r = (String) map.get(key);
        if (r.toLowerCase().equals("true"))
            return true;
        else
            return false;
    }
    
    public boolean contains(String key) {
    	if (map.containsKey(key))
    		return true;
    	return false;    		
    }
    
    public HashMap<String, Object> getHashMap() {
    	return map;
    }
    
    public HashMap<String, Object> getListStartWith(String prefix) {
    	HashMap<String, Object> subset = new HashMap<String, Object>();
    	for (String key : map.keySet()) {
    		if (key.startsWith(prefix)) {
    			subset.put(key, map.get(key));
    		}
    	}
    	return subset;
    }
    
    public static void addFileLogger(Logger logger, String filename) throws IOException {
    	String layout = "[%d{HH:mm:ss.SSS}] %-5p %C{1}.%M (%F:%L) - %m%n";
    	String logfilename = filename;
    	PatternLayout patternlayout = new PatternLayout(layout);
    	FileAppender appender = new FileAppender(patternlayout, logfilename, true);
    	logger.addAppender(appender);
    	logger.setLevel(Level.INFO); 
    }
}
