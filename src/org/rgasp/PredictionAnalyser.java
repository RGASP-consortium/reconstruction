package org.rgasp;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jargs.gnu.*;

import org.rgasp.GTFRecord.RecordClass;
import org.rgasp.GTFRecord.Strand;

public class PredictionAnalyser {
	private static String reference_file;
	private static String prediction_file;
	private static String output_file;
	private static BufferedWriter gene_writer;
	private static BufferedWriter transcript_writer;
	private static BufferedWriter feature_writer;
	private static BufferedWriter feature_writer_unique;
	private static BufferedWriter statistics_writer;

	private static Map<String, Map<String, GTFTranscript>> transcripts_reference_new;
	private static Map<String, Map<String, GTFTranscript>> transcripts_prediction_new;
	private static Map<String, Map<String, GTFGene>> genes_reference;
	private static Map<String, Map<String, GTFGene>> genes_prediction;
	private static Map<String,List<GTFRecord>> exons_reference;
	private static Map<String,List<GTFRecord>> exons_prediction;
	public static Map<String, Map<Strand, Set<Integer>>> start_positions;
	public static Map<String, Map<Strand, Set<Integer>>> end_positions;
	public static Map<String, Map<Strand, Map<Integer, List<GTFRecord>>>> start_positions_map;
	public static Map<String, Map <Strand, Map<Integer, List<GTFRecord>>>> end_positions_map;
	
	private static Map<String, Statistics> stats_map;
	
	private static GFFParser input_parser;
	public enum Feature {EXON, CDS};
	public enum Level {FIXED, FLEXIBLE};

	private static Feature[] features = {Feature.CDS};
	private static Feature feature;
	private static Level[] levels = {Level.FIXED};
	private static Level level;

	public static void printUsage() {
		System.err.println("Usage: PredictionAnalyser");
		System.err.println("-l, --level\tone of cds, exon, both");
		System.err.println("-m, --mode\tone of fixed, flexible, both");
		System.err.println("-a, --annotation\tannotation_file");
		System.err.println("-p, --prediction\tprediction_file");
		System.err.println("-o, output\toutput_directory");
	}
	
	public static void main(String[] args) {

		CmdLineParser parser = new CmdLineParser();
		CmdLineParser.Option level_opt = parser.addStringOption('l', "level");
		CmdLineParser.Option mode = parser.addStringOption('m', "mode");
		CmdLineParser.Option annotation = parser.addStringOption('a', "annotation");
		CmdLineParser.Option prediction = parser.addStringOption('p', "prediction");
		CmdLineParser.Option output = parser.addStringOption('o', "output");
		
		try {
			parser.parse(args);
		}
		catch (CmdLineParser.OptionException e) {
			System.err.println(e.getMessage());
			printUsage();
			System.exit(2);
		}
		

		String option_level = ((String)parser.getOptionValue(level_opt));
		if (option_level != null && option_level.equals("cds")) features = new Feature[] {Feature.CDS};
		else if (option_level != null && option_level.equals("cds")) features = new Feature[] {Feature.EXON};
		else if (option_level != null) features = new Feature[] {Feature.CDS, Feature.EXON};
		reference_file = (String)parser.getOptionValue(annotation);
		prediction_file = (String)parser.getOptionValue(prediction);
		String option_mode = ((String)parser.getOptionValue(mode));
		if (option_mode != null && option_mode.equals("fixed")) levels = new Level[] {Level.FIXED};
		else if (option_mode != null && option_mode.equals("flexible")) levels = new Level[] {Level.FLEXIBLE};
		else if (option_mode != null) levels = new Level[] {Level.FIXED, Level.FLEXIBLE};
		output_file = (String)parser.getOptionValue(output);
		
		if (option_level == null || reference_file == null || prediction_file == null || option_mode == null || output_file == null) {
			System.err.println("Mising arguments ...");
			printUsage();
			System.exit(2);
		}
		else if (!(option_level.equals("cds") || option_level.equals("exon") || option_level.equals("both"))) {
			System.err.println("Invalid level argument ...");
			printUsage();
			System.exit(2);
		}
		else if (!(option_mode.equals("fixed") || option_mode.equals("flexible") || option_mode.equals("both"))) {
			System.err.println("Invalid mode argument ...");
			printUsage();
			System.exit(2);
		}
		
		for (int fi = 0; fi < features.length; fi++) {
			feature = features[fi];
			transcripts_reference_new = new HashMap<String, Map<String, GTFTranscript>>();
			exons_reference = new HashMap<String, List<GTFRecord>>();
			genes_reference = new HashMap<String, Map<String, GTFGene>>();
			stats_map = new HashMap<String, Statistics>();

			GTFTranscript.existing_transcripts = new HashMap<String, GTFTranscript>();
			try {
				input_parser = new GFFParser(reference_file, feature);
				input_parser.getTranscriptList(exons_reference, transcripts_reference_new, genes_reference);
			}
			catch(IOException e) {
				System.err.println("Problem while parsing reference file ...");
				e.printStackTrace();
			}
			
			for (int l = 0; l < levels.length; l++) {
				level = levels[l];
				try {
					stats_map = new HashMap<String, Statistics>();
					gene_writer = new BufferedWriter(new FileWriter(output_file + feature + "_" + level + "_gene_tmp.tbl"));
					transcript_writer = new BufferedWriter(new FileWriter(output_file + feature + "_" + level + "_transcript_tmp.tbl"));
					feature_writer = new BufferedWriter(new FileWriter(output_file + feature + "_" + level + "_exon_tmp.tbl"));
					feature_writer_unique = new BufferedWriter(new FileWriter(output_file + feature + "_" + level + "_exon_unique_tmp.tbl"));
					statistics_writer = new BufferedWriter(new FileWriter(output_file + feature + "_" + level + "_statistics.txt"));

					gene_writer.append("submission");
					transcript_writer.append("submission");
					feature_writer.append("submission");
					feature_writer_unique.append("submission");

					boolean first = true;
					StringBuffer secondLine = new StringBuffer();
					Map<String, Map<Integer, Map<Integer, Map<Strand, GTFRecord>>>> unique_exons = new HashMap<String, Map<Integer, Map<Integer, Map<Strand, GTFRecord>>>>();
					for (String sequence : transcripts_reference_new.keySet()) {
						secondLine.append("Sequence");
						if (!stats_map.containsKey(sequence)) stats_map.put(sequence, new Statistics());
						if (first) {
							statistics_writer.append("Submission\tSequence\t" + stats_map.get(sequence).header());
							statistics_writer.newLine();
							first = false;
						}
//						stats_map.get(sequence).m_Num_Exons_Ref += exons_reference.get(sequence).size();
						stats_map.get(sequence).m_Num_Transcr_Ref += transcripts_reference_new.get(sequence).size();
						stats_map.get(sequence).m_Num_Genes_Ref += genes_reference.get(sequence).size();
						for (String gene_id : genes_reference.get(sequence).keySet()) {
							gene_writer.append("\t" + gene_id);
						}
						gene_writer.newLine();
						for (String transcript_id : transcripts_reference_new.get(sequence).keySet()) {
							transcript_writer.append("\t" + transcript_id);
							List<GTFRecord> tmp = transcripts_reference_new.get(sequence).get(transcript_id).getExons();

//							if (transcripts_reference_new.get(sequence).get(transcript_id).strand == Strand.Strand_Positive) stats_map.get(sequence).m_Num_Transcr_Ref_Forw++;
//							else stats_map.get(sequence).m_Num_Transcr_Ref_Rev++;

							for (int i = 0; i < tmp.size(); i++) {
								GTFRecord rec = tmp.get(i);
								secondLine.append("\t" + sequence);
								feature_writer.append("\t" + transcript_id + "-" + (rec.strand == Strand.Strand_Positive ? i + 1 : tmp.size() - i));
								if (!unique_exons.containsKey(rec.sequence)) unique_exons.put(rec.sequence, new HashMap<Integer, Map<Integer, Map<Strand, GTFRecord>>>());
								if (!unique_exons.get(rec.sequence).containsKey(rec.start_pos)) unique_exons.get(rec.sequence).put(rec.start_pos, new HashMap<Integer, Map<Strand, GTFRecord>>());
								if (!unique_exons.get(rec.sequence).get(rec.start_pos).containsKey(rec.end_pos)) unique_exons.get(rec.sequence).get(rec.start_pos).put(rec.end_pos, new HashMap<Strand, GTFRecord>());
								if (!unique_exons.get(rec.sequence).get(rec.start_pos).get(rec.end_pos).containsKey(rec.strand)) {
									unique_exons.get(rec.sequence).get(rec.start_pos).get(rec.end_pos).put(rec.strand, rec);
									Collections.sort(rec.container);
									String tr_id = rec.container.get(0);
									List<GTFRecord> exon_list = transcripts_reference_new.get(sequence).get(tr_id).exons;
									int index = exon_list.indexOf(rec);
									int size = exon_list.size();
									feature_writer_unique.append("\t" + tr_id + "-" + (rec.strand == Strand.Strand_Positive ? index + 1 : size - index));

//									if (rec.strand == Strand.Strand_Positive) stats_map.get(sequence).m_Num_Exons_Ref_Forw++;
//									else stats_map.get(sequence).m_Num_Exons_Ref_Rev++;
									boolean single=true;
									boolean first_exon=true;
									boolean last=true;
									for (String t_id : rec.container) {
										List<GTFRecord> exons = transcripts_reference_new.get(rec.sequence).get(t_id).exons;
										if (exons.size() > 1) {
											single=false;
										}
										if ((rec.strand == Strand.Strand_Positive && exons.indexOf(rec) > 0) ||
												(rec.strand == Strand.Strand_Negative && exons.indexOf(rec) < exons.size() - 1)) {
											first_exon=false;
										}
										if ((rec.strand == Strand.Strand_Positive && exons.indexOf(rec) < exons.size() - 1) ||
												(rec.strand == Strand.Strand_Negative && exons.indexOf(rec) > 0)) {
											last=false;
										}
									}
									stats_map.get(rec.sequence).m_Num_Exons_Ref++;
									if (single) {
										rec.record_class = RecordClass.SINGLE;
										stats_map.get(rec.sequence).m_Num_Exons_Ref_SINGLE++;
									}
									else if (first_exon) {
										stats_map.get(rec.sequence).m_Num_Exons_Ref_FIRST++;
										rec.record_class = RecordClass.FIRST;
									}
									else if (last) {
										rec.record_class = RecordClass.TERMINAL;
										stats_map.get(rec.sequence).m_Num_Exons_Ref_TERMINAL++;
									}
									else {
										rec.record_class = RecordClass.INTERNAL;
										stats_map.get(sequence).m_Num_Exons_Ref_INTERNAL++;
									}
								}
							}
						}
					}
					transcript_writer.newLine();
					feature_writer.newLine();
					feature_writer_unique.newLine();
					statistics_writer.newLine();
				} catch (IOException e) {
					System.out.println("Problem writing output files ...");
					e.printStackTrace();
				}
				File p_file = new File(prediction_file);
				if (p_file.isDirectory()) {
					File[] files = p_file.listFiles();
					System.out.println("Directory contains " + files.length + " many files");
					for (File f : files) {
						prediction_file = f.getAbsolutePath();
						if (f.length() == 0) continue;
						GTFTranscript.existing_transcripts = new HashMap<String, GTFTranscript>();
						start_positions = new HashMap<String, Map<Strand, Set<Integer>>>();
						start_positions_map = new HashMap<String, Map<Strand, Map<Integer, List<GTFRecord>>>>();
						end_positions = new HashMap<String, Map<Strand, Set<Integer>>>();
						end_positions_map = new HashMap<String, Map<Strand, Map<Integer, List<GTFRecord>>>>();
						GTFTranscript.existing_transcripts = new HashMap<String, GTFTranscript>();
						System.gc();
						analysePrediction();
						writeOutput();
						for (String sequence : exons_reference.keySet()) {
							stats_map.get(sequence).reset();
							for (GTFRecord rec : exons_reference.get(sequence)) {
								rec.overlap = 0;
								rec.matchingExons.clear();
								rec.overlapingExons.clear();
							}
						}
					}
				}
				else {
					if (p_file.length() == 0) continue;
					prediction_file = p_file.getAbsolutePath();
					GTFTranscript.existing_transcripts = new HashMap<String, GTFTranscript>();
					start_positions = new HashMap<String, Map<Strand, Set<Integer>>>();
					start_positions_map = new HashMap<String, Map<Strand, Map<Integer, List<GTFRecord>>>>();
					end_positions = new HashMap<String, Map<Strand, Set<Integer>>>();
					end_positions_map = new HashMap<String, Map<Strand, Map<Integer, List<GTFRecord>>>>();
					GTFTranscript.existing_transcripts = new HashMap<String, GTFTranscript>();
					System.gc();
					analysePrediction();
					writeOutput();
					for (String sequence : exons_reference.keySet()) {
						stats_map.get(sequence).reset();
						for (GTFRecord rec : exons_reference.get(sequence)) {
							rec.overlap = 0;
							rec.matchingExons.clear();
							rec.overlapingExons.clear();
						}
					}
				}
				try {
					gene_writer.close();
					transcript_writer.close();
					feature_writer.close();
					feature_writer_unique.close();
					statistics_writer.close();
					
					BufferedReader tmp_reader = new BufferedReader(new InputStreamReader(PredictionAnalyser.class.getClassLoader().getResourceAsStream("transpose.sh")));
					File tmp = File.createTempFile("transpose", ".sh");
					tmp.deleteOnExit();
					tmp.setExecutable(true);
					
					BufferedWriter tmp_writer = new BufferedWriter(new FileWriter(tmp.getAbsolutePath()));
					while (tmp_reader.ready()) {
						tmp_writer.append(tmp_reader.readLine());
						tmp_writer.append("\n");
					}
					tmp_writer.close();
					tmp_reader.close();
					
					Process p = Runtime.getRuntime().exec(new String[]{tmp.getAbsolutePath(), (output_file + feature + "_" + level + "_gene_tmp.tbl"), (output_file + feature + "_" + level + "_gene.tbl")," && rm run.sh./transpose.sh"});
					p.waitFor();
					Runtime.getRuntime().exec("rm " + output_file + feature + "_" + level + "_gene_tmp.tbl");
					p = Runtime.getRuntime().exec(new String[]{tmp.getAbsolutePath(), (output_file + feature + "_" + level + "_exon_tmp.tbl"), (output_file + feature + "_" + level + "_exon.tbl")});
					p.waitFor();
					Runtime.getRuntime().exec("rm " + output_file + feature + "_" + level + "_exon_tmp.tbl");
					p = Runtime.getRuntime().exec(new String[]{tmp.getAbsolutePath(), (output_file + feature + "_" + level + "_exon_unique_tmp.tbl"), (output_file + feature + "_" + level + "_exon_unique.tbl")});
					p.waitFor();
					Runtime.getRuntime().exec("rm " + output_file + feature + "_" + level + "_exon_unique_tmp.tbl");
					p = Runtime.getRuntime().exec(new String[]{tmp.getAbsolutePath(), (output_file + feature + "_" + level + "_transcript_tmp.tbl"), (output_file + feature + "_" + level + "_transcript.tbl")});
					p.waitFor();
					Runtime.getRuntime().exec("rm " + output_file + feature + "_" + level + "_transcript_tmp.tbl");
				}
				catch(IOException e) {
					e.printStackTrace();
				} catch (InterruptedException e) {
					e.printStackTrace();
				} 
			}
		}
	}
	
	public static void writeOutput() {
		Map<String, Map<Integer, Map<Integer, Map<Strand, GTFRecord>>>> unique_exons = new HashMap<String, Map<Integer, Map<Integer, Map<Strand, GTFRecord>>>>();
		Map<String, Map<Integer, Map<Integer, Map<Strand, GTFRecord>>>> unique_exons_pred = new HashMap<String, Map<Integer, Map<Integer, Map<Strand, GTFRecord>>>>();
		Map<String, Map<Integer, Map<Integer, Map<Strand, GTFRecord>>>> unique_correct_exons = new HashMap<String, Map<Integer, Map<Integer, Map<Strand, GTFRecord>>>>();
		List<String> unique_transcripts_pred = new ArrayList<String>();
		List<String> unique_genes_pred = new ArrayList<String>();
		try {
			String submission_name = prediction_file.substring(prediction_file.lastIndexOf("/") + 1);
			submission_name = submission_name.replace(".rgasp_submission", "");
			gene_writer.append(submission_name);
			transcript_writer.append(submission_name);
			feature_writer.append(submission_name);
			feature_writer_unique.append(submission_name);
			List<GTFRecord> transcript;
			for (String sequence : transcripts_reference_new.keySet()) {
				if (genes_prediction.containsKey(sequence)) stats_map.get(sequence).m_Num_Genes_Pred += genes_prediction.get(sequence).size();
				for (String gene_id : genes_reference.get(sequence).keySet()) {
					if (genes_reference.get(sequence).get(gene_id).isPredicted(level).size() > 0) {
						stats_map.get(sequence).m_TP_Genes++;
						for (String g_id : genes_reference.get(sequence).get(gene_id).isPredicted(level)) {
							if (!unique_genes_pred.contains(g_id)) {
								unique_genes_pred.add(g_id);
								stats_map.get(sequence).m_TP_Genes_Pred++;
							}
						}
					}
					gene_writer.append("\t" + ((genes_reference.get(sequence).get(gene_id).isPredicted(level).size() > 0) ? "1" : "0"));
				}
				if (transcripts_prediction_new.containsKey(sequence)) {
					stats_map.get(sequence).m_Num_Transcr_Pred += transcripts_prediction_new.get(sequence).keySet().size();
					for (GTFTranscript t : transcripts_prediction_new.get(sequence).values()) {
//						if (t.strand == Strand.Strand_Positive) stats_map.get(sequence).m_Num_Transcr_Pred_Forw++;
//						else stats_map.get(sequence).m_Num_Transcr_Pred_Rev++;
						for (GTFRecord rec : t.getExons()) {
							if (!unique_exons_pred.containsKey(rec.sequence)) unique_exons_pred.put(rec.sequence, new HashMap<Integer, Map<Integer, Map<Strand, GTFRecord>>>());
							if (!unique_exons_pred.get(rec.sequence).containsKey(rec.start_pos)) unique_exons_pred.get(rec.sequence).put(rec.start_pos, new HashMap<Integer, Map<Strand, GTFRecord>>());
							if (!unique_exons_pred.get(rec.sequence).get(rec.start_pos).containsKey(rec.end_pos)) unique_exons_pred.get(rec.sequence).get(rec.start_pos).put(rec.end_pos, new HashMap<Strand, GTFRecord>());
							if (!unique_exons_pred.get(rec.sequence).get(rec.start_pos).get(rec.end_pos).containsKey(rec.strand)) {
								unique_exons_pred.get(rec.sequence).get(rec.start_pos).get(rec.end_pos).put(rec.strand, rec);
								stats_map.get(sequence).m_Num_Exons_Pred++;
							}
						}
					}
				}
				for (String transcript_id : transcripts_reference_new.get(sequence).keySet()) {
					transcript =  transcripts_reference_new.get(sequence).get(transcript_id).getExons();
					if (transcripts_reference_new.get(sequence).get(transcript_id).isPredicted(level).size() > 0) {
						stats_map.get(sequence).m_TP_Transcr++;
						for (String tr_id : transcripts_reference_new.get(sequence).get(transcript_id).isPredicted(level)) {
							if (!unique_transcripts_pred.contains(tr_id)) {
								unique_transcripts_pred.add(tr_id);
								stats_map.get(sequence).m_TP_Transcr_Pred++;
							}
						}
					}
					for (GTFRecord rec : transcript) {
						if (!unique_exons.containsKey(rec.sequence)) unique_exons.put(rec.sequence, new HashMap<Integer, Map<Integer, Map<Strand, GTFRecord>>>());
						if (!unique_exons.get(rec.sequence).containsKey(rec.start_pos)) unique_exons.get(rec.sequence).put(rec.start_pos, new HashMap<Integer, Map<Strand, GTFRecord>>());
						if (!unique_exons.get(rec.sequence).get(rec.start_pos).containsKey(rec.end_pos)) unique_exons.get(rec.sequence).get(rec.start_pos).put(rec.end_pos, new HashMap<Strand, GTFRecord>());
						if (!unique_exons.get(rec.sequence).get(rec.start_pos).get(rec.end_pos).containsKey(rec.strand)) {
							unique_exons.get(rec.sequence).get(rec.start_pos).get(rec.end_pos).put(rec.strand, rec);
							
							if (rec.matchingExons.size() > 0 || (level == Level.FLEXIBLE && transcript.size() == 1 && rec.overlap >= 0.9)  || ((level == Level.FLEXIBLE && transcript.indexOf(rec) == 0 && end_positions.containsKey(sequence) &&
									end_positions.get(sequence).containsKey(rec.strand) && end_positions.get(sequence).get(rec.strand).contains(rec.end_pos)) ||
											(level == Level.FLEXIBLE && transcript.indexOf(rec) == transcript.size() - 1 && start_positions.containsKey(sequence) &&
											start_positions.get(sequence).containsKey(rec.strand) && start_positions.get(sequence).get(rec.strand).contains(rec.start_pos)))) {
								if (rec.matchingExons.size() > 0) {
									for (GTFRecord gr : rec.matchingExons) {
										if (!unique_correct_exons.containsKey(gr.sequence)) unique_correct_exons.put(gr.sequence, new HashMap<Integer, Map<Integer, Map<Strand, GTFRecord>>>());
										if (!unique_correct_exons.get(gr.sequence).containsKey(gr.start_pos)) unique_correct_exons.get(gr.sequence).put(gr.start_pos, new HashMap<Integer, Map<Strand, GTFRecord>>());
										if (!unique_correct_exons.get(gr.sequence).get(gr.start_pos).containsKey(gr.end_pos)) unique_correct_exons.get(gr.sequence).get(gr.start_pos).put(gr.end_pos, new HashMap<Strand, GTFRecord>());
										if (!unique_correct_exons.get(gr.sequence).get(gr.start_pos).get(gr.end_pos).containsKey(gr.strand)) {
											unique_correct_exons.get(gr.sequence).get(gr.start_pos).get(gr.end_pos).put(gr.strand, gr);
											stats_map.get(sequence).m_TP_Exons_Pred++;
											break;
										}
									}
								}
								if (transcript.size() == 1 && rec.overlap >= 0.6) {
									for (GTFRecord gr : rec.overlapingExons) {
										if (!unique_correct_exons.containsKey(gr.sequence)) unique_correct_exons.put(gr.sequence, new HashMap<Integer, Map<Integer, Map<Strand, GTFRecord>>>());
										if (!unique_correct_exons.get(gr.sequence).containsKey(gr.start_pos)) unique_correct_exons.get(gr.sequence).put(gr.start_pos, new HashMap<Integer, Map<Strand, GTFRecord>>());
										if (!unique_correct_exons.get(gr.sequence).get(gr.start_pos).containsKey(gr.end_pos)) unique_correct_exons.get(gr.sequence).get(gr.start_pos).put(gr.end_pos, new HashMap<Strand, GTFRecord>());
										if (!unique_correct_exons.get(gr.sequence).get(gr.start_pos).get(gr.end_pos).containsKey(gr.strand)) {
											unique_correct_exons.get(gr.sequence).get(gr.start_pos).get(gr.end_pos).put(gr.strand, gr);
											stats_map.get(sequence).m_TP_Exons_Pred++;
											break;
										}
									}
								}
								else if (transcript.indexOf(rec) == 0 && end_positions.containsKey(sequence) &&
										end_positions.get(sequence).containsKey(rec.strand) && end_positions.get(sequence).get(rec.strand).contains(rec.end_pos)) {
									for (GTFRecord gr : end_positions_map.get(sequence).get(rec.strand).get(rec.end_pos)) {
										if (!unique_correct_exons.containsKey(gr.sequence)) unique_correct_exons.put(gr.sequence, new HashMap<Integer, Map<Integer, Map<Strand, GTFRecord>>>());
										if (!unique_correct_exons.get(gr.sequence).containsKey(gr.start_pos)) unique_correct_exons.get(gr.sequence).put(gr.start_pos, new HashMap<Integer, Map<Strand, GTFRecord>>());
										if (!unique_correct_exons.get(gr.sequence).get(gr.start_pos).containsKey(gr.end_pos)) unique_correct_exons.get(gr.sequence).get(gr.start_pos).put(gr.end_pos, new HashMap<Strand, GTFRecord>());
										if (!unique_correct_exons.get(gr.sequence).get(gr.start_pos).get(gr.end_pos).containsKey(gr.strand)) {
											unique_correct_exons.get(gr.sequence).get(gr.start_pos).get(gr.end_pos).put(gr.strand, gr);
											stats_map.get(sequence).m_TP_Exons_Pred++;
											break;
										}
									}
								}
								else if (transcript.indexOf(rec) == transcript.size() - 1 && start_positions.containsKey(sequence) &&
										start_positions.get(sequence).containsKey(rec.strand) && start_positions.get(sequence).get(rec.strand).contains(rec.start_pos)) {
									for (GTFRecord gr : start_positions_map.get(sequence).get(rec.strand).get(rec.start_pos)) {
										if (!unique_correct_exons.containsKey(gr.sequence)) unique_correct_exons.put(gr.sequence, new HashMap<Integer, Map<Integer, Map<Strand, GTFRecord>>>());
										if (!unique_correct_exons.get(gr.sequence).containsKey(gr.start_pos)) unique_correct_exons.get(gr.sequence).put(gr.start_pos, new HashMap<Integer, Map<Strand, GTFRecord>>());
										if (!unique_correct_exons.get(gr.sequence).get(gr.start_pos).containsKey(gr.end_pos)) unique_correct_exons.get(gr.sequence).get(gr.start_pos).put(gr.end_pos, new HashMap<Strand, GTFRecord>());
										if (!unique_correct_exons.get(gr.sequence).get(gr.start_pos).get(gr.end_pos).containsKey(gr.strand)) {
											unique_correct_exons.get(gr.sequence).get(gr.start_pos).get(gr.end_pos).put(gr.strand, gr);
											stats_map.get(sequence).m_TP_Exons_Pred++;
											break;
										}
									}
								}
								feature_writer_unique.append("\t1");

								stats_map.get(rec.sequence).m_TP_Exons++;
								if (rec.record_class == RecordClass.SINGLE) stats_map.get(rec.sequence).m_TP_Exons_SINGLE++;
								else if (rec.record_class == RecordClass.FIRST) stats_map.get(rec.sequence).m_TP_Exons_FIRST++;
								else if (rec.record_class == RecordClass.TERMINAL) stats_map.get(rec.sequence).m_TP_Exons_TERMINAL++;
								else stats_map.get(sequence).m_TP_Exons_INTERNAL++;
							}
							else {
								feature_writer_unique.append("\t0");
							}
						}
					}
					transcript_writer.append("\t" + (transcripts_reference_new.get(sequence).get(transcript_id).isPredicted(level).size() > 0 ? "1" : "0"));
					for (int i = 0; i < transcript.size(); i++) {
						if (level == Level.FIXED || (i > 0 && i < (transcript.size() - 1))) {
							feature_writer.append("\t" + (transcript.get(i).matchingExons.size() > 0 ? "1" : "0"));
////						chr_feature_writer.append("\t" + (transcript.get(i).overlap == 1 ? "1" : "0"));
						}
						else if (level == Level.FLEXIBLE && transcript.size() == 1) {
							feature_writer.append("\t" + (transcript.get(0).overlap >= 0.9 ? "1" : "0"));
						}
						else if (level == Level.FLEXIBLE && i == 0) {
							feature_writer.append("\t" + (end_positions.containsKey(sequence) &&
								end_positions.get(sequence).containsKey(transcript.get(i).strand) && end_positions.get(sequence).get(transcript.get(i).strand).contains(transcript.get(i).end_pos) ? "1" : "0")); 
						}
						else if (level == Level.FLEXIBLE && i == (transcript.size() - 1)) {
							feature_writer.append("\t" + (start_positions.containsKey(sequence) &&
									start_positions.get(sequence).containsKey(transcript.get(i).strand) && start_positions.get(sequence).get(transcript.get(i).strand).contains(transcript.get(i).start_pos) ? "1" : "0"));
						}
					}
				}
				statistics_writer.append(submission_name);
				statistics_writer.append("\t" + sequence + "\t" + stats_map.get(sequence).data());
				statistics_writer.newLine();
			}
			gene_writer.newLine();
			transcript_writer.newLine();
			feature_writer.newLine();
			feature_writer_unique.newLine();
		}
		catch(IOException e) {
			System.out.println("Problem writing output files ...");
			e.printStackTrace();
		}
	}
	
	private static void analysePrediction() {
		try {
			input_parser = new GFFParser(prediction_file, feature);
			transcripts_prediction_new = new HashMap<String, Map<String, GTFTranscript>>();
			exons_prediction = new HashMap<String, List<GTFRecord>>();
			genes_prediction = new HashMap<String, Map<String, GTFGene>>();
			input_parser.getTranscriptList(exons_prediction, transcripts_prediction_new, genes_prediction);
			System.gc();
			analyseExonOverlap();
		} catch (IOException e) {
			System.out.println("Problem parsing prediction file " + prediction_file);
			e.printStackTrace();
		}
	}
	
	private static void analyseExonOverlap() {
		for (String sequence : exons_prediction.keySet()) {
			if (exons_prediction.get(sequence).size() == 0 || !exons_reference.containsKey(sequence)) continue;
			int ref_index = 0;
			int pred_index = 0;
			GTFRecord reference_rec = (ref_index < exons_reference.get(sequence).size()) ? exons_reference.get(sequence).get(ref_index) : null,
					prediction_rec = (pred_index < exons_prediction.get(sequence).size()) ? exons_prediction.get(sequence).get(pred_index) : null;
			if (reference_rec == null || prediction_rec == null) {
				continue;
			}
			while(true) {

				int start_pos = prediction_rec.start_pos;
				int end_pos = prediction_rec.end_pos;
				Strand strand = prediction_rec.strand;
				
				if (!start_positions.containsKey(sequence)) start_positions.put(sequence, new HashMap<Strand, Set<Integer>>());
				if (!start_positions.get(sequence).containsKey(strand)) start_positions.get(sequence).put(strand, new HashSet<Integer>());
				if (!start_positions.get(sequence).get(strand).contains(start_pos)) start_positions.get(sequence).get(strand).add(start_pos);
				
				if (!end_positions.containsKey(sequence)) end_positions.put(sequence, new HashMap<Strand, Set<Integer>>());
				if (!end_positions.get(sequence).containsKey(strand)) end_positions.get(sequence).put(strand, new HashSet<Integer>());
				if (!end_positions.get(sequence).get(strand).contains(end_pos)) end_positions.get(sequence).get(strand).add(end_pos);
				
				if (!start_positions_map.containsKey(sequence)) start_positions_map.put(sequence, new HashMap<Strand, Map<Integer, List<GTFRecord>>>());
				if (!start_positions_map.get(sequence).containsKey(strand)) start_positions_map.get(sequence).put(strand, new HashMap<Integer, List<GTFRecord>>());
				if (!start_positions_map.get(sequence).get(strand).containsKey(start_pos)) start_positions_map.get(sequence).get(strand).put(start_pos, new ArrayList<GTFRecord>());
				if (!start_positions_map.get(sequence).get(strand).get(start_pos).contains(prediction_rec)) start_positions_map.get(sequence).get(strand).get(start_pos).add(prediction_rec);
				
				if (!end_positions_map.containsKey(sequence)) end_positions_map.put(sequence, new HashMap<Strand, Map<Integer, List<GTFRecord>>>());
				if (!end_positions_map.get(sequence).containsKey(strand)) end_positions_map.get(sequence).put(strand, new HashMap<Integer, List<GTFRecord>>());
				if (!end_positions_map.get(sequence).get(strand).containsKey(end_pos)) end_positions_map.get(sequence).get(strand).put(end_pos, new ArrayList<GTFRecord>());
				if (!end_positions_map.get(sequence).get(strand).get(end_pos).contains(prediction_rec)) end_positions_map.get(sequence).get(strand).get(end_pos).add(prediction_rec);

				float overlap = 0.f;
				
				if (reference_rec.compareTo(prediction_rec) == 0) {
					reference_rec.matchingExons.add(prediction_rec);
					overlap = 1.f;
					reference_rec.overlap = overlap;
				}
				if (reference_rec.record_class == RecordClass.SINGLE && level == Level.FLEXIBLE && overlap != 1) {
//					System.out.println("Single exon transcript");
					int pred_index2 = pred_index;
					while (pred_index2 < exons_prediction.get(sequence).size() && prediction_rec.sequence.equals(reference_rec.sequence) &&
							prediction_rec.strand == reference_rec.strand && prediction_rec.start_pos <= reference_rec.end_pos) {
						overlap = 0;
						if (reference_rec.compareTo(prediction_rec) < 0) {
							if (reference_rec.getEnd_pos() >= prediction_rec.start_pos &&
									reference_rec.strand.equals(prediction_rec.strand)) {
								if (reference_rec.getEnd_pos() > prediction_rec.end_pos) {
									overlap = ((float)(prediction_rec.end_pos - prediction_rec.start_pos + 1)) / ((float)reference_rec.getLength());
								}
								else if (reference_rec.getEnd_pos() < prediction_rec.end_pos) {
									if (reference_rec.start_pos == prediction_rec.start_pos) overlap = ((float)(prediction_rec.end_pos - prediction_rec.start_pos + 1)) / ((float)reference_rec.getLength());
									else overlap = ((float)(reference_rec.getEnd_pos() - prediction_rec.start_pos + 1)) / ((float)reference_rec.getLength());
								}
								else if (reference_rec.getEnd_pos() == prediction_rec.end_pos) {
									overlap = ((float)(prediction_rec.end_pos - prediction_rec.start_pos + 1)) / ((float)reference_rec.getLength());
								}
							}
						}
						else if (reference_rec.compareTo(prediction_rec) > 0) {
							if (prediction_rec.end_pos >= reference_rec.getStart_pos() &&
									reference_rec.strand.equals(prediction_rec.strand)) {
								if (prediction_rec.end_pos > reference_rec.getEnd_pos()) {
									overlap = ((float)(prediction_rec.end_pos - prediction_rec.start_pos + 1)) / ((float)reference_rec.getLength());
								}
								else if (prediction_rec.end_pos < reference_rec.getEnd_pos()) {
									overlap = ((float)(prediction_rec.end_pos - reference_rec.getStart_pos() + 1)) / ((float)reference_rec.getLength());
								}
								else if (prediction_rec.end_pos == reference_rec.getEnd_pos()) {
									overlap = ((float)(prediction_rec.end_pos - prediction_rec.start_pos + 1)) / ((float)reference_rec.getLength());
								}
							}
						}
//						System.out.println("Reference inside " + reference_rec.start_pos + " " + reference_rec.end_pos);
//						System.out.println("Prediction inside " + prediction_rec.start_pos + " " + prediction_rec.end_pos);
//						System.out.println("Overlap " + overlap);
						if (overlap >= 0.9) {
							if (reference_rec.overlap < overlap) reference_rec.overlap = overlap;
//							System.out.println("Detecting overlapping exon");
							reference_rec.overlapingExons.add(prediction_rec);
						}
						pred_index2++;
						if (pred_index2 < exons_prediction.get(sequence).size()) prediction_rec = exons_prediction.get(sequence).get(pred_index2);
						else break;
					}
					prediction_rec = exons_prediction.get(sequence).get(pred_index);
				}

				if (overlap == 1) {
					if (ref_index + 1 < exons_reference.get(sequence).size()) {
						ref_index++;
						reference_rec = exons_reference.get(sequence).get(ref_index);
					}
					else break;
					if (pred_index + 1 < exons_prediction.get(sequence).size()) {
						pred_index++;
						prediction_rec = exons_prediction.get(sequence).get(pred_index);
					}
					else break;
				}
				else if (reference_rec.compareTo(prediction_rec) < 0) {
					if (ref_index + 1 < exons_reference.get(sequence).size()) {
						ref_index++;
						reference_rec = exons_reference.get(sequence).get(ref_index);
					}
					else break;
				}
				else if (reference_rec.compareTo(prediction_rec) > 0) {
					if (pred_index + 1 < exons_prediction.get(sequence).size()) {
						pred_index++;
						prediction_rec = exons_prediction.get(sequence).get(pred_index);
					}
					else break;
				}
				else {
					System.out.println("Something went wrong here!");
				}
			}
		}
	}
}
