package org.rgasp;

import java.util.ArrayList;
import java.util.List;

import org.rgasp.PredictionAnalyser.Level;

public class GTFGene extends GTFRecord {
	
	protected List<GTFTranscript> transcripts;
	protected String gene_id;
	
	public GTFGene(GTFTranscript transcript) {
		super(transcript.sequence, RecordType.GENE.toString(), transcript.start_pos, transcript.end_pos, transcript.strand.toString(), null);
		transcripts = new ArrayList<GTFTranscript>();
		transcripts.add(transcript);
	}
	
	public void addTranscript(GTFTranscript transcript) {
		transcripts.add(transcript);
		if (transcript.end_pos > this.end_pos) {
			this.end_pos = transcript.end_pos;
		}
		if (transcript.start_pos < this.start_pos) {
			this.start_pos = transcript.start_pos;
		}
	}
	
	//
	public List<GTFTranscript> getTranscripts() {
		return this.transcripts;
	}
	
	public List<String> isPredicted(Level level) {
		List<String> matching_pred_genes = new ArrayList<String>();
		for (GTFTranscript t : transcripts) {
			if (t.isPredicted(level).size() > 0) {
				for(String tr_id : t.isPredicted(level)) {
					for (String g_id : GTFTranscript.existing_transcripts.get(tr_id).container) {
						if (!matching_pred_genes.contains(g_id)) matching_pred_genes.add(g_id);
					}
				}
			}
		}
		return matching_pred_genes;
	}
	
	public int numTranscriptsPred(Level level) {
		int numTranscriptsPred = 0;
		for (GTFTranscript t : transcripts) {
			if (t.isPredicted(level).size() > 0) numTranscriptsPred++;
		}
		return numTranscriptsPred;
	}
	
}
