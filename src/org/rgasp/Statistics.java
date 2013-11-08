package org.rgasp;

import java.lang.reflect.Field;

public class Statistics {
	
	// Reference don't include in reset
	public int m_Num_Genes_Ref;
	public int m_Num_Transcr_Ref;
//	public int m_Num_Transcr_Ref_Forw;
//	public int m_Num_Transcr_Ref_Rev;
	public int m_Num_Exons_Ref;
//	public int m_Num_Exons_Ref_Forw;
//	public int m_Num_Exons_Ref_Rev;
	public int m_Num_Exons_Ref_FIRST;
	public int m_Num_Exons_Ref_INTERNAL;
	public int m_Num_Exons_Ref_TERMINAL;
	public int m_Num_Exons_Ref_SINGLE;
	
	// reference features that have a match in the prediction
	// in contrast to predicted features that match a reference feature
	public int m_TP_Genes;
	public int m_TP_Transcr;
	//public int m_Num_Transcr_Ref_Forw;
	//public int m_Num_Transcr_Ref_Rev;
	public int m_TP_Exons;
	//public int m_Num_Exons_Forw;
	//public int m_Num_Exons_Rev;
	public int m_TP_Exons_FIRST;
	public int m_TP_Exons_INTERNAL;
	public int m_TP_Exons_TERMINAL;
	public int m_TP_Exons_SINGLE;
	
	//gene prediction
	public int m_Num_Genes_Pred;
	public int m_TP_Genes_Pred;
	
	// transcript wise
	// number of transcripts predicted
	public int m_Num_Transcr_Pred;
	//public int m_Num_Transcr_Pred_Forw;
	//public int m_Num_Transcr_Pred_Rev;
	
	// number of validated transcript
	public int m_TP_Transcr_Pred;
	//public int m_TP_Transcr_Pred_Forw;
	//public int m_TP_Transcr_Pred_Rev;
	
	// exon wise
	// number of exons predicted
	public int m_Num_Exons_Pred;
	//public int m_Num_Exons_Pred_Forw;
	//public int m_Num_Exons_Pred_Rev;
	//public int m_Num_Exons_Pred_FIRST;
	//public int m_Num_Exons_Pred_INTERNAL;
	//public int m_Num_Exons_Pred_TERMINAL;
	//public int m_Num_Exons_Pred_SINGLE;
	
	// validated exons
	public int m_TP_Exons_Pred;
	//public int m_TP_Exons_Pred_Forw;
	//public int m_TP_Exons_Pred_Rev;
	//public int m_TP_Exons_Pred_FIRST;
	//public int m_TP_Exons_Pred_INTERNAL;
	//public int m_TP_Exons_Pred_TERMINAL;
	//public int m_TP_Exons_Pred_SINGLE;
	
	public void reset() {
		Class<Statistics> c = Statistics.class;
		for (Field f : c.getFields()) {
			if (!f.getName().contains("Ref"))
				try {
					f.setInt(this, 0);
				} catch (IllegalArgumentException e) {
					e.printStackTrace();
				} catch (IllegalAccessException e) {
					e.printStackTrace();
				}
		}
	}
	
	public String header() {
		StringBuffer text = new StringBuffer();
		Class<Statistics> c = Statistics.class;
		for (Field f : c.getFields()) {
			text.append("\t" + f.getName());
		}
		return text.toString();
	}
	
	public String data() {
		StringBuffer text = new StringBuffer();
		Class<Statistics> c = Statistics.class;
		for (Field f : c.getFields()) {
			try {
				text.append("\t" + f.getInt(this));
			} catch (IllegalArgumentException e) {
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			}
		}
		return text.toString();
	}
}
