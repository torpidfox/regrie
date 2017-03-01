package utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import objects.BasePairs;
import objects.DNAsequence;
import objects.PFM;
import objects.TargetSitesGroupLogic;
import utils.Constants;
import utils.Utils;

/**
 * this class contains static fields and methods that are used for the biology part of the application
 * @author created: n.r.zabet@gen.cam.ac.uk, modified: dmitrav@inbox.ru
 *
 */
public class CellUtils {

	public static final BasePairs bps =new BasePairs();
	public static final TargetSitesGroupLogic tgsl =new TargetSitesGroupLogic();
	
	public static final byte bpANYID = bps.getANYID();
	/**
	 * generate random sequence of nucleotides ( can be used for DBD of TFs or DNA sequence the probability of nucleotides are uniformly distributed
	 * @param generator the random number generator
	 * @param length the length in nucleotides
	 * @return a string of bytes represnting the nucleotide sequence
	 */
	public static byte[] generateRandomDNASequence(Random generator,int length, double proportionOfA, double proportionOfT, double proportionOfC, double proportionOfG){
		byte[] seq = new byte[length];
		double buffer;
		
		double quotaOfA = proportionOfA;
		double quotaOfT = quotaOfA+proportionOfT;
		double quotaOfC = quotaOfT+proportionOfC;
		double quotaOfG = quotaOfC+proportionOfG;
		
		byte current;
		
		for(int i=0; i<length;i++){
			buffer = Utils.generateNextDouble(generator, 0, quotaOfG);
			
			if(buffer<quotaOfA){
				current = bps.bpsID.get("A");
			} else if(buffer<quotaOfT){
				current = bps.bpsID.get("T");
			} else if(buffer<quotaOfC){
				current = bps.bpsID.get("C");
			} else {
				current = bps.bpsID.get("G");
			}
			
			seq[i] =current;
		}
		
		return seq;
	} 
	
	
	
	

	
	
	/**
	 * comptes the affinities between a TF and DNA
	 * @param strand the DNA strand
	 * @param TFseq the recognise DNA sequence
	 * @param sizeLeft The size on the left of the DBD that the TF occupies on the DNA.
	 * @param sizeTotal The total number of bp  that the TF occupies on the DNA.
	 * @param es the specific energy
	 * @return the affinity vector
	 */
	public static double[] computeTFAffinities(Random generator, byte[] strand, byte[] TFseq, int sizeLeft, int sizeTotal, double es, int direction, double roughness){
		double[] affinities = new double[strand.length];
		
		if(TFseq!=null && TFseq.length>0){
			if(direction == 0){
				for(int i=0; i<strand.length-sizeTotal;i++){
					affinities[i] = computeTFAffinityLR(strand, i, TFseq, sizeLeft, es);
				}
			}
			if(direction == 1){
				for(int i=0; i<strand.length-sizeTotal;i++){
					affinities[i] = computeTFAffinityRL(strand, i, TFseq, sizeLeft, es);
				}
			}
		} else{
			for(int i=0; i<strand.length-sizeTotal;i++){
				affinities[i] = Utils.generateNextNormalDistributedDouble(generator, es, roughness, 0);
			}
		}
		
		for(int i=strand.length-sizeTotal; i< strand.length;i++){
			affinities[i] = Constants.NONE;
		}
		
		return affinities;
	} 
	
	
	/**
	 * computes the affinity between a TF and the DNA at a specific position using the Gerland 2002 two ways from 5' to 3' and from 3' to 5'
	 * @param DNAseq the DNA sequence
	 * @param DNApos the index of the starting bp of the position on the DNA where the affinity is computed
	 * @param TFseq the recognise DNA sequence
	 * @param ens the non-specific energy
	 * @param es the specific energy
	 * @return the value of the affinity
	 */
	public static double computeTFAffinity2Way(byte[] DNAseq, int DNApos, byte[] TFseq, int sizeLeft, double es){
		double sumLR = 0;
		double sumRL = 0;

		byte[] revSeq = getReversedComplementSequences(TFseq);
		
		for(int i=0;i<TFseq.length;i++){
			if(TFseq[i]!= bpANYID && DNAseq[DNApos+i+sizeLeft] != TFseq[i]){
				sumLR += es;
			}
			
			if(revSeq[i]!= bpANYID && DNAseq[DNApos+i+sizeLeft] != revSeq[i]){
				sumRL += es;
			} 
		}
		 return sumLR>sumRL?sumRL:sumLR;
	}
	
	
	
	/**
	 * computes the affinity between a TF and the DNA at a specific position using the Gerland 2002 one way from 5' to 3'
	 * @param DNAseq the DNA sequence
	 * @param DNApos the index of the starting bp of the position on the DNA where the affinity is computed
	 * @param TFseq the recognise DNA sequence
	 * @param ens the non-specific energy
	 * @param es the specific energy
	 * @return the value of the affinity
	 */
	public static double computeTFAffinityLR(byte[] DNAseq, int DNApos, byte[] TFseq, int sizeLeft, double es){
		double sumLR = 0;

	
		for(int i=0;i<TFseq.length;i++){
			if(TFseq[i]!= bpANYID && DNAseq[DNApos+i+sizeLeft] != TFseq[i]){
				sumLR += es;
			}
			

		}
		 return sumLR;
	}
	
	
	/**
	 * computes the affinity between a TF and the DNA at a specific position using the Gerland 2002 two ways from 3' to 5'
	 * @param DNAseq the DNA sequence
	 * @param DNApos the index of the starting bp of the position on the DNA where the affinity is computed
	 * @param TFseq the recognise DNA sequence
	 * @param ens the non-specific energy
	 * @param es the specific energy
	 * @return the value of the affinity
	 */
	public static double computeTFAffinityRL(byte[] DNAseq, int DNApos, byte[] TFseq, int sizeLeft, double es){
		double sumRL = 0;

		byte[] revSeq = getReversedComplementSequences(TFseq);
		
		for(int i=0;i<TFseq.length;i++){
			
			if(revSeq[i]!= bpANYID && DNAseq[DNApos+i+sizeLeft] != revSeq[i]){
				sumRL += es;
			} 
		}
		 return sumRL;
	}
	
	
	/**
	 * computes the affinity between a TF and the DNA at a specific position using a known sequence affinity 3' -> 5'
	 * @param DNAseq the DNA sequence
	 * @param DNApos the index of the starting bp of the position on the DNA where the affinity is computed
	 * @param seqsAffinities the list of sequences and their affinties
	 * @param sizeLeft the left size
	 * @param sizeMotif the size of the motif
	 * @param defaultAffinity the default affinity if not found in the list
	 * @return the value of the affinity
	 */
	public static double computeTFAffinityLR(byte[] DNAseq, int DNApos, HashMap<DNAsequence,Double> seqsAffinities, int sizeLeft, int sizeMotif, double defaultAffinity){
		double sumLR = defaultAffinity;
	
		byte[] seq = new byte[sizeMotif-sizeLeft];
		
		for(int i=0;i<sizeMotif;i++){
			seq[i]=DNAseq[i+DNApos+sizeLeft];
		}

		DNAsequence buffer = new DNAsequence(seq);

		
		if(seqsAffinities.containsKey(buffer)){
			sumLR=seqsAffinities.get(buffer);
		}
		
		 return sumLR;
	}
	

	/**
	 * computes the affinity between a TF and the DNA at a specific position using a known sequence affinity 5' -> 3'
	 * @param DNAseq the DNA sequence
	 * @param DNApos the index of the starting bp of the position on the DNA where the affinity is computed
	 * @param seqsAffinities the list of sequences and their affinties
	 * @param sizeLeft the left size
	 * @param sizeMotif the size of the motif
	 * @param defaultAffinity the default affinity if not found in the list
	 * @return the value of the affinity
	 */
	public static double computeTFAffinityRL(byte[] DNAseq, int DNApos, HashMap<DNAsequence,Double> seqsAffinities, int sizeLeft, int sizeMotif, double defaultAffinity){
		double sumLR = defaultAffinity;
	
		byte[] seq = new byte[sizeMotif-sizeLeft];
		
		for(int i=0;i<sizeMotif;i++){
			seq[i]=DNAseq[DNApos+sizeLeft+ (sizeMotif-1) - i];
		}

		DNAsequence buffer = new DNAsequence(seq);

		if(seqsAffinities.containsKey(buffer)){
			sumLR=seqsAffinities.get(buffer);
		}
		
		 return sumLR;
	}
	
	
	/**
	 * comptes the affinities between a TF and DNA
	 * @param strand the DNA strand
	 * @param TFseq the recognise DNA sequence
	 * @param sizeLeft The size on the left of the DBD that the TF occupies on the DNA.
	 * @param sizeTotal The total number of bp  that the TF occupies on the DNA.
	 * @param es the specific energy
	 * @return the affinity vector
	 */
	public static double[] computeTFAffinities(Random generator, byte[] strand, PFM pfm, int sizeLeft, int sizeTotal, double es, int direction, double roughness, double[] bpFreq){
		double[] affinities = new double[strand.length];
		
		if(pfm!=null && pfm.isCorrect && pfm.motifSize>0){
			if(direction == 0){
				for(int i=0; i<strand.length-sizeTotal;i++){
					affinities[i] = computeTFAffinityLR(strand, i, pfm, sizeLeft, es, bpFreq);
				}
			}
			if(direction == 1){
				for(int i=0; i<strand.length-sizeTotal;i++){
					affinities[i] = computeTFAffinityRL(strand, i, pfm, sizeLeft, es, bpFreq);
				}
			}
		} else{
			for(int i=0; i<strand.length-sizeTotal;i++){
				affinities[i] = Utils.generateNextNormalDistributedDouble(generator, es, roughness, 0);
				
			}
		}
		
		for(int i=strand.length-sizeTotal; i< strand.length;i++){
			affinities[i] = Constants.NONE;
		}
		
		return affinities;
	} 
	

	/**
	 * comptes the affinities between a TF and DNA
	 * @param strand the DNA strand
	 * @param TFseq the recognise DNA sequence
	 * @param sizeLeft The size on the left of the DBD that the TF occupies on the DNA.
	 * @param sizeTotal The total number of bp  that the TF occupies on the DNA.
	 * @param es the specific energy
	 * @return the affinity vector
	 */
	public static double[] computeTFAffinities(byte[] strand, HashMap<DNAsequence,Double> seqsAffinities, double defaultAffinity, int sizeLeft, int sizeMotif, int sizeTotal, int direction){
		double[] affinities = new double[strand.length];
		
		if(seqsAffinities!=null && seqsAffinities.size()>0){
			if(direction == 0){
				for(int i=0; i<strand.length-sizeTotal;i++){
					//before
					//affinities[i] = computeTFAffinityLR(strand, i, seqsAffinities, sizeLeft, sizeMotif, defaultAffinity);
					//after
					//computing affinities according to the thermodynamic approach
					affinities[i] = Math.exp(computeTFAffinityLR(strand, i, seqsAffinities, sizeLeft, sizeMotif, defaultAffinity));

				}
			}
			if(direction == 1){
				for(int i=0; i<strand.length-sizeTotal;i++){
					//before
					//affinities[i] = computeTFAffinityRL(strand,  i, seqsAffinities, sizeLeft, sizeMotif, defaultAffinity);
					//after
					//computing affinities according to the thermodynamic approach
					affinities[i] = computeTFAffinityRL(strand,  i, seqsAffinities, sizeLeft, sizeMotif, defaultAffinity);
				}
			}
		} else{
			for(int i=0; i<strand.length-sizeTotal;i++){
				affinities[i] =defaultAffinity;
				
			}
		}
		
		for(int i=strand.length-sizeTotal; i< strand.length;i++){
			affinities[i] = Constants.NONE;
		}
		
		return affinities;
	} 
	
	/**
	 * Edited by AD
	 * computes the affinity between a TF and the DNA at a specific position using the Gerland 2002 one way from 5' to 3'
	 * @param DNAseq the DNA sequence
	 * @param DNApos the index of the starting bp of the position on the DNA where the affinity is computed
	 * @param TFseq the recognise DNA sequence
	 * @param ens the non-specific energy
	 * @param es the specific energy
	 * @return the value of the affinity - value of max affinity
	 */
	public static double computeTFAffinityLR(byte[] DNAseq, int DNApos, PFM pfm, int sizeLeft, double es, double[] bpFreq){
		double sumLR = 0;
		double sumMax = 0;
		for(int i=0;i<pfm.motifSize;i++){
			//before
			//sumLR += es*pfm.getScorePFM(DNAseq[DNApos+i+sizeLeft], i, bpFreq[DNAseq[DNApos+i+sizeLeft]]);
			//after
			sumLR += pfm.getScorePFM(DNAseq[DNApos+i+sizeLeft], i, bpFreq[DNAseq[DNApos+i+sizeLeft]]);
			sumMax += pfm.getMaxScorePFM(i);
		}
		 return sumLR - sumMax;

	}

	/**
	 * AD
	 * for repression
	 * computes the affinity between a TF and the DNA at a specific position
	 * @param DNAseq the DNA sequence
	 * @param DNApos the index of the starting bp of the position on the DNA where the affinity is computed
	 * @param TFseq the recognise DNA sequence
	 * @param ens the non-specific energy
	 * @param es the specific energy
	 * @return the value of the affinity - value of max affinity
	 */
	public static double computeTFAffinity(byte[] DNAseq, int DNApos, PFM pfm, int sizeLeft){
		double sumLR = 0;
		double sumMax = 0;
		for(int i=0;i<pfm.motifSize;i++){
			//before
			//sumLR += es*pfm.getScorePFM(DNAseq[DNApos+i+sizeLeft], i, bpFreq[DNAseq[DNApos+i+sizeLeft]]);
			//after

			//this is a kind of duct tape, unfortunately
			if (DNApos+i == DNAseq.length)
			{
				break;
			} else {
				sumLR += pfm.getScorePFM(DNAseq[DNApos+i+sizeLeft], i);
				sumMax += pfm.getMaxScorePFM(i);
			}

		}
		return Math.exp(sumLR - sumMax);

	}

	
	/**
	 * Edited by AD
	 * computes the affinity between a TF and the DNA at a specific position using the Gerland 2002 two ways from 3' to 5'
	 * @param DNAseq the DNA sequence
	 * @param DNApos the index of the starting bp of the position on the DNA where the affinity is computed
	 * @param TFseq the recognise DNA sequence
	 * @param ens the non-specific energy
	 * @param es the specific energy
	 * @return the value of the affinity - value of max affinity
	 */
	public static double computeTFAffinityRL(byte[] DNAseq, int DNApos, PFM pfm, int sizeLeft, double es, double[] bpFreq){
		double sumRL = 0;	
		double sumMax = 0;
		byte[] revComplement = getReversedComplementSequences(DNAseq,DNApos+sizeLeft, pfm.motifSize);
		for(int i=0;i<pfm.motifSize;i++){
			//sumRL += es*pfm.getComplementScorePFM(DNAseq[DNApos+i+sizeLeft] , i, bpFreq[DNAseq[DNApos+i+sizeLeft]]);
			//before
			//sumRL += es*pfm.getScorePFM(revComplement[i] , i, bpFreq[revComplement[i]]);
			//after
			sumRL += pfm.getScorePFM(revComplement[i], i, bpFreq[revComplement[i]]);
			sumMax += pfm.getMaxScorePFM(i);
		}
		 return sumRL - sumMax;
	}
	
	
	/**
	 * computes the rate at which a bound TF will take a decision whether to move or not from current position.
	 * @param nonSpecificWaitingTime
	 * @param affinity
	 * @return
	 */
	public static double computeAvgMoveRate(double specificWaitingTime, double bindingEnergy){
		return (double)1.0/(specificWaitingTime*Math.exp(-bindingEnergy));
	} 
	
	/**
	 * converts waiting time into binding energy
	 * @param specificWaitingTime
	 * @param TFaffinity
	 * @return
	 */
	public static double computeBindingEnergy(double specificWaitingTime, double TFwaitingTime){
		return -Math.log(TFwaitingTime/specificWaitingTime);
	} 
	
	/**
	 * counts the poly(A) 
	 * @param strand the DNA strand where we count poly A content
	 * @return an array with the polyA content
	 */
	public static int[] getPolyA(byte[] strand){
		//initialise vector
		int[] polyA = new int[strand.length];
		for(int i=0;i<polyA.length;i++){
			polyA[i] = 0;
		}
		
		//count the appeareaces of consecutive As
		byte aValue = bps.bpsID.get("A");
		int count = 0;
		boolean foundA=false;
		for(int i=0; i<strand.length;i++){
			if(strand[i]==aValue){
				if(!foundA){
					foundA=true;
					count = 0;
				} else{
					count++;
				}
			} else{
				if(foundA && count>1){
					polyA[count]++;
					foundA=false;
				} 
			}
		}


		
		return polyA;
	}
	
	
	/**
	 * returns the C+G content of a DNA strand
	 * @param strand
	 * @return
	 */
	public static double computeGCcontent(byte[] strand){
		double result = 0;
		if(strand.length>0){
			byte cValue = bps.bpsID.get("C");
			byte gValue = bps.bpsID.get("G");
			for(int i=0;i<strand.length;i++){
				if(strand[i] == cValue || strand[i] == gValue){
					result++;
				}
			}
			
			result = (double) result/strand.length;			
		}
		
		return result;
	}
	
	/**
	 * computes the freqency of appearance of a sequence in a strand.
	 * @param strand the long DNA strand
	 * @param seq the short sequence
	 * @return
	 */
	public static double getFrequency(byte[] strand, ArrayList<Byte> seq, boolean doubleWay){
		double freq=0;
		
		if(seq.size()>0 && seq.size()<strand.length){
			boolean match;
			for(int i=0;i<strand.length - seq.size();i++){
				match =true;
				for (int j=0;j<seq.size() && match;j++){
					if(strand[i + j]!=seq.get(j)){
						match =false;
					}
				} 
				
				if(match){
					freq++;
				}
			}
			
			if(doubleWay){
				ArrayList<Byte> revSeq = getReversedComplementSequences(seq);
				for(int i=0;i<strand.length - revSeq.size();i++){
					match =true;
					for (int j=0;j<revSeq.size() && match;j++){
						if(strand[i + j]!=revSeq.get(j)){
							match =false;
						}
					} 
					
					if(match){
						freq++;
					}
				}
				
			}
			freq=(double) freq/(strand.length - seq.size());
		}
		return freq;
	}
	

	
	/**
	 * converts a string of DNA sequence into the corresponding vector of bytes
	 * @param strand the text
	 * @return
	 */
	public static ArrayList<Byte> getSequenceIDs(String strand){
		
		ArrayList<Byte> buffer = new ArrayList<Byte>();
		
		strand = strand.trim();
		strand = strand.toUpperCase();
		
		String letter;
		
		for(int i=0;i<strand.length();i++){
			letter= strand.substring(i, i+1);
			if(bps.bpsID.containsKey(letter)){
				buffer.add(bps.bpsID.get(letter));
			}
		}
		
				
		return buffer;
	}
	
	/**
	 * converts a string of DNA sequence into the corresponding vector of bytes
	 * @param strand the text
	 * @return
	 */
	public static byte[] getSeqIDs(String strand){
		
		
		byte[] buffer = new byte[strand.length()];
		
		strand = strand.trim();
		strand = strand.toUpperCase();
		
		String letter;
		
		for(int i=0;i<strand.length();i++){
			buffer[i] =Constants.NONE;
			letter= strand.substring(i, i+1);
			if(bps.bpsID.containsKey(letter)){
				buffer[i] = bps.bpsID.get(letter);
			}
		}
		
				
		return buffer;
	}
	
	/**
	 * returns the literal string of the byte specified sequence
	 * @param seq
	 * @return
	 */
	public static String sequenceToString(byte[] seq){
		String result="";
		
		if(seq!=null && seq.length > 0){
			for(byte b:seq){
				if(b>=0 && b<bps.bps.length){
					result+=bps.bps[b];
				}
			}
		}
		return result;
		
	}
	
	
	/**
	 * returns the revered complement of a sequence
	 * @param seq
	 * @return
	 */
	public static ArrayList<Byte> getReversedComplementSequences(ArrayList<Byte> seq){
		ArrayList<Byte> revSeq = new ArrayList<Byte>();
		for(int i=seq.size()-1; i>=0;i--){
			revSeq.add(bps.getComplement(seq.get(i)));
		}
		return revSeq;
	}
	
	
	/**
	 * returns the revered complement of a sequence
	 * @param seq
	 * @return
	 */
	public static byte[] getReversedComplementSequences(byte[] seq){
		byte[] revSeq = new byte[seq.length];
		int j=0;
		for(int i=seq.length-1; i>=0;i--){
			revSeq[j]=bps.getComplement(seq[i]);
			j++;
		}
		return revSeq;
	}
	
	
	/**
	 * returns the revered complement of a sub sequence
	 * @param seq the DNA sequence
	 * @param start the start position
	 * @param length the length of the sequence
	 * @return
	 */
	public static byte[] getReversedComplementSequences(byte[] seq, int start, int length){
		if(length <= 0 || ((start + length) >=seq.length)){
			length=seq.length;
			start=0;
		}
		byte[] revSeq = new byte[length];
		
		int j=0;
		
		for(int i=(length+start)-1; i>=start;i--){
			revSeq[j]=bps.getComplement(seq[i]);
			j++;
		}
		return revSeq;
	}
	
	
	
	/**
	 * computes the PWM score for 1 nucleotide
	 * @param freqInMotif
	 * @param freqInDNA
	 * @return
	 */
	public static double computePWMscore(double freqInMotif, double freqInDNA){
		return Utils.log2(freqInMotif/freqInDNA);
	}
	
	
	/**
	 * compare two sequences 
	 * @param seq1
	 * @param seq2
	 * @return
	 */
	public static boolean areSequencesEqual(byte[] seq1, byte[] seq2){
		boolean result=false;
		
		if(seq1.length==seq2.length && seq1.length>0){
			result=true;
			for(int i=0;i<seq1.length&& result;i++){
				if(seq1[i]!=seq2[i] && seq1[i] !=bpANYID && seq2[i]!=bpANYID){
					result=false;
				}
			}
			
		}		
		return result;
	}
	

	/**
	 * print a DNA sequence to a file
	 * @param path the path were to save the file
	 * @param filename the filename
	 * @param description the description
	 * @param seq the DNA sequence
	 */
	public static void printSequence(String path, String filename, String description, byte[] seq){
		try {
		    BufferedWriter out = new BufferedWriter(new FileWriter(new File(path,filename)));
		    out.write(">");
		    out.write(description);
		    out.newLine();
			int sectorSize =100;
			int sectors = (int) Math.ceil((double)seq.length/sectorSize);
			String buffer;
			for(int i=0;i<sectors;i++){
				buffer="";
				
				for(int j=i*sectorSize;j< Math.min((i+1)*sectorSize, seq.length);j++){
					buffer+=CellUtils.bps.bps[seq[j]];
				}
				out.write(buffer);
				out.newLine();
			}
		    
		    out.close();
		} catch (IOException e) {
		}
	}
	
	
	/**
	 * concatenates two sequences and returns the result
	 * @param seq1
	 * @param seq2
	 * @return
	 */
	public static byte[] concatenateDNAseq(byte[] seq1, byte[] seq2){
		byte[] result;
		result = new byte[seq1.length+seq2.length];
		for(int i=0;i<seq1.length;i++){
			result[i]= seq1[i];
		}
		for(int i=0;i<seq2.length;i++){
			result[i+seq1.length]= seq2[i];
		}
		return result;
	}
	
	
	/**
	 * creates a copy of a DNA sequence and returns it
	 * @param seq
	 * @return
	 */
	public static byte[] copySequence(byte[] seq){
		byte[] result=null;
		
		if(seq !=null && seq.length>0){
			result= new byte[seq.length];
			for(int i=0;i<seq.length;i++){
				result[i] = seq[i];
			}
		}
		
		return result;
	}
	
	/**
	 * replaces the DNA sequence in the strand from position pos with the subsequence subseq
	 * @param strand
	 * @param subSeq
	 * @param pos
	 * @return
	 */
	public static byte[] replaceDNAseq(byte[] strand, byte[] subSeq, int pos){
		byte[] result = copySequence(strand);
		
		if(pos > 0 && pos < strand.length){
			for(int i=pos;i<Math.min(pos+subSeq.length, strand.length);i++){
				if(subSeq[i-pos]!=bpANYID){
					result[i] = subSeq[i-pos];
				} else{
					//System.out.println("test");
					
				}
			}
		}
		return result;
	}
	
	
	/**
	 * generates a DNA strand containing only N
	 * @param length
	 * @return
	 */
	public static byte[] generateEmptyDNAStrand(int length){
		byte[] result = new byte[length];
		
		for(int i=0;i<length;i++){
			result[i]=bpANYID;
		}
		
		return result;
		
	}
	
	/**
	 * generates a DNA seq that contains only one bp
	 * @param length
	 * @param bp
	 * @return
	 */
	public static byte[] generateDNAStrand(int length, String bp){
		byte[] result = new byte[length];
		byte bpID=bps.bpsID.get(bp);
		for(int i=0;i<length;i++){
			result[i]=bpID;
		}
		
		return result;
	}
	
	/**
	 * extracts a subsequence form a sequence
	 * @param seq the original sequence
	 * @param start the start 
	 * @param end the end (exclusive)
	 * @return the subsequence
	 */
	public static byte[] extractSubSequence(byte[] seq, int start, int end){
		int s = Math.max(0, start);
		int e = Math.min(seq.length, end);
		
		byte[] result = null;
		
		if(e>s){
			result = new byte[e-s];
			for(int i=s;i<e;i++){
				result[i-s]=seq[i];
			}
		}
		
		return result;
		
	}
	
}
