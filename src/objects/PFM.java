package objects;

import java.io.Serializable;
import java.util.ArrayList;

import environment.Cell;
import utils.CellUtils;
import utils.Constants;
import utils.Utils;

/**
 * class that imports a Position Frequency Matrix 
 * @author created: n.r.zabet@gen.cam.ac.uk, modified: dmitrav@inbox.ru
 *
 */
public class PFM  implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1372756026390102207L;
	public ArrayList<ArrayList<Double>> pfm;
	public ArrayList<ArrayList<Double>> normPFM;
	public ArrayList<ArrayList<Double>> energyPFM;
	
	public int[] nucleotidePosition;
	public int motifSize;
	public double normSum;
	public boolean isCorrect;
	
	public boolean isInfo;
	public boolean isPWM;
	public double correction;
	/**
	 * class constructor
	 * an example of a correct PFM specification is A=[3, 1, 5, 7, 3, 6, 4, 7, 1, 9, 8, 5, 4, 2]; C=[1, 1, 2, 0, 0, 0, 1, 0, 8, 0, 0, 0, 0, 3]; G=[4, 1, 1, 1, 0, 0, 4, 1, 0, 0, 1, 3, 0, 2]; T=[1, 6, 1, 1, 6, 3, 0, 1, 0, 0, 0, 1, 5, 2]
	 * @param str
	 */
	public PFM(String str, Cell n){
		
		if(str.startsWith(Constants.DBD_TYPE_PWM)){
			str = str.replaceAll(Constants.DBD_TYPE_PWM, "").trim();
			isPWM = true;
			isInfo = true;
		}
		
		if(str.startsWith(Constants.DBD_TYPE_PFM)){
			str = str.replaceAll(Constants.DBD_TYPE_PFM, "").trim();
			isPWM = false;
		}
		
		pfm = new ArrayList<ArrayList<Double>>();
		motifSize = Constants.NONE;
		ArrayList<Double> bufferPFM;
		isCorrect = true;
		nucleotidePosition = new int[CellUtils.bps.numberOfBP];
		for(int i=0;i<nucleotidePosition.length;i++){
			nucleotidePosition[i]=Constants.NONE;
		}
		
		if(str.startsWith(Constants.DBD_TYPE_PFM_INFO)){
			str=str.replaceAll(Constants.DBD_TYPE_PFM_INFO, "");
			this.isInfo=true;
		}
		
		if(str.startsWith(Constants.DBD_TYPE_PFM_ENERGY)){
			str=str.replaceAll(Constants.DBD_TYPE_PFM_ENERGY, "");
			this.isInfo=false;
		}
				
		//correction so the log is not -infinity
		correction = Constants.DBD_TYPE_PFM_ENERGY_CORRECTION_DEFAULT;
		if(str.contains(Constants.DBD_TYPE_SEPARATOR)){
			String[] bufferStr = str.split(Constants.DBD_TYPE_SEPARATOR);
			if(bufferStr.length==2){
				correction = Utils.parseDouble(bufferStr[0], correction);
				str = bufferStr[1];
			}
		}
		

		if(str.contains(Constants.PFM_NUCLEOTIDE_SEPARATOR)){
			String[] bufferNucleotide, bufferNucleotideContainer;
			int nucleotideID, matrixID;
			String nucleotideKey, nucleotidePFM;
			
			bufferNucleotide = str.split(Constants.PFM_NUCLEOTIDE_SEPARATOR);
			
			for(String buffer:bufferNucleotide){
				if(buffer.contains(Constants.PFM_NUCLEOTIDE_ASSIGNMENT)){
					bufferNucleotideContainer = buffer.split(Constants.PFM_NUCLEOTIDE_ASSIGNMENT);				
					if(bufferNucleotideContainer.length == 2){
						nucleotideKey = bufferNucleotideContainer[0].trim();
						if(nucleotideKey.length()==1 && CellUtils.bps.bpsID.containsKey(nucleotideKey)){
							nucleotideID = CellUtils.bps.bpsID.get(nucleotideKey);

							nucleotidePFM = bufferNucleotideContainer[1].trim();
							nucleotidePFM = nucleotidePFM.substring(1);
							nucleotidePFM = nucleotidePFM.substring(0,nucleotidePFM.length()-1);
						
							if(nucleotidePFM.contains(Constants.PFM_NUCLEOTIDE_CONTAINER_SEPRATOR)){
								bufferPFM = new ArrayList<Double>();
								for(String bufferNucleotidePosPFM:nucleotidePFM.split(Constants.PFM_NUCLEOTIDE_CONTAINER_SEPRATOR)){
									bufferPFM.add(Utils.parseDouble(bufferNucleotidePosPFM.trim(), 0));
								}
								if(motifSize==Constants.NONE){
									motifSize=bufferPFM.size();
								}
								
								if(bufferPFM.size()!=motifSize){
									isCorrect =false;
									n.printDebugInfo("pfm " +str+ " various size rows pfm");
								}
								
								pfm.add(bufferPFM);
								matrixID = pfm.size()-1;
								this.nucleotidePosition[nucleotideID] = matrixID;
							}
						}
					}
					
				}
			}
		}
		//not all fiilll with zero
		if(isCorrect && pfm.size()!=CellUtils.bps.numberOfBP){
			int id;
			for(String buffer: BasePairs.bps){
				if(CellUtils.bps.bpsID.containsKey(buffer)){
					id = CellUtils.bps.bpsID.get(buffer);
					if(id!=CellUtils.bps.getANYID() && nucleotidePosition[id]==Constants.NONE){
						bufferPFM = new ArrayList<Double>();
						for(int i=0;i<motifSize;i++){
							bufferPFM.add(0.0);
						}
						pfm.add(bufferPFM);
						nucleotidePosition[id] =pfm.size()-1;
					}
				}
			}
		}	
		
		if(!this.isPWM){
			//check for correctness
			if(isCorrect){
				normSum = Constants.NONE;
				double bufferSum;
				
				for (int j=0;j<motifSize; j++){
					bufferSum=0;
					for(int i=0;i<this.pfm.size();i++){
						bufferSum+=this.pfm.get(i).get(j);
					}
					
					if(normSum==Constants.NONE){
						normSum=bufferSum;
					}
				
					if(Math.floor(normSum*100.0)!=Math.floor(bufferSum*100.0)){
						isCorrect=false;
						n.printDebugInfo("PFM "+str+" does not add to 1.");
					}
				}
			}
			


			//normalise
			/*
			if(isCorrect){
				double bufferMax;
				this.normPFM = new ArrayList<ArrayList<Double>>();
				ArrayList<Double> bufferNormPFM;
				this.energyPFM = new ArrayList<ArrayList<Double>>();
				ArrayList<Double> bufferEnergyPFM;
				
				for(int i=0;i<this.pfm.size();i++){
					bufferNormPFM = new ArrayList<Double>();
					bufferEnergyPFM = new ArrayList<Double>();
					bufferMax = Utils.getMax(this.pfm.get(i));
					for (int j=0;j<motifSize; j++){
						bufferNormPFM.add((double)(this.pfm.get(i).get(j)+correction)/(normSum+4*correction));
						bufferEnergyPFM.add(Math.log((double) (bufferMax+correction)/(this.pfm.get(i).get(j)+correction)));
					}
				
					this.normPFM.add(bufferNormPFM);
					this.energyPFM.add(bufferEnergyPFM);
				}
			}
			*/

		} else{
			this.normPFM = this.pfm;
		}

	}
	
	/**
	 * normalizes PFM with background frequencies 
	 * @param bpFreq
	 */
	public void normalizePFM(double[] bpFreq){
		if(!this.isPWM){
			//normalise
			if(isCorrect){
				this.normPFM = new ArrayList<ArrayList<Double>>();
				ArrayList<Double> bufferNormPFM;
				
				for(int i=0;i<this.pfm.size();i++){
					bufferNormPFM = new ArrayList<Double>();
					for (int j=0;j<motifSize; j++){
						bufferNormPFM.add(Math.log(((this.pfm.get(i).get(j)+correction*bpFreq[i]) /(normSum+correction)) /bpFreq[i]));
					}
					//System.out.println(bufferNormPFM);
					this.normPFM.add(bufferNormPFM);
				}
			}
		}
		
	}
	
	/**
	 * reutrns the normalised frequency of the nucleotide
	 * @param nucleotide current nucleotide
	 * @param position in the motif
	 * @return
	 */
	public double getNormPFM(byte nucleotide, int position){
		return normPFM.get(this.nucleotidePosition[nucleotide]).get(position);
	}
	
	
	/**
	 * reutrns the normalised frequency of the nucleotide
	 * @param nucleotide current nucleotide
	 * @param position in the motif
	 * @return
	 */
	public double getComplementNormPFM(byte nucleotide, int position){
		return normPFM.get(this.nucleotidePosition[CellUtils.bps.getComplement(nucleotide)]).get(position);
	}
	
	

	/**
	 * reutrns the normalised frequency of the nucleotide
	 * @param nucleotide current nucleotide
	 * @param position in the motif
	 * @return
	 */
	public double getEnergyPFM(byte nucleotide, int position){
		return energyPFM.get(this.nucleotidePosition[nucleotide]).get(position);
	}
	
	
	/**
	 * reutrns the normalised frequency of the nucleotide
	 * @param nucleotide current nucleotide
	 * @param position in the motif
	 * @return
	 */
	public double getComplementEnergyPFM(byte nucleotide, int position){
		return energyPFM.get(this.nucleotidePosition[CellUtils.bps.getComplement(nucleotide)]).get(position);
	}
	
	


	/**
	 * reutrns the normalised frequency of the nucleotide
	 * @param nucleotide current nucleotide
	 * @param position in the motif
	 * @return
	 */
	public double getScorePFM(byte nucleotide, int position, double freqInDNA){

		return normPFM.get(this.nucleotidePosition[nucleotide]).get(position);

		/*
		if(isInfo){
			return -normPFM.get(this.nucleotidePosition[nucleotide]).get(position);
		}
		return energyPFM.get(this.nucleotidePosition[nucleotide]).get(position);
		*/
	}

	/**
	 * AD
	 * for repression
	 * returns score
	 * @param nucleotide current nucleotide
	 * @param position in the motif
	 * @return
	 */
	public double getScorePFM(byte nucleotide, int position){

		return normPFM.get(this.nucleotidePosition[nucleotide]).get(position);

		/*
		if(isInfo){
			return -normPFM.get(this.nucleotidePosition[nucleotide]).get(position);
		}
		return energyPFM.get(this.nucleotidePosition[nucleotide]).get(position);
		*/
	}

	/**
	 * AD
	 * returns max PWM score for a position on the DNA
	 * @param position
	 * @return
     */
	public double getMaxScorePFM(int position){
		double max = 0;
		for (int i = 0; i < this.nucleotidePosition.length; i++) {
			if (normPFM.get(this.nucleotidePosition[i]).get(position) > max){
				max = normPFM.get(this.nucleotidePosition[i]).get(position);
			}
		}
		return max;
	}
	
	
	/**
	 * reutrns the normalised frequency of the nucleotide
	 * @param nucleotide current nucleotide
	 * @param position in the motif
	 * @return
	 */
	public double getComplementScorePFM(byte nucleotide, int position, double freqInDNA){
		if(isInfo){
			return -normPFM.get(this.nucleotidePosition[CellUtils.bps.getComplement(nucleotide)]).get(position);
		}
		return energyPFM.get(this.nucleotidePosition[CellUtils.bps.getComplement(nucleotide)]).get(position);
	}
	
	
	/**
	 * generates a string description of the class
	 * @return
	 */
	public String toString(double[] bpFreq){
		StringBuffer str= new StringBuffer();
		
		
		
		/*if(this.isPWM){
			str.append(Constants.DBD_TYPE_PWM);
		} else{
			str.append(Constants.DBD_TYPE_PFM);
			if(this.isInfo){
				str.append(Constants.DBD_TYPE_PFM_INFO);
			} else{
				str.append(Constants.DBD_TYPE_PFM_ENERGY);
			}
			str.append(correction);
			str.append(Constants.DBD_TYPE_SEPARATOR);
		}*/
		str.append(Constants.DBD_TYPE_PWM);
		
		if(this.isPWM){	
			for(int i=0;i<this.pfm.size();i++){	
				str.append(BasePairs.bps[i]);
				str.append("=[");
				for (int j=0;j<motifSize; j++){
					str.append(this.normPFM.get(i).get(j));			
					if(j<motifSize-1){
						str.append(", ");
					}
				}
				if(i<this.pfm.size()-1){
					str.append("]; ");
				} else{
					str.append("]");
				}
			}
		} else{
			
			for(int i=0;i<this.pfm.size();i++){	
				str.append(BasePairs.bps[i]);
				str.append("=[");
				for (int j=0;j<motifSize; j++){
					
					str.append(-getScorePFM((byte)i,j,bpFreq[i]));	
					
					
					if(j<motifSize-1){
						str.append(", ");
					}
				}
				if(i<this.pfm.size()-1){
					str.append("]; ");
				} else{
					str.append("]");
				}
			}
		}
		
		
		return str.toString();
	}
	
	
}
