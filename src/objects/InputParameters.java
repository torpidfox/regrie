package objects;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Serializable;
import java.util.ArrayList;

import objects.Parameter;

import utils.Constants;
import utils.Utils;

/**
 * class that contains all the parameters that are read from parameter files
 * @author created: n.r.zabet@gen.cam.ac.uk, modified: dmitrav@inbox.ru
 *
 */
public class InputParameters  implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = -2308322873169482087L;
	
	public String outputDirectory;
	
	//parameters set externally
	//SIMULATION PARAMATERS
	public Parameter<Double> STOP_TIME;
	public Parameter<Integer> ENSEMBLE_SIZE;
	public Parameter<Integer> RANDOM_SEED;
	public Parameter<Integer> COMPUTED_AFFINITY_PRECISION;
	public Parameter<Integer> DNA_SECTOR_SIZE;
	public Parameter<Integer> EVENT_LIST_SUBGROUP_SIZE;
	public Parameter<Boolean> EVENT_LIST_USES_FR;
	
	
	//SIMULATION-OUTPUT PARAMATERS
	public Parameter<String> OUTPUT_FOLDER;
	public Parameter<String> OUTPUT_FILENAME;
	public Parameter<Double> PRINT_INTERMEDIARY_RESULTS_AFTER;
	public Parameter<Boolean> PRINT_FINAL_OCCUPANCY;
	public Parameter<Boolean> DEBUG_MODE;
	public Parameter<String> OUTPUT_TF;
	public Parameter<Integer> OUTPUT_TF_POINTS;
	public Parameter<Boolean> FOLLOW_TS;
	public Parameter<Boolean> OUTPUT_AFFINITY_LANDSCAPE;
	public Parameter<Boolean> OUTPUT_BINDING_ENERGY;
	public Parameter<Boolean> OUTPUT_DNA_OCCUPANCY;
	public Parameter<Boolean> DNA_OCCUPANCY_FULL_MOLECULE_SIZE;
	public Parameter<Boolean> OUTPUT_SLIDING_LENGTHS;
	public Parameter<Integer> WIG_STEP;
	public Parameter<Double> WIG_THRESHOLD;

	//TF PARAMETERS
	public Parameter<String> TF_FILE;
	public Parameter<String> TS_FILE;
	public Parameter<String> TF_COOPERATIVITY_FILE;

	
	
	//TF_RANDOM PARAMETERS
	public Parameter<Integer> TF_DBD_LENGTH_MIN;
	public Parameter<Integer> TF_DBD_LENGTH_MAX;
	public Parameter<Integer> TF_SPECIES_COUNT;
	public Parameter<Integer> 	TF_COPY_NUMBER_MIN;
	public Parameter<Integer> TF_COPY_NUMBER_MAX;
	public Parameter<Double> TF_ES;
	public Parameter<Integer> TF_SIZE_LEFT;
	public Parameter<Integer> TF_SIZE_RIGHT;
	public Parameter<Double> TF_ASSOC_RATE;
	public Parameter<Boolean> TF_READ_IN_BOTH_DIRECTIONS;
	public Parameter<Double>  TF_PREBOUND_PROPORTION;
	public Parameter<Boolean>  TF_PREBOUND_TO_HIGHEST_AFFINITY;
	public Parameter<Boolean> SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE;

	//TF_REPRESSION PARAMETERS
	public Parameter<Boolean> REPRESSOR;
	public Parameter<Integer> TF_REPLENLEFT;
	public Parameter<Integer> TF_REPLENRIGHT;
	public Parameter<Double> PWM_REP_THRESHOLD;
	public Parameter<Double> REPRESSION_PROBABILITY;



	//DNA PARAMETERS
	public Parameter<String> DNA_SEQUENCE_FILE;
	public Parameter<String> DNA_BTRACK_FILE;


	//DNA_RANDOM PARAMETERS
	public Parameter<Integer> DNA_LENGTH;
	public Parameter<Double> DNA_PROPORTION_OF_A;
	public Parameter<Double> DNA_PROPORTION_OF_T;
	public Parameter<Double> DNA_PROPORTION_OF_C;
	public Parameter<Double> DNA_PROPORTION_OF_G;
	public Parameter<String> DNA_BOUNDARY_CONDITION;

	//TF RANDOM WALK PARAMETERS
	public Parameter<Boolean>  TF_IS_IMMOBILE;
	public Parameter<Double> TF_UNBINDING_PROBABILITY;
	public Parameter<Double> TF_SLIDE_LEFT_PROBABILITY;
	public Parameter<Double> TF_SLIDE_RIGHT_PROBABILITY;
	public Parameter<Double> TF_JUMPING_PROBABILITY;
	public Parameter<Double> TF_HOP_STD_DISPLACEMENT;
	public Parameter<Double> TF_SPECIFIC_WAITING_TIME;
	public Parameter<Integer> TF_STEP_LEFT_SIZE;
	public Parameter<Integer> TF_STEP_RIGHT_SIZE;
	public Parameter<Integer> TF_UNCORRELATED_DISPLACEMENT_SIZE;
	public Parameter<Boolean> TF_STALLS_IF_BLOCKED;
	public Parameter<Double> TF_COLLISION_UNBIND_PROBABILITY;
	public Parameter<Double> TF_AFFINITY_LANDSCAPE_ROUGHNESS;
	public Parameter<Boolean> CHECK_OCCUPANCY_ON_BINDING;
	public Parameter<Boolean> CHECK_OCCUPANCY_ON_SLIDING;
	public Parameter<Boolean> CHECK_OCCUPANCY_ON_REBINDING;
	public Parameter<Boolean> IS_BIASED_RANDOM_WALK;
	public Parameter<Boolean> IS_TWO_STATE_RANDOM_WALK;

	public InputParameters(String parametersFile){
		BufferedReader br ;
		InputStream in;
		InputStreamReader is;
		
		String default_filename="";
		//initialise the internal parameters
		initialiseInputParameters();
		
		if(System.getProperty("os.name").toUpperCase().contains("WINDOWS")){
			default_filename = Constants.DEFAULT_PARAMS_FILE_WIN;
		} else{
			default_filename = Constants.DEFAULT_PARAMS_FILE;
		}
	    try {

			if(!default_filename.isEmpty()){
				File defaultFile = new File(default_filename);
				if(!defaultFile.exists()){
					in = this.getClass().getClassLoader().getResourceAsStream(default_filename);
					//System.out.println("from resources");
					//System.out.println(in.available());
					if(in.available() >0){
						is = new InputStreamReader(in);
						br = new BufferedReader(is);
						loadInitialInputParameters(br);
					} else{
						default_filename = "";
						default_filename = Constants.DEFAULT_PARAMS_FILE_EMPTY;
						defaultFile = new File(default_filename);
						if(!defaultFile.exists()){
							default_filename = "";
							System.err.println("system_empty.ini file is missing from the current directory!");
							System.exit(0);
						} else{
								br = new BufferedReader(new FileReader(default_filename));
								loadInitialInputParameters(br);
						}
					}
				} else{
		            br = new BufferedReader(new FileReader(default_filename));
					loadInitialInputParameters(br);
				}
			}
			
			if(parametersFile!=null && !parametersFile.isEmpty()){
				br = new BufferedReader(new FileReader(parametersFile));
				loadParameters(br);
			}
	    } catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	private void initialiseInputParameters(){
		//SIMULATION PARAMATERS
		this.STOP_TIME = new Parameter<Double>();
		this.ENSEMBLE_SIZE = new Parameter<Integer>();
		this.RANDOM_SEED = new Parameter<Integer>();
		this.COMPUTED_AFFINITY_PRECISION = new Parameter<Integer>();
		this.DNA_SECTOR_SIZE = new Parameter<Integer>();
		this.EVENT_LIST_SUBGROUP_SIZE = new Parameter<Integer>();
		
		//SIMULATION-OUTPUT PARAMATERS
		this.OUTPUT_FOLDER = new Parameter<String>();
		this.OUTPUT_FILENAME = new Parameter<String>();
		this.PRINT_INTERMEDIARY_RESULTS_AFTER = new Parameter<Double>();
		this.PRINT_FINAL_OCCUPANCY = new Parameter<Boolean>();	
		this.DEBUG_MODE = new Parameter<Boolean>();
		this.OUTPUT_TF = new Parameter<String>();
		this.OUTPUT_TF_POINTS = new Parameter<Integer>();
		this.FOLLOW_TS = new Parameter<Boolean>();
		this.EVENT_LIST_USES_FR = new Parameter<Boolean>();	
		this.OUTPUT_AFFINITY_LANDSCAPE = new Parameter<Boolean>();
		this.OUTPUT_BINDING_ENERGY = new Parameter<Boolean>();
		this.OUTPUT_DNA_OCCUPANCY = new Parameter<Boolean>();
		this.DNA_OCCUPANCY_FULL_MOLECULE_SIZE = new Parameter<Boolean>();
		this.OUTPUT_SLIDING_LENGTHS = new Parameter<Boolean>();
		this.WIG_STEP = new Parameter<Integer>();
		this.WIG_THRESHOLD = new Parameter<Double>();




		
		//TF PARAMETERS
		this.TF_FILE= new Parameter<String>();
		this.TF_COOPERATIVITY_FILE= new Parameter<String>();
		this.TS_FILE= new Parameter<String>();

		
		
		
		//TF_RANDOM PARAMETERS
		this.TF_DBD_LENGTH_MIN = new Parameter<Integer>();
		this.TF_DBD_LENGTH_MAX = new Parameter<Integer>();
		this.TF_SPECIES_COUNT = new Parameter<Integer>();
		this.TF_COPY_NUMBER_MIN = new Parameter<Integer>();
		this.TF_COPY_NUMBER_MAX = new Parameter<Integer>();
		this.TF_ES = new Parameter<Double>();
		this.TF_SIZE_LEFT = new Parameter<Integer>();
		this.TF_SIZE_RIGHT = new Parameter<Integer>();
		this.TF_ASSOC_RATE = new Parameter<Double>();
		this.TF_READ_IN_BOTH_DIRECTIONS = new Parameter<Boolean>();
		this.TF_PREBOUND_PROPORTION = new Parameter<Double>();
		this.TF_PREBOUND_TO_HIGHEST_AFFINITY= new Parameter<Boolean>();
		this.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE = new Parameter<Boolean>();
		
		//DNA PARAMETERS
		this.DNA_SEQUENCE_FILE= new Parameter<String>();
		this.DNA_BTRACK_FILE= new Parameter<String>();

		
		//DNA_RANDOM PARAMETERS
		this.DNA_LENGTH = new Parameter<Integer>();
		this.DNA_PROPORTION_OF_A = new Parameter<Double>();
		this.DNA_PROPORTION_OF_T = new Parameter<Double>();
		this.DNA_PROPORTION_OF_C = new Parameter<Double>();
		this.DNA_PROPORTION_OF_G = new Parameter<Double>();
		this.DNA_BOUNDARY_CONDITION = new Parameter<String>();


		
		//TF RANDOM WALK PARAMETERS
		this.TF_IS_IMMOBILE= new Parameter<Boolean>();
		this.TF_UNBINDING_PROBABILITY = new Parameter<Double>();
		this.TF_SLIDE_LEFT_PROBABILITY = new Parameter<Double>();
		this.TF_SLIDE_RIGHT_PROBABILITY = new Parameter<Double>();
		this.TF_JUMPING_PROBABILITY = new Parameter<Double>();
		this.TF_HOP_STD_DISPLACEMENT = new Parameter<Double>();
		this.TF_SPECIFIC_WAITING_TIME = new Parameter<Double>();
		this.TF_STEP_LEFT_SIZE = new Parameter<Integer>();
		this.TF_STEP_RIGHT_SIZE = new Parameter<Integer>();
		this.TF_UNCORRELATED_DISPLACEMENT_SIZE = new Parameter<Integer>();
		this.TF_STALLS_IF_BLOCKED = new Parameter<Boolean>();
		this.TF_COLLISION_UNBIND_PROBABILITY = new Parameter<Double>();
		this.TF_AFFINITY_LANDSCAPE_ROUGHNESS = new Parameter<Double>();
		this.CHECK_OCCUPANCY_ON_BINDING = new Parameter<Boolean>();
		this.CHECK_OCCUPANCY_ON_SLIDING = new Parameter<Boolean>();
		this.CHECK_OCCUPANCY_ON_REBINDING = new Parameter<Boolean>();
		this.IS_BIASED_RANDOM_WALK = new Parameter<Boolean>();
		this.IS_TWO_STATE_RANDOM_WALK = new Parameter<Boolean>();

		//TF_REPRESSION PARAMETERS
		this.REPRESSOR = new Parameter<>();
		this.TF_REPLENLEFT = new Parameter<>();
		this.TF_REPLENRIGHT = new Parameter<>();
		this.PWM_REP_THRESHOLD = new Parameter<>();
		this.REPRESSION_PROBABILITY = new Parameter<>();
	}

	/**
	 * sets the default values for the parameters provided from outside 
	 */
	private void loadInitialInputParameters(BufferedReader br){
		
		ArrayList<ArrayList<String>> params = getParamArray(br, Constants.PARAMS_FILE_COMMENT_CHAR, Constants.PARAMS_FILE_LINE_ENDING, Constants.PARAMS_FILE_ASSIGNMENT_CHAR);
		String name = "", label = "", description = "", category = "", value = "";
		
		
		String lastLoadParam = "";
		for(int i=0; i< params.size();i++){
			if(params.get(i).get(0).equals("name")){
				if(!name.isEmpty()){
					setParameter(name, label, description, category, value);
					lastLoadParam = name;
				}
				name = params.get(i).get(1);
			} else if(params.get(i).get(0).equals("label")){
				label = params.get(i).get(1);
			} else if(params.get(i).get(0).equals("description")){
				description = params.get(i).get(1);
			} else if(params.get(i).get(0).equals("category")){
				category = params.get(i).get(1);
			} else if(params.get(i).get(0).equals("value")){
				value = params.get(i).get(1);
			} 

		}
		
		// load the last one
		if(!name.isEmpty() && !lastLoadParam.equals(name)){
			setParameter(name, label, description, category, value);
		}

		
	}
	
	/**
	 * loads parameter values from a file
	 * @param filename
	 */
	public void loadParameters(BufferedReader br){
		
		ArrayList<ArrayList<String>> params = getParamArray(br, Constants.PARAMS_FILE_COMMENT_CHAR, Constants.PARAMS_FILE_LINE_ENDING, Constants.PARAMS_FILE_ASSIGNMENT_CHAR);
		String name = "", value = "";
		
		for(int i=0; i< params.size();i++){
			name = params.get(i).get(0);
			value = params.get(i).get(1);
				if(!name.isEmpty()){
					setParameter(name, "", "", "", value);
				}
				
		}
	}
	

	


	
	/**
	 * copys the parameter file to a new file
	 * @param parameterFilename
	 */
	public File exportParameterFile(String parameterFilename){
		File result=null;
		String filename="params_";
		boolean fileCreated=false;
		
		
		File directory;
		
		
		directory= new File(".");
		this.outputDirectory = directory.getPath()+File.separator;

		
		//if directory doesn't exists then create it.
		if(!this.OUTPUT_FOLDER.value.trim().isEmpty()){
			directory = new File(this.OUTPUT_FOLDER.value.trim());
			if(!directory.exists()){
				if(!directory.mkdir()){
					this.OUTPUT_FOLDER.value = "";
				}
			}
			this.outputDirectory = directory.getPath()+File.separator;

		}
		
		
		
		    try {
		    	
		    	if(parameterFilename!=""){
		    		filename = parameterFilename;
		    		result = new File(filename);
		    		fileCreated = true;
		    	} else if( !this.OUTPUT_FILENAME.value.trim().isEmpty()){
		    		if(this.OUTPUT_FILENAME.value.endsWith(Constants.PARAMETR_FILE_EXTENSION)){
		    			filename= this.OUTPUT_FILENAME.value.replaceAll(Constants.PARAMETR_FILE_EXTENSION, "_params_");
		    		} else{
		    			filename= this.OUTPUT_FILENAME.value + "_params_";
		    		}
		    		//result = new File(this.OUTPUT_FOLDER.value, filename);
		    		
		    		if(this.OUTPUT_FOLDER.value.isEmpty()){
		    			directory = new File(".");
		    		} else{
			    	  	directory = new File(this.OUTPUT_FOLDER.value);
		    		}
				this.outputDirectory = directory.getPath()+File.separator;

		    	  	
			    	result = File.createTempFile(filename, Constants.PARAMETR_FILE_EXTENSION,directory);
			    	filename = result.getName();
		    		
		    		fileCreated = true;
		    		
		    		/*if(result.exists()){
		    			fileCreated=false;
		    			filename = this.OUTPUT_FILENAME.value + "_params_";
		    		}*/
		    	}
		    	
		    	if(!fileCreated){
		    		  // Create temporary file.
		    		if(this.OUTPUT_FOLDER.value.isEmpty()){
		    			directory = new File(".");
		    		} else{
			    	  	directory = new File(this.OUTPUT_FOLDER.value);
		    		}
				
		    		this.outputDirectory = directory.getPath()+File.separator;

		    		result = File.createTempFile(filename, Constants.PARAMETR_FILE_EXTENSION,directory);
			    	filename = result.getName();
		    	}

		    		// Write to temp file
		        BufferedWriter out = new BufferedWriter(new FileWriter(result));

		        
		      //SIMULATION PARAMATERS		
		      out.write("#SIMULATION PARAMATERS\n\n");
		      out.write("#"+this.STOP_TIME.description+"\n");
		      out.write("STOP_TIME = "+this.STOP_TIME.value+";\n\n");
		      out.write("#"+this.ENSEMBLE_SIZE.description+"\n");
		      out.write("ENSAMBLE_SIZE = "+this.ENSEMBLE_SIZE.value+";\n\n");	      
		      out.write("#"+this.RANDOM_SEED.description+"\n");
		      out.write("RANDOM_SEED = "+this.RANDOM_SEED.value+";\n\n");
		      out.write("#"+this.OUTPUT_FOLDER.description+"\n");
		      out.write("#"+this.COMPUTED_AFFINITY_PRECISION.description+"\n");
		      out.write("COMPUTED_AFFINITY_PRECISION = "+this.COMPUTED_AFFINITY_PRECISION.value+";\n\n"); 
		      out.write("#"+this.DNA_SECTOR_SIZE.description+"\n");
		      out.write("DNA_SECTOR_SIZE = "+this.DNA_SECTOR_SIZE.value+";\n\n"); 
		      out.write("#"+this.EVENT_LIST_SUBGROUP_SIZE.description+"\n");
		      out.write("EVENT_LIST_SUBGROUP_SIZE = "+this.EVENT_LIST_SUBGROUP_SIZE.value+";\n\n"); 
		      out.write("#"+this.EVENT_LIST_USES_FR.description+"\n");
		      out.write("EVENT_LIST_USES_FR = "+this.EVENT_LIST_USES_FR.value+";\n\n");
		      
		      //SIMULATION-OUTPUT PARAMATERS		
		      out.write("OUTPUT_FOLDER = \""+this.OUTPUT_FOLDER.value+"\";\n\n");
		      out.write("#"+this.OUTPUT_FILENAME.description+"\n");
		      out.write("OUTPUT_FILENAME = \""+this.OUTPUT_FILENAME.value+"\";\n\n");
		      
		      out.write("#"+this.PRINT_INTERMEDIARY_RESULTS_AFTER.description+"\n");
		      out.write("PRINT_INTERMEDIARY_RESULTS_AFTER = \""+this.PRINT_INTERMEDIARY_RESULTS_AFTER.value+"\";\n\n");		      
		      out.write("#"+this.PRINT_FINAL_OCCUPANCY.description+"\n");
		      out.write("PRINT_FINAL_OCCUPANCY = \""+this.PRINT_FINAL_OCCUPANCY.value+"\";\n\n");
		      out.write("#"+this.DEBUG_MODE.description+"\n");
		      out.write("DEBUG_MODE = "+this.DEBUG_MODE.value+";\n\n");
		      out.write("#"+this.OUTPUT_TF.description+"\n");
		      out.write("OUTPUT_TF = \""+this.OUTPUT_TF.value+"\";\n\n");
		      out.write("#"+this.OUTPUT_TF_POINTS.description+"\n");
		      out.write("OUTPUT_TF_POINTS = "+this.OUTPUT_TF_POINTS.value+";\n\n"); 
		      out.write("#"+this.FOLLOW_TS.description+"\n");
		      out.write("FOLLOW_TS = "+this.FOLLOW_TS.value+";\n\n");
		      out.write("#"+this.OUTPUT_AFFINITY_LANDSCAPE.description+"\n");
		      out.write("OUTPUT_AFFINITY_LANDSCAPE = "+this.OUTPUT_AFFINITY_LANDSCAPE.value+";\n\n"); 
		      out.write("#"+this.OUTPUT_BINDING_ENERGY.description+"\n");
		      out.write("OUTPUT_BINDING_ENERGY = "+this.OUTPUT_BINDING_ENERGY.value+";\n\n"); 
		      out.write("#"+this.OUTPUT_DNA_OCCUPANCY.description+"\n");
		      out.write("OUTPUT_DNA_OCCUPANCY = "+this.OUTPUT_DNA_OCCUPANCY.value+";\n\n"); 
		      out.write("#"+this.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.description+"\n");
		      out.write("DNA_OCCUPANCY_FULL_MOLECULE_SIZE = "+this.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.value+";\n\n"); 
		      out.write("#"+this.OUTPUT_SLIDING_LENGTHS.description+"\n");
		      out.write("OUTPUT_SLIDING_LENGTHS = "+this.OUTPUT_SLIDING_LENGTHS.value+";\n\n"); 
		      out.write("#"+this.WIG_STEP.description+"\n");
		      out.write("WIG_STEP = "+this.WIG_STEP.value+";\n\n"); 
		      out.write("#"+this.WIG_THRESHOLD.description+"\n");
		      out.write("WIG_THRESHOLD = "+this.WIG_THRESHOLD.value+";\n\n"); 


		      
		      //TF PARAMETERS
		      out.write("#TF PARAMATERS\n\n");
		      out.write("#"+this.TF_FILE.description+"\n");
		      out.write("TF_FILE = \""+this.TF_FILE.value+"\";\n\n");
		      out.write("#"+this.TF_COOPERATIVITY_FILE.description+"\n");
		      out.write("TF_COOPERATIVITY_FILE = \""+this.TF_COOPERATIVITY_FILE.value+"\";\n\n");  	
		      out.write("#"+this.TS_FILE.description+"\n");
		      out.write("TS_FILE = \""+this.TS_FILE.value+"\";\n\n");  	
		      
		      
		      
			  //TF_RANDOM PARAMETERS
		      out.write("#TF_RANDOM PARAMATERS\n\n");
		      out.write("#"+this.TF_DBD_LENGTH_MIN.description+"\n");
		      out.write("TF_DBD_LENGTH_MIN = "+this.TF_DBD_LENGTH_MIN.value+";\n\n");
		      out.write("#"+this.TF_DBD_LENGTH_MAX.description+"\n");
		      out.write("TF_DBD_LENGTH_MAX = "+this.TF_DBD_LENGTH_MAX.value+";\n\n");
		      out.write("#"+this.TF_SPECIES_COUNT.description+"\n");
		      out.write("TF_SPECIES_COUNT = "+this.TF_SPECIES_COUNT.value+";\n\n");
		      out.write("#"+this.TF_COPY_NUMBER_MIN.description+"\n");
		      out.write("TF_COPY_NUMBER_MIN = "+this.TF_COPY_NUMBER_MIN.value+";\n\n");
		      out.write("#"+this.TF_COPY_NUMBER_MAX.description+"\n");
		      out.write("TF_COPY_NUMBER_MAX = "+this.TF_COPY_NUMBER_MAX.value+";\n\n");
		      out.write("#"+this.TF_ES.description+"\n");
		      out.write("TF_ES = "+this.TF_ES.value+";\n\n");
		      out.write("#"+this.TF_SIZE_LEFT.description+"\n");
		      out.write("TF_SIZE_LEFT = "+this.TF_SIZE_LEFT.value+";\n\n");
		      out.write("#"+this.TF_SIZE_RIGHT.description+"\n");
		      out.write("TF_SIZE_RIGHT = "+this.TF_SIZE_RIGHT.value+";\n\n");
		      out.write("#"+this.TF_ASSOC_RATE.description+"\n");
		      out.write("TF_ASSOC_RATE = "+this.TF_ASSOC_RATE.value+";\n\n");  
		      out.write("#"+this.TF_READ_IN_BOTH_DIRECTIONS.description+"\n");
		      out.write("TF_READ_IN_BOTH_DIRECTIONS = "+this.TF_READ_IN_BOTH_DIRECTIONS.value+";\n\n");  
		      out.write("#"+this.TF_PREBOUND_PROPORTION.description+"\n");
		      out.write("TF_PREBOUND_PROPORTION = "+this.TF_PREBOUND_PROPORTION.value+";\n\n");		     	
		      out.write("#"+this.TF_PREBOUND_TO_HIGHEST_AFFINITY.description+"\n");
		      out.write("TF_PREBOUND_TO_HIGHEST_AFFINITY = "+this.TF_PREBOUND_TO_HIGHEST_AFFINITY.value+";\n\n");		
		      out.write("#"+this.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.description+"\n");
		      out.write("SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE = "+this.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.value+";\n\n");

				//TF REPRESSION PARAMETERS
				out.write("#TF_REPRESSION PARAMETERS\n\n");
				out.write("#"+this.REPRESSOR.description+"\n");
				out.write("TF_REPRESSES = "+this.REPRESSOR.value+";\n\n");
				out.write("#"+this.TF_REPLENLEFT.description+"\n");
				out.write("TF_REPLENLEFT = "+this.TF_REPLENLEFT.value+";\n\n");
				out.write("#"+this.TF_REPLENRIGHT.description+"\n");
				out.write("TF_REPLENRIGHT = "+this.TF_REPLENRIGHT.value+";\n\n");
				out.write("#"+this.PWM_REP_THRESHOLD.description+"\n");
				out.write("PWM_REP_THRESHOLD = "+this.PWM_REP_THRESHOLD.value+";\n\n");
				out.write("#"+this.REPRESSION_PROBABILITY.description+"\n");
				out.write("REPRESSION_PROBABILITY = "+this.REPRESSION_PROBABILITY.value+";\n\n");

			  //DNA PARAMETERS
		      out.write("#DNA PARAMETERS\n\n");
		      out.write("#"+this.DNA_SEQUENCE_FILE.description+"\n");
		      out.write("DNA_SEQUENCE_FILE = \""+this.DNA_SEQUENCE_FILE.value+"\";\n\n");
		      out.write("#"+this.DNA_BTRACK_FILE.description+"\n");
		      out.write("DNA_BTRACK_FILE = \""+this.DNA_BTRACK_FILE.value+"\";\n\n");


		      
				//DNA_RANDOM PARAMETERS
		      out.write("#DNA_RANDOM PARAMATERS\n\n");
		      out.write("#"+this.DNA_LENGTH.description+"\n");
		      out.write("DNA_LENGTH = "+this.DNA_LENGTH.value+";\n\n");
		      out.write("#"+this.DNA_PROPORTION_OF_A.description+"\n");
		      out.write("DNA_PROPORTION_OF_A = "+this.DNA_PROPORTION_OF_A.value+";\n\n");
		      out.write("#"+this.DNA_PROPORTION_OF_T.description+"\n");
		      out.write("DNA_PROPORTION_OF_T = "+this.DNA_PROPORTION_OF_T.value+";\n\n");
		      out.write("#"+this.DNA_PROPORTION_OF_C.description+"\n");
		      out.write("DNA_PROPORTION_OF_C = "+this.DNA_PROPORTION_OF_C.value+";\n\n");
		      out.write("#"+this.DNA_PROPORTION_OF_G.description+"\n");
		      out.write("DNA_PROPORTION_OF_G = "+this.DNA_PROPORTION_OF_G.value+";\n\n");
		      out.write("#"+this.DNA_BOUNDARY_CONDITION.description+"\n");
		      out.write("DNA_BOUNDARY_CONDITION = "+this.DNA_BOUNDARY_CONDITION.value+";\n\n");
		      

			//TF RANDOM WALK PARAMETERS
		      out.write("#"+this.TF_IS_IMMOBILE.description+"\n");
		      out.write("TF_IS_IMMOBILE = "+this.TF_IS_IMMOBILE.value+";\n\n");		
		      out.write("#TF_RANDOM_WALK  PARAMATERS\n\n");
		      out.write("#"+this.TF_UNBINDING_PROBABILITY.description+"\n");
		      out.write("TF_UNBINDING_PROBABILITY = "+this.TF_UNBINDING_PROBABILITY.value+";\n\n");
		      out.write("#"+this.TF_SLIDE_LEFT_PROBABILITY.description+"\n");
		      out.write("TF_SLIDE_LEFT_PROBABILITY = "+this.TF_SLIDE_LEFT_PROBABILITY.value+";\n\n");
		      out.write("#"+this.TF_SLIDE_RIGHT_PROBABILITY.description+"\n");
		      out.write("TF_SLIDE_RIGHT_PROBABILITY = "+this.TF_SLIDE_RIGHT_PROBABILITY.value+";\n\n");
		      out.write("#"+this.TF_JUMPING_PROBABILITY.description+"\n");
		      out.write("TF_JUMPING_PROBABILITY = "+this.TF_JUMPING_PROBABILITY.value+";\n\n");
		      out.write("#"+this.TF_HOP_STD_DISPLACEMENT.description+"\n");
		      out.write("TF_HOP_STD_DISPLACEMENT = "+this.TF_HOP_STD_DISPLACEMENT.value+";\n\n");
		      out.write("#"+this.TF_SPECIFIC_WAITING_TIME.description+"\n");
		      out.write("TF_SPECIFIC_WAITING_TIME = "+this.TF_SPECIFIC_WAITING_TIME.value+";\n\n");      
		      out.write("#"+this.TF_STEP_LEFT_SIZE.description+"\n");
		      out.write("TF_STEP_LEFT_SIZE = "+this.TF_STEP_LEFT_SIZE.value+";\n\n");      
		      out.write("#"+this.TF_STEP_RIGHT_SIZE.description+"\n");
		      out.write("TF_STEP_RIGHT_SIZE = "+this.TF_STEP_RIGHT_SIZE.value+";\n\n");      
		      out.write("#"+this.TF_UNCORRELATED_DISPLACEMENT_SIZE.description+"\n");
		      out.write("TF_UNCORRELATED_DISPLACEMENT_SIZE = "+this.TF_UNCORRELATED_DISPLACEMENT_SIZE.value+";\n\n");      
		      out.write("#"+this.TF_STALLS_IF_BLOCKED.description+"\n");
		      out.write("TF_STALLS_IF_BLOCKED = "+this.TF_STALLS_IF_BLOCKED.value+";\n\n");      
		      out.write("#"+this.TF_COLLISION_UNBIND_PROBABILITY.description+"\n");
		      out.write("TF_COLLISION_UNBIND_PROBABILITY = "+this.TF_COLLISION_UNBIND_PROBABILITY.value+";\n\n");      
		      out.write("#"+this.TF_AFFINITY_LANDSCAPE_ROUGHNESS.description+"\n");
		      out.write("TF_AFFINITY_LANDSCAPE_ROUGHNESS = "+this.TF_AFFINITY_LANDSCAPE_ROUGHNESS.value+";\n\n");    
		      out.write("#"+this.CHECK_OCCUPANCY_ON_BINDING.description+"\n");
		      out.write("CHECK_OCCUPANCY_ON_BINDING = "+this.CHECK_OCCUPANCY_ON_BINDING.value+";\n\n");   
		      out.write("#"+this.CHECK_OCCUPANCY_ON_SLIDING.description+"\n");
		      out.write("CHECK_OCCUPANCY_ON_SLIDING = "+this.CHECK_OCCUPANCY_ON_SLIDING.value+";\n\n");   
		      out.write("#"+this.CHECK_OCCUPANCY_ON_REBINDING.description+"\n");
		      out.write("CHECK_OCCUPANCY_ON_REBINDING = "+this.CHECK_OCCUPANCY_ON_REBINDING.value+";\n\n");   
		      out.write("#"+this.IS_BIASED_RANDOM_WALK.description+"\n");
		      out.write("IS_BIASED_RANDOM_WALK = "+this.IS_BIASED_RANDOM_WALK.value+";\n\n");   
		      out.write("#"+this.IS_TWO_STATE_RANDOM_WALK.description+"\n");
		      out.write("IS_TWO_STATE_RANDOM_WALK = "+this.IS_TWO_STATE_RANDOM_WALK.value+";\n\n");   
		      
		      
		      out.close();

		    } catch (IOException e) {
		    	e.printStackTrace();
		    }
		    

		    
		    return result;
	}
	
	/**
	 * extracts an array list for parameter value from a file 
	 * @param filename the parameters file
	 * @param comment comment lines start with
	 * @param lineEnding the line ends with
	 * @param assignment the assignment sign 
	 * @return
	 */
	private ArrayList<ArrayList<String>> getParamArray(BufferedReader reader, String comment, String lineEnding, String assignment){
		ArrayList<ArrayList<String>> result = new ArrayList<ArrayList<String>> ();
		String[] buffer = new String[2];
		ArrayList<String> param; 
		    try{
	            String text = null;
	            while ((text = reader.readLine()) != null)
	            {
	            	text=text.trim();
	            	if(!text.isEmpty() && !text.startsWith(comment)){
	            		if(text.endsWith(lineEnding)){
	            			buffer = Utils.extractParameterFromCommandLine(text, assignment);
	            			if(buffer.length==2 && !buffer[0].isEmpty()){
	            				param = new ArrayList<String>();
	            				param.add(buffer[0]);
	            				param.add(buffer[1]);
	            				result.add(param);
	            			}
	            		} else{
	            			System.out.println("Error interpretting input line:"+text);
	            		}
	            	} 
	            }
	        } catch (Exception e) {
	        	System.out.println(e.toString());
            } finally {
                try{
                    if (reader != null){
                        reader.close();
                    }
                } catch (IOException e){
                	System.out.println(e.toString());
                }
           	}  
		return result;
	}
	
	
	/**
	 * sets a parameter
	 * @param name
	 * @param label
	 * @param description
	 * @param category
	 * @param value
	 * @return
	 */
	public boolean setParameter(String name, String label, String description, String category,String value){
		boolean found = false;
		//SIMULATION PARAMATERS
		if(name.equals("STOP_TIME")){
			this.STOP_TIME.value = Utils.parseDouble(value, Constants.NONE);
			if(!label.isEmpty()){this.STOP_TIME.label = label;}
			if(!description.isEmpty()){this.STOP_TIME.description = description;}
			if(!category.isEmpty()){this.STOP_TIME.category = category;}
			found = true;
		} else if(name.equals("ENSAMBLE_SIZE")){
			this.ENSEMBLE_SIZE.value =  Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.ENSEMBLE_SIZE.label = label;}
			if(!description.isEmpty()){this.ENSEMBLE_SIZE.description = description;}
			if(!category.isEmpty()){this.ENSEMBLE_SIZE.category = category;}
			found = true;
		} else if(name.equals("RANDOM_SEED")){
			this.RANDOM_SEED.value =  Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.RANDOM_SEED.label = label;}
			if(!description.isEmpty()){this.RANDOM_SEED.description = description;}
			if(!category.isEmpty()){this.RANDOM_SEED.category = category;}
			found = true;
		} else if(name.equals("COMPUTED_AFFINITY_PRECISION")){
			this.COMPUTED_AFFINITY_PRECISION.value = Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.COMPUTED_AFFINITY_PRECISION.label = label;}
			if(!description.isEmpty()){this.COMPUTED_AFFINITY_PRECISION.description = description;}
			if(!category.isEmpty()){this.COMPUTED_AFFINITY_PRECISION.category = category;}
			found = true;	
		} else if(name.equals("DNA_SECTOR_SIZE")){
			this.DNA_SECTOR_SIZE.value = Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.DNA_SECTOR_SIZE.label = label;}
			if(!description.isEmpty()){this.DNA_SECTOR_SIZE.description = description;}
			if(!category.isEmpty()){this.DNA_SECTOR_SIZE.category = category;}
			found = true;	
		} else if(name.equals("EVENT_LIST_SUBGROUP_SIZE")){
			this.EVENT_LIST_SUBGROUP_SIZE.value = Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.EVENT_LIST_SUBGROUP_SIZE.label = label;}
			if(!description.isEmpty()){this.EVENT_LIST_SUBGROUP_SIZE.description = description;}
			if(!category.isEmpty()){this.EVENT_LIST_SUBGROUP_SIZE.category = category;}
			found = true;	
		} else if(name.equals("EVENT_LIST_USES_FR")){
			this.EVENT_LIST_USES_FR.value = Utils.parseBoolean(value,false);
			if(!label.isEmpty()){this.EVENT_LIST_USES_FR.label = label;}
			if(!description.isEmpty()){this.EVENT_LIST_USES_FR.description = description;}
			if(!category.isEmpty()){this.EVENT_LIST_USES_FR.category = category;}
			found = true;
		} 
		//SIMULATION-OUTPUT PARAMATERS
		else if(name.equals("OUTPUT_FOLDER")){
			this.OUTPUT_FOLDER.value = value;
			if(!label.isEmpty()){this.OUTPUT_FOLDER.label = label;}
			if(!description.isEmpty()){this.OUTPUT_FOLDER.description = description;}
			if(!category.isEmpty()){this.OUTPUT_FOLDER.category = category;}
			found = true;
		} else if(name.equals("OUTPUT_FILENAME")){
			this.OUTPUT_FILENAME.value = value;
			if(!label.isEmpty()){this.OUTPUT_FILENAME.label = label;}
			if(!description.isEmpty()){this.OUTPUT_FILENAME.description = description;}
			if(!category.isEmpty()){this.OUTPUT_FILENAME.category = category;}
			found = true;
		} else if(name.equals("PRINT_INTERMEDIARY_RESULTS_AFTER")){
			this.PRINT_INTERMEDIARY_RESULTS_AFTER.value = Utils.parseDouble(value, Constants.NONE);
			if(!label.isEmpty()){this.PRINT_INTERMEDIARY_RESULTS_AFTER.label = label;}
			if(!description.isEmpty()){this.PRINT_INTERMEDIARY_RESULTS_AFTER.description = description;}
			if(!category.isEmpty()){this.PRINT_INTERMEDIARY_RESULTS_AFTER.category = category;}
			found = true;
		} else if(name.equals("PRINT_FINAL_OCCUPANCY")){
			this.PRINT_FINAL_OCCUPANCY.value =  Utils.parseBoolean(value,false);
			if(!label.isEmpty()){this.PRINT_FINAL_OCCUPANCY.label = label;}
			if(!description.isEmpty()){this.PRINT_FINAL_OCCUPANCY.description = description;}
			if(!category.isEmpty()){this.PRINT_FINAL_OCCUPANCY.category = category;}
			found = true;	      
		} else if(name.equals("DEBUG_MODE")){
			this.DEBUG_MODE.value = Utils.parseBoolean(value,false);
			if(!label.isEmpty()){this.DEBUG_MODE.label = label;}
			if(!description.isEmpty()){this.DEBUG_MODE.description = description;}
			if(!category.isEmpty()){this.DEBUG_MODE.category = category;}
			found = true;
		} else if(name.equals("OUTPUT_TF")){
			this.OUTPUT_TF.value = value;
			if(!label.isEmpty()){this.OUTPUT_TF.label = label;}
			if(!description.isEmpty()){this.OUTPUT_TF.description = description;}
			if(!category.isEmpty()){this.OUTPUT_TF.category = category;}
			found = true;		
		} else if(name.equals("OUTPUT_TF_POINTS")){
			this.OUTPUT_TF_POINTS.value = Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.OUTPUT_TF_POINTS.label = label;}
			if(!description.isEmpty()){this.OUTPUT_TF_POINTS.description = description;}
			if(!category.isEmpty()){this.OUTPUT_TF_POINTS.category = category;}
			found = true;	
		} else if(name.equals("FOLLOW_TS")){
			this.FOLLOW_TS.value = Utils.parseBoolean(value,false);
			if(!label.isEmpty()){this.FOLLOW_TS.label = label;}
			if(!description.isEmpty()){this.FOLLOW_TS.description = description;}
			if(!category.isEmpty()){this.FOLLOW_TS.category = category;}
			found = true;			
		} else if(name.equals("OUTPUT_AFFINITY_LANDSCAPE")){
			this.OUTPUT_AFFINITY_LANDSCAPE.value = Utils.parseBoolean(value,false);
			if(!label.isEmpty()){this.OUTPUT_AFFINITY_LANDSCAPE.label = label;}
			if(!description.isEmpty()){this.OUTPUT_AFFINITY_LANDSCAPE.description = description;}
			if(!category.isEmpty()){this.OUTPUT_AFFINITY_LANDSCAPE.category = category;}
			found = true;
		} else if(name.equals("OUTPUT_BINDING_ENERGY")){
			this.OUTPUT_BINDING_ENERGY.value = Utils.parseBoolean(value,false);
			if(!label.isEmpty()){this.OUTPUT_BINDING_ENERGY.label = label;}
			if(!description.isEmpty()){this.OUTPUT_BINDING_ENERGY.description = description;}
			if(!category.isEmpty()){this.OUTPUT_BINDING_ENERGY.category = category;}
			found = true;
		} else if(name.equals("OUTPUT_DNA_OCCUPANCY")){
			this.OUTPUT_DNA_OCCUPANCY.value = Utils.parseBoolean(value,false);
			if(!label.isEmpty()){this.OUTPUT_DNA_OCCUPANCY.label = label;}
			if(!description.isEmpty()){this.OUTPUT_DNA_OCCUPANCY.description = description;}
			if(!category.isEmpty()){this.OUTPUT_DNA_OCCUPANCY.category = category;}
			found = true;
		} else if(name.equals("DNA_OCCUPANCY_FULL_MOLECULE_SIZE")){
			this.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.value = Utils.parseBoolean(value,false);
			if(!label.isEmpty()){this.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.label = label;}
			if(!description.isEmpty()){this.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.description = description;}
			if(!category.isEmpty()){this.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.category = category;}
			found = true;
		} else if(name.equals("OUTPUT_SLIDING_LENGTHS")){
			this.OUTPUT_SLIDING_LENGTHS.value = Utils.parseBoolean(value,false);
			if(!label.isEmpty()){this.OUTPUT_SLIDING_LENGTHS.label = label;}
			if(!description.isEmpty()){this.OUTPUT_SLIDING_LENGTHS.description = description;}
			if(!category.isEmpty()){this.OUTPUT_SLIDING_LENGTHS.category = category;}
			found = true;
		} else if(name.equals("WIG_STEP")){
			this.WIG_STEP.value = Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.WIG_STEP.label = label;}
			if(!description.isEmpty()){this.WIG_STEP.description = description;}
			if(!category.isEmpty()){this.WIG_STEP.category = category;}
			found = true;
		} else if(name.equals("WIG_THRESHOLD")){
			this.WIG_THRESHOLD.value = Utils.parseDouble(value, Constants.NONE);
			if(!label.isEmpty()){this.WIG_THRESHOLD.label = label;}
			if(!description.isEmpty()){this.WIG_THRESHOLD.description = description;}
			if(!category.isEmpty()){this.WIG_THRESHOLD.category = category;}
			found = true;	
		}       	
		//TF PARAMETERS
		else if(name.equals("TF_FILE")){
			this.TF_FILE.value = value;
			if(!label.isEmpty()){this.TF_FILE.label = label;}
			if(!description.isEmpty()){this.TF_FILE.description = description;}
			if(!category.isEmpty()){this.TF_FILE.category = category;}
			found = true;	
		}  else if(name.equals("TF_COOPERATIVITY_FILE")){
			this.TF_COOPERATIVITY_FILE.value = value;
			if(!label.isEmpty()){this.TF_COOPERATIVITY_FILE.label = label;}
			if(!description.isEmpty()){this.TF_COOPERATIVITY_FILE.description = description;}
			if(!category.isEmpty()){this.TF_COOPERATIVITY_FILE.category = category;}
			found = true;
		}  else if(name.equals("TS_FILE")){
			this.TS_FILE.value = value;
			if(!label.isEmpty()){this.TS_FILE.label = label;}
			if(!description.isEmpty()){this.TS_FILE.description = description;}
			if(!category.isEmpty()){this.TS_FILE.category = category;}
			found = true;
		}
		//TF_RANDOM PARAMETERS
		else if(name.equals("TF_DBD_LENGTH_MIN")){
			this.TF_DBD_LENGTH_MIN.value = Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.TF_DBD_LENGTH_MIN.label = label;}
			if(!description.isEmpty()){this.TF_DBD_LENGTH_MIN.description = description;}
			if(!category.isEmpty()){this.TF_DBD_LENGTH_MIN.category = category;}
			found = true;
		} else if(name.equals("TF_DBD_LENGTH_MAX")){
			this.TF_DBD_LENGTH_MAX.value =  Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.TF_DBD_LENGTH_MAX.label = label;}
			if(!description.isEmpty()){this.TF_DBD_LENGTH_MAX.description = description;}
			if(!category.isEmpty()){this.TF_DBD_LENGTH_MAX.category = category;}
			found = true;
		} else if(name.equals("TF_SPECIES_COUNT")){
			this.TF_SPECIES_COUNT.value =  Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.TF_SPECIES_COUNT.label = label;}
			if(!description.isEmpty()){this.TF_SPECIES_COUNT.description = description;}
			if(!category.isEmpty()){this.TF_SPECIES_COUNT.category = category;}
			found = true;
		} else if(name.equals("TF_COPY_NUMBER_MIN")){
			this.TF_COPY_NUMBER_MIN.value = Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.TF_COPY_NUMBER_MIN.label = label;}
			if(!description.isEmpty()){this.TF_COPY_NUMBER_MIN.description = description;}
			if(!category.isEmpty()){this.TF_COPY_NUMBER_MIN.category = category;}
			found = true;
		} else if(name.equals("TF_COPY_NUMBER_MAX")){
			this.TF_COPY_NUMBER_MAX.value = Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.TF_COPY_NUMBER_MAX.label = label;}
			if(!description.isEmpty()){this.TF_COPY_NUMBER_MAX.description = description;}
			if(!category.isEmpty()){this.TF_COPY_NUMBER_MAX.category = category;}
			found = true;
		} else if(name.equals("TF_ES")){
			this.TF_ES.value = Utils.parseDouble(value, Constants.NONE);
			if(!label.isEmpty()){this.TF_ES.label = label;}
			if(!description.isEmpty()){this.TF_ES.description = description;}
			if(!category.isEmpty()){this.TF_ES.category = category;}
			found = true;
		} else if(name.equals("TF_SIZE_LEFT")){
			this.TF_SIZE_LEFT.value =  Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.TF_SIZE_LEFT.label = label;}
			if(!description.isEmpty()){this.TF_SIZE_LEFT.description = description;}
			if(!category.isEmpty()){this.TF_SIZE_LEFT.category = category;}
			found = true;
		} else if(name.equals("TF_SIZE_RIGHT")){
			this.TF_SIZE_RIGHT.value =  Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.TF_SIZE_RIGHT.label = label;}
			if(!description.isEmpty()){this.TF_SIZE_RIGHT.description = description;}
			if(!category.isEmpty()){this.TF_SIZE_RIGHT.category = category;}
			found = true;
		} else if(name.equals("TF_ASSOC_RATE")){
			this.TF_ASSOC_RATE.value =  Utils.parseDouble(value, Constants.NONE);
			if(!label.isEmpty()){this.TF_ASSOC_RATE.label = label;}
			if(!description.isEmpty()){this.TF_ASSOC_RATE.description = description;}
			if(!category.isEmpty()){this.TF_ASSOC_RATE.category = category;}
			found = true;
		} else if(name.equals("TF_READ_IN_BOTH_DIRECTIONS")){
			this.TF_READ_IN_BOTH_DIRECTIONS.value =  Utils.parseBoolean(value, true);
			if(!label.isEmpty()){this.TF_READ_IN_BOTH_DIRECTIONS.label = label;}
			if(!description.isEmpty()){this.TF_READ_IN_BOTH_DIRECTIONS.description = description;}
			if(!category.isEmpty()){this.TF_READ_IN_BOTH_DIRECTIONS.category = category;}
			found = true;			
		} else if(name.equals("TF_PREBOUND_PROPORTION")){
			this.TF_PREBOUND_PROPORTION.value =  Utils.parseDouble(value, Constants.NONE);
			if(!label.isEmpty()){this.TF_PREBOUND_PROPORTION.label = label;}
			if(!description.isEmpty()){this.TF_PREBOUND_PROPORTION.description = description;}
			if(!category.isEmpty()){this.TF_PREBOUND_PROPORTION.category = category;}
			found = true;
		} else if(name.equals("TF_PREBOUND_TO_HIGHEST_AFFINITY")){
			this.TF_PREBOUND_TO_HIGHEST_AFFINITY.value =  Utils.parseBoolean(value, true);
			if(!label.isEmpty()){this.TF_PREBOUND_TO_HIGHEST_AFFINITY.label = label;}
			if(!description.isEmpty()){this.TF_PREBOUND_TO_HIGHEST_AFFINITY.description = description;}
			if(!category.isEmpty()){this.TF_PREBOUND_TO_HIGHEST_AFFINITY.category = category;}
			found = true;
		} else if(name.equals("SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE")){
			this.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.value =  Utils.parseBoolean(value, true);
			if(!label.isEmpty()){this.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.label = label;}
			if(!description.isEmpty()){this.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.description = description;}
			if(!category.isEmpty()){this.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.category = category;}
			found = true;
		}

		//TF_REPRESSION PARAMETERS
		else if(name.equals("REPRESSOR")){
			this.REPRESSOR.value =  Utils.parseBoolean(value, false);
			if(!label.isEmpty()){this.REPRESSOR.label = label;}
			if(!description.isEmpty()){this.REPRESSOR.description = description;}
			if(!category.isEmpty()){this.REPRESSOR.category = category;}
			found = true;
		}
		else if(name.equals("TF_REPLENLEFT")){
			this.TF_REPLENLEFT.value =  Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.TF_REPLENLEFT.label = label;}
			if(!description.isEmpty()){this.TF_REPLENLEFT.description = description;}
			if(!category.isEmpty()){this.TF_REPLENLEFT.category = category;}
			found = true;
		}
		else if(name.equals("TF_REPLENRIGHT")){
			this.TF_REPLENRIGHT.value =  Utils.parseInteger(value, Constants.NONE);
			if(!label.isEmpty()){this.TF_REPLENRIGHT.label = label;}
			if(!description.isEmpty()){this.TF_REPLENRIGHT.description = description;}
			if(!category.isEmpty()){this.TF_REPLENRIGHT.category = category;}
			found = true;
		}
		else if(name.equals("PWM_REP_THRESHOLD")){
			this.PWM_REP_THRESHOLD.value =  Utils.parseDouble(value, Constants.NONE);
			if(!label.isEmpty()){this.PWM_REP_THRESHOLD.label = label;}
			if(!description.isEmpty()){this.PWM_REP_THRESHOLD.description = description;}
			if(!category.isEmpty()){this.PWM_REP_THRESHOLD.category = category;}
			found = true;
		}
		else if(name.equals("REPRESSION_PROBABILITY")){
			this.REPRESSION_PROBABILITY.value =  Utils.parseDouble(value, Constants.NONE);
			if(!label.isEmpty()){this.REPRESSION_PROBABILITY.label = label;}
			if(!description.isEmpty()){this.REPRESSION_PROBABILITY.description = description;}
			if(!category.isEmpty()){this.REPRESSION_PROBABILITY.category = category;}
			found = true;
		}

		//DNA PARAMETERS
		else if(name.equals("DNA_SEQUENCE_FILE")){
				this.DNA_SEQUENCE_FILE.value = value;
				if(!label.isEmpty()){this.DNA_SEQUENCE_FILE.label = label;}
				if(!description.isEmpty()){this.DNA_SEQUENCE_FILE.description = description;}
				if(!category.isEmpty()){this.DNA_SEQUENCE_FILE.category = category;}
				found = true;
		}
		//BTRACK FILE
		else if(name.equals("DNA_BTRACK_FILE")){
				this.DNA_BTRACK_FILE.value = value;
				if(!label.isEmpty()){this.DNA_BTRACK_FILE.label = label;}
				if(!description.isEmpty()){this.DNA_BTRACK_FILE.description = description;}
				if(!category.isEmpty()){this.DNA_BTRACK_FILE.category = category;}
				found = true;
		}
		//DNA_RANDOM PARAMETERS
		else if(name.equals("DNA_LENGTH")){
				this.DNA_LENGTH.value = Utils.parseInteger(value, Constants.NONE);
				if(!label.isEmpty()){this.DNA_LENGTH.label = label;}
				if(!description.isEmpty()){this.DNA_LENGTH.description = description;}
				if(!category.isEmpty()){this.DNA_LENGTH.category = category;}
				found = true;
			} else if(name.equals("DNA_PROPORTION_OF_A")){
				this.DNA_PROPORTION_OF_A.value =   Utils.parseDouble(value, Constants.NONE);
				if(!label.isEmpty()){this.DNA_PROPORTION_OF_A.label = label;}
				if(!description.isEmpty()){this.DNA_PROPORTION_OF_A.description = description;}
				if(!category.isEmpty()){this.DNA_PROPORTION_OF_A.category = category;}
				found = true;
			} else if(name.equals("DNA_PROPORTION_OF_T")){
				this.DNA_PROPORTION_OF_T.value =  Utils.parseDouble(value, Constants.NONE);
				if(!label.isEmpty()){this.DNA_PROPORTION_OF_T.label = label;}
				if(!description.isEmpty()){this.DNA_PROPORTION_OF_T.description = description;}
				if(!category.isEmpty()){this.DNA_PROPORTION_OF_T.category = category;}
				found = true;
			} else if(name.equals("DNA_PROPORTION_OF_C")){
				this.DNA_PROPORTION_OF_C.value = Utils.parseDouble(value, Constants.NONE);
				if(!label.isEmpty()){this.DNA_PROPORTION_OF_C.label = label;}
				if(!description.isEmpty()){this.DNA_PROPORTION_OF_C.description = description;}
				if(!category.isEmpty()){this.DNA_PROPORTION_OF_C.category = category;}
				found = true;
			} else if(name.equals("DNA_PROPORTION_OF_G")){
				this.DNA_PROPORTION_OF_G.value = Utils.parseDouble(value, Constants.NONE);
				if(!label.isEmpty()){this.DNA_PROPORTION_OF_G.label = label;}
				if(!description.isEmpty()){this.DNA_PROPORTION_OF_G.description = description;}
				if(!category.isEmpty()){this.DNA_PROPORTION_OF_G.category = category;}
				found = true;
			} else if(name.equals("DNA_BOUNDARY_CONDITION")){
				this.DNA_BOUNDARY_CONDITION.value = value;
				if(!label.isEmpty()){this.DNA_BOUNDARY_CONDITION.label = label;}
				if(!description.isEmpty()){this.DNA_BOUNDARY_CONDITION.description = description;}
				if(!category.isEmpty()){this.DNA_BOUNDARY_CONDITION.category = category;}
				found = true;			
		} //TF RANDOM WALK PARAMETERS
			else if(name.equals("TF_IS_IMMOBILE")){
				this.TF_IS_IMMOBILE.value =  Utils.parseBoolean(value, true);
				if(!label.isEmpty()){this.TF_IS_IMMOBILE.label = label;}
				if(!description.isEmpty()){this.TF_IS_IMMOBILE.description = description;}
				if(!category.isEmpty()){this.TF_IS_IMMOBILE.category = category;}
				found = true;
			} else if(name.equals("TF_UNBINDING_PROBABILITY")){
				this.TF_UNBINDING_PROBABILITY.value =  Utils.parseDouble(value, Constants.NONE);
				if(!label.isEmpty()){this.TF_UNBINDING_PROBABILITY.label = label;}
				if(!description.isEmpty()){this.TF_UNBINDING_PROBABILITY.description = description;}
				if(!category.isEmpty()){this.TF_UNBINDING_PROBABILITY.category = category;}
				found = true;
			} else if(name.equals("TF_SLIDE_LEFT_PROBABILITY")){
				this.TF_SLIDE_LEFT_PROBABILITY.value =   Utils.parseDouble(value, Constants.NONE);
				if(!label.isEmpty()){this.TF_SLIDE_LEFT_PROBABILITY.label = label;}
				if(!description.isEmpty()){this.TF_SLIDE_LEFT_PROBABILITY.description = description;}
				if(!category.isEmpty()){this.TF_SLIDE_LEFT_PROBABILITY.category = category;}
				found = true;
			} else if(name.equals("TF_SLIDE_RIGHT_PROBABILITY")){
				this.TF_SLIDE_RIGHT_PROBABILITY.value =  Utils.parseDouble(value, Constants.NONE);
				if(!label.isEmpty()){this.TF_SLIDE_RIGHT_PROBABILITY.label = label;}
				if(!description.isEmpty()){this.TF_SLIDE_RIGHT_PROBABILITY.description = description;}
				if(!category.isEmpty()){this.TF_SLIDE_RIGHT_PROBABILITY.category = category;}
				found = true;
			} else if(name.equals("TF_JUMPING_PROBABILITY")){
				this.TF_JUMPING_PROBABILITY.value = Utils.parseDouble(value, Constants.NONE);
				if(!label.isEmpty()){this.TF_JUMPING_PROBABILITY.label = label;}
				if(!description.isEmpty()){this.TF_JUMPING_PROBABILITY.description = description;}
				if(!category.isEmpty()){this.TF_JUMPING_PROBABILITY.category = category;}
				found = true;
			} else if(name.equals("TF_HOP_STD_DISPLACEMENT")){
				this.TF_HOP_STD_DISPLACEMENT.value = Utils.parseDouble(value, Constants.NONE);
				if(!label.isEmpty()){this.TF_HOP_STD_DISPLACEMENT.label = label;}
				if(!description.isEmpty()){this.TF_HOP_STD_DISPLACEMENT.description = description;}
				if(!category.isEmpty()){this.TF_HOP_STD_DISPLACEMENT.category = category;}
				found = true;
			} else if(name.equals("TF_SPECIFIC_WAITING_TIME")){
				this.TF_SPECIFIC_WAITING_TIME.value = Utils.parseDouble(value, Constants.NONE);
				if(!label.isEmpty()){this.TF_SPECIFIC_WAITING_TIME.label = label;}
				if(!description.isEmpty()){this.TF_SPECIFIC_WAITING_TIME.description = description;}
				if(!category.isEmpty()){this.TF_SPECIFIC_WAITING_TIME.category = category;}
				found = true;
			} else if(name.equals("TF_STEP_LEFT_SIZE")){
				this.TF_STEP_LEFT_SIZE.value = Utils.parseInteger(value, Constants.NONE);
				if(!label.isEmpty()){this.TF_STEP_LEFT_SIZE.label = label;}
				if(!description.isEmpty()){this.TF_STEP_LEFT_SIZE.description = description;}
				if(!category.isEmpty()){this.TF_STEP_LEFT_SIZE.category = category;}
				found = true;
			} else if(name.equals("TF_STEP_RIGHT_SIZE")){
				this.TF_STEP_RIGHT_SIZE.value = Utils.parseInteger(value, Constants.NONE);
				if(!label.isEmpty()){this.TF_STEP_RIGHT_SIZE.label = label;}
				if(!description.isEmpty()){this.TF_STEP_RIGHT_SIZE.description = description;}
				if(!category.isEmpty()){this.TF_STEP_RIGHT_SIZE.category = category;}
				found = true;	
			} else if(name.equals("TF_UNCORRELATED_DISPLACEMENT_SIZE")){
				this.TF_UNCORRELATED_DISPLACEMENT_SIZE.value = Utils.parseInteger(value, Constants.NONE);
				if(!label.isEmpty()){this.TF_UNCORRELATED_DISPLACEMENT_SIZE.label = label;}
				if(!description.isEmpty()){this.TF_UNCORRELATED_DISPLACEMENT_SIZE.description = description;}
				if(!category.isEmpty()){this.TF_UNCORRELATED_DISPLACEMENT_SIZE.category = category;}
				found = true;
			} else if(name.equals("TF_STALLS_IF_BLOCKED")){
				this.TF_STALLS_IF_BLOCKED.value = Utils.parseBoolean(value, false);
				if(!label.isEmpty()){this.TF_STALLS_IF_BLOCKED.label = label;}
				if(!description.isEmpty()){this.TF_STALLS_IF_BLOCKED.description = description;}
				if(!category.isEmpty()){this.TF_STALLS_IF_BLOCKED.category = category;}
				found = true;
			} else if(name.equals("TF_COLLISION_UNBIND_PROBABILITY")){
				this.TF_COLLISION_UNBIND_PROBABILITY.value = Utils.parseDouble(value, Constants.NONE);
				if(!label.isEmpty()){this.TF_COLLISION_UNBIND_PROBABILITY.label = label;}
				if(!description.isEmpty()){this.TF_COLLISION_UNBIND_PROBABILITY.description = description;}
				if(!category.isEmpty()){this.TF_COLLISION_UNBIND_PROBABILITY.category = category;}
				found = true;
			} else if(name.equals("TF_AFFINITY_LANDSCAPE_ROUGHNESS")){
				this.TF_AFFINITY_LANDSCAPE_ROUGHNESS.value = Utils.parseDouble(value, Constants.NONE);
				if(!label.isEmpty()){this.TF_AFFINITY_LANDSCAPE_ROUGHNESS.label = label;}
				if(!description.isEmpty()){this.TF_AFFINITY_LANDSCAPE_ROUGHNESS.description = description;}
				if(!category.isEmpty()){this.TF_AFFINITY_LANDSCAPE_ROUGHNESS.category = category;}
				found = true;		
			} else if(name.equals("CHECK_OCCUPANCY_ON_BINDING")){
				this.CHECK_OCCUPANCY_ON_BINDING.value = Utils.parseBoolean(value, false);
				if(!label.isEmpty()){this.CHECK_OCCUPANCY_ON_BINDING.label = label;}
				if(!description.isEmpty()){this.CHECK_OCCUPANCY_ON_BINDING.description = description;}
				if(!category.isEmpty()){this.CHECK_OCCUPANCY_ON_BINDING.category = category;}
				found = true;
			} else if(name.equals("CHECK_OCCUPANCY_ON_SLIDING")){
				this.CHECK_OCCUPANCY_ON_SLIDING.value = Utils.parseBoolean(value, false);
				if(!label.isEmpty()){this.CHECK_OCCUPANCY_ON_SLIDING.label = label;}
				if(!description.isEmpty()){this.CHECK_OCCUPANCY_ON_SLIDING.description = description;}
				if(!category.isEmpty()){this.CHECK_OCCUPANCY_ON_SLIDING.category = category;}
				found = true;
			} else if(name.equals("CHECK_OCCUPANCY_ON_REBINDING")){
				this.CHECK_OCCUPANCY_ON_REBINDING.value = Utils.parseBoolean(value, false);
				if(!label.isEmpty()){this.CHECK_OCCUPANCY_ON_REBINDING.label = label;}
				if(!description.isEmpty()){this.CHECK_OCCUPANCY_ON_REBINDING.description = description;}
				if(!category.isEmpty()){this.CHECK_OCCUPANCY_ON_REBINDING.category = category;}
				found = true;
			} else if(name.equals("IS_BIASED_RANDOM_WALK")){
				this.IS_BIASED_RANDOM_WALK.value = Utils.parseBoolean(value, false);
				if(!label.isEmpty()){this.IS_BIASED_RANDOM_WALK.label = label;}
				if(!description.isEmpty()){this.IS_BIASED_RANDOM_WALK.description = description;}
				if(!category.isEmpty()){this.IS_BIASED_RANDOM_WALK.category = category;}
				found = true;		
			} else if(name.equals("IS_TWO_STATE_RANDOM_WALK")){
				this.IS_TWO_STATE_RANDOM_WALK.value = Utils.parseBoolean(value, false);
				if(!label.isEmpty()){this.IS_TWO_STATE_RANDOM_WALK.label = label;}
				if(!description.isEmpty()){this.IS_TWO_STATE_RANDOM_WALK.description = description;}
				if(!category.isEmpty()){this.IS_TWO_STATE_RANDOM_WALK.category = category;}
				found = true;		
			} else{
			System.out.println(name+" not recognised.");
		}		
		



		return found;
	}
}
