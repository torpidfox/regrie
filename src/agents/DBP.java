package agents;

import java.io.Serializable;

import environment.Cell;
import event.ProteinEvent;
import utils.Constants;
import utils.Utils;

/**
 * this class specifies all DNA Binding Proteins
 * @author created: n.r.zabet@gen.cam.ac.uk, modified: dmitrav@inbox.ru
 *
 */
public abstract class DBP  implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 121399550800359214L;
	public int ID;  // the ID of the TF
	protected int 	position; // the position on the DNA  (-1 free)
	protected int lastPosition; // the last position of this TF
	protected double timeOfLastPositionChange; // the last time the position was changed
	public int size;
	public int state;
	public int speciesID; //the ID of the species to which this TF belongs
	public boolean hasDNAbasedCooperativity;
	public boolean hasDirectCooperativity;
	public boolean repressionState; //the state of TF repressor ()
	
	public ProteinEvent pe;
	
	protected int direction;// the direction on which it is bound 0 for 5' -> 3' or 1 for 3' -> 5' 
	
	protected double timeBound;
	protected boolean wasBound;

	public double moveRate;
	
	public int leftNeighbour;
	public int rightNeighbour;

	public int stickToLeft;
	public int stickToRight;

	public int leftMostPosition;
	public int rightMostPosition;
	public int observedLeftMostPosition;
	public int observedRightMostPosition;
	
	public int slidingEvents;
	

	public int boundToDNA;
	
	
	/**
	 * empty class constructor
	 */
	public DBP(){
		ID = Constants.NONE;
		position = Constants.NONE;
		lastPosition = Constants.NONE;
		timeOfLastPositionChange = 0;
		size=0;
		timeBound = 0;
		moveRate =Constants.NONE;
		repressionState = false;
		
		leftNeighbour=Constants.NONE;
		rightNeighbour=Constants.NONE;	
		pe=null;
		stickToLeft = Constants.NONE;
		stickToRight = Constants.NONE;
		
		leftMostPosition = Constants.NONE;
		rightMostPosition = Constants.NONE;
		observedLeftMostPosition = Constants.NONE;
		observedRightMostPosition = Constants.NONE;
		slidingEvents = 0;
		wasBound = false;
		boundToDNA = Constants.NONE;
	}
	
	/**
	 * class constructor
	 * @param ID the ID of the TF
	 * @param position  the position on the DNA  (-1 free)
	 * @param lastPosition  the last position of this TF
	 * @param timeOfLastPositionChange the last time the position was changed
	 */
	public DBP(int ID, int position, int lastPosition, int timeOfLastPositionChange, int size){
		this.ID = ID;
		this.position = position;
		this.lastPosition = lastPosition;
		this.timeOfLastPositionChange = timeOfLastPositionChange;
		this.size = size;
		timeBound = 0;
		moveRate =Constants.NONE;
		repressionState = false;
		
		leftNeighbour=Constants.NONE;
		rightNeighbour=Constants.NONE;
		pe = null;
		
		stickToLeft = Constants.NONE;
		stickToRight = Constants.NONE;
		
		
		leftMostPosition = Constants.NONE;
		rightMostPosition = Constants.NONE;
		observedLeftMostPosition = Constants.NONE;
		observedRightMostPosition = Constants.NONE;
		wasBound = false;
		boundToDNA = Constants.NONE;
	}
	
	public abstract void act(Cell n, ProteinEvent pe);
	
	public abstract void updateBoundTime(Cell n, double timeBound, int direction, int position);

	/**
	 * returns current position
	 * @return
	 */
	public int getPosition(){
		return this.position;
	}

	public abstract void setPosition( int newPosition, double timeOfLastPositionChange);
	
	public abstract void clearPosition(double timeOfLastPositionChange);

	
	/**
	 * sets the sliding extremes
	 */
	public void initSlidingExtremes(){
		this.leftMostPosition = position;
		this.rightMostPosition = position;
		this.slidingEvents = 0;
	}
	
	
	/**
	 * sets the sliding extremes when sliding to left
	 */
	public void setLeftSlidingExtreme(){
		if(position < this.leftMostPosition){  
			this.leftMostPosition = position;
		}
		slidingEvents++;
	}
		
	/**
	 * sets the sliding extremes when sliding to right
	 */
	public void setRightSlidingExtreme(){
		if(position > this.rightMostPosition){  
			this.rightMostPosition = position;
		}
		slidingEvents++;
	}
	
	/**
	 * increments sliding events
	 */
	public void incrementSlidingEvents(){
		slidingEvents++;
	}
	
	/**
	 * returns the current slided length
	 * @return
	 */
	public int getSlidingLength(){
		return this.rightMostPosition-this.leftMostPosition;
	}
	
	/**
	 * returns the sliding events
	 * @return
	 */
	public int getSlidingEvents(){
		return slidingEvents;
	}
	
	/**
	 * sets the sliding extremes
	 */
	public void initObservedSlidingExtremes(){
		this.observedLeftMostPosition = position;
		this.observedRightMostPosition = position;
	}
	
	/**
	 * sets the sliding extremes
	 */
	public void initObservedSlidingExtremes(int left, int right){
		this.observedLeftMostPosition = left;
		this.observedRightMostPosition = right;
	}
	
	/**
	 * sets the sliding extremes when sliding to left
	 */
	public void setObservedLeftSlidingExtreme(){
		if(position < this.observedLeftMostPosition){  
			this.observedLeftMostPosition = position;
		}
	}
		
	/**
	 * sets the sliding extremes when sliding to right
	 */
	public void setObservedRightSlidingExtreme(){
		if(position > this.observedRightMostPosition){  
			this.observedRightMostPosition = position;
		}
	}

	/**
	 * returns the observed slided length (with hopping)
	 * @return
	 */
	public int getObservedSlidingLength(){
		return  this.observedRightMostPosition-this.observedLeftMostPosition;
	}
	
	/**
	 * returns current direction
	 * @return
	 */
	public int getDirection(){
		return this.direction;
	}

	/**
	 * sets the direction of the molecule
	 * @param direction
	 */
	public void setDirection(int direction){
		this.direction = direction;
	}

	
	/**
	 * returns current move Rate
	 * @return
	 */
	public double getMoveRate(){
		return this.moveRate;
	}

	/**
	 * sets the move rate of the molecule
	 * @param direction
	 */
	public void setMoveRate(Cell n){
		moveRate = n.dna.TFavgMoveRate[speciesID][position][direction];
	}
	
	
	/**
	 * returns the last occupied position
	 * @return
	 */
	public int getLastPosition(){
		return lastPosition;
	}

	/**
	 * AD
	 * Done for the purpose of setting position -1 to the TF which cannot find a position to bind
	 * @param position
	 */
	public void setLastPosition(int position){ lastPosition = position; }

	/**
	 * gets the time of last position change
	 * @param timeOfLastPositionChange
	 */
	public double getTimeOfLastPositionChange(){
		return this.timeOfLastPositionChange;
	}

	
	/**
	 * updates the time of last position change
	 * @param timeOfLastPositionChange
	 */
	public void setTimeOfLastPositionChange(double timeOfLastPositionChange){
		this.timeOfLastPositionChange = timeOfLastPositionChange;
	}

	/**
	 * returns the time this molecule was bound
	 * @return
	 */
	public double getTimeBound(){
		return this.timeBound;
	}
	
	public abstract void setCooperativityArea(Cell n, double time);

	public abstract void resetCooperativityArea(Cell n, double time);

	
	
	public abstract void setCooperativityLeft(Cell n, double time);
	public abstract void resetCooperativityLeft(Cell n, double time);
	public abstract void setCooperativityRight(Cell n, double time);
	public abstract void resetCooperativityRight(Cell n, double time);
	
	
	
	
	
	/**
	 * this methods sets the neighbours of the current molecule once it bound to the DNA
	 * @param n pointer to the environment
	 */
	public void setNeighboursOnBinding(Cell n){
		
		//set left neighbour	
		this.leftNeighbour = n.dna.getLeftNeighbour(this.position);
		if(this.leftNeighbour!=Constants.NONE ){
			n.dbp[this.leftNeighbour].rightNeighbour=this.ID;
		}
		
		
		//set right neighbour
		this.rightNeighbour = n.dna.getRightNeighbour(this.position+this.size-1);
		if(this.rightNeighbour!=Constants.NONE){
			n.dbp[this.rightNeighbour].leftNeighbour = this.ID;
		}
		
	}
	
	
	
	

	/**
	 * resets TF neighbours when current TF unbinds
	 * @param moleculeID
	 */
	public void setNeighboursOnUnbinding(Cell n){
		//reset left neighbour
		if(this.leftNeighbour!=Constants.NONE){
			n.dbp[this.leftNeighbour].rightNeighbour=Constants.NONE;
			this.leftNeighbour=Constants.NONE;
		}
		
		//reset right neighbour
		if(this.rightNeighbour!=Constants.NONE){
			n.dbp[this.rightNeighbour].leftNeighbour=Constants.NONE;
			this.rightNeighbour=Constants.NONE;

		}		
	}
	
	
	/**
	 * resets TF neighbours on slide left
	 * @param moleculeID
	 */
	public void setNeighboursOnSlideLeft(Cell n){
		//set left neiBoundProteinghbour
		this.leftNeighbour = n.dna.getLeftNeighbour(this.position);
		if(this.leftNeighbour!=Constants.NONE ){
			n.dbp[this.leftNeighbour].rightNeighbour=this.ID;
		}
		
		//reset right neighbour
		if(this.rightNeighbour!=Constants.NONE){
			n.dbp[this.rightNeighbour].leftNeighbour=Constants.NONE;
			this.rightNeighbour=Constants.NONE;

		}		
	}
	
	/**
	 * resets TF neighbours on slide right
	 * @param moleculeID
	 */
	public void setNeighboursOnSlideRight(Cell n){
		//reset left neighbour
		if(this.leftNeighbour!=Constants.NONE){
			n.dbp[this.leftNeighbour].rightNeighbour=Constants.NONE;
			this.leftNeighbour=Constants.NONE;
		}
		
		//set right neighbour
		this.rightNeighbour = n.dna.getRightNeighbour(this.position+this.size-1);
		if(this.rightNeighbour!=Constants.NONE){
			n.dbp[this.rightNeighbour].leftNeighbour = this.ID;
		}
	}
	
	
	public abstract boolean unbindMolecule(Cell n, double time);
	public abstract int bindMolecule(Cell n, double time, int newPosition);
	public abstract int bindMolecule(Cell n, double time, int newPosition, int direction);
	public abstract int hopMolecule(Cell n, double time, int newPosition);
	public abstract int slideLeftMolecule(Cell n, double time, int newPosition, boolean isHopEvent);
	public abstract int slideRightMolecule(Cell n, double time, int newPosition, boolean isHopEvent);

	/**
	 * changes the current direction
	 * @param n
	 */
	public void changeDirection(Cell n, double time){
		if(n.TFreadingDirection>1){
			int newDirection = Utils.generateNextInteger(n.randomGenerator, 0, n.TFreadingDirection);
			if(newDirection!=direction){
				setDirection(newDirection);//draw a random direction
				setMoveRate(n); //get new move rate				
				//n.eventQueue.TFRandomWalkEventQueue.updateNextTFRandomWalkEvent(n, this.ID, time); //updates the waiting time
			}
		}
	} 
	
	/**
	 * returns true if this molecule previously bound to the DNA
	 * @return
	 */
	public boolean wasBound(){
		return wasBound;
	}
	
}
