package event;

import java.util.ArrayList;
import java.util.PriorityQueue;

import utils.CellUtils;
import utils.Constants;
import utils.Gillespie;
import utils.Utils;

import environment.Cell;

/**
 * random walk event class using the clustered First Reaction method
 * @author created: n.r.zabet@gen.cam.ac.uk, modified: dmitrav@inbox.ru
 *
 */
public class TFRandomWalkEventQueueFRopt  extends TFRandomWalkEventQueue{

	/**
	 * 
	 */
	private static final long serialVersionUID = -8520471978546008987L;
	public ArrayList<PriorityQueue<ProteinEvent>> randomWalkEvents;
	public PriorityQueue<ProteinEvent> randomWalkEventsMin;
	
	public int size;
	
	public int groupsNo;
	public int groupSize;
	public int startID;
	public int TFcount;
	
	

	/**
	 * class constructor. Initialkises the event list
	 */
	public TFRandomWalkEventQueueFRopt(Cell n){
		randomWalkEvents = new ArrayList<PriorityQueue<ProteinEvent>>();
	
		this.startID = 0;
		this.TFcount = n.dbp.length;
		
		
		this.groupSize = n.ip.EVENT_LIST_SUBGROUP_SIZE.value;
		if(this.groupSize <= 0){
			this.groupSize = (int) Math.floor(Math.sqrt(TFcount));
		} else if(groupSize > TFcount){
			this.groupSize = TFcount;
		}
		
		this.groupsNo  = (int)  (int) Math.ceil((double)TFcount/groupSize);
			
		for(int i=0; i< groupsNo; i++){
			randomWalkEvents.add(new PriorityQueue<ProteinEvent>());
		}
		randomWalkEventsMin = new PriorityQueue<ProteinEvent>();

		size = 0;
	}
	
	
	/**
	 * adds a new event to the list
	 * @param re the new event
	 */
	public void add(ProteinEvent newEvent){
		ProteinEvent oldEvent = null;
		int groupID= this.getGroupID(newEvent.proteinID);
		if(!randomWalkEvents.get(groupID).isEmpty()){
			oldEvent = randomWalkEvents.get(groupID).peek();
		}
		
		// add the event to the specific codon event list
		randomWalkEvents.get(groupID).add(newEvent);

		size++;
		
		// update the min list per codon if necessary
		updateElongationEventQueueBindingCodonMinList(oldEvent,randomWalkEvents.get(groupID).peek());
				
	}
	
	/**
	 * updates the list of soonest event per codon when a new event is added to the list
	 * @param oldEvent
	 * @param newEvent
	 * @return
	 */
	private boolean updateElongationEventQueueBindingCodonMinList(ProteinEvent oldEvent, ProteinEvent newEvent){
		boolean result=false;
		
		if(oldEvent!=null){
			randomWalkEventsMin.remove(oldEvent);
		}
		if(newEvent!=null){
			randomWalkEventsMin.add(newEvent);
			result = true;
		}
		
		return result;
	}
	
	
	
	/**
	 * polls the soonest event
	 * @return
	 */
	public ProteinEvent pop(){
		ProteinEvent result=null;
		if(!randomWalkEventsMin.isEmpty()){
			result = randomWalkEventsMin.poll();
			int groupID= this.getGroupID(result.proteinID);
			randomWalkEvents.get(groupID).poll();
			if(!randomWalkEvents.get(groupID).isEmpty()){
				randomWalkEventsMin.add(randomWalkEvents.get(groupID).peek());
			}
			size--;
		}
		
		
		return result;
		
	} 
	
	/**
	 * peeks the soonest event
	 * @return
	 */
	public ProteinEvent peek(){
		return randomWalkEventsMin.peek();
	} 
	

	/**
	 * returns true if the list of events is empty or false otherwise
	 */
	public boolean isEmpty(){
		return size==0;
	}
	
	/**
	 * returns the number of events in the list
	 */
	public int size(){
		return size;
	}
	
	/**
	 * gets the groupt to which a molecule is assigned
	 * @param moleculeID
	 * @return
	 */
	private int getGroupID(int moleculeID){
		return (moleculeID-startID)/groupSize;
	}
	
	
	
	/**
	 * schedules the next non-cognate TF random walk event
	 */
	public void scheduleNextTFRandomWalkEvent(Cell n, int moleculeID, double time){
		
	
		if(n.dbp[moleculeID].getPosition()!=Constants.NONE){
			
			int position = n.dbp[moleculeID].getPosition();
			int newPosition=position;
			int nextAction = Constants.NONE;
			//double avgMoveRate;
			int direction = n.dbp[moleculeID].getDirection();

			int speciesID = n.dbp[moleculeID].speciesID;
			boolean isHoppingEvent = false;
			boolean isLeftSlidingEvent = false;
			boolean isRightSlidingEvent = false;

			double affinity = direction == 1 ?
					CellUtils.computeTFAffinityLR(n.dna.strand, position, n.TFspecies[speciesID].pfm, 0, true) :
					CellUtils.computeTFAffinityRL(n.dna.strand, position, n.TFspecies[speciesID].pfm, 0, true);

			boolean isTimeZero = affinity < n.specificBindingThres;

			//compute the time the TF stays stucked 
			double nextTime = isTimeZero ? 0:
					Gillespie.computeNextReactionTime( n.dbp[moleculeID].getMoveRate(), n.randomGenerator);

			
			double randomNumber=n.randomGenerator.nextDouble()*n.TFspecies[speciesID].slideRightNo;
			if(randomNumber < n.TFspecies[speciesID].jumpNo){
				nextAction = Constants.EVENT_TF_RANDOM_WALK_JUMP;
				newPosition = Constants.NONE;
			} else if(randomNumber < n.TFspecies[speciesID].hopNo){
				isHoppingEvent = true;
				nextAction = Constants.EVENT_TF_RANDOM_WALK_HOP;
				newPosition  = Utils.generateNextNormalDistributedInteger(n.randomGenerator, position, n.TFspecies[speciesID].hopSTDdisplacement);
				
				
				if(newPosition < 0){
					newPosition = 0;
				} else if( newPosition >= n.dna.strand.length-n.dbp[moleculeID].size){
					newPosition =  n.dna.strand.length-n.dbp[moleculeID].size;
				} else if (Math.abs(newPosition-position) > n.TFspecies[speciesID].uncorrelatedDisplacementSize){
					//the hop size is to big and the hop becomes a jump 
					nextAction = Constants.EVENT_TF_RANDOM_WALK_JUMP;
					newPosition = Constants.NONE;
				} else if(newPosition > position && n.dbp[moleculeID].size > (newPosition-position)){
					// the size it too small and the molecule slides to  right
					isRightSlidingEvent = true;
					nextAction = Constants.EVENT_TF_RANDOM_WALK_SLIDE_RIGHT;
				} else if(newPosition < position && n.dbp[moleculeID].size > (position - newPosition)){
					// the size it too small and the molecule slides to  left
					isLeftSlidingEvent = true;
					nextAction = Constants.EVENT_TF_RANDOM_WALK_SLIDE_LEFT;
				}
			} else if(randomNumber < n.dna.TFSlideLeftNo[speciesID][position][direction]){
				isLeftSlidingEvent = true;
				nextAction = Constants.EVENT_TF_RANDOM_WALK_SLIDE_LEFT;
				newPosition= position-n.TFspecies[speciesID].stepLeftSize;
			} else if(randomNumber <  n.dna.TFSlideRightNo[speciesID][position][direction]){
				isRightSlidingEvent = true;
				nextAction = Constants.EVENT_TF_RANDOM_WALK_SLIDE_RIGHT;
				newPosition= position+n.TFspecies[speciesID].stepRightSize;
			}			
			
			ProteinEvent e = new ProteinEvent(time+nextTime, moleculeID, newPosition, true, nextAction, isHoppingEvent, isLeftSlidingEvent, isRightSlidingEvent);
			n.dbp[moleculeID].pe = e;
			this.add(e);
		}
	}
	
	/**
	 * updates the event of a bound molecule. Because the DM implementation always checks for the latest move rate and updates it in the list and then redraws a new event we can just schedule the event list
	 */
	public void updateNextTFRandomWalkEvent(Cell n, int moleculeID, double time){
		boolean removed = this.randomWalkEvents.remove(n.dbp[moleculeID].pe);
		if(removed){
			this.scheduleNextTFRandomWalkEvent(n, moleculeID, time);
		}
	}

	/**
	 * returns the time of the next event for a species 
	 */
	public double getNextTFRandomWalkEventTime(int moleculeID){
		
		
		return -1;
	}
	
}
