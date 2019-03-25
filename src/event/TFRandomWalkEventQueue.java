package event;

import java.io.Serializable;

import environment.Cell;

/**
 * interface for random walk event classes
 * @author created: n.r.zabet@gen.cam.ac.uk, modified: dmitrav@inbox.ru
 *
 */
public abstract class TFRandomWalkEventQueue  implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -7798070719703888097L;

	public abstract void add(ProteinEvent pe);
	
	public abstract ProteinEvent peek();
	
	public abstract ProteinEvent pop();
	
	
	public abstract int  size();
	public abstract boolean isEmpty();
	
	public abstract void scheduleNextTFRandomWalkEvent(Cell n, int moleculeID, double time);

	
	public abstract void updateNextTFRandomWalkEvent(Cell n, int moleculeID, double time);

	public abstract double getNextTFRandomWalkEventTime(int moleculeID);

	
}
