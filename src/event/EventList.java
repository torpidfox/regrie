package event;

import java.io.Serializable;

import utils.Constants;
import utils.Gillespie;
import environment.Cell;

/**
 * class that contains the event list
 *
 * @author created: n.r.zabet@gen.cam.ac.uk, modified: dmitrav@inbox.ru
 */
public class EventList implements Serializable {
    /**
     *
     */
    private static final long serialVersionUID = 8623732338279526421L;
    /**
     *
     */
    public TFBindingEventQueue TFBindingEventQueue;
    public TFRandomWalkEventQueue TFRandomWalkEventQueue;


    public EventList(Cell n) {


        TFBindingEventQueue = new TFBindingEventQueueDNAOccupancy(n);


        // TF random walk event list 1D diffusion
        if (n.ip.EVENT_LIST_USES_FR.value) {
            if (n.ip.EVENT_LIST_SUBGROUP_SIZE.value >= 0 && n.ip.EVENT_LIST_SUBGROUP_SIZE.value < n.dbp.length) {
                TFRandomWalkEventQueue = new TFRandomWalkEventQueueFRopt(n);
            } else {
                TFRandomWalkEventQueue = new TFRandomWalkEventQueueFR(n);
            }
        } else {
            TFRandomWalkEventQueue = new TFRandomWalkEventQueueDM(n);
        }
    }


    /**
     * returns the next TF binding event and removes it from the list
     *
     * @return
     */
    public ProteinEvent popNextTFBindingEvent() {
        return TFBindingEventQueue.pop();
    }


    /**
     * returns the next TF random walk event and removes it from the list
     *
     * @return
     */
    public ProteinEvent popNextTFRandomWalkEvent() {
        return TFRandomWalkEventQueue.pop();
    }


    /**
     * returns a number which encodes whether the next event is TF binding ...?
     *
     * @return
     */
    public int getNextEventType() {
        int result = Constants.NEXT_EVENT_IS_NONE;
        double nextEventTime = Double.MAX_VALUE;

        if (!TFBindingEventQueue.isEmpty() && nextEventTime > TFBindingEventQueue.peek().time) {
            nextEventTime = TFBindingEventQueue.peek().time;
            result = Constants.NEXT_EVENT_IS_TF_BINDING;
        }


        if (!TFRandomWalkEventQueue.isEmpty() && nextEventTime > TFRandomWalkEventQueue.peek().time) {
            nextEventTime = TFRandomWalkEventQueue.peek().time;
            result = Constants.NEXT_EVENT_IS_TF_RANDOM_WALK;
        }

        return result;
    }

    /**
     * returns the soones event
     *
     * @return
     */
    public Event getNextEvent() {
        Event e = null;
        int nextEventType = getNextEventType();


        switch (nextEventType) {
            case Constants.NEXT_EVENT_IS_NONE:
                break;
            case Constants.NEXT_EVENT_IS_TF_BINDING:
                e = this.popNextTFBindingEvent();
                break;
            case Constants.NEXT_EVENT_IS_TF_RANDOM_WALK:
                e = this.popNextTFRandomWalkEvent();
                break;

            default:
                e = null;
        }

        return e;
    }

    /**
     * checks whether there is any event left in the entire list
     *
     * @return
     */
    public boolean isEmpty() {
        return (TFBindingEventQueue == null || TFBindingEventQueue.isEmpty()) && (TFRandomWalkEventQueue == null || TFRandomWalkEventQueue.isEmpty());
    }


    /**
     * computes the total number of events in the lists
     *
     * @return
     */
    public int size() {
        int result = 0;


        if (TFBindingEventQueue != null && !TFBindingEventQueue.isEmpty()) {
            result++;
        }

        if (TFRandomWalkEventQueue != null) {
            result += TFRandomWalkEventQueue.size();
        }

        return result;
    }

    /**
     * schedules the next TF binding event
     */
    public void scheduleNextTFBindingEvent(Cell n, double time) {
        // if no TF then halt
        if (n.freeTFmoleculesTotal > 0 && this.TFBindingEventQueue.proteinBindingPropensitySum > 0) {

            //generate next reaction time
            double nextTime = Gillespie.computeNextReactionTime(this.TFBindingEventQueue.proteinBindingPropensitySum,
					n.randomGenerator);

            //find next reaction
            int nextTFspecies =
					Gillespie.getNextReaction(this.TFBindingEventQueue.proteinBindingPropensitySum * n.randomGenerator.nextDouble(), this.TFBindingEventQueue.proteinBindingPropensity);

            if (nextTFspecies > Constants.NONE && nextTFspecies < n.TFspecies.length) {
                int TFID = n.getFreeTFmolecule(nextTFspecies);
                if (TFID != Constants.NONE) {
                    int position = Constants.NONE;// Gillespie.getNextReaction(n.randomGenerator.nextDouble()*n.dna.effectiveTFAffinitiesSum[nextTFspecies], n.dna.effectiveTFaffinities[nextTFspecies]);
                    this.TFBindingEventQueue.add(new ProteinEvent(time + nextTime, TFID, position, true,
                            Constants.EVENT_TF_BINDING, false, false, false));
                }
            }
        }
    }


    /**
     * schedules the next TF random walk event
     */
    public void scheduleNextTFRandomWalkEvent(Cell n, int moleculeID, double time) {
        this.TFRandomWalkEventQueue.scheduleNextTFRandomWalkEvent(n, moleculeID, time);
    }
}
