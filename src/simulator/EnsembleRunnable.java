package simulator;

import environment.Cell;
import objects.InputParameters;

import java.io.*;

public class EnsembleRunnable implements Runnable {
    private int ensemble_index;
    private int steps = 10;
    public Cell cell;

    EnsembleRunnable(int ensemble_index, InputParameters ip) {
        this.ensemble_index = ensemble_index;


        try {
            this.cell = new Cell(ip, null, true);
            this.cell.resetOutputDir("set" + ensemble_index);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

    }

    public void run() {
        double stopTime = cell.getTotalStopTime();
        double timeStep = stopTime / steps;
        int i = 0;
        double time = 0;
        double elapsedTime = 0;

        System.out.println("thread " + this.ensemble_index);

        //run intervals
        if (cell.totalStopTime > 0) {
            while (i < steps) {
                //System.out.println("stopAfterBackup="+stopAfterBackup +"; wasSaved="+wasSaved);
                time += timeStep;
                try {
                    elapsedTime += cell.runInterval(time, elapsedTime);
                }
                catch (FileNotFoundException e) {
                    e.printStackTrace();
                }

                System.out.println("to perform step " + (i % steps) + " of set " + steps + " for ensemble " + ensemble_index);

                if (i % steps == (steps - 1)) {
                    time = stopTime;
                    System.out.println("last step of the set to finish for thread " + this.ensemble_index);
                }

                i++;
            }

        } else {
            try {
                cell.runUntilTSReached();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }

        System.out.println("thread " + ensemble_index + " finished");
    }

}
