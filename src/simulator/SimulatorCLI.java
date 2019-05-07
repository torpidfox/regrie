package simulator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;

import objects.InputParameters;
import objects.TFspecies;
import objects.TargetSitesGroup;
import utils.BtrackFileParser;
import utils.Constants;
import utils.TFfileParser;
import utils.Utils;

import environment.Cell;

/**
 * runnable class from command line
 *
 * @author created: n.r.zabet@gen.cam.ac.uk, modified: dmitrav@inbox.ru
 */
public class SimulatorCLI {


    public static void oldMain(String[] args) throws IOException {
        boolean stopAfterBackup = false;
        boolean wasSaved = false;

        Cell cell = null;

        String parametersFilename = "";


        int steps = 10, ensembleSteps;
        double backupAfter;

        if (args.length > 0) {
            parametersFilename = args[0];
        }

        if (args.length > 1) {
            steps = Utils.parseInteger(args[1], 1);
        }

        backupAfter = Constants.MIN_INTERMEDIARY_STATE;

        if (args.length > 2) {
            backupAfter = Utils.parseDouble(args[2], backupAfter);
        }

        if (args.length > 3) {
            stopAfterBackup = Utils.parseBoolean(args[3], stopAfterBackup);
        }


        boolean backupRestarted = false, canRestore;

        int start = 0;
        double time = 0;
        double elapsedTime = 0, estimatedTime;
        int ensemble = 0;
        String lastIntermediaryFilename = "", backupFile, tmpBackupFile;


        //restart a backup file
        File dir = new File(System.getProperty("user.dir"));
        String[] files = dir.list();

        tmpBackupFile = "";
        if (files != null && files.length > 0) {

            int tmp_ensemble = 0, tmp_start = 0, tmp_steps = 0;
            for (int i = 0; i < files.length; i++) {
                //System.out.println("file "+files[i]);
                //located the first backup file
                if (files[i].startsWith("backup_") && files[i].endsWith(Constants.INFO_FILE_EXTENSION)) {
                    ArrayList<String> buffer = new ArrayList<String>();
                    String thisLine = "";
                    BufferedReader infoFile = null;
                    canRestore = false;
                    try {
                        infoFile = new BufferedReader(new FileReader(files[i]));
                        while ((thisLine = infoFile.readLine()) != null) {
                            buffer.add(thisLine);
                        }

                        if (buffer.size() >= 4) {
                            //Construct the BufferedWriter object
                            backupFile = buffer.get(0);
                            canRestore = true;
                            if (backupFile != null && !backupFile.isEmpty()) {
                                File f = new File(backupFile);
                                if (!f.exists()) {
                                    f = new File("old_" + backupFile);
                                    if (f.exists()) {
                                        backupFile = "old_" + backupFile;
                                    } else {
                                        canRestore = false;
                                        System.out.println("cannot restore file " + "old_" + backupFile);
                                    }
                                }
                            } else {
                                canRestore = false;
                                System.out.println("cannot restore file " + files[i]);
                            }

                            if (canRestore) {
                                tmp_ensemble = Utils.parseInteger(buffer.get(1), 0);
                            }

                            if (canRestore) {
                                tmp_start = Utils.parseInteger(buffer.get(2), 0);
                            }

                            if (canRestore) {
                                tmp_steps = Utils.parseInteger(buffer.get(3), steps);
                                //System.out.println("buffer.get(3) is " + buffer.get(3));
                            }

                            if (canRestore) {
                                tmpBackupFile = backupFile;
                            }

                            if (canRestore && (tmp_start) > (start)) {
                                tmpBackupFile = backupFile;
                                lastIntermediaryFilename = backupFile;
                                ensemble = tmp_ensemble;
                                start = tmp_start;
                                steps = tmp_steps;
                            } else {
                                System.out.println("did not restore file " + backupFile);
                            }

                        }
                    } catch (FileNotFoundException ex) {
                        ex.printStackTrace();
                    } catch (IOException ex) {
                        ex.printStackTrace();
                    }
                }
            }
        }

        if (tmpBackupFile != null && !tmpBackupFile.isEmpty()) {
            FileInputStream fis = null;
            ObjectInputStream in = null;
            try {
                fis = new FileInputStream(tmpBackupFile);
                in = new ObjectInputStream(fis);
                cell = (Cell) in.readObject();
                backupRestarted = true;
                time = cell.cellTime;
                cell.ensemble = ensemble;

                int actualStart = (int) Math.round(cell.cellTime / (cell.totalStopTime / steps));
                System.out.println("start=" + start + "; actualStart = " + actualStart + "");

                if (actualStart != start) {
                    start = actualStart;
                }

                in.close();
                System.out.println("restore file " + tmpBackupFile);
                backupRestarted = true;
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        if (!backupRestarted) {
            cell = new Cell(parametersFilename, null, true);
        }

        double stopTime = cell.getTotalStopTime();
        double timeStep = stopTime / steps;
        int i = start;

        double lastSave = elapsedTime;

        System.out.println("reGRiE v1.2");

        ensembleSteps = steps * cell.ip.ENSEMBLE_SIZE.value;

        String intermediaryFilename = cell.outputIntermediaryBackupFile;

        //debug info
		/*
		for (TFspecies tf:cell.TFspecies) {
			System.out.println(tf.name + " " + tf.repressor);
		}
		System.out.println(cell.ip.TF_FILE.value);
		*/

        //run intervals
        if (cell.totalStopTime > 0) {
            while (i < ensembleSteps && (!stopAfterBackup || !wasSaved)) {
                //System.out.println("stopAfterBackup="+stopAfterBackup +"; wasSaved="+wasSaved);
                ensemble = cell.ensemble;
                time += timeStep;
                System.out.println("to perform step " + (i % steps) + " of set " + (steps));
                if (i % steps == (steps - 1) || i == ensembleSteps - 1) {
                    time = stopTime;
                    System.out.println("last step of the set to finish");
                }

                elapsedTime += cell.runInterval(time, elapsedTime);
                estimatedTime = steps * elapsedTime / (i + 1) * cell.ip.ENSEMBLE_SIZE.value;
                System.out.println((i + 1) + "/" + ensembleSteps + " elapsed time: " + elapsedTime + "s; estimated " +
						"total time:" + estimatedTime + "s");
                i++;
                if (ensemble != cell.ensemble) {
                    time = 0;
                    System.out.println("set finished");
                    //added
                    cell.resetOutputDir("set" + (i / steps));
                    cell.initialiseInternalParameters();

                    for (TargetSitesGroup group : cell.tsg.tsg) {
                        group.isAvailable = group.isAvailable(cell);
                    }
                }

                //if it has been more than 1h from last save then save the state
                if (backupAfter > 0 && elapsedTime - lastSave > backupAfter) {
                    System.out.println("backup started");
                    double currentTime = System.currentTimeMillis();

                    FileOutputStream fos = null;
                    ObjectOutputStream out = null;
                    try {

                        intermediaryFilename = cell.outputIntermediaryBackupFile;

                        fos = new FileOutputStream("tmp_" + intermediaryFilename);
                        out = new ObjectOutputStream(fos);
                        out.writeObject(cell);
                        out.close();

                        if (!lastIntermediaryFilename.isEmpty()) {
                            File f = new File(lastIntermediaryFilename);
                            if (f.exists() && f.canWrite()) {
                                f.renameTo(new File("old_" + lastIntermediaryFilename));
                            }
                        }

                        // Rename file (or directory)
                        File file = new File("tmp_" + intermediaryFilename);
                        file.renameTo(new File(intermediaryFilename));
                        System.out.println("backup file created " + intermediaryFilename + "; ensemble: " + ensemble + "; step:" + i + "; total:" + steps);


                        BufferedWriter infoFile = null;
                        try {
                            //Construct the BufferedWriter object
                            infoFile = new BufferedWriter(new FileWriter(cell.outputIntermediaryInfoFile));
                            infoFile.write(intermediaryFilename);
                            infoFile.newLine();
                            infoFile.write(ensemble + "");
                            infoFile.newLine();
                            infoFile.write(i + "");
                            infoFile.newLine();
                            infoFile.write(steps + "");
                            infoFile.newLine();
                            infoFile.flush();
                            infoFile.close();
                        } catch (FileNotFoundException ex) {
                            ex.printStackTrace();
                        } catch (IOException ex) {
                            ex.printStackTrace();
                        }

                        if (!lastIntermediaryFilename.isEmpty()) {
                            File f = new File("old_" + lastIntermediaryFilename);
                            if (f.exists() && f.canWrite()) {
                                f.delete();
                            }
                        }

                        lastIntermediaryFilename = intermediaryFilename;
                        wasSaved = true;
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                    lastSave = elapsedTime + cell.computeElapsedTime(currentTime);
                }

            }
            //if simulations reached the end the backup files can be deleted
            if (i >= ensembleSteps) {
                if (!lastIntermediaryFilename.isEmpty()) {
                    File f = new File(lastIntermediaryFilename);
                    if (f.exists() && f.canWrite()) {
                        f.delete();
                    }
                }
                File f = new File(cell.outputIntermediaryInfoFile);
                if (f.exists() && f.canWrite()) {
                    f.delete();
                }
            }
        } else {
            cell.runUntilTSReached();
        }
    }

    public static void threadableMain(String params) {
        int ensemble_count = 1;
        InputParameters ip = new InputParameters(params);
        Thread[] threads = new Thread[ensemble_count];

        for (int i = 0; i < ensemble_count; i++) {
            threads[i] = new Thread(new EnsembleRunnable(i, ip));
            threads[i].start();
        }
    }
    /**
     * @param args
     * @throws FileNotFoundException
     */
    public static void main(String[] args) throws IOException {
        oldMain(args);
    }
}
