package utils;

import objects.DNA;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * class to parse experimentally revealed regions of open chromatin
 * Created by AD on 14.05.16.
 */
public class BtrackFileParser {

    /**
     *
     * @param filename is a path to the file containing info of open chromatin in the considerable DNA region!
     * @param dna
     * @throws IOException
     */
    public static void fileParser(String filename, DNA dna) throws IOException {

        ArrayList<Byte> bufferChromatin = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String buffer;

        while ((buffer = reader.readLine()) != null){
            bufferChromatin.add(Byte.parseByte(buffer));
        }
        reader.close();


        if (bufferChromatin.size() == dna.strand.length){
            dna.chromatin = new boolean[bufferChromatin.size()];
            for (int i = 0; i < bufferChromatin.size(); i++) {
                dna.chromatin[i] = bufferChromatin.get(i) == 1;
            }
        }
        else {
            System.out.println("Error: dna.chromatin and dna.strand.length have different lengths!");
        }

    }



}
