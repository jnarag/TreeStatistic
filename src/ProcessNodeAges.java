import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.RootedTree;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import org.apache.commons.math3.stat.*;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;


/**
 * Created with IntelliJ IDEA.
 * User: jayna
 * Date: 12/06/2013
 * Time: 11:34
 * To change this template use File | Settings | File Templates.
 */

public class ProcessNodeAges {

    String treeFileNamesList1;
    String treeFileNamesList2;

    public ProcessNodeAges() {

        treeFileNamesList1 = "tree1_20140909.txt";
        treeFileNamesList2 = "tree2_20140909.txt";
    }

    public void listFilesForFolder(File folder, List<File> filelist) {

        if(folder.isDirectory()) {
            for (File fileEntry : folder.listFiles()) {
                if (fileEntry.isDirectory()) {
                    listFilesForFolder(fileEntry, filelist);
                } else {
                    //System.out.println(fileEntry.getName());
                    if(fileEntry.getName().contains("tree")) {
                        filelist.add(fileEntry);
                    }
                }
            }
        }

    }

    public void getRootAges() {

        String[] psi_parts = new String[]{"psi_0.0", "psi_0.5", "psi_1.0"};

        File folder = new File("/Users/jayna/Documents/duke_reassortment_BACKUP/2015-04-15");

        try {
            List<File> tree1list = new ArrayList<File>();
            List<File> tree2list = new ArrayList<File>();

            List<File> filelist = new ArrayList<File>();

            listFilesForFolder(folder, filelist);

            FileWriter treeFile1 = new FileWriter(new File(treeFileNamesList1));
            FileWriter treeFile2 = new FileWriter(new File(treeFileNamesList2));


            for (File f : filelist) {

                if (f.getName().contains("tree_100N")) {
                    if (f.getName().contains("segment_1.tre")) {

                        tree1list.add(f);
                        treeFile1.write(f.getAbsolutePath() + "\n");
                        treeFile1.flush();


                    }
                    if (f.getName().contains("segment_2.tre")) {

                        tree2list.add(f);
                        treeFile2.write(f.getAbsolutePath() + "\n");
                        treeFile2.flush();


                    }
                }

            }
            treeFile1.close();
            treeFile2.close();

            BufferedReader reader1 = new BufferedReader(new FileReader(treeFileNamesList1));
            BufferedReader reader2 = new BufferedReader(new FileReader(treeFileNamesList2));

            int i = 1;
            int j = 1;
            int counter = 0;
            System.out.println("Season,Segment1,Segment2,Diff");
            StringBuilder sBuilder1 = new StringBuilder();
            StringBuilder sBuilder2 = new StringBuilder();
            StringBuilder sBuilder3 = new StringBuilder();
            StringBuilder sBuilder4 = new StringBuilder();

            String seg1 = "";
            String seg2 = "";

            StringBuilder sBuilder_diff = new StringBuilder();
            while (reader1.ready() && reader2.ready()) {

                String tree1 = reader1.readLine();
                String tree2 = reader2.readLine();

                //System.out.println(tree1);

                Reader readerTree1 = new FileReader(tree1);
                Reader readerTree2 = new FileReader(tree2);

                NexusImporter importer1 = new NexusImporter(readerTree1);
                NexusImporter importer2 = new NexusImporter(readerTree2);

                RootedTree t1 = (RootedTree) importer1.importNextTree();
                RootedTree t2 = (RootedTree) importer2.importNextTree();

                double[] fitnessStats = getFitnessStats(t1);


                //System.out.println(i+","+getRootAge(tree1)+","+getRootAge(tree2)+","+(getRootAge(tree1)-getRootAge(tree2)));

                sBuilder1.append(getRootAge(t1)).append(",");
                sBuilder2.append(getRootAge(t2)).append(",");
                sBuilder_diff.append(getRootAge(t1) - getRootAge(t2)).append(",");
                sBuilder3.append(fitnessStats[0]).append(",");
                sBuilder4.append(fitnessStats[1]).append(",");


                if (i % 40 == 0) {
                    sBuilder1.append("\n");
                    sBuilder2.append("\n");
                    sBuilder_diff.append("\n");
                    sBuilder3.append("\n");
                    sBuilder4.append("\n");


                }
                //System.out.println(i+","+treeFile1);

                Double rootAge_tree1 = getRootAge(t1);
                Double rootAge_tree2 = getRootAge(t2);
                Double diff = rootAge_tree1 - rootAge_tree2;

//                if(counter == 19) {
//                    //System.out.println(treeFile1);
//                    seg1 = "segment 1,"+j+sBuilder1.toString()+"\n";
//                    seg2 = "segment 2,"+j+sBuilder2.toString()+"\n";
//
//                    System.out.println(seg1);
//                    System.out.println(seg2);
//
//                    sBuilder1 = new StringBuilder();
//
//                    sBuilder2 = new StringBuilder();
//                    j++;
//                    counter = 0;
//
//
//                }
//                else{
//
//                    sBuilder1.append(","+rootAge_tree1);
//                    sBuilder2.append(","+rootAge_tree2);
//                    counter++;
//                }


                i++;



            }System.out.println(sBuilder1.toString());
            System.out.println();
            System.out.println(sBuilder2.toString());
            System.out.println();
            System.out.println(sBuilder_diff.toString());
            System.out.println();
            System.out.println(sBuilder3.toString());
            System.out.println();
            System.out.println(sBuilder4.toString());

            System.out.println(folder.getAbsolutePath());

            FileWriter writer = new FileWriter(new File(folder.getAbsoluteFile()+"/reassortment_results_20150415.csv"));

            writer.write("\n"+sBuilder1.toString()+"\n\n");
            writer.write(sBuilder2.toString()+"\n\n");
            writer.write(sBuilder3.toString()+"\n\n");
            writer.write(sBuilder4.toString()+"\n");

        }catch(FileNotFoundException e){
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }catch(IOException e){
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }catch(ImportException e){
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }



    }




    public void getBranches(RootedTree tree) {

        Node rootNode = tree.getRootNode();
        Set<Node> nodes = tree.getNodes();

        for(Node n: nodes) {

            if(!tree.isRoot(n)) {
                System.out.println(tree.getHeight(n));
            }
        }
    }

    public double[] getFitnessStats(RootedTree tree) {

        Set<Node> nodes = tree.getExternalNodes();

        double[] results = new double[2];
        double[] fitness = new double[nodes.size()];
        int i = 0;
        for(Node n: nodes) {

            fitness[i] = (Double)n.getAttribute("fitness");
            i++;


        }

        DescriptiveStatistics statistics = new DescriptiveStatistics(fitness);

        double mean = statistics.getMean();
        double var = statistics.getVariance();

        return new double[]{mean, var};


    }

    private double getRootAge(RootedTree tree) {

        Node rootNode = tree.getRootNode();
        return tree.getHeight(rootNode);
    }
}
