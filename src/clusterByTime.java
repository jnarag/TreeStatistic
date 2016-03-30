import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.RootedTree;

import java.io.*;
import java.util.*;

/**
 * Created by jayna on 27/05/15.
 */
public class clusterByTime {

    String treeFile;
    double TMRCA;
    String outputFile;
    double max_h;
    double min_h;
    public clusterByTime(String treeFile, double TMRCA) {

        this.treeFile = treeFile;
        this.TMRCA = TMRCA;
        this.outputFile = "output_clustersByTime_" + TMRCA + "yrs.csv";
        this.max_h = 0;
        this.min_h = 0;

    }

    public void findClusters() {

        try {
            NexusImporter importer = new NexusImporter(new FileReader(treeFile));

            BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
            int tree_no = 1;
            while(importer.hasTree()) {

                RootedTree tree = (RootedTree) importer.importNextTree();

                Set<Node> nodes = tree.getInternalNodes();


                int cluster_no = 1;
                for (Node node : nodes) {

                    List<Double> heights = new ArrayList<Double>();
                    Set<Node> clade = getTips(tree, node, heights);

                    if ((max_h - min_h) < TMRCA && (max_h - min_h) > 0.0 && clade.size()>1) {

                        for(Node n: clade) {
                            writer.write("tree" + tree_no + ",cluster" + cluster_no + "," +cluster_no + "," + tree.getTaxon(n).getName()+"\n");
                            writer.flush();

                        }
                        cluster_no++;
                    }
                }

                tree_no++;


            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public Set<Node> getTips(RootedTree tree, Node node, List<Double> heights) {

        Set<Node> tips = new HashSet<Node>();

        heights.add(tree.getHeight(node));
        if (tree.isExternal(node)) {


            heights.add(tree.getHeight(node));

            tips.add(node);


        } else {
            Node child1 = tree.getChildren(node).get(0);

            Set<Node> tips1 = getTips(tree, child1, heights);
            Node child2 = tree.getChildren(node).get(1);
            Set<Node> tips2 = getTips(tree, child2, heights);
            tips.addAll(tips1);
            tips.addAll(tips2);


        }


        //if(tips.size()>1) {
        max_h = Collections.max(heights);
        min_h = Collections.min(heights);


        if ((max_h - min_h) < TMRCA && (max_h - min_h) > 0.0 && tips.size()>1) {

            System.out.println("Cluster < TMRCA "+TMRCA+": " + (max_h) + " " + "," + min_h + "," + (max_h - min_h));

            for (Node n : tips) {


                System.out.println(tree.getTaxon(n).getName());
            }
            System.out.println();


        }

        return tips;

    }


    public static void main(String[] args){

        clusterByTime tree = new clusterByTime(args[0], Double.parseDouble(args[1]));
        tree.findClusters();

    }

}
