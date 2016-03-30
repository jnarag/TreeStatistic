import com.sun.javafx.collections.transformation.SortedList;
import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.RootedTree;

import java.io.*;
import java.util.*;

/**
 * Created by jayna on 01/12/14.
 */
public class TraitValues {

    String treeFile = "";
    RootedTree tree = null;
    double mrsd = 0;
    String trait_name = "height";

    private  Comparator<Node>  comparator = new Comparator<Node>() {
        @Override
        public int compare(Node o1, Node o2) {

            double diff = tree.getHeight(o1)-tree.getHeight(o2);
            if(diff == 0) {

                if(countChildren(tree, o1) == countChildren(tree, o2)) {
                    return 0;
                }
                else if(countChildren(tree, o1) > countChildren(tree, o2)) {
                    return 1;

                }
                else{
                    return -1;
                }
            }
            else if(diff > 0) {
                return 1;
            }
            else{
                return -1;
            }
        }
    };

    public TraitValues(String treeFile) {

        this.treeFile = treeFile;
        Reader reader = null;
        try {
            reader = new FileReader(treeFile);

            NexusImporter importer = new NexusImporter(reader);

            tree = (RootedTree) importer.importNextTree();

        } catch (IOException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        }

    }

    public void getTraitMatrix() {

        Set<Node> int_nodes = tree.getInternalNodes();
        Set<Node> lead_nodes = tree.getExternalNodes();

        List<Double> nodesX = new ArrayList<Double>();
        List<Double> nodesY = new ArrayList<Double>();
        List<Double> nodesXp = new ArrayList<Double>();
        List<Double> nodesYp = new ArrayList<Double>();

        List<Double> leavesX = new ArrayList<Double>();
        List<Double> leavesY = new ArrayList<Double>();

//        List<double[]> X = new ArrayList<double[]>();
//        List<double[]> Y = new ArrayList<double[]>();

        Set<Node> nodes = tree.getNodes();
        nodes.remove(tree.getRootNode());

        List<Node> sortedNodes = new ArrayList<Node>();
        sortedNodes.addAll(nodes);

        Collections.sort(sortedNodes, comparator);

//        for(int i=0; i<sortedNodes.size(); i++) {
//
//            for(int j=0; j < sortedNodes.size(); j++) {
//
//                Node node_i = sortedNodes.get(i);
//                Node node_j = sortedNodes.get(j);
//
//                if(tree.getHeight(sortedNodes.get(i)) > tree.getHeight(sortedNodes.get(j))) {
//
//                    sortedNodes.set(i, node_j);
//                    sortedNodes.set(j, node_i);
//
//                }
//
//
//            }
//        }





//        for(int i=0; i < sortedNodes.size(); i++) {
//
//            List<Node> children = tree.getChildren(sortedNodes.get(i));
//
//
//            Node leftChild = children.get(0);
//            Node rightChild = children.get(1);
//
//
//            int r_index = sortedNodes.indexOf(rightChild);
//            int l_index = sortedNodes.indexOf(leftChild);
//
//            Node newLeft = null;
//
//            if (tree.isExternal(leftChild) && tree.isExternal(rightChild)) {
//
//            } else if (tree.isExternal(leftChild) && !tree.isExternal(rightChild)) {
//
//            } else if (!tree.isExternal(leftChild) && tree.isExternal(rightChild)) {
//
//                newLeft = rightChild;
//                rightChild = leftChild;
//                leftChild = newLeft;
//
//                sortedNodes.set(r_index, rightChild);
//                sortedNodes.set(l_index, leftChild);
//            } else if (!tree.isExternal(leftChild) && !tree.isExternal(rightChild)) {
//
//                if (countChildren(tree, leftChild) > countChildren(tree, rightChild)) {
//                    newLeft = rightChild;
//                    rightChild = leftChild;
//                    leftChild = newLeft;
//
//                    sortedNodes.set(r_index, rightChild);
//                    sortedNodes.set(l_index, leftChild);
//                } else if (countChildren(tree, leftChild) < countChildren(tree, rightChild)) {
//
//                } else {
//
//                    if (tree.getHeight(leftChild) > tree.getHeight(rightChild)) {
//                        newLeft = rightChild;
//                        rightChild = leftChild;
//                        leftChild = newLeft;
//
//                        sortedNodes.set(r_index, rightChild);
//                        sortedNodes.set(l_index, leftChild);
//                    }
//                }
//            }
//        }


        //}


        //}

        //sortedNodes.remove(tree.getRootNode());

        for(Node n: sortedNodes) {

            //System.out.println(tree.getHeight(n));

            Node p = tree.getParent(n);


            double height = mrsd-tree.getHeight(n);
            double trait = (Double)n.getAttribute(trait_name);

            double height_p = mrsd-tree.getHeight(p);
            double trait_p;

            if(tree.isRoot(p)) {

                trait_p = trait;

            }
            else{
                trait_p = (Double)p.getAttribute(trait_name);
            }



            nodesX.add(height);
            nodesXp.add(height_p);
            nodesY.add(trait);
            nodesYp.add(trait_p);


//            else{
//                nodesX.add(new double[]{height, height_p});
//                nodesY.add(new double[]{trait, trait_p});
//            }


            if(tree.isExternal(n)) {
                leavesX.add(height);
                leavesY.add(trait);
            }


            //X.addAll(leavesX);

            //Y.addAll(leavesY);

        }

        System.out.println(leavesX.size());
//        System.out.println(Y);


        try {
            BufferedWriter writer1 = new BufferedWriter(new FileWriter(new File("/Users/jayna/Dropbox/duke_reassortment/2014-11-23/1/X.csv")));
            BufferedWriter writer2 = new BufferedWriter(new FileWriter(new File("/Users/jayna/Dropbox/duke_reassortment/2014-11-23/1/Y.csv")));
            BufferedWriter writer3 = new BufferedWriter(new FileWriter(new File("/Users/jayna/Dropbox/duke_reassortment/2014-11-23/1/leaves_X.csv")));
            BufferedWriter writer4 = new BufferedWriter(new FileWriter(new File("/Users/jayna/Dropbox/duke_reassortment/2014-11-23/1/leaves_Y.csv")));
            BufferedWriter writer5 = new BufferedWriter(new FileWriter(new File("/Users/jayna/Dropbox/duke_reassortment/2014-11-23/1/Xp.csv")));
            BufferedWriter writer6 = new BufferedWriter(new FileWriter(new File("/Users/jayna/Dropbox/duke_reassortment/2014-11-23/1/Yp.csv")));


            for(double d: nodesX) {
                writer1.write(d+"\n");
                writer1.flush();
            }

            for(double d: nodesY) {
                writer2.write(d+"\n");
                writer2.flush();
            }

            for(double d: leavesX) {
                writer3.write(d+"\n");
                writer3.flush();
            }

            for(double d: leavesY) {
                writer4.write(d+"\n");
                writer4.flush();
            }
            for(double d: nodesXp) {
                writer5.write(d+"\n");
                writer5.flush();
            }

            for(double d: nodesYp) {
                writer6.write(d+"\n");
                writer6.flush();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }



//        for(double[] d: Y) {
//            System.out.println(d[0]+","+d[1]);
//        }

    }
    private int countChildren (RootedTree tree, Node node) {

        if(node == null) {
            return 0;
        }
        if( tree.isExternal(node) ) {
            return 1;
        }

        if( tree.isExternal(tree.getChildren(node).get(0)) && tree.isExternal(tree.getChildren(node).get(1)) ) {
            return 2;
        } else {
            return countChildren(tree,  tree.getChildren(node).get(0)) + countChildren(tree,  tree.getChildren(node).get(1));
        }
    }

    public static void main(String [] args) {

        //TraitValues tv = new TraitValues("/Users/jayna/Dropbox/InfluenzaB/2014-09-03/Vic_mds.FigTreeFormat.mcc.tre");
        //TraitValues tv = new TraitValues("/Users/jayna/Dropbox/duke_reassortment/2014-11-23/1/psi_0.0/antigenicMu_0.001_s_0.0125/tree_300N_3000000_antigenicMu_0.001_s_0.0125_psi_0.0_psij_0.5_samplePeriod_20.0_to_60.0_simNo_0_segment_1.tre");
        //TraitValues tv = new TraitValues("/Users/jayna/Dropbox/duke_reassortment/2014-11-23/1/psi_1.0/antigenicMu_0.000001_s_0.107/tree_300N_3000000_antigenicMu_1.0E-6_s_0.107_psi_1.0_psij_0.5_samplePeriod_20.0_to_60.0_simNo_0_segment_1.tre");
        TraitValues tv = new TraitValues("/Users/jayna/Dropbox/duke_reassortment/2014-11-23/1/psi_1.0/small_and_large/tree_300N_3000000_antigenicMu_0.001_s_0.0125_psi_1.0_psij_0.5_samplePeriod_20.0_to_60.0_simNo_0_segment_1.tre");
        tv.mrsd = 60;
        tv.trait_name = "fitness";
        tv.getTraitMatrix();
    }


}
