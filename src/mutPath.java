import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NewickExporter;
import jebl.evolution.io.NexusExporter;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.trees.SimpleRootedTree;

import java.io.*;
import java.util.*;

/**
 * Created by jayna on 11/08/2014.
 */
public class mutPath {

    String treeFileName;
    String key;
    int no_of_sites;
    double mrsd;
    int cutoff;

    public mutPath() {

        treeFileName = "/Users/jayna/Documents/Projects/fluB/data/trees_from_pinky/mutational_mapping/aa/2014-10-30/yam/combined_sanger_bedford_YAM_18K_JTT_const_strict.trees";
        treeFileName= "/Users/jayna/Documents/Projects/fluB/data/trees_from_pinky/mutational_mapping/aa/2014-10-30/vic/combined_vic_4runs_JTT_aa.trees";
        //treeFileName = "combined_yam_const_aa_18K_MCC.tre";
        //treeFileName = "combined_sanger_bedford_YAM_18K_JTT_const_strict_MCC.tre";
        //treeFileName = "/Users/jayna/Dropbox/fluB/data/trees_from_pinky/mutational_mapping/aa/2014-09-10/yam/const_strict/combined_yam_const_aa_18K_MCC.tre";
        //treeFileName = "combined_run_2_3_vic_aa.MCC.tre";
        //treeFileName = "combined_run_2_3_vic_aa.MCC.tre";
        //treeFileName = "/Users/jayna/Documents/Projects/fluB/data/trees_from_pinky/mutational_mapping/aa/2014-09-10/vic/const_strict/combined_run_2_3_vic_aa.MCC.tre";
        //treeFileName = "B_NA_subset_aa_JTT_const_20M.MCC.tre";//jayna/Documents/Projects/fluB/data/trees_from_pinky/mutational_mapping/aa/2014-10-17/NA/1/B_NA_subset_aa_JTT_const_20M.MCC.tre";
        //treeFileName = "combined_vic_4runs_JTT_aa_MCC.tre";
        //treeFileName = "test_mutpath_vic.tre";
        //key = "subset_yam_HA_fullseq_JTT_BSP20_aa_20M";
        //key = "subset_vic_HA_JTT_strict_BSP20_50M";
        //key = "B_NA_subset_aa";
        key = "states";
        no_of_sites = 585;
        mrsd = 2014.081;
        cutoff=363;
    }

    public RootedTree readTrees() {

        RootedTree tree = null;
        try{

            Reader reader = new FileReader(treeFileName);
            NexusImporter importer = new NexusImporter(reader);
            tree = (RootedTree) importer.importNextTree();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return tree;
    }

    public void iterateTrees(RootedTree tree) {

        Map<String, Double> mutmapDates = new HashMap<String, Double>();

        mapTrunkMutations(tree, mutmapDates);
        mapAminoAcidSimilarity(tree, key);

        try {
            NexusExporter exporter = new NexusExporter(new BufferedWriter(new FileWriter("test_mutpath_vic_sanger_bedford_aaMapped.tre")));
            exporter.exportTree(tree);
        } catch (IOException e) {
            e.printStackTrace();
        }

        int[] timeArray = timeArray(1988,2015);
        int[] mutation_per_time = new int[timeArray.length];
        Arrays.fill(mutation_per_time,0);


        for(String s: mutmapDates.keySet()) {


            double time = mrsd-mutmapDates.get(s);

            for(int i = 0; i < timeArray.length; i++) {

                //System.out.println(time+","+(int)Math.floor(Math.round(time)));
                if((int)Math.floor(Math.round(time))==timeArray[i]) {

                    mutation_per_time[i]+=1;
                }
            }
            //System.out.println((mrsd-mutmapDates.get(s))+","+s);
        }


        for(int t=0; t<timeArray.length; t++) {

            System.out.println(timeArray[t]+","+mutation_per_time[t]);
        }

    }

    public void mapTrunkMutations(RootedTree tree, Map<String, Double> mutmapDates) {
        Set<Node> nodes = tree.getInternalNodes();
        int i = 1;
        int[] counts = new int[no_of_sites];
        Arrays.fill(counts, 0);


        for (Node n : nodes) {

            Node p = tree.getParent(n);
            //System.out.println(i+","+countChildren(tree, p));
            if (!tree.isRoot(n)) {

                int n_children = countChildren(tree, n);

                i++;

                if (n_children > 5) {

                    Map<String, Object> attributeMap = n.getAttributeMap();
                    if (attributeMap.containsKey(key)) {

                        String node_state = (String) attributeMap.get(key);
                        String parent_state = (String) p.getAttribute(key);

                        StringBuilder sBuilder = new StringBuilder();

                        sBuilder.append("time," + p.getAttribute("height_median"));

                        //System.out.println(">"+node_state.length());
                        //for each amino acid...
                        for (int a = 0; a < cutoff; a++) {



                            if(node_state.length() > no_of_sites) {
                                break;
                            }
                            //I may not need this if I do not ignore gaps/ambiguous AA states
                            int correction = 0;
//                                if (a > 177) {
//                                    correction++;
//
//                                }

                            char n_state = node_state.charAt(a);
                            char p_state = parent_state.charAt(a);


//                        if (n.getAttribute("posterior") != null && (Double) n.getAttribute("posterior") > 0.8) {


                            if (n_state != p_state) {
                                n.setAttribute("state_" + (a + 1 + correction), String.valueOf(n_state));
                                p.setAttribute("state_" + (a + 1 + correction), String.valueOf(p_state));
                                counts[a + correction]++;
                                //System.out.println((a+1)+": "+n_state + (a+1) + p_state);
                                sBuilder.append("," + n_state + (a + 1 + correction) + p_state);

                                StringBuilder s = new StringBuilder();
                                s.append(n_state + "," + String.valueOf(a + 1) + "," + p_state);
                                mutmapDates.put(s.toString(), tree.getHeight(n));

                            }
//                        }

                        }
                        System.out.println(sBuilder.toString());

                        i++;

                    }

                    //}
                }

            }
        }
        for(int c=0; c<cutoff; c++) {

            if(counts[c] > 0) {
                System.out.println((c + 1) + ", " + counts[c]);
            }
        }
    }

    public void mapAminoAcidSimilarity(RootedTree tree, String key) {

        Node root = tree.getRootNode();
        Map<String, Object> attributeMap = root.getAttributeMap();
        String root_seq = "";
        Double root_time = (Double)attributeMap.get("height");
        if(attributeMap.containsKey(key)) {

            root_seq = (String)attributeMap.get(key);
            System.out.println(root_seq.length());
        }

        Set<Node> nodes = tree.getNodes();
        nodes.remove(root);

//        List<Node> children = tree.getChildren(root);
//
//        String child1 = (String)children.get(0).getAttribute(key);
//        String child2 = (String)children.get(1).getAttribute(key);
//
//        System.out.println(root_seq);
//        System.out.println(child1);
//        System.out.println(child2);


        List<Node> weird_nodes = new ArrayList<Node>();
        for(Node n: nodes) {

            Map<String, Object> map = n.getAttributeMap();
            if(map.containsKey(key)) {
                String n_seq = (String) n.getAttribute(key);

                if (n_seq.length() > no_of_sites) {


                    weird_nodes.add(n);
                }
            }
        }
        nodes.removeAll(weird_nodes);


        for(Node n: nodes) {

            Map<String, Object> n_attributeMap = n.getAttributeMap();
            String n_seq = "";

            if(n_attributeMap.containsKey(key)) {

                n_seq = (String)n_attributeMap.get(key);

                if(n_seq.length() > no_of_sites) {
                    break;
                }
            }

            int count_r = 0; //count from root to node
            int count_p = 0; //count from parent to child node
            for(int i=0; i < cutoff; i++ ) {


                String base_r = root_seq.substring(i,i+1);
                String base_n = n_seq.substring(i,i+1);


                if(base_r.compareTo(base_n)!=0) {

                    count_r++;
                }


            }

            Node p = tree.getParent(n);
            String p_seq = (String)p.getAttribute(key);

            for(int i=0; i < cutoff; i++) {


                String base_p = p_seq.substring(i,i+1);
                String base_n = n_seq.substring(i,i+1);


                if(base_p.compareTo(base_n)!=0) {

                    count_p++;
                }
            }

//            Node temp_n = n;
//            int total_change_between_root_and_n = 0;
//            Double total_time = 0.0;
//            while(!tree.isRoot(p)) {
//
//
//                n_seq = (String)temp_n.getAttribute(key);
//                if(n_seq.length() > no_of_sites) {
//                    break;
//                }
//
//                for (int i = 0; i < n_seq.length(); i++) {
//
//                    String base_r = root_seq.substring(i,i+1);
//                    String base_n = n_seq.substring(i, i + 1);
//
//
//                    if (base_r.compareTo(base_n) != 0) {
//
//                        total_change_between_root_and_n++;
//
//                    }
//
//                }
//                total_time += (Double)n.getAttribute("height");
//                temp_n = p;
//                p = tree.getParent(temp_n);
//
//                n_seq = (String)temp_n.getAttribute(key);
//            }

            Double n_time = (Double)n_attributeMap.get("height");
            Double pwdiff = (double)count_r/root_seq.length();
            Double corrected_dist = (-19.0/20.0)*(Math.log10(1-((20/19)*pwdiff)));

            Double rate = pwdiff/(root_time-n_time);
            Double corrected_rate = corrected_dist/(root_time-n_time);

            Double pwdiff_p = (double)count_p/n_seq.length();

            //Double total_rate = ((double)total_change_between_root_and_n/(n_seq.length()))/(root_time-n_time);



            n.setAttribute("pw_aa_diff", pwdiff);
            n.setAttribute("raw_aa_rate", rate);
            //n.setAttribute("total_aa_rate", total_rate);
            n.setAttribute("corrected_aa_rate", corrected_rate);
            n.setAttribute("pw_aa_diff_p", pwdiff_p);

            Map<String, Object> map = n.getAttributeMap();

        }


    }
    public void mapHighlySupportedMutations(RootedTree tree, Map <String, Double> mutmapDates) {

        Set<Node> nodes = tree.getNodes();

        int i = 1;
        int[] counts = new int[no_of_sites];
        Arrays.fill(counts, 0);

        for(Node n: nodes) {

            if(!tree.isRoot(n)) {

                Map<String, Object> attributeMap = n.getAttributeMap();
                if(attributeMap.containsKey(key)) {

                    String  node_state = (String)attributeMap.get(key);
                    Node p = tree.getParent(n);
                    String parent_state = (String) p.getAttribute(key);


                    //System.out.println(i + ": " +tree.getHeight(n));

                    StringBuilder sBuilder = new StringBuilder();

                    sBuilder.append("time,"+(mrsd-tree.getHeight(n)));

                    for(int a=0; a<node_state.length(); a++) {

                        //I may not need this if I do not ignore gaps/ambiguous AA states
                        int correction = 0;
                        if(a>177) {
                            correction++;

                        }

                        char n_state = node_state.charAt(a);
                        char p_state = parent_state.charAt(a);

                        n.setAttribute("state_"+(a+1 + correction), String.valueOf(n_state));
                        p.setAttribute("state_"+(a+1 + correction), String.valueOf(p_state));


                        if( n.getAttribute("posterior")!=null && (Double)n.getAttribute("posterior")>0.8 ) {
                            if (n_state != p_state) {



                                counts[a+correction]++;
                                //System.out.println((a+1)+": "+n_state + (a+1) + p_state);
                                sBuilder.append("," + n_state + (a + 1+correction) + p_state);

                                StringBuilder s = new StringBuilder();
                                s.append(n_state +","+ String.valueOf(a + 1) +","+ p_state);
                                mutmapDates.put(s.toString(), tree.getHeight(n));

                            }
                        }

                    }
                    System.out.println(sBuilder.toString());

                    i++;

                }

            }

        }

        for(int c=0; c<no_of_sites; c++) {

            if(counts[c] == 1) {
                System.out.println((c + 1) + ", " + counts[c]);
            }
        }


    }

    public void getTraits(String filename, String edgeLengthFile) {


        try {
            List<Double> edgeLengths = new ArrayList<Double>();
            BufferedReader reader = new BufferedReader(new FileReader(edgeLengthFile));

            while(reader.ready()) {

                String line = reader.readLine();
                String[] parts = line.split(",");
                edgeLengths.add(Double.parseDouble(parts[1].trim()));


            }

            NexusImporter importer = new NexusImporter(new FileReader(filename));

            RootedTree tree = (RootedTree)importer.importNextTree();

            Set<Node> nodes = tree.getNodes();

            List<Double> traitvalues = new ArrayList<Double>();


            //System.out.println(tree.getLength(n));
            Node chosen_n = null;
            for(int i=0; i<edgeLengths.size(); i++) {


                for(Node n: nodes) {

                    if(Math.round(edgeLengths.get(i)-(Double.valueOf(tree.getLength(n))))==0) {
                        //System.out.println(edgeLengths.get(i)+","+tree.getLength(n));
                        //System.out.println(i+","+edgeLengths.get(i)+","+(Double)n.getAttribute("corrected_aa_rate"));
                        System.out.println((Double)n.getAttribute("corrected_aa_rate"));
                        traitvalues.add((Double)n.getAttribute("corrected_aa_rate"));
                        break;
                    }

                }


            }

//            for(int i=0; i<traitvalues.size(); i++){
//
//                System.out.println(traitvalues.get(i));
//            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        }

    }

    public void removeAttribute() {

        List<RootedTree> trees = new ArrayList<RootedTree>();

        RootedTree tree = null;
        try{

            Reader reader = new FileReader(treeFileName);
            NexusImporter importer = new NexusImporter(reader);


            int count = 1;
            while(importer.hasTree()) {

                tree = (RootedTree) importer.importNextTree();
                Set<Node> nodes = tree.getNodes();

                for(Node n: nodes) {

                    n.removeAttribute(key);
                }
                trees.add(tree);
                System.out.println(count);
                count++;

            }

            NexusExporter exporter = new NexusExporter(new BufferedWriter(new FileWriter(treeFileName.replace(".trees", "_stripped.trees"))));
            exporter.exportTrees(trees);

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }



    public int[] timeArray(int start, int end) {

        int[] array = new int[end-start];

        for(int i=0; i<array.length; i++) {

            array[i] = start;
            start++;
        }

        return array;
    }

//    public void getTranslateBlock(String filename) {
//
//
//        try {
//            NexusImporter importer = new NexusImporter(new FileReader(filename));
//            RootedTree tree = (RootedTree) importer.importNextTree();
//            Set<Node> nodes = tree.getNodes();
//            Set<Taxon> taxa = tree.getTaxa();
//            int count = 1;
//            Map<Integer, Taxon> translateMap = new HashMap<Integer, Taxon>();
//
////            for(Taxon t: taxa) {
////
////                Taxon new_taxon = new Taxon(String.valueOf(count));
////                System.out.println(count+"\t"+t.getName());
////
////                tree.renameTaxa(t,new_taxon);
////                count++;
////                    //n.setAttribute("Names", count);
////                //t.setAttribute("name", count);
////            }
//
//            for(Node n: nodes) {
//
//                n.removeAttribute("subset_vic_HA_JTT_strict_BSP20_50M.set");
//                n.removeAttribute("subset_vic_HA_JTT_strict_BSP20_50M");
//                n.removeAttribute("subset_vic_HA_JTT_strict_BSP20_50M.set.prob");
//
//
//            }
//
//            NexusExporter exporter = new NexusExporter(new BufferedWriter(new FileWriter("translated_"+filename)));
//            exporter.exportTree(tree);
//        } catch (IOException e) {
//            e.printStackTrace();
//        } catch (ImportException e) {
//            e.printStackTrace();
//        }
//
//
//    }


    private BitSet addClades(RootedTree tree, Node node, Map<Node, BitSet> cladeMap, Map<Taxon, Integer> taxaMap) {


        BitSet bits = new BitSet();
        if (tree.isExternal(node)) {

            int index = Integer.valueOf(taxaMap.get(tree.getTaxon(node)));
            bits.set(index);


        } else {


            for (Node child : tree.getChildren(node)) {
                bits.or(addClades(tree, child, cladeMap, taxaMap));
            }
        }
        cladeMap.put(node, bits);
        return bits;
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

        mutPath mutPath = new mutPath();

        mutPath.removeAttribute();
//        RootedTree tree = mutPath.readTrees();
//        mutPath.iterateTrees(tree);

        //mutPath.getTraits("test_mutpath_vic_aaMapped.tre","vic_edgelengths.csv");
    }



}
