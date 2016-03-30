/**
 * Created by jayna on 21/05/2014.
 */

import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusExporter;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.trees.SubtreeRootedTree;

import java.io.*;
import java.util.*;

public class dropTips {


    List<SubtreeRootedTree> subtrees;
    String outputTreeFileName;
    List<String> subtreeList;
    String subtreelistFileName;

    public dropTips(String subtreelistFileName) {

        this.subtreelistFileName = subtreelistFileName;
        outputTreeFileName = subtreelistFileName.replaceAll(".txt", "_output.trees");
        subtrees = new ArrayList<SubtreeRootedTree>();

    }

    public dropTips() {

        outputTreeFileName = "out.trees";
        subtrees = new ArrayList<SubtreeRootedTree>();
    }

    public void writeSubtrees(String inputTreeFileName, String substring, boolean excludedTaxa) {

        try {
            NexusImporter importer = new NexusImporter(new FileReader(inputTreeFileName));



            RootedTree firstTree = (RootedTree) importer.importNextTree();

            Set<Taxon> taxa = firstTree.getTaxa();
            Set<Taxon> subtreeTaxa = new HashSet<Taxon>();

            List<String> subtreeTaxa_names = new ArrayList<String>();
            List<String> taxa_names = new ArrayList<String>();


            for(Taxon t: taxa) {
                    taxa_names.add(t.getName());

            }

            List<String> included_taxaNames = new ArrayList<String>();

            for(String n: taxa_names) {

                if(excludedTaxa) {

                    if(!n.contains(substring)) {

                        included_taxaNames.add(n);
                    }
                }
                else{
                    if(n.contains(substring)) {

                        included_taxaNames.add(n);
                    }
                }
            }

            for(Taxon t: taxa) {

                for(String s: included_taxaNames) {

                    if(t.getName().compareToIgnoreCase(s) == 0) {

                        subtreeTaxa.add(t);
                    }
                }
            }



            SubtreeRootedTree firstSubtree = new SubtreeRootedTree(firstTree, subtreeTaxa);
            subtrees.add(firstSubtree);

            System.out.println("included taxa = "+included_taxaNames.size());
            int count = 1;
            while(importer.hasTree())  {

                count++;
                System.out.println("tree "+ count);
                RootedTree tree = (RootedTree) importer.importNextTree();
                SubtreeRootedTree subtree = new SubtreeRootedTree(tree, subtreeTaxa);

                subtrees.add(subtree);

            }

            NexusExporter exporter = new NexusExporter(new BufferedWriter(new FileWriter(outputTreeFileName)));
            exporter.exportTrees(subtrees);

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }
    public void writeSubtrees(String inputTreeFileName, boolean excludedTaxa) {

        try {
            NexusImporter importer = new NexusImporter(new FileReader(inputTreeFileName));
            BufferedReader reader = new BufferedReader(new FileReader(subtreelistFileName));



            RootedTree firstTree = (RootedTree) importer.importNextTree();

            Set<Taxon> taxa = firstTree.getTaxa();
            Set<Taxon> subtreeTaxa = new HashSet<Taxon>();

            List<String> subtreeTaxa_names = new ArrayList<String>();
            List<String> taxa_names = new ArrayList<String>();
            while(reader.ready()) {

                String taxon_name = reader.readLine().trim();
                subtreeTaxa_names.add(taxon_name);
                //System.out.println(taxon_name);

            }

            for(Taxon t: taxa) {

                taxa_names.add(t.getName());
            }

            List<String> included_taxaNames = new ArrayList<String>();
            if(excludedTaxa) {

                taxa_names.removeAll(subtreeTaxa_names);
                included_taxaNames.addAll(taxa_names);

            }
            else{

                included_taxaNames.addAll(subtreeTaxa_names);
            }

            for(Taxon t: taxa) {

                for(String s: included_taxaNames) {

                    if(t.getName().compareToIgnoreCase(s) == 0) {

                        subtreeTaxa.add(t);
                    }
                }
            }



            SubtreeRootedTree firstSubtree = new SubtreeRootedTree(firstTree, subtreeTaxa);
            subtrees.add(firstSubtree);

            System.out.println("included taxa = "+included_taxaNames.size());
            int count = 1;
            while(importer.hasTree())  {

                count++;
                System.out.println("tree "+ count);
                RootedTree tree = (RootedTree) importer.importNextTree();
                SubtreeRootedTree subtree = new SubtreeRootedTree(tree, subtreeTaxa);

                subtrees.add(subtree);

            }

            NexusExporter exporter = new NexusExporter(new BufferedWriter(new FileWriter(outputTreeFileName)));
            exporter.exportTrees(subtrees);

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public List<Double> getMeanRate(List<SubtreeRootedTree> trees, String rate_identifier){
        double totalMeanRate = 0;
        List<Double> meanRates = new ArrayList<Double>();

        String logFileName = subtreelistFileName.replaceAll(".txt", "_meanRate.log");
        BufferedWriter writer = null;

        try {
            writer = new BufferedWriter(new FileWriter(logFileName));
            writer.write("state,rate");

        } catch (IOException e) {
            e.printStackTrace();
        }


        int i = 0;
        for(RootedTree t: trees) {

            Set<Node> nodes = t.getNodes();

            double totalBranchLength = 0;
            double totalTime = 0;
            double totalbranchRate = 0;
            double totalSubstitutions = 0;

            for(Node n: nodes) {

                if(!t.isRoot(n)) {
                    double branchLength = t.getLength(n);
                    totalBranchLength += branchLength;

                    Map<String, Object> attributeMap = n.getAttributeMap();
                    if(attributeMap.containsKey(rate_identifier)) {
                        double rate = (Double) attributeMap.get(rate_identifier);
                        totalbranchRate += rate;

                        totalSubstitutions += (rate * branchLength);
                    }

                }


            }

            if (writer != null) {

                try {
                    writer.write(i+","+ (totalbranchRate/nodes.size())+"\n");
                    writer.flush();
                    i++;
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            meanRates.add(totalbranchRate / nodes.size());
            totalMeanRate += (totalbranchRate/nodes.size());
//            System.out.println("Total branch rate "+ (totalbranchRate));
//            System.out.println("Total branch length (time) "+ totalBranchLength);
//            System.out.println("Total branch substitutions "+ totalSubstitutions);
//            System.out.println(totalbranchRate/nodes.size());
        }
        try {
            assert writer != null;
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        totalMeanRate/=trees.size();
        System.out.println("Average mean rate: "+ totalMeanRate);



        return meanRates;

    }



    public static void main(String [] args) {

        dropTips dropTips = new dropTips();
        //System.out.println(Boolean.valueOf(args[2]));
        dropTips.writeSubtrees("sanger_ivr_gisaid.HA.raxml.tre", "YAM88", false);
        //dropTips.getMeanRate(dropTips.subtrees, args[3]);

    }



}
