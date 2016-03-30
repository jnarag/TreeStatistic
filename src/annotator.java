import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusExporter;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;

import java.io.*;
import java.util.*;

/**
 * Created by jayna on 17/07/15.
 */
public class annotator {



    String treeFile;
    String outputFile;
    String trait;
    int position_in_string;

    public annotator(String treeFile, String trait, int position_in_string){

        this.trait = trait;
        this.treeFile = treeFile;
        this.position_in_string = position_in_string;
        outputFile = treeFile.replace(".tre","_annotatedByTime.tre");

    }

    public void labelBranch(boolean BEAST) {

        try {
            NexusImporter importer = new NexusImporter(new FileReader(treeFile));
            NexusExporter exporter = new NexusExporter(new BufferedWriter(new FileWriter(outputFile)));

            SimpleRootedTree tree = null;
            Set<String> timepoints = new HashSet<String>();

            if(importer.hasTree()) {
                tree = (SimpleRootedTree)importer.importNextTree();

                Set<Taxon> taxa = tree.getTaxa();

                for(Taxon t: taxa) {

                    String name = t.getName();
                    System.out.println(name);
                    String [] parts = name.split("_");
                    String trait = parts[position_in_string-1];

                    tree.getNode(t).setAttribute(this.trait, trait);
                 }
            }

            //String outputTreeFile = treeFile.replace(".tre", "_annotated.tre");

            exporter.exportTree(tree);

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    public void labelBranch() {

        String tree_string = null;

        try {
            BufferedReader reader = new BufferedReader(new FileReader(new File("/Users/jayna/Dropbox/liver_hcv/gismondi_liver_HCV/tree.txt")));
            tree_string = reader.readLine();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }


        //System.out.println(tree_string.length());
        Map<String, String> trait_map = new HashMap<String, String>();
        try {
            NexusImporter importer = new NexusImporter(new FileReader(treeFile));
            //NexusExporter exporter = new NexusExporter(new BufferedWriter(new FileWriter(outputFile)));

            SimpleRootedTree tree = null;
            Set<String> trait_list = new HashSet<String>();

            while (importer.hasTree()) {

                tree = (SimpleRootedTree) importer.importNextTree();

                Set<Taxon> taxa = tree.getTaxa();

                int i = 1;
                //for(Node n: nodes){

                for (Taxon t : taxa) {

                    Node n = tree.getNode(t);

                    if (tree.isExternal(n)) {
                        String name = t.getName();

                        String[] parts = name.split("_");

                        String tr = parts[position_in_string-1].trim();

                        trait_list.add(tr);

                        trait_map.put(name, tr);

//                        n.setAttribute("time", tr);
//                        if(tree.getLength(n)<=0.0) {
//                            System.out.println(i+","+name+","+tree.getLength(n));
//
//                            tree.setLength(n, 1E-27);
//                            i++;
//                        }


                    }
                }


            }

            System.out.println(tree_string);

            importer = new NexusImporter(new FileReader(treeFile));
            List<Taxon> taxaList = importer.parseTaxaBlock();
            for(Taxon t: taxaList) {


                Integer i = taxaList.indexOf(t)+1;
                String t_name = t.getName();
                //System.out.println(t_name);
                String new_string = "('"+t_name+"'[&"+trait+"=\""+trait_map.get(t_name)+"\",";
                System.out.println(new_string);
                tree_string = tree_string.replace("('"+t_name+"'[&", new_string);
                new_string = ",'"+t_name+"'[&"+trait+"=\""+trait_map.get(t_name)+"\",";
                tree_string = tree_string.replace(",'"+t_name+"'[&", new_string);



            }

            System.out.println(tree_string);

            //exporter.exportTree(tree);

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        }



    }

    public static void main(String [] args) {

        annotator treeAnnotator = new annotator("/Users/jayna/Dropbox/liver_hcv/gismondi_liver_HCV/p4/gismondi_HCV_patient4_BSP10_SDR06_UCLN_10M.MCC.tre", "location", 2);
        treeAnnotator.labelBranch();
    }


}
