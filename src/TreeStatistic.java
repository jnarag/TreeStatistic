import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.RootedTree;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: jayna
 * Date: 12/06/2013
 * Time: 11:22
 * To change this template use File | Settings | File Templates.
 */
public class TreeStatistic {


    RootedTree tree = null;
    public TreeStatistic() {




    }

    public void getBranches(String treeFileName) {

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

        Set<Node> nodes = tree.getNodes();

        for(Node n: nodes) {

            if(!tree.isRoot(n)) {
                System.out.println(tree.getLength(n));
            }
        }
    }

    public static void main(String [] args) {


        ProcessNodeAges pr = new ProcessNodeAges();
        pr.getRootAges();

//        TreeStatistic t = new TreeStatistic();
//
//        //t.getBranches("/Users/jayna/Dropbox/duke_reassortment/2014-09-21/psi_0.0/tree_40yrs/antigenicMu_0.000001_s_0.107/tree_500N_3000000_antigenicMu_1.0E-6_s_0.107_psi_0.0_psij_0.5_samplePeriod_20.0_to_60.0_simNo_0_segment_1.tre");
//        t.getBranches("/Users/jayna/Dropbox/duke_reassortment/2014-09-21/psi_1.0/tree_40yrs/antigenicMu_0.000001_s_0.107/tree_500N_3000000_antigenicMu_1.0E-6_s_0.107_psi_1.0_psij_0.5_samplePeriod_20.0_to_60.0_simNo_0_segment_1.tre");
        //t.getBranches("/Users/jayna/Dropbox/duke_reassortment/2014-09-21/psi_0.0/tree_40yrs/antigenicMu_0.000001_s_0.107/tree_500N_3000000_antigenicMu_1.0E-6_s_0.107_psi_0.0_psij_0.5_samplePeriod_20.0_to_60.0_simNo_0_segment_1.tre");

    }

}
