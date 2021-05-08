package gr.dgk.utils.genes;

import java.io.File;
import java.util.List;

public class Main {
    private static String leftPrimer;
    private static String rightPrimer;
    private static List<String> pairOfPrimersInCsvFormat;
//
//    public static void main(String[] args) {
//	leftPrimer = "TCTGAAACTAGGCGGCAGAG";
//	rightPrimer = "GCTCCAGAGGTGCAGTTCTTT";
//	pairOfPrimersInCsvFormat = new ArrayList<String>();
//	pairOfPrimersInCsvFormat.add("TCTGAAACTAGGCGGCAGAG|GCTCCAGAGGTGCAGTTCTTT");
//	pairOfPrimersInCsvFormat.add("TCTGAAACTAGGCGGCAGAG|TCCAGAGGTGCAGTTCTTTTT");
//	pairOfPrimersInCsvFormat.add("TCTGAAACTAGGCGGCAGAG|TTTGGCCGGAGTAAGCTG");
//	pairOfPrimersInCsvFormat.add("TCTGAAACTAGGCGGCAGAG|TCCAGAGGTGCAGTTCTTTTT");
//	pairOfPrimersInCsvFormat.add("TTCTGAAACTAGGCGGCAGA|GCCGGAGTAAGCTGACAAAA");
//
//	String result02 = PrimerUtils.findPairsBetweenPrimersListOfCsv(pairOfPrimersInCsvFormat);
//	System.out.println(result02);
//    }

    public static void main(String[] args) throws Exception {
	PrimerUtils.findPairsBetweenPrimersFromCSVFile("primers_excel_input.txt", File.class);
    }
}
