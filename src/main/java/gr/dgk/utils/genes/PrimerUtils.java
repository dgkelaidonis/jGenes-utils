package gr.dgk.utils.genes;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import com.sun.istack.internal.NotNull;

public class PrimerUtils {

	private static BufferedReader reader;

	/**
	 * This method checks a sequence of ATCG letters for potential loops. The given
	 * string must be comprised by even number of letters. Otherwise, it is not
	 * checked.
	 * 
	 * @param primer
	 * @return A structured string as text, that presents analytics info about the
	 *         given input.
	 */
	public static String checkPrimerForLoops(@NotNull String primer) {
		/* evaluate the given value */
		if (!isValidSequence(primer)) {
			return null;
		}
		/* start analysis */
		int loops = 0;
		/* groom the primer string */
		primer = primer.toUpperCase().replaceAll("\\s", "");

		String cArray = "";
		String pArray = "";
		int p = 1;
		for (char s : primer.toCharArray()) {
			cArray = cArray + ((p < 10) ? ("| " + s) : ("| " + s));
			pArray = pArray + ((p < 10) ? ("| " + p) : ("|" + p));
			p++;
		}
		String resultText = "\n-----------------------------------------------------------------------------------------------\n"
				+ cArray + "\n" + pArray
				+ "\n-----------------------------------------------------------------------------------------------\n";
		/*
		 * set the size of array based on the middle position minus 2 for the last left
		 * and the last right position (e.g. for 20 -> mid is 10 and table size should
		 * be 10 numbers, that tart from 0, to 9).
		 */
		int primerLength = primer.length();
		/* we check only primers with even number of letters */
		if (primerLength % 2 == 0) {
			int midIndex = (primerLength / 2);
			// int endOfLength = middleOfLength - 2;
			String leftPartStr = primer.substring(0, midIndex).toUpperCase().replaceAll("\\s", "");
			String rightPartStr = primer.substring(midIndex, primerLength).toUpperCase().replaceAll("\\s", "");
			resultText = resultText.concat("Check two half parts: <" + leftPartStr + "> | <" + rightPartStr + ">\n");
			/* get the substring as leftPart <--middle--> rightPart */
			char[] leftPart = leftPartStr.toCharArray();
			char[] rightPart = reverseString(rightPartStr)
					.toCharArray();/* reverse the right part for bringing the last position to first */

			/*
			 * Iterate on char arrays that they have the same size and the right array
			 * should be reversed
			 */
			for (int i = 0; i < midIndex; i++) {
				if ((leftPart[i] == 'A' && rightPart[i] == 'T') || (leftPart[i] == 'T' && rightPart[i] == 'A')) {
					resultText = resultText.concat(leftPart[i] + "-" + rightPart[i] + " at [Left=" + (i + 1)
							+ ", Right=" + (primerLength - i) + "]\n");
					loops++;
				} else if ((leftPart[i] == 'G' && rightPart[i] == 'C') || (leftPart[i] == 'C' && rightPart[i] == 'G')) {
					resultText = resultText.concat(leftPart[i] + "-" + rightPart[i] + " at [Left=" + (i + 1)
							+ ", Right=" + (primerLength - i) + "]\n");
					loops++;
				}
			}
			/* set numbers */
			resultText = resultText.concat("----------------\nTotal number=" + loops + "\n");
		} else {
			resultText = resultText.concat("Loops check, not applicable on primers with no even number of  letters.\n");
		}
		return resultText;
	}

	/**
	 * This method gets as input a CSV file and checks its contents, that
	 * corresponds to a sequence of ATCG letters, for potential loops. The given
	 * string must be comprised by even number of letters. Otherwise, it is not
	 * checked.
	 * 
	 * @param csvFilepath
	 * @return A structured string as text, that presents analytics info about the
	 *         given input. The text result is written into an output file.
	 * @throws IOException
	 */
	public static File checkPrimerForLoopsFromCSV(String csvFilepath) throws IOException {
		File file = new File(csvFilepath);
		reader = new BufferedReader(new FileReader(file));
		String csvLine, resultTxt = "";
		int sequenceNumber = 0;
		/* Build the output file path and setup the writer stream */
		String resultFilepath = csvFilepath.concat("_result_" + System.currentTimeMillis() + ".txt");
		BufferedWriter f_writer = new BufferedWriter(new FileWriter(resultFilepath));

		/* source the file contents line by line */
		while ((csvLine = reader.readLine()) != null) {
			/* get the comma separated values that corresponds to sequences */
			for (String primer : csvLine.split(",")) {
				resultTxt = resultTxt.concat(PrimerUtils.checkPrimerForLoops(primer));
				sequenceNumber++;
				/* Write the string content into the file */
				f_writer.write(resultTxt);
			}
		}
		/* Write the string content into the file */
		f_writer.write(
				resultTxt.concat("\n***********************\nTotal number of analyzed sequences = " + sequenceNumber));
		f_writer.close();

		return new File(resultFilepath);
	}

	/**
	 * This is a wrapper method for the specific symmetric/asymmetric methods, that
	 * are used to check pairs between primers that are either symmetric or
	 * asymmetric.
	 * 
	 * @param leftPrimer  Sequence of ATGC gene encoding.
	 * @param rightPrimer Sequence of ATGC gene encoding.
	 * @return A structured string as text, that presents analytics info about the
	 *         given input.
	 */
	public static String findPairsBetweenPrimers(@NotNull String leftPrimer, @NotNull String rightPrimer) {
		/* evaluate the given value */
		if (!isValidSequence(leftPrimer) || !isValidSequence(rightPrimer)) {
			return null;
		}

		/* groom primers since they are valid */
		leftPrimer = leftPrimer.toUpperCase().replace("\\s", "");
		rightPrimer = rightPrimer.toUpperCase().replace("\\s", "");

		/* check the number of letters' difference */
		int difference = (leftPrimer.length() - rightPrimer.length());
		difference = (difference < 0) ? difference * -1 : difference;

		/* select the appropriate function for checking the pairs */
		if (difference < 3) {
			return findPairsBetweenSymmetricPrimers(leftPrimer, rightPrimer);
		} else {
			return findPairsBetweenAsymmetricPrimers(leftPrimer, rightPrimer, 6);
		}
	}

	/**
	 * This method gets 2 symmetric primers (difference between letters < 3) and it
	 * checks whether they can create pairs with each other. For such length
	 * difference, do explicit comparison by leaving outside the 1 or 2 redundant
	 * letters.
	 * 
	 * @param leftPrimer  Sting with ATCG letters combination.
	 * @param rightPrimer Sting with ATCG letters combination.
	 * @return A structured string as text, that presents analytics info about the
	 *         given input.
	 */
	private static String findPairsBetweenAsymmetricPrimers(@NotNull String leftPrimer, @NotNull String rightPrimer,
			@NotNull int checkedPositions) {
		String resultText = "";
		int pairs = 0;
		resultText = resultText
				.concat("-----------------------------------------------------------------------------------\n");
		resultText = resultText.concat("Asymmetric Primers<" + leftPrimer + "(" + leftPrimer.length() + "), "
				+ rightPrimer + "(" + rightPrimer.length() + ")>\n");

		/*
		 * Create primers chunks by cutting them for checking the 6 last letters from
		 * the left primer and the 6 first from the right primers
		 */
		leftPrimer = leftPrimer.substring(leftPrimer.length() - checkedPositions, leftPrimer.length());
		rightPrimer = reverseString(rightPrimer.substring(0, checkedPositions));
		resultText = resultText.concat("Check between chunks <" + leftPrimer + " (" + leftPrimer.length()
				+ " last letters)," + rightPrimer + " (" + rightPrimer.length() + " first letters reversed)>\n");
		resultText = resultText
				.concat("-----------------------------------------------------------------------------------\n");

		/* break to char array */
		char[] leftPrimerCharArray = new char[leftPrimer.length()];
		char[] rightPrimerCharArray = new char[rightPrimer.length()];
		leftPrimerCharArray = leftPrimer.toCharArray();
		rightPrimerCharArray = rightPrimer.toCharArray();

		/*
		 * Iterate on available sequence's letters and find potential pairs between A-T,
		 * C-G
		 */
		resultText = resultText.concat("Available Pairs:\n");

		for (int i = 0; i < 6; i++) {
			if ((leftPrimerCharArray[i] == 'A' && rightPrimerCharArray[i] == 'T')
					|| (leftPrimerCharArray[i] == 'T' && rightPrimerCharArray[i] == 'A')) {
				resultText = resultText.concat(
						leftPrimerCharArray[i] + " - " + rightPrimerCharArray[i] + ": [at position=" + (i + 1) + "]\n");
				pairs++;
			} else if ((leftPrimerCharArray[i] == 'C' && rightPrimerCharArray[i] == 'G')
					|| (leftPrimerCharArray[i] == 'G' && rightPrimerCharArray[i] == 'C')) {
				resultText = resultText.concat(
						leftPrimerCharArray[i] + " - " + rightPrimerCharArray[i] + ": [at position=" + (i + 1) + "]\n");
				pairs++;
			}
		}

		resultText = resultText.concat("----------------\nTotal number=" + pairs + "\n\n");
		return resultText;
	}

	/**
	 * This method checks asymmetric primers. Specifically, this method should be
	 * used if the primers' length difference is >= 3 then do comparison between
	 * Left-5-last & Right-5-first.
	 * 
	 * @param leftPrimer  Sequence of ATGC gene encoding.
	 * @param rightPrimer Sequence of ATGC gene encoding.
	 * @return A structured string as text, that presents analytics info about the
	 *         given input.
	 */
	private static String findPairsBetweenSymmetricPrimers(@NotNull String leftPrimer, @NotNull String rightPrimer) {
		String resultText = "";
		int pairs = 0;

		/* break to char array */
		char[] leftPrimerCharArray = new char[leftPrimer.length()];
		char[] B = new char[rightPrimer.length()];
		resultText = resultText
				.concat("-----------------------------------------------------------------------------------\n");
		resultText = resultText.concat("Symmetric Primers<" + leftPrimer + "(" + leftPrimer.length() + "), "
				+ rightPrimer + "(" + rightPrimer.length() + ")>\n");
		resultText = resultText
				.concat("-----------------------------------------------------------------------------------\n");
		leftPrimerCharArray = leftPrimer.toCharArray();
		B = rightPrimer.toCharArray();
		/*  */
		int maxIndex = (leftPrimerCharArray.length >= B.length) ? B.length : leftPrimerCharArray.length;
		/*
		 * Iterate on available sequence's letters and find potential pairs between A-T,
		 * C-G
		 */
		resultText = resultText.concat("Available Pairs:\n");
		for (int i = 0; i < maxIndex; i++) {
			if ((leftPrimerCharArray[i] == 'A' && B[i] == 'T') || (leftPrimerCharArray[i] == 'T' && B[i] == 'A')) {
				resultText = resultText
						.concat(leftPrimerCharArray[i] + " - " + B[i] + ": [at position=" + (i + 1) + "]\n");
				pairs++;
			} else if ((leftPrimerCharArray[i] == 'C' && B[i] == 'G')
					|| (leftPrimerCharArray[i] == 'G' && B[i] == 'C')) {
				resultText = resultText
						.concat(leftPrimerCharArray[i] + " - " + B[i] + ": [at position=" + (i + 1) + "]\n");
				pairs++;
			}
		}
		resultText = resultText.concat("----------------\nTotal number=" + pairs + "\n\n");
		return resultText;
	}

	/**
	 * It gets a list of strings in form <"left_primer,right_primer"> and it
	 * analyses the overall input by printing the results.
	 * 
	 * @param pairOfPrimersInCsvFormat correspond to a List of Strings, that each
	 *                                 String equals to format
	 *                                 "left_primer,right_primer".
	 * @return A structured string as text, that presents analytics info about the
	 *         given input.
	 */
	public static String findPairsBetweenPrimersListWithCsvFormat(@NotNull List<String> pairOfPrimersInCsvFormat) {
		String result = "";
		for (String pairOfPrimersCsv : pairOfPrimersInCsvFormat) {
			/* 2 primers per record in the List */
			if (pairOfPrimersCsv.split(",").length == 2) {
				result = result
						.concat(findPairsBetweenPrimers(pairOfPrimersCsv.split(",")[0], pairOfPrimersCsv.split(",")[1]))
						.concat("\n");
			}
		}
		return result;
	}

	/**
	 * This method gets as input a CSV text file that includes in each line a comma
	 * separated set of primers, in form of left_primer, right_primer. It source the
	 * data from the text, using input stream technique and it grooms the data of
	 * each line, by trimming the spaces and converting any input to upper case, so
	 * as to be compatible with genes sequences format.
	 * 
	 * @param inputCsvFilepath: path on the FS of the CSV txt file with the primers'
	 *        pairs.
	 * @param clazz: File or String, as the given Object type.
	 * @return A structured string as text, that presents analytics info about the
	 *         given input. Depending on the given clazz value, the output is
	 *         written either into a File, or it is returned as string.
	 *
	 * @throws Exception
	 */
	public static <T> T findPairsBetweenPrimersFromCSVFile(@NotNull String inputCsvFilepath, Class<T> clazz)
			throws Exception {
		File file = new File(inputCsvFilepath);
		reader = new BufferedReader(new FileReader(file));
		String csvLine, resultTxt = "";
		int sequenceNumber = 0;
		while ((csvLine = reader.readLine()) != null) {
			resultTxt = resultTxt.concat(
					PrimerUtils.findPairsBetweenPrimers(csvLine.split(",")[0].toUpperCase().replaceAll("\\s", ""),
							csvLine.split(",")[1].toUpperCase().replaceAll("\\s", "")));
			sequenceNumber++;
		}
		resultTxt = resultTxt
				.concat("\n***********************\nTotal number of analyzed sequences = " + sequenceNumber);
		if (clazz.getTypeName().equals(File.class.getTypeName())) {
			if (inputCsvFilepath.isEmpty()) {
				throw new Exception(
						"Invalid filepath. Please give a valid path, that will lead to the CSV text file, location.");
			}
			/* Build the output file path and setup the writer stream */
			String resultFilepath = inputCsvFilepath.concat("_result_" + System.currentTimeMillis() + ".txt");
			BufferedWriter f_writer = new BufferedWriter(new FileWriter(resultFilepath));
			/* Write the string content into the file */
			f_writer.write(resultTxt);
			f_writer.close();
			return (T) new File(resultFilepath);
		} else if (clazz.getTypeName().equals(String.class.getTypeName())) {
			return (T) resultTxt;
		} else {
			throw new Exception(
					"Invalid data output format. Please select between the two (2) of the followings: File, String.");
		}
	}

	/**
	 * This method is used so as to evaluate a string sequence of nucleotides. If
	 * the string includes a letter, different than ATCG, or the string is
	 * null/empty, the evaluation returns false. Otherwise, it returns true.
	 * 
	 * @param primer A sequence of ATCG letters.
	 * @return True or False depending on the result of the validation process.
	 */
	public static boolean isValidSequence(@NotNull String primer) {
		String pattern = "ATCG";
		String groomedPrimer = primer.toUpperCase().replaceAll("\\s", "");
		if (primer.isEmpty()) {
			return false;
		}
		/*
		 * iterate on given string chars and ensure that it does not includes any
		 * invalid char
		 */
		for (char c : groomedPrimer.toCharArray()) {
			if (pattern.indexOf(c) == -1) {
				return false;
			}
		}
		return true;
	}

	/**
	 * This method gets a string, at maximum length of 100 characters, and it
	 * reverse it.
	 * 
	 * @param str The string that should be reversed.
	 * @return The given string in the reversed view (e.g. Hello --> olleH).
	 */
	private static String reverseString(@NotNull String str) {
		if (!str.isEmpty() && str.length() > 100) {
			return null;
		}
		String reversedStr = "";
		int rI = str.length() - 1;
		for (int i = 0; i < str.length(); i++) {
			reversedStr = reversedStr.concat(Character.toString(str.toCharArray()[rI - i]));
		}
		return reversedStr;
	}

}
