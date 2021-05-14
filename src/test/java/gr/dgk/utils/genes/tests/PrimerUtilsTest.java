package gr.dgk.utils.genes.tests;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Before;
import org.junit.Test;

import gr.dgk.utils.genes.PrimerUtils;

public class PrimerUtilsTest {
	private String loopCheckInput;
	private String leftPrimerInput;
	private String rightPrimerInput;
	private String invalidPrimerInput;
	private String validPrimerInput;
	private int loopCheckOutputLength;
	private int pairsCheckOutputLength;

	@Before
	public void init() {
		loopCheckInput = "TAGCGCTTGGCATACTGTATCTATAT";
		leftPrimerInput = "TAAACAGATCAGCAATTTCTTAACCAATGAATAC";
		rightPrimerInput = "ATCGTTAACCAATGAATACTAAACAGATCAGCAA";
		invalidPrimerInput = "AMAGAGATCGTCATATAGTATCTATAG";
		validPrimerInput = "ATAGAGATCGTCATATAGTATCTATAG";
		loopCheckOutputLength = 492;
		pairsCheckOutputLength = 481;
	}

	@Test
	public void test_checkPrimerForLoops() {
		String result = PrimerUtils.checkPrimerForLoops(loopCheckInput);
		System.out.println(result);
		assertEquals("Test <checkPrimerForLoops()> method", loopCheckOutputLength, result.length());
	}

	@Test
	public void test_findPairsBetweenPrimers() {
		String result = PrimerUtils.findPairsBetweenPrimers(leftPrimerInput, rightPrimerInput);
		System.out.println(result);
		assertEquals("Test <findPairsBetweenPrimers()> method - with symmetric primers", pairsCheckOutputLength, result.length());
	}

	@Test
	public void test_isValidSequence() {
		assertFalse(PrimerUtils.isValidSequence(invalidPrimerInput));
		assertTrue(PrimerUtils.isValidSequence(validPrimerInput));
	}

}
