# jGenes-utils ACGT analyser
Java code utils for genes decoding and bioinformatics applications for genomics. The ACGT constitutes the four types of bases found in a DNA molecule.
## Where it is applicable
The given open-source utils, may be used for processing sequences of ‘DNA atcg’.
ACGT is an acronym for the four types of bases found in a DNA molecule: adenine (A), cytosine (C), guanine (G), and thymine (T). A DNA molecule consists of two strands wound around each other, with each strand held together by bonds between the bases. Adenine pairs with thymine, and cytosine pairs with guanine. For the genomics applications, they are used millions of ACGT combinations for performing experiments towards the genomic therapies. 
The deep analysis of the given combinations, constitutes a very important step for the bio-scientists, to ensure the success of the experiment. Invalid pairs, loops’ on sequences, asymmetric/symmetric sequences, etc., constitutes some of parts that should be exam before an experiment start. This java package, provides such basic operations, to help biologists to efficiently process their sequences developments.

## Package overview
The open-source version of this utilities package, includes the very basic operations that are used by biologists and bioinformatics applications. In particular,
- Find loops in a given sequence; it is used to identifying potential loops between the bases in a given ACGT sequence, with even number of bases.
- Find potential pairs among two difference sequences; supports both symmetric and asymmetric sequences. As Symmetric they are considered the sequences with delta, less than 3 bases (ACGT). Otherwise, the sequences with delta greater than 3 bases, are considered as asymmetric.
- Check validity of the given sequence. It gets a sequence and it checks its contents in terms of the bases letters. The given sequence must include only the four DNA bases decoded letters, i.e. ACGT. Otherwise, it is considered as invalid.
The provided methods, are based on streaming data-sourcing, to allow the real-time processing of huge amount of sequences. A sequence of bases, may construct tenth of Gigabytes of text files, that equals to millions of ACGT letters’ combinations. Thus, we should use input/output streaming mechanisms for handling and processing the data, instead to consume memory on the host.

## Indicative examples
Finding potential loops on sequences with even number of letters. Such sequences are used as primers for experiments. A loop among the bases may destroy the overall experiment and the operation of a primer. Thus the avoidance of the loops during the design phase, constitutes an important factor for the success. An indicative output of the check on a sequence, is shown below. The method 'checkPrimerForLoops(str)', gets as input the and it return the followings:

-----------------------------------------------------------------------------------------------
| T| A| G| C| G| C| T| T| G| G| C| A| T| A| C| T| G| T| A| T| C| T| A| T| A| T
| 1| 2| 3| 4| 5| 6| 7| 8| 9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26
-----------------------------------------------------------------------------------------------
Check two half parts: <TAGCGCTTGGCAT> | <ACTGTATCTATAT>
T-A at [Left=8, Right=19]
T-A at [Left=13, Right=14]
----------------
Total number=2

## Reusability options
The utility class will be used in the backend of a wider online system. It will be based on the microservices architecture and it will be accessible over an REST-API. The business logic of the microservices it will be implemented by this open-source package.
