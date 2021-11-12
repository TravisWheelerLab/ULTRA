# ULTRA

**U**LTRA **L**ocates **T**andemly **R**epetitive **A**reas

ULTRA is a tool to find and annotate tandem repeats inside genomic sequence. It
is able to find repeats of any length and of any period (up to a maximum period
of 4000). It can find highly decayed repeats missed by other software, and it
will also be able to find very large repeats in highly repetitive sequence,
regardless of the size of sequence or length of repeats. ULTRA offers meaningful
annotation scores and can produce annotation P-values at user request.

## Install

ULTRA uses the CMake build system (version 3.0 or greater). This makes building
the application pretty straightforward.

  1. Generate the Makefile: `cmake .`
  2. Build the `ultra` executable: `make`

Note that on some systems you may need to install CMake 3 as a separate package
and run `cmake3` specifically.

## Run

After building ULTRA you can use it by running:

```
./ultra [options] <input file>
```

ULTRA will produce high quality annotations out of the box, and so we can
analyze one of our example files by simply running.
 
```
ultra examples/ex1.fa
```

If you want to direct the output to a file you can do so with the -f flag:

```
ultra -f ex1_output.json examples/ex1.fa
```

Now your current working directory will contain an `ex1_output.json` file. At
the top of this file will be a complete description of the parameters used to
run ULTRA, and then there will be a list of repeats found, and their associated
annotation information. Output is in JSON format.

By default, ULTRA has a maximum detectable repeat period of 15. The sequence in
`examples/ex_large.fa` has a period 1000 repeat with a 10% substitution rate. To
detect larger period repeats, we will use the -p argument.

```
ultra -p 1000 -f ex_large.json examples/ex_large.fa
```

It's important to note that ULTRA's runtime grows linearly with maximum period
size. It’s also important to note that as a repeat increases in period, it
becomes harder to pin down the exact repeat period. 

Notice how ULTRA labels the period 1000 repeat inside `examples/ex_large.fa` as
being a period 999 repeat with 3 insertions. Period 999 ends up being the most
probable explanation of the `examples/ex_large.fa` repeat given the high
substitution rate, even though we have artificially created the repeat and
therefore know that it is in fact a period 1000 repeat with no indels.

ULTRA's default parameters are robust against different levels of AT richness.
For example, running ULTRA on ex_at80.fa (random sequence with 80% AT bias and a
single tandem repeat in the middle) with default parameters (run: `ultra
examples/ex_at80.fa`) yields a correct repeat annotation with one single tandem
repeat found.

However, running ULTRA with default parameters on `examples/ex_at90.fa` (random
sequence with 90% AT bias and a single tandem repeat in the middle) will result
in many false positive repeat annotations. We can clean up our annotation by
adjusting ULTRA's expected nucleotide distribution using the `-at` flag:

```
ultra -at 0.9 examples/ex_at90.fa
```

Using `-at 0.9`, ULTRA will correctly find the only tandem repeat inside the
sequence without reporting false positives. 

To see a full list of ULTRA arguments you can use ‘ultra -h’, or see the list of
arguments below.

```
------------------
0.99.17ultra <arguments> <input sequence path>
------------------
-flag	name	description [default value / behavior]
------------------
-f	Output File Path:                   Path to output csv. [STDOUT]
-s	Score Threshold:                    Minimum score necessary for a repeat to be recorded. [-10000.000000]
-ml	Region Length Trheshold:            Minimum total repeat length necessary for a repeat to be recorded. [10]
-mu	Repeat Unit Threshold:              Minimum number of repeat units necessary for a repeat to be recorded. [3]
-at	AT Richness:                        Frequency of As and Ts in the input sequence. [0.600000]
-atcg	ATCG Distribution:                  Frequency of As Ts Cs and Gs in the input sequence. [A=0.300000,T=0.300000,C=0.200000,G=0.200000]
-m	Match Probability:                  Probability of two nucleotides in consecutive repeat units being the same. [0.800000]

-n	Number of Threads:                  Number of threads used by ultra to analyze the input sequence. [1]
-p	Maximum Period:                     Largest repeat period detectable by Ultra. [15]
-mi	Maximum Consecutive Insertions:     Maximum number of insertions that can occur in tandem. [8]
-md	Maximum Consecutive Deletions:      Maximum number of deletions that can occur in tandem. [7]

-nr	Repeat Probability:                 Probability of transitioning from the nonrepetitive state to a repetitive state. [0.010000]
-rn	Repeat Stop Probability:            Probability of transitioning from a repetitive state to the nonrepetitive state. [0.050000]
-pd	Repeat probability decay:           **reminder to write a good description later. [0.850000]
-ri	Insertion Probability:              Probability of transitioning from a repetitive state to an insertion state. [0.020000]
-rd	Deletion Probability:               Probability of transitioning from a repetitive state to a deletion state. [0.020000]
-ii	Consecutive Insertion Probability:  The probability of transitioning from an insertion state to another insertion state. [0.020000]
-dd	Consecutive Deletion Probability:   The probability of transitioning from a deletion state to another deletion state. [0.020000]

-sr	***NOT FUNCTIONAL, DON'T USE***:    Split repeats such as ATATATGCGCGCGC into subrepeats ATATAT and GCGCGCGC. [0]
-sd	Split Depth:                        Number of repeat units to consider when splitting repeat. [5]
-sc	Split Cutoff:                       Cutoff value used during splitting (smaller is more conservative). [2]
-msp	Split Max Period:                   Maximum repeat period that will be considered for splitting. [6]

-hs	Hide Repeat Sequence:               Don't show the repetitive sequences in the results CSV. [Show Repeat Sequence]
-ss	Show score Deltas:                  Output the score change per residue. [False]
-st	Show traceback:                     Output the Viterbi traceback in the results CSV. [Hide traceback]
-wid	Show Window ID:                     Display the windowID corresponding to a repeat in the results CSV. [Hide Window ID]
-sl	Show logo numbers:                  Output the corresponding logo annotation for a given repeat. [Hide logo numbers]
-json	Read JSON file:                     Process all passes in JSON file. [False]
-jpass	Process passes in JSON file:        Process selected passes in JSON file. [None]
-pid	Pass ID:                            Assigns a custom pass ID. [Smallest unused positive pass ID]
-R	Completely Read File:               Read the entire input file during initialization.. [Do not read whole file]
-ws	Window Size:                        The number of nucleotides per sequence window. [8192]
-os	Window Overlap:                     The number of nucleotides overlaped between two consecutive windows. [100]
-wn	Number of Windows:                  Maximum number of windows stored in memory at once. [1024]
-v	Ultra Version:                      Shows Ultra's version. []

-doc	Debug overlap correction:           Report overlap correction in repeat information. [0]
```

## Development

ULTRA is written in C++ 17. Source files can be found in `src/`.

Run `cmake .` to generate the Makefile, then the following targets should work:

  * `make format` - format the code using clang-format
  * `make examples` - run all examples

To run the code formatter (`clang-format`) use `./tool/run-format.sh`.

There are two Docker containers, one for building the software (used in CI) and
another for running the software (this is the principle means by which we
distribute ULTRA).

Build and push the builder image: `./tool/build-builder-image.sh` and
`./tool/push-builder-image.sh`. Build and push the runner image:
`./tool/build-runner-image.sh` and `./tool/push-runner-image.sh`.

Please make a GitHub pull request to propose changes to the code.

## License

BSD 3-clause license. See LICENSE file.

## Authors

See `AUTHORS` file for up-to-date authors list.

  - Daniel Olson <daniel.olson@umontana.edu>
  - Travis Wheeler <travis.wheeler@umontana.edu>
  - George Lesica <george.lesica@umontana.edu>

