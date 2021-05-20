# ULTRA

**U**LTRA **L**ocates **T**andemly **R**epetitive **A**reas

**TODO(daniel): Describe what ULTRA is for and who might want to use it

## Install

ULTRA uses the CMake build system (version 3.0 or greater). This makes building
the application pretty straightforward.

  1. Generate the Makefile: `cmake .`
  2. Build the `ultra` executable: `make`

Note that on some systems you may need to install CMake 3 as a separate package
and run `cmake3` specifically.

## Run

TODO(daniel): Describe command line options and show examples

```
Usage: ultra <arguments> <input sequence path>

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

-sr	Split Divergent Repeats:            Split repeats such as ATATATGCGCGCGC into subrepeats. [0]
-sd	Split Depth:                        Number of repeat units to consider when splitting repeat. [5]
-sc	Split Cutoff:                       Cutoff value used during splitting (smaller is more conservative). [2]
-msp	Split Max Period:                   Maximum repeat period that will be considered for splitting. [6]

-hs	Hide Repeat Sequence:               Don't show the repetitive sequences in the results CSV. [Show Repeat Sequence]
-ss	Show scores:                        Output the score change per residue. [False]
-st	Show traceback:                     Output the Viterbi traceback in the results CSV. [Hide traceback]
-wid	Show Window ID:                     Display the windowID corresponding to a repeat in the results CSV. [Hide Window ID]
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

### Examples

There are some example data files in the `examples/` directory. See below for
some example invocations of the `ultra` command.

It is also possible to run all of the examples in this directory at once (which
happens in CI as a smoke test) using `make examples`.

## Development

ULTRA is written in C++ 17. Source files can be found in `src/`.

Run `cmake .` to generate the Makefile, then there are a couple custom targets
beyond what is necessary to build the software.

  * `make format` - format the code using clang-format
  * `make check-format` - verify that formatting is correct
  * `make container-build` - build the Docker container used in CI
  * `make container-push` - push the container to Docker Hub
  * `make examples` - run all examples

### Packaging

It is possible to create Debian and RPM packages using the local operating
system or a Docker container.

```
# Build the Debian builder image
make deb-container-image

# Build the Debian package in the builder container
make deb-container

# Build the Debian package with the local operating system
make deb
```

Substitute `deb` for `rpm` in the commands above to create an RPM package
instead.

## License

BSD 3-clause license. See LICENSE file.

## Authors

See `AUTHORS` file for up-to-date authors list.

  - Daniel Olson <daniel.olson@umontana.edu>
  - Travis Wheeler <travis.wheeler@umontana.edu>
  - George Lesica <george.lesica@umontana.edu>

