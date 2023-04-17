# ULTRA

**U**LTRA **L**ocates **T**andemly **R**epetitive **A**reas

(Note: The current stable release is 1.0.0 beta. As we move out of beta we will update the README to be more informative.)

ULTRA is a tool to find and annotate tandem repeats inside genomic sequence. It
is able to find repeats of any length and of any period (up to a maximum period
of 4000). It can find highly decayed repeats missed by other software, and it
will also be able to find very large repeats in highly repetitive sequence,
regardless of the size of sequence or length of repeats. ULTRA offers meaningful
annotation scores and can produce annotation P-values at user request.

## Install

ULTRA uses the CMake build system (version 3.12 or greater). This makes building
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

If you want to direct the output to a file you can do so with the -o flag:

```
ultra -o ex1_output.bed examples/ex1.fa
```

Now your current working directory will contain an `ex1_output.bed` file. At
the top of this file will be a complete description of the parameters used to
run ULTRA, and then there will be a list of repeats found, and their associated
annotation information. Output is in JSON format.

By default, ULTRA has a maximum detectable repeat period of 15. The sequence in
`examples/ex_large.fa` has a period 1000 repeat with a 10% substitution rate. To
detect larger period repeats, we will use the -p argument.

```
ultra -p 1000 --json -o ex_large.json examples/ex_large.fa
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

To see a full list of ULTRA arguments you can use ‘ultra -h’.

## License

BSD 3-clause license. See LICENSE file.

## Authors

See `AUTHORS` file for up-to-date authors list.

  - Daniel Olson <daniel.olson@umontana.edu>
  - Travis Wheeler <travis.wheeler@umontana.edu>
  - George Lesica <george.lesica@umontana.edu>

