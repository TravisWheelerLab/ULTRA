# ULTRA

**U**LTRA **L**ocates **T**andemly **R**epetitive **A**reas
## About
ULTRA is a tool for finding and annotating tandem repeats within genomic sequence. Model details and evaluation can be found in our release paper, [ULTRA-Effective Labeling of Tandem Repeats in Genomic Sequence](https://www.biorxiv.org/content/10.1101/2024.06.03.597269v1)

## Building
ULTRA requires a compiler supporting C++11 or higher and CMake 3.12 or higher.  To download and build ULTRA run the following commands:
```
git clone https://github.com/TravisWheelerLab/ULTRA
cd ULTRA
cmake .
make
```
## Basic usage
A list of all flags and options can be seen with `ultra -h`. To annotate tandem repeats with ULTRA use `ultra [options] <path to FASTA file>`.  The following examples demonstrate common use cases.
  
<details>
<summary>Example 1 - Default settings</summary>

`examples/example_1.fa` contains randomly generated sequence with three inserted tandem repeats. We can use ULTRA to annotate the file by running:
```
ultra examples/example_1.fa
```

Running the above command will cause ULTRA to output (to stdout) the settings being used followed by a TSV formatted annotation of the repeats found within `examples/example_1.fa`.
ULTRA annotations can be directed to an output file using the `-o <output file path>` option. For example:

```
ultra -o examples/example_1_ultra.tsv examples/at_repeat.fa
```

Running this command will result in two files, `examples/example_1_ultra.tsv` (containing the repeat annotations) and `examples/example_1_ultra.tsv.settings` (containing the ULTRA settings used). ULTRA's default TSV format consists of a header row followed by one row per annotated tandem repeat. Looking inside `examples/example_1_ultra.tsv` we find:
```
SeqID Start End Period  Score Consensus #Subrepeats SubrepeatStarts SubrepeatConsensi
Rand_seq_1  447 592 2 175.659958  AT  2 0,81  AT,CG
Rand_seq_1  859 920 4 49.729370 AGGC  1 0 AGGC
Rand_seq_1  1104  1167  3 71.850250 AGT 1 0 AGT

```
Columns 1-6 (`SeqID`, `Start`, `End`, `Period`, `Score`, and `Consensus`) describe the overall repetitive region. Columns 7-9 (`#Subrepeats`, `SubrepeatStarts`, and `SubrepeatConsensi`) describe changes in repetitive pattern. Running ULTRA with the `--show_seq` flag will result in an additional column, `Sequence`, that contains the complete repetitive sequence.

The second and third repeats (AGGC and AGT respectively) have no changes in repetitive pattern and so their `#Subrepeats` field is "1". The first repeat contains two unique repetitive patterns though, and so its `#Subrepeats` field is "2". The `SubrepeatStarts` field contains a comma separated list of start locations (relative to the repeat `Start` field) and the `SubrepeatConsensi` field contains a comma separated list of subrepeat consensus patterns.  When a letter cannot be resolved for a subrepeat consensus pattern the letter will be displayed as "\*".

ULTRA can be used to mask repeats with the `--mask <masked output path>` option. For example:

```
ultra --mask examples/example_1_masked.fa -o examples/example_1_ultra.tsv examples/at_repeat.fa
```
By default ULTRA will use lower-case masking, although ULTRA will use N-masking if supplied with the `--nmask` flag.

</details>

<details>
<summary>Example 2 - Large period repeats</summary>
`examples/example_2.fa` contains a period 1000 repeat, which is larger than ULTRA's default maximum detectable repetitive period (100). To find the large period repeat we must adjust ULTRA's maximum detectable repetitive period using the `-p <max repeat period>` option.

`ultra -o examples/example_2_ultra.tsv -p 1000 examples/example_2.fa`

  After running the above command `examples/example_2_ultra.tsv` will contain the following output:

```
SeqID Start End Period  Score Consensus #Subrepeats SubrepeatStarts SubrepeatConsensi
period_1000_repeat  0 17999 1000  22938.433594  AAAATCCATACCGCTCATTCACCAGGCTGCGAAGCCTACACTGGTATATGAATCCGAGCTGGAGCAGGGCCCCTAAAATTCGGAGTCGTAGGTGCTCAATACTCCAATCGGTTTTCTCGTGCACCACCGCGGGTGGCTGACAGGGGTTTGACATCGAGAGGCAAGGCAGTTCCGGGCTGAAAGTAGCGCCGGGTAAGGTACGCGCCCGGTATGGCGGGGCCATGGGGCCAATACAGAGGCTGCGCCCTCACTCGGGTGGACGGAAACGCAGAACTATGGTTACTCCTTGGATACGTGAAACGTGTCCCACGGTAGCCCAAGGACTCGGGAGTCTATCACCCCTAGGGCCCATTCCCGGATATAGACGCCAGGTTGAATCCGCATTTGGAGGTACGGTGGATCAGTCTGGGTGGGGCGCGCCCCACTTATACCCTGCGCAGGCTGGACCGAGGGCCGCAAGATGCGACGGTGCACAAGTAGTTGACGACAGACCGTCGTGTTTTCATTACGGTACCAGGATCTTCGGGCCGAGTCAATCAAGCTCGGATTGCGGTGTCTACCGTTAGATCATACCCAACGCCGCAGAGGTGACACGGCGCCGATGGGTACCGGACTTTGGGTCGGCCGCAGTTCGGCAGGGGAGAGGCCCTGCGGCGCGCTTCACTCTGTATGTGCAACGTGCCCAAGTGGCGCCAGGCAGGTCTCAGCCGGTTCCTGCGTCAGCTCGAGGCTGGGCGCGGGAGCTGATCGAACATGGGCCGGGGGCCTCGAACCGTCGAGGACCCCATAGTACCCCGGGCACCAAGTAGGGCAGCCTATAGCCTGAAGCAGTACCGTTTCAGGGGGGGAGCCCTCATGGTCTCCTCTACTGATGACTCAACACGCCAGGGGCGTGAAGCCGGTTCCTTCGGTGGTTATAGATCAAAGGCTCAGAGTGCGGTCTGGAGCGCCCATCTAGCGGCACGCGTCTCGATTGCTCGGTCGCCCTTCACACTCCGCG  1 0 AAAATCCATACCGCTCATTCACCAGGCTGCGAAGCCTACACTGGTATATGAATCCGAGCTGGAGCAGGGCCCCTAAAATTCGGAGTCGTAGGTGCTCAATACTCCAATCGGTTTTCTCGTGCACCACCGCGGGTGGCTGACAGGGGTTTGACATCGAGAGGCAAGGCAGTTCCGGGCTGAAAGTAGCGCCGGGTAAGGTACGCGCCCGGTATGGCGGGGCCATGGGGCCAATACAGAGGCTGCGCCCTCACTCGGGTGGACGGAAACGCAGAACTATGGTTACTCCTTGGATACGTGAAACGTGTCCCACGGTAGCCCAAGGACTCGGGAGTCTATCACCCCTAGGGCCCATTCCCGGATATAGACGCCAGGTTGAATCCGCATTTGGAGGTACGGTGGATCAGTCTGGGTGGGGCGCGCCCCACTTATACCCTGCGCAGGCTGGACCGAGGGCCGCAAGATGCGACGGTGCACAAGTAGTTGACGACAGACCGTCGTGTTTTCATTACGGTACCAGGATCTTCGGGCCGAGTCAATCAAGCTCGGATTGCGGTGTCTACCGTTAGATCATACCCAACGCCGCAGAGGTGACACGGCGCCGATGGGTACCGGACTTTGGGTCGGCCGCAGTTCGGCAGGGGAGAGGCCCTGCGGCGCGCTTCACTCTGTATGTGCAACGTGCCCAAGTGGCGCCAGGCAGGTCTCAGCCGGTTCCTGCGTCAGCTCGAGGCTGGGCGCGGGAGCTGATCGAACATGGGCCGGGGGCCTCGAACCGTCGAGGACCCCATAGTACCCCGGGCACCAAGTAGGGCAGCCTATAGCCTGAAGCAGTACCGTTTCAGGGGGGGAGCCCTCATGGTCTCCTCTACTGATGACTCAACACGCCAGGGGCGTGAAGCCGGTTCCTTCGGTGGTTATAGATCAAAGGCTCAGAGTGCGGTCTGGAGCGCCCATCTAGCGGCACGCGTCTCGATTGCTCGGTCGCCCTTCACACTCCGCG
```
The large consensus can be cumbersome and some users may prefer ULTRA's output to not include consensus patterns for large repeats. Consensus patterns can be limited to smaller period repeats using the `--max_consensus <repetitive period>` option. For example:

```
ultra -o examples/example_2_ultra.tsv -p 1000 --max_consensus 10 examples/example_2.fa
```

This command will result in the following output:
```
SeqID Start End Period  Score Consensus #Subrepeats SubrepeatStarts SubrepeatConsensi
period_1000_repeat  0 17999 1000  22938.433594  . 1 0 .
```
</details>
<details>
<summary>Example 3 - Tuning and FDR</summary>
`examples/example_3.fa` contains randomly generated 80% AT rich sequence along with two inserted tandem repeats (an "AAAGC" repeat and an "AAAATAC" repeat). The large AT bias is far outside ULTRA's default expectation, and as a result ULTRA will have a high false discovery rate, as seen by running:
```
ultra --fdr -o examples/example_3_def.tsv examples/example_3.fa
```
After this command runs, ULTRA will print to standard out: `Estimated false discovery rate: 0.576698`, meaning that ULTRA expects 58% of the repeats annotated in `examples/example_3_def.tsv` to be a result of sequence bias.  ULTRA uses random sequence shuffling to estimate false discovery rate (see the tuning section of [our paper](https://www.biorxiv.org/content/10.1101/2024.06.03.597269v1)), and as a result the exact number you see will be different. 

We can improve our results by automatically tuning ULTRA's parameters using the `--tune` flag:
```
ultra --tune --tune_indel --fdr -o examples/example_3_tuned.tsv examples/example_3.fa
```
`--tune` causes ULTRA to perform a parameter search that maximizes coverage while keeping the estimated false discovery rate under some threshold (by default 0.05  and adjustable via `--tune_fdr <fdr threshold>`). By default ULTRA will tune itself using a lighter-weight repeat model (one without indel states) in an attempt to decrease tuning time. The `--tune_indel` flag will enable the full ULTRA model and result in the highest quality tuning results. ULTRA's tuned annotation, `examples/example_3_tuned.tsv`, contains only the inserted "AAAGC" and "AAAATAC" repeats:
```
SeqID Start End Period  Score Consensus #Subrepeats SubrepeatStarts SubrepeatConsensi
80_AT 4404  4535  5 132.823410  AAAGC 1 0 AAAGC
80_AT 8406  9134  7 716.489990  AAAATAC 1 0 AAAATAC
```
</details>

## Output formats and tuning guide
<details>
<summary>ULTRA TSV format</summary>

ULTRA's default output format is a tab-separated-values format that includes a descriptive header. Default settings will result in the following columns:
```
SeqID Start End Period  Score Consensus #Subrepeats SubrepeatStarts SubrepeatConsensi
```
`SeqID` describes the SequenceID, `Start` describes the repeat starting location (starting with 0), `End` describes the ending location (total repeat length is `End` - `Start`), `Period` describes the repetitive period, `Score` gives the ULTRA score, `Consensus` gives the repetitive consensus, `#Subrepeats` describes the number of repetitive patterns found within the region, `SubrepeatStarts` contains a comma separated list of where subrepeats begin (the first will always be 0), `SubrepeatConsensi` contains a comma separated list of subrepeat consensus patterns. 

Using `--pval` will cause the `Score` column to change to "`Score,PValue`".
Using `--max_split -1` will remove the columns `#Subrepeats SubrepeatStarts SubrepeatConsensi` columns.
Using `--max_consensus -1` will remove the columns `Consensus SubrepeatConsensi`.
Using `--show_seq` will include an additional `Sequence` column containing the complete repetitive sequence.
</details>
<details>
<summary>JSON format</summary>
Using ULTRA with `--json` will result in JSON formatted output. The JSON output will contain an objects with a `Repeats` array. Each object in the `Repeats` array will contain descriptive fields such as "Start", "Length",  "Period", "Score",  "Substitutions" (the number of mismatches), "Insertions" (the number of insertions), "Deletions" (the number of deletions), "Consensus", and additional fields depending on the specific settings being used. 

Repeats that contain subrepeats will have a "Subrepeats" array, each object in the array containing a "Start" field (that describes the subrepeat's starting location relative to the overall repeat), and a "Consensus" field.
</details>
<details>
<summary>BED format</summary>
Using ULTRA with `--bed`will result in a BED file with four columns (`Sequence ID`, `Start`, `End`,  `Consensus`).
</details>
<details>
<summary>Multiformat output</summary>
When using the `-o <file path>` option for saving output it is possible to provide ULTRA with multiple output formats. For example:
```
ultra  --tsv --json --bed -o examples/example_1_multi examples/at_repeat.fa
```
Running the above command will result in four output files: 
```
examples/example_1_multi.settings
examples/example_1_multi.tsv
examples/example_1_multi.json
examples/example_1_multi.bed
```
Note that when using multiple output formats ULTRA will automatically choose the extension based on the output format type.
</details>

<details>
  <summary>Tuning guide</summary>
ULTRA's automatic parameter tuning (via `--tune`) can greatly improve annotation performance by testing several parameter sets and then using the parameter set that achieved the highest coverage within some threshold of estimated false discovery rate (by default 0.05, and adjusted with `--tune_fdr <fdr threshold>`) . By default `--tune` will test 18 parameter sets, `--tune_medium` will test 40 parameter sets, and `--tune_large` will test 252 parameter sets.  Each parameter contains different emission probabilities and transition probabilities, but does not affect the repeat period. By default tuning will disable indel states in order to decrease the overall runtime. Best tuning performance will be achieved by tuning with indel states, using the `--tune_indel` flag. 

Users can also specify their own parameter sets to tune ULTRA against by using the `--tune_file <tune file path>` option. Each line in the tune file should contain arguments to run ULTRA against. An example can be seen in `examples/tune_file` which contains:
```
-p 1
-p 4
-p 16
-p 64
-p 256
-p 10 -m 0.5 --at 0.3
-p 10 -m 0.55 --at 0.2
-p 10 -m 0.5 --at 0.1 -i 2 -d 2
```
Here `examples/tune_file` tests different repeat periods (`-p <repeat period>`), different repeat emission probabilities (`-m <repeat emission probability>`),  different background AT frequencies (`--at <AT rate>`), and different indel states (`-i <max insertions> -d <max deletions>`).
</details>

## Citing
```
@article {Olson2024ultra,
  author = {Olson, Daniel R. and Wheeler, Travis J.},
  title = {ULTRA-Effective Labeling of Repetitive Genomic Sequence},
  elocation-id = {2024.06.03.597269},
  year = {2024},
  doi = {10.1101/2024.06.03.597269},
  publisher = {Cold Spring Harbor Laboratory},
  URL = {https://www.biorxiv.org/content/early/2024/06/04/2024.06.03.597269},
  eprint = {https://www.biorxiv.org/content/early/2024/06/04/2024.06.03.597269.full.pdf},
  journal = {bioRxiv}
}
```
