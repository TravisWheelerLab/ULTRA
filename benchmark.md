# M1 Benchmarks

Benchmarks on an first-gen M1 MacBook Air.

ARM:

```
ULTRA on  prerelease [?]
➜ make

ULTRA on  prerelease [?] took 48s
➜ file ultra
ultra: Mach-O 64-bit executable arm64

ULTRA on  prerelease [?]
➜ time ./ultra -ss chr21.fa > ultra_data.txt
./ultra -ss chr21.fa > ultra_data.txt  46.42s user 0.15s system 99% cpu 46.600 total
```

Intel emulation (Rosetta 2):

```
ULTRA on  prerelease [?] took 47s
➜ make clean && arch -x86_64 make

ULTRA on  prerelease [?] took 8s
➜ file ultra
ultra: Mach-O 64-bit executable x86_64

ULTRA on  prerelease [?]
➜ time ./ultra -ss chr21.fa > ultra_data.txt
./ultra -ss chr21.fa > ultra_data.txt  74.10s user 0.25s system 99% cpu 1:14.91 total
```

Intel native on 2020 MBA:

```
ULTRA on  benchmarking [$]
➜ make

ULTRA on  benchmarking [$] took 11s
➜ file ultra
ultra: Mach-O 64-bit executable x86_64

ULTRA on  benchmarking [$]
➜ time ./ultra -ss chr21.fa > ultra_data.txt
./ultra -ss chr21.fa > ultra_data.txt  81.20s user 0.14s system 99% cpu 1:22.09 total
```

