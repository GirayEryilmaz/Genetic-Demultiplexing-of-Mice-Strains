#!/usr/bin/env python

# MIT License

# Copyright (c) 2024 drneavin

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
import demuxalot
from demuxalot import Demultiplexer, BarcodeHandler, ProbabilisticGenotypes, count_snps
from demuxalot.cellranger_specific import parse_read
import pandas as pd
import argparse
import subprocess
import os

__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

print("loaded packages")

parser = argparse.ArgumentParser(
    description="wrapper for demuxalot for demultiplexing of multiplexed single-cell data.")
parser.add_argument("-b", "--barcodes", required = True, help = "Input barcodes file")
parser.add_argument("-a", "--bamfile", required = True, help = "Input bam file")
parser.add_argument("-n", "--indiv_file", required = True, default = None, help = "File containing individuals in pool separated by line. No header.")
parser.add_argument("-v", "--vcf", required = True, default = None, help = "The vcf file that has the SNP genotypes for each donor in the pool.")
parser.add_argument("-o", "--outdir", required = True, default = None, help = "The output directory.")
parser.add_argument("-r", "--refine", required = True, default = True, help = "Whether to run genotype refinement.")
parser.add_argument("-c", "--celltag", required = False, default = "CB", help = "optional: SAM tag used for cell barcodes; default: `CB`. ")
parser.add_argument("-u", "--umitag", required = False, default = "UB", help = "optional: SAM tag used for UMIs; default: `UB`.")
parser.add_argument("-p", "--nproc", required = False, type = int, default = 1, help = "(max.) number of parallel jobs; -1: number of CPU cores/threads; default: 1.")
args = parser.parse_args()

print("read arguments")


### Load in file with individual IDs as list
with open(args.indiv_file) as f:
    names = f.read().splitlines()

names_ordered = np.sort(names)

# Load genotypes
print("loading genotypes")
genotypes = ProbabilisticGenotypes(genotype_names=names_ordered)
genotypes.add_vcf(args.vcf)


# Load barcodes
print("loading barcodes")
barcode_handler = BarcodeHandler.from_file(args.barcodes, tag=args.celltag)

parse_read_custom = lambda read: parse_read(read, umi_tag = args.umitag)
snps = count_snps(
    bamfile_location = args.bamfile,
    chromosome2positions=genotypes.get_chromosome2positions(),
    barcode_handler=barcode_handler, 
    joblib_n_jobs=args.nproc,
    parse_read=parse_read_custom,
)

print("estimating posterior probs and likelihoods")
# returns two dataframes with likelihoods and posterior probabilities 
likelihoods, posterior_probabilities = Demultiplexer.predict_posteriors(
    snps,
    genotypes=genotypes,
    barcode_handler=barcode_handler,
)

print("writing unrefined results")
# Write output
likelihoods.to_csv(args.outdir + "/likelihoods.tsv.gz", sep = "\t", index = True)
posterior_probabilities.to_csv(args.outdir + "/posterior_probabilities.tsv.gz", sep = "\t", index = True)
posterior_probabilities.idxmax(axis = 1).to_csv(args.outdir + "/assignments.tsv.gz", sep = "\t", index = True)

if args.refine:

    print("inferring refined genotypes")
    # Infer refined genotypes 
    refined_genotypes, _posterior_probabilities = Demultiplexer.learn_genotypes(snps, genotypes=genotypes, n_iterations=5, barcode_handler=barcode_handler)

    print("estimating likelihoods and probabilities")
    # Use learnt genotypes for demultiplexing
    likelihoods, posterior_probabilities = Demultiplexer.predict_posteriors(
        snps,
        genotypes=refined_genotypes,
        barcode_handler=barcode_handler,
    )

    print("writing second set of output")
    # Write output
    likelihoods.to_csv(args.outdir + "/likelihoods_refined.tsv.gz", sep = "\t", index = True)
    posterior_probabilities.to_csv(args.outdir + "/posterior_probabilities_refined.tsv.gz", sep = "\t", index = True)
    posterior_probabilities.idxmax(axis = 1).to_csv(args.outdir + "/assignments_refined.tsv.gz", sep = "\t", index = True)

print("Done!")