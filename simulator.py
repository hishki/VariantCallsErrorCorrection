 #! /usr/bin/python
# -*- coding: utf-8 -*-

'''
Created on Mon Jul  6 11:17:22 2015
@author: Peter Edge
'''

# imports
import math
import numpy as np
import os
import random

# the core function of the simulation program
# if insert_size set > 0, then read_length refers to total length of both mates

def simulate_haplotypes(frag_file, vcf_file, read_length, coverage, ref_length=int(2.5e8), miscall_rate=0.02, missing_rate=0, per_frag_sw_err=0, snp_rate=0.0008, read_length_by_variants=False, span_by_variants=False, coverage_by_variants=False, hic_span = 0, ref_name='chr20'):

    print("SIMULATING HAPLOTYPES...".format(coverage))

    # if "use_variants" is True, then read_length and insert_size refer to read length in variants,
    # also, reads will be created until the average coverage PER VARIANT is >= coverage,
    # whereas the normal behavior is that coverage refers to average coverage of reads over base pairs
    read_length_bak = read_length

    read_length = int(read_length) if not read_length_by_variants else int((read_length / (1-missing_rate) / (1-2/3*miscall_rate)) / snp_rate)

    mate_length = read_length

    if hic_span > 0:
        read_length = int(hic_span) if not span_by_variants else int(hic_span / snp_rate)
        num_reads_for_read_coverage = int(coverage * ref_length / ((2*mate_length)*(1-missing_rate)*(1-2/3*miscall_rate)))
    else:
        # calculate number of reads based on coverage
        # account for the fact that we'll have missing calls from missing and miscall rate
        num_reads_for_read_coverage = int(coverage * ref_length / ((read_length)*(1-missing_rate)*(1-2/3*miscall_rate)))


    # determine positions of SNPs
    snp_logical = np.less(np.random.random_sample(ref_length), snp_rate)   # logical array specifying occurence of SNPs or not for length of reference
    snp_ix = np.where(snp_logical)[0]                                      # indices of SNPs in reference
    num_snps = np.sum(snp_logical)                                         # total number of SNPs

    # create random haplotypes of size num_snps
    hap1 = np.less(np.array(np.random.random_sample(num_snps)), 0.5)   # first haplotype is a random 50/50 binary array
    hap2 = np.logical_not(hap1)                                        # second hap is just first hap's complement
    haps = [hap1, hap2]

    q = '~' if miscall_rate < 5.011872336272714e-10 else chr(int(33-10*math.log10(miscall_rate)))

    total_read_variants = 0
    actual_read_count = 0
    insert_spanning_read_count = 0
    total_insert_snps = 0
    reads_produced = 0
    done = False

    with open(frag_file, 'w') as ff:

        # loop and generate fragments
        while True:
            if done:
                break

            if hic_span > 0:
                #curr_read_length = read_length+1
                #while curr_read_length > read_length:
                #    curr_read_length = np.random.exponential(0.25*read_length)
                curr_read_length = random.randrange(2*mate_length+2, read_length)
            else:
                curr_read_length = read_length

            hap_ix = random.randrange(0,2)                                 # which haplotype was fragment generated by
            start_index = random.randrange(0, ref_length-curr_read_length)      # get a random start spot
            frag_snp = np.searchsorted(snp_ix, start_index)                # find corresponding spot in snp array

            fs = []

            read_end_gix = start_index + curr_read_length
            mate1_end_gix   = start_index + mate_length
            mate2_start_gix = read_end_gix - mate_length
            spans_mate = False
            insert_snps = 0

            while(frag_snp < num_snps) and (snp_ix[frag_snp] <= read_end_gix):
                if hic_span == 0 or (snp_ix[frag_snp] < mate1_end_gix or snp_ix[frag_snp] >= mate2_start_gix):
                    a = (frag_snp,str(int(haps[hap_ix][frag_snp])),q)
                    fs.append(a)

                if snp_ix[frag_snp] >= mate1_end_gix and snp_ix[frag_snp] < mate2_start_gix:
                    insert_snps += 1

                if snp_ix[frag_snp] >= mate2_start_gix:
                    spans_mate = True
                frag_snp += 1

            reads_produced += 1

            if not coverage_by_variants:

                interval = int(num_reads_for_read_coverage/10)

                if reads_produced % interval == 0:
                    frac = int(reads_produced / interval)*10
                    print("{}%...".format(frac))

                if reads_produced > num_reads_for_read_coverage:
                    done = True

            if len(fs) < 2:       # noninformative fragment, no point in continuing
                continue

            if per_frag_sw_err > 0:
                new_fs = []
                # add switch errors to the fragment
                alpha = per_frag_sw_err / (len(fs)-1) # per base switch_errors
                switched = False # currently in a "switched" state or not
                for pos, call, qual in fs:
                    if random.random() < alpha:
                        switched = not switched

                    if switched:
                        if call == '1':
                            call = '0'
                        elif call == '0':
                            call = '1'

                    new_fs.append((pos,call,qual))
                fs = new_fs

            if missing_rate > 0:
                new_fs = []
                for a in fs:
                    if random.random() >= missing_rate:
                        new_fs.append(a)
                fs = new_fs

            # randomly add errors to fragment based on miscall rate
            # if base is miscalled to non-ref and non-alternate, convert to '-' (thrown away)
            if miscall_rate > 0:

                new_fs = []

                for pos,call,qual in fs:
                    # miscalled to ref or variant
                    if random.random() < miscall_rate:
                        if random.randrange(0,3) == 0: # 1 in 3 chance of switching to the other base that is a ref or variant
                            if call == '1':
                                new_fs.append((pos,'0',qual))
                            elif call == '0':
                                new_fs.append((pos,'1',qual))
                    else:
                        new_fs.append((pos,call,qual))

                fs = new_fs

            # don't print an empty fragment
            if len(fs) < 2:
                continue

            if spans_mate:
                insert_spanning_read_count += 1
                total_insert_snps += insert_snps

            name = "FRAG{}".format(reads_produced)
            actual_read_count += 1
            print_fragment(name,fs,ff)
            total_read_variants += len(fs)

            if coverage_by_variants and total_read_variants / num_snps >= coverage:

                break

    with open (vcf_file, 'w') as vcf:
        header ='''##fileformat=VCFv4.1
##contig=<ID={},length={}>
##INFO=<ID=mt,Number=1,Type=String,Description="Variant Type: SUBSTITUTE/INSERT/DELETE">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=FP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SIM_INDIVIDUAL'''.format(ref_name,ref_length)

        # print header to vcf file
        print(header, file=vcf)

        # iterate over SNPs and print one line per SNP
        bases = ['A','T','G','C']
        for snp, genome_ix in np.ndenumerate(snp_ix):
            snp = snp[0]
            ref_ix = random.randrange(0,4)
            alt_ix = ref_ix
            while(alt_ix == ref_ix):
                alt_ix = random.randrange(0,4)

            genotype_field = '{}|{}'.format(int(hap1[snp]),int(hap2[snp]))

            ID = '.'
            ref_snp = bases[ref_ix]
            alt_snp = bases[alt_ix]
            qual = 100
            fltr = 'PASS'
            info = 'mt=SUBSTITUTE'
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tGT\t{}'.format(ref_name, genome_ix, ID, ref_snp, alt_snp, qual, fltr, info, genotype_field), file=vcf)

    print("DONE.")
    readspan = read_length if hic_span == 0 else 2*mate_length
    estimated_coverage = reads_produced * (readspan) * (1-missing_rate) * (1-2/3*miscall_rate) / ref_length

    print("ESTIMATED READ COVERAGE:      {}".format(estimated_coverage))
    if hic_span > 0:
        print("TARGET MAX MATE BASE SPAN:    {}".format(hic_span))
        print("MATE LENGTH:                  {}".format(mate_length))
    else:
        print("READ LENGTH:                  {}".format(read_length))

   #print("TARGET COVERAGE (in SNPs):    {}"
    if coverage_by_variants:
        print("TARGET COVERAGE PER SNP:      {}".format(coverage))
    print("MEAN COVERAGE PER SNP:        {}".format(total_read_variants/num_snps))

    if read_length_by_variants:
        print("TARGET SNPs PER READ:         {}".format(read_length_bak))
    print("MEAN SNPs PER READ:           {}".format(total_read_variants/actual_read_count))

    if hic_span > 0:

        if span_by_variants:
            print("TARGET MAX SNPs PER INSERT:   {}".format(hic_span))
        print("MEAN SNPs PER INSERT:         {}".format(total_insert_snps/insert_spanning_read_count))

# make directory, create path if necessary
def mkdir_p(path):
    if not os.path.exists(path):
        os.makedirs(path)

# print out fragment object to fragment matrix file
# frag.seq should be in
def print_fragment(name, fs, ff):
    # build line content

    fragstr = ''
    num_pairs = 0
    prev_snp_ix = -2
    qual = ' '
    for snp_ix, allele, q_char in fs:

        diff = snp_ix - prev_snp_ix

        if diff == 1:
            fragstr += allele
        else:
            num_pairs += 1
            fragstr += ' {} {}'.format(snp_ix+1, allele)

        prev_snp_ix = snp_ix
        qual += q_char

    fragstr += qual

    prefix = '{} {}'.format(num_pairs,name)
    fragstr = prefix + fragstr

    # print line to file
    print(fragstr, file=ff)


if __name__ == '__main__':
    simulate_haplotypes(frag_file='sample.fragments', vcf_file='sample.vcf', read_length=10000, coverage=10)
