from utility import *

MAX_NUM_READS = 100 * 1000
FROM_POS = 100 * 1000
TO_POS = 20 * 1000 * 1000  # 200 * MAX_NUM_READS
INTERSECT_FILENAME = 'data/intersect.vcf'
ALL_VARIANTS_FILENAME = 'data/4.0.final_genotypes.vcf'
# ALL_VARIANTS_FILENAME = 'data/sample.vcf'
LONGSHOT_OUTPUT_FILENAME = 'data/out.vcf'
FRAGMENT_FILENAME = 'data/fragments.txt'
# FRAGMENT_FILENAME = 'data/sample.fragments'
EDGES_LIMIT = 100
GROUND_TRUTH_FILENAME = 'data/chr_20_HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf'
# GROUND_TRUTH_FILENAME = 'data/sample.vcf'
MIN_HETERO_READS = 3
LOCAL_SERACH_THRESHOLD = 1
ITERATION_NUM = 50
ITERATION_RANDOM_CHANGE = .1
ITERATION_RANDOM_SWAP = 0.1
GREEDY_THRESHOLD_VALUE = -100
METHOD = 'fisher'  # fisher or barnard
MIN_SUM_EDGE = 8
SIGNIFICANT_THRESHOLD = .10


vcf_dict, flist = read_fragment_matrix(FRAGMENT_FILENAME, ALL_VARIANTS_FILENAME, MAX_NUM_READS, FROM_POS, TO_POS)
# vcf_dict, flist = keep_bed_intersect(vcf_dict, flist, INTERSECT_FILENAME, FROM_POS, TO_POS)
vcf_dict, flist = keep_output_intersect(vcf_dict, flist, LONGSHOT_OUTPUT_FILENAME, FROM_POS, TO_POS)
variant_count = get_variant_count(flist)
homo_variants_pos, v_set = vertex_filter(vcf_dict, variant_count, MIN_HETERO_READS)
with open('variant_count.txt', 'w') as file:
    file.write(str(variant_count))

edges = build_edges(flist, v_set, EDGES_LIMIT)
merged_edges = merge_edge(edges, METHOD, MIN_SUM_EDGE, SIGNIFICANT_THRESHOLD)
# with open('edges.txt', 'w') as file:
#     file.write(str(edges))

# incorrect_variants, output_true_variants_set = greedy_algorithm(merged_edges.copy(), GREEDY_THRESHOLD_VALUE)
output_variants_pos = iterative_local_search(vcf_dict, merged_edges.copy(), LOCAL_SERACH_THRESHOLD, ITERATION_NUM, ITERATION_RANDOM_CHANGE, ITERATION_RANDOM_SWAP, None)
# with open('incorrect_variants_pos.txt', 'w') as file:
#     file.write(str(incorrect_variants_pos))

hetero_truth_pos = dict_to_set(read_vcf_file(GROUND_TRUTH_FILENAME, ['chr20'], from_pos=FROM_POS, to_pos=TO_POS, return_homo=False))
homo_truth_pos = dict_to_set(read_vcf_file(GROUND_TRUTH_FILENAME, ['chr20'], from_pos=FROM_POS, to_pos=TO_POS, return_hetero=False))
known_pos = dict_to_set(read_vcf_file(INTERSECT_FILENAME, ['chr20'], from_pos=FROM_POS, to_pos=TO_POS))

accuracy(merged_edges, output_variants_pos, homo_variants_pos, known_pos, hetero_truth_pos, homo_truth_pos, vcf_dict)
