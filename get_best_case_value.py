from utility import *

MAX_NUM_READS = 100 * 1000
INTERSECT_FILENAME = 'data/intersect.vcf'
ALL_VARIANTS_FILENAME = 'data/4.0.final_genotypes.vcf'
LONGSHOT_OUTPUT_FILENAME = 'data/out.vcf'
# FRAGMENT_FILENAME = 'data/fragments.txt'
FRAGMENT_FILENAME = 'sample.fragments'
EDGES_LIMIT = 1000
# GROUND_TRUTH_FILENAME = 'data/chr_20_HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf'
GROUND_TRUTH_FILENAME = 'data/sample.vcf'
MIN_HETERO_READS = 3
LOCAL_SERACH_THRESHOLD = 1
ITERATION_NUM = 5
ITERATION_RANDOM_CHANGE = 10
GREEDY_THRESHOLD_VALUE = -100

vcf_dict, flist = read_fragment_matrix(FRAGMENT_FILENAME, ALL_VARIANTS_FILENAME, MAX_NUM_READS)
# vcf_dict, flist = keep_bed_intersect(vcf_dict, flist, INTERSECT_FILENAME)
# vcf_dict, flist = keep_output_intersect(vcf_dict, flist, LONGSHOT_OUTPUT_FILENAME)
variant_count = get_variant_count(flist)
homo_variants, v_set = vertex_filter(variant_count, MIN_HETERO_READS)
# with open('variant_count.txt', 'w') as file:
#     file.write(str(variant_count))

edges = build_edges(flist, v_set, EDGES_LIMIT)
merged_edges = merge_edge(edges)
# with open('edges.txt', 'w') as file:
#     file.write(str(edges))

# incorrect_variants, output_true_variants_set = greedy_algorithm(merged_edges.copy(), GREEDY_THRESHOLD_VALUE)
incorrect_variants, output_true_variants_set = iterative_local_search(merged_edges.copy(), LOCAL_SERACH_THRESHOLD, ITERATION_NUM, ITERATION_RANDOM_CHANGE, None)
incorrect_variants_pos = set()
for idx in incorrect_variants:
    incorrect_variants_pos.add(vcf_dict[idx])
# with open('incorrect_variants_pos.txt', 'w') as file:
#     file.write(str(incorrect_variants_pos))
accuracy(GROUND_TRUTH_FILENAME, vcf_dict, incorrect_variants_pos, merged_edges, homo_variants, edges, output_true_variants_set)
