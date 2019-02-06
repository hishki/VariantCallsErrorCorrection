from utility import *

MAX_NUM_READS = 8 * 1000
INTERSECT_FILENAME = 'data/intersect.vcf'
ALL_VARIANTS_FILENAME = 'data/4.0.final_genotypes.vcf'
LONGSHOT_OUTPUT_FILENAME = 'data/out.vcf'
FRAGMENT_FILENAME = 'data/fragments.txt'
EDGES_LIMIT = 50
GROUND_TRUTH_FILENAME = 'data/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf'
MIN_HETERO_READS = 2

vcf_dict, flist = read_fragment_matrix(FRAGMENT_FILENAME, ALL_VARIANTS_FILENAME, MAX_NUM_READS)
vcf_dict, flist = keep_bed_intersect(vcf_dict, flist, INTERSECT_FILENAME)
vcf_dict, flist = keep_output_intersect(vcf_dict, flist, LONGSHOT_OUTPUT_FILENAME)
variant_count = get_variant_count(flist)
homo_variants, v_set = vertex_filter(variant_count, MIN_HETERO_READS)
with open('variant_count.txt', 'w') as file:
    file.write(str(variant_count))
incorrect_variants = set()
edges = build_edges(flist, v_set, EDGES_LIMIT)
merged_edges = merge_edge(edges)
with open('edges.txt', 'w') as file:
    file.write(str(edges))

incorrect_variants = greedy_algorithm(incorrect_variants, merged_edges.copy())
incorrect_variants_pos = set()
for idx in incorrect_variants:
    incorrect_variants_pos.add(vcf_dict[idx])
with open('incorrect_variants_pos.txt', 'w') as file:
    file.write(str(incorrect_variants_pos))
accuracy(GROUND_TRUTH_FILENAME, vcf_dict, incorrect_variants_pos, merged_edges, homo_variants, edges)