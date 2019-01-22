def read_fragment_matrix(frag_matrix, vcf_file):

    snp_ix = 0
    vcf_dict = dict()
    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            genomic_pos = int(el[1])-1
            vcf_dict[snp_ix] = genomic_pos
            snp_ix += 1

    flist = []

    with open(frag_matrix,"r") as fm:
        cnt = 0
        cnt_quality = {}
        for line in fm:
            cnt += 1
            if cnt > 8 * 1000:
                break
            if len(line) < 2:
                continue

            el = line.strip().split()

            num_blks = int(el[0])
            name = el[1]

            call_list = el[2:(2+2*num_blks)]              # extract base call part of line
            call_list = zip(*[iter(call_list)]*2)             # list -> tuple list conversion: credit to http://stackoverflow.com/questions/23286254/convert-list-to-a-list-of-tuples-python
            call_list = [(int(a)-1, b) for a,b in call_list]  # convert index to 0-based integer
            call_list2 = []

            for ix, blk in call_list:
                curr_ix = ix
                for a in blk:
                    call_list2.append((curr_ix, a))
                    curr_ix += 1

            qlist = el[-1]
            for q in qlist:
                x = (ord(q) - 33)
                if x not in cnt_quality:
                    cnt_quality[x] = 0
                cnt_quality[x] += 1
            qlist = [10**((ord(q) - 33) * -0.1) for q in qlist]
            alist = [(a,b,c) for ((a,b),c) in zip(call_list2,qlist)]
            flist.append(alist)
    print(cnt_quality)
    print(flist[-1])
    return vcf_dict, flist


def build_edges(flist, v_set):
    edges = {}
    cnt = 0
    for fragment in flist:
        cnt += 1
        print(cnt, len(fragment))
        for i in range(len(fragment)):
            v = fragment[i][0]
            if v not in v_set:
                continue
            if edges.get(v) is None:
                edges[v] = {}
            for j in range(len(fragment)):
                u = fragment[j][0]
                if u not in v_set:
                    continue
                if u == v:
                    continue
                if edges[v].get(u) is None:
                    edges[v][u] = [0, 0, 0, 0]
                idx = int(fragment[i][1])*2 + int(fragment[j][1])
                edges[v][u][idx] += 1
    return edges


def get_variant_count(flist):
    print(flist[-1][-1][0])
    variant_count = [[0, 0] for _ in range(flist[-1][-1][0]+1000)]
    for fragment in flist:
        for i in range(len(fragment)):
            print(fragment[i][0], int(fragment[i][1]))
            variant_count[fragment[i][0]][int(fragment[i][1])] += 1
    return variant_count


def vertex_filter(variant_count):
    filtered_v = set()
    incorrect_variants = set()
    for i in range(len(variant_count)):
        if variant_count[i][0] > 2 and variant_count[i][1] > 2:
            filtered_v.add(i)
        elif variant_count[i][0] > 2:
            incorrect_variants.add(i)
    return incorrect_variants, filtered_v


# def make_haplotype(variant_count, edges):
#     haplotype = [{}, {}]
#     for v in edges:
#         state
#         for u in edges[v]:
#             if u >= v:
#                 break
#             if haplotype[0][u]
#             if edges[v][u]

def merge_edge(edges):
    merged_edges = {}
    for v in edges:
        merged_edges[v] = {}
        for u in edges[v]:
            t = edges[v][u]
            # t.sort()
            # merged_edges[v][u] = (t[0]+t[1], t[2]+t[3])
            merged_edges[v][u] = (t[0]*t[3] + t[1]*t[2], t[0]*t[1] + t[2]*t[3])
    return merged_edges


def greedy_algorithm(incorrect_variants, merged_edges):
    variant_set = set()
    variant_score = {}
    for v in merged_edges:
        good_edges, bad_edges = 0, 0
        for u in merged_edges[v]:
            good_edges += merged_edges[v][u][0]
            bad_edges += merged_edges[v][u][1]
        variant_set.add((good_edges - bad_edges, v))
        variant_score[v] = good_edges - bad_edges

    while len(variant_set) > 0:
        v = min(variant_set)
        variant_set.remove(v)
        print(v, len(variant_set))
        if v[0] < -100:
            incorrect_variants.add(v[1])
            for u in merged_edges[v[1]]:
                variant_set.remove((variant_score[u], u))
                variant_score[u] -= merged_edges[u][v[1]][0] - merged_edges[u][v[1]][1]
                variant_set.add((variant_score[u], u))
                merged_edges[u].pop(v[1])
    return incorrect_variants


def accuracy(vcf_file, all_variants_pos, incorrect_variants_pos, merged_edges, edges, all_variants_pos_to_idx):
    vcf_pos = set()
    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue
            if el[0] != '20':
                continue
            genomic_pos = int(el[1])-1
            vcf_pos.add(genomic_pos)
    vcf_pos = set([x for x in vcf_pos if x <= max(incorrect_variants_pos)])
    incorrect_intersect = incorrect_variants_pos.intersection(vcf_pos)
    all_variants_intersect = all_variants_pos.intersection(vcf_pos)
    for pos in incorrect_intersect:
        v = all_variants_pos_to_idx[pos]
        good_edges, bad_edges = 0, 0
        for u in merged_edges[v]:
            good_edges += merged_edges[v][u][0]
            bad_edges += merged_edges[v][u][1]
        print(v, good_edges, bad_edges)
    print(len(incorrect_intersect), len(all_variants_intersect), len(vcf_pos), len(incorrect_variants_pos), len(all_variants_pos))


vcf_dict, flist = read_fragment_matrix('data/fragments.txt', 'data/out.vcf')
variant_count = get_variant_count(flist)
incorrect_variants, v_set = vertex_filter(variant_count)
with open('variant_count.txt', 'w') as file:
    file.write(str(variant_count))
incorrect_variants = set()
edges = build_edges(flist, v_set)
merged_edges = merge_edge(edges)
# edges = filtering_edges(flist)
with open('edges.txt', 'w') as file:
    file.write(str(edges))

incorrect_variants = greedy_algorithm(incorrect_variants, merged_edges.copy())
incorrect_variants_pos = set()
for idx in incorrect_variants:
    incorrect_variants_pos.add(vcf_dict[idx])
with open('incorrect_variants_pos.txt', 'w') as file:
    file.write(str(incorrect_variants_pos))
all_variants_pos = set([x for x in vcf_dict.values() if x <= max(incorrect_variants_pos)])
all_variants_pos_to_idx = {pos: idx for idx, pos in vcf_dict.items()}
accuracy('data/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf', all_variants_pos, incorrect_variants_pos, merged_edges, edges, all_variants_pos_to_idx)
#
# def parse_vcf_phase(vcf_file, CHROM, indels = False):
#
#     #block = []
#     PS_index = None
#     blocks = defaultdict(list)
#
#     with open(vcf_file, 'r') as vcf:
#
#         snp_ix = 0
#
#         for line in vcf:
#             if line[0] == '#':
#                 continue
#
#             el = line.strip().split('\t')
#             if len(el) < 10:
#                 continue
#             if len(el) != 10:
#                 print("VCF file must be single-sample.")
#                 exit(1)
#
#             consider = True
#
#             phase_data = el[9]
#
#             a0 = el[3]
#             a1 = el[4]
#             a2 = None
#
#             if ',' in a1:
#                 alt_lst = a1.split(',')
#                 if len(alt_lst) == 2:
#                     a1,a2 = alt_lst
#                 else:
#                     consider = False
#
#             # get the index where the PS information is
#             for i,f in enumerate(el[8].split(':')):
#                 if i == 0:
#                     assert(f == 'GT')
#                 if f == 'PS':
#                     if snp_ix == 0:
#                         PS_index = i
#                     else:
#                         assert(PS_index == i)
#                     break
#
#             dat = el[9].split(':')
#             genotype = dat[0]
#
#             if not (len(genotype) == 3 and genotype[0] in ['0','1','2'] and
#                     genotype[1] in ['|'] and genotype[2] in ['0','1','2']):
#                 consider = False
#
#             if genotype[0] == genotype[2]:
#                 consider = False
#
#             if consider and (not indels) and (('0' in genotype and len(a0) != 1) or
#                 ('1' in genotype and len(a1) != 1) or ('2' in genotype and len(a2) != 1)):
#                 consider = False
#
#             ps = None
#             if consider and PS_index != None and len(dat) > PS_index:
#                 ps = dat[PS_index]
#                 if ps == '.':
#                     consider = False
#             elif PS_index == None:
#                 ps = "1" # just put everytthing in one block
#             else:
#                 ps = None
#
#             chrom = el[0]
#
#             if chrom != CHROM:
#                 print("ERROR: Chromosome in reference haplotype VCF doesn't match chromosome in VCF used for phasing")
#                 print("reference haplotype VCF: " + vcf_file)
#                 print("{} != {}".format(CHROM, chrom))
#                 exit(1)
#
#             pos = int(el[1])-1
#             if ps != None and consider and phase_data[1] == '|':
#                 blocks[ps].append((snp_ix, pos, phase_data[0:1], phase_data[2:3], a0, a1, a2))
#
#             snp_ix += 1
#
#     return [v for k,v in sorted(list(blocks.items())) if len(v) > 1]
