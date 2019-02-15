import numpy as np
from sortedcontainers import SortedList


def read_vcf_file(filename):
    snp_ix = 0
    vcf_dict = dict()
    with open(filename, 'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            genomic_pos = int(el[1]) - 1
            vcf_dict[snp_ix] = genomic_pos
            snp_ix += 1
    return vcf_dict


def read_fragment_matrix(frag_matrix, vcf_file, max_num_reads):
    vcf_dict = read_vcf_file(vcf_file)
    flist = []

    with open(frag_matrix,"r") as fm:
        cnt = 0
        cnt_quality = {}
        for line in fm:
            cnt += 1
            if cnt > max_num_reads:
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
    # print(cnt_quality)
    # print(flist[-1])
    return vcf_dict, flist


def keep_bed_intersect(vcf_dict, flist, intersect_file):
    intersect_pos = set()
    snp_ix = 0
    intersect_dict = dict()
    reverse_intersect_dict = dict()
    with open(intersect_file, 'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue
            genomic_pos = int(el[1])-1
            intersect_dict[snp_ix] = genomic_pos
            reverse_intersect_dict[genomic_pos] = snp_ix
            intersect_pos.add(genomic_pos)
            snp_ix += 1

    new_flist = []
    for fragment in flist:
        new_fragment = []
        for snp in fragment:
            if vcf_dict[snp[0]] in intersect_pos:
                snp = (reverse_intersect_dict[vcf_dict[snp[0]]], snp[1])
                new_fragment.append(snp)
        if len(new_fragment) > 0:
            new_flist.append(new_fragment)

    return intersect_dict, new_flist


def keep_output_intersect(vcf_dict, flist, out_file):
    intersect_pos = set()
    snp_ix = 0
    with open(out_file, 'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue
            genomic_pos = int(el[1])-1
            intersect_pos.add(genomic_pos)
            snp_ix += 1

    new_flist = []
    res_vcf_dict = {}
    for fragment in flist:
        new_fragment = []
        for snp in fragment:
            if vcf_dict[snp[0]] in intersect_pos:
                new_fragment.append(snp)
                res_vcf_dict[snp[0]] = vcf_dict[snp[0]]
        if len(new_fragment) > 0:
            new_flist.append(new_fragment)
    return res_vcf_dict, new_flist


def build_edges(flist, v_set, edges_limit):
    edges = {}
    cnt = 0
    for fragment in flist:
        cnt += 1
        # print(cnt, len(fragment))
        for i in range(len(fragment)):
            v = fragment[i][0]
            if v not in v_set:
                continue
            if edges.get(v) is None:
                edges[v] = {}
            for j in range(max(0, i-edges_limit), min(i+edges_limit+1, len(fragment))):
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
    # print(flist[-1][-1][0])
    variant_count = {}
    for fragment in flist:
        for snp in fragment:
            if snp[0] not in variant_count.keys():
                variant_count[snp[0]] = [0, 0]
            variant_count[snp[0]][int(snp[1])] += 1
    return variant_count


def vertex_filter(variant_count, min_hetero_reads):
    filtered_v = set()
    homo_variants = [set(), set()]
    for snp_idx, snp_count in variant_count.items():
        if snp_count[0] > min_hetero_reads and snp_count[1] > min_hetero_reads:
            filtered_v.add(snp_idx)
        else:
            is_variant = snp_count[1] > snp_count[0]
            homo_variants[is_variant].add(snp_idx)
    return homo_variants, filtered_v


def merge_edge(edges):
    merged_edges = {}
    for v in edges:
        merged_edges[v] = {}
        for u in edges[v]:
            t = edges[v][u]
            # t.sort()
            # merged_edges[v][u] = (t[0]+t[1], t[2]+t[3])
            merged_edges[v][u] = (t[0]*t[3] + t[1]*t[2], (t[0]*t[1] + t[2]*t[3] + t[0]*t[2] + t[1]*t[3])*2)
    return merged_edges


def greedy_algorithm(merged_edges, limit_threshold):
    incorrect_variants = set()
    correct_variants = set()
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
        # print(v, len(variant_set))
        if v[0] < limit_threshold:
            incorrect_variants.add(v[1])
            for u in merged_edges[v[1]]:
                variant_set.remove((variant_score[u], u))
                variant_score[u] -= merged_edges[u][v[1]][0] - merged_edges[u][v[1]][1]
                variant_set.add((variant_score[u], u))
                merged_edges[u].pop(v[1])
        else:
            correct_variants.add(v[1])

    return incorrect_variants, correct_variants


def iterative_local_search(merged_edges, threshold, iteration_num, random_change_factor, init_false_group_set):
    group = {}  # [idx] = 0 -> error, 1 -> True
    group_size = [0, 0]
    # initialized
    for v in merged_edges:
        edge_score = 0
        for u in merged_edges[v]:
            edge_score += merged_edges[v][u][0] - merged_edges[v][u][1]
        if init_false_group_set is not None:
            group[v] = v not in init_false_group_set
        else:
            group[v] = int(edge_score > 0)
            # group[v] = np.random.randint(2)
        group_size[group[v]] += 1

    objective_value, moves_score, sorted_list = get_objective_value(group, merged_edges, group_size)
    # local_search
    for i in range(iteration_num):
        print('Iteration #{i}'.format(i=i))
        random_group_change(merged_edges, group, group_size, random_change_factor)
        objective_value, moves_score, sorted_list = get_objective_value(group, merged_edges, group_size)
        local_search_step(objective_value, sorted_list, moves_score, group_size, threshold, group, merged_edges)
        get_objective_value(group, merged_edges, group_size)

    get_objective_value(group, merged_edges, group_size)

    variants = [set(), set()]
    for idx, is_true in group.items():
        variants[is_true].add(idx)
    return variants[0], variants[1]


# will mutate arguments
def random_group_change(merged_edges, group, group_size, random_change_factor):
    for v in merged_edges:
        if np.random.randint(random_change_factor) == 0:
            group_size[group[v]] -= 1
            group[v] = 1 - group[v]
            group_size[group[v]] += 1
    return group, group_size


# will mutate arguments
def local_search_step(objective_value, sorted_list, moves_score, group_size, threshold, group, merged_edges):
    while True:
        move_score, v = sorted_list.pop(-1)
        if move_score < threshold:
            break
        group_size[group[v]] -= 1
        group[v] = 1 - group[v]
        group_size[group[v]] += 1
        moves_score[v] = -move_score
        sorted_list.add((moves_score[v], v))
        for u in merged_edges[v]:
            edge_value = merged_edges[v][u][0] - merged_edges[v][u][1]
            sorted_list.remove((moves_score[u], u))
            if group[v] == group[u]:
                objective_value[group[v]] += edge_value
                objective_value['middle'] -= edge_value
                moves_score[u] -= 2 * edge_value
            else:
                objective_value['middle'] += edge_value
                objective_value[group[u]] -= edge_value
                moves_score[u] += 2 * edge_value
            sorted_list.add((moves_score[u], u))


def get_objective_value(group, merged_edges, group_size):
    objective_value = {0: 0, 1: 0, 'middle': 0}
    sorted_list = SortedList()  # (move_score, idx)
    moves_score = {}  # [idx] = -move_score
    for v in merged_edges:
        move_score = 0
        for u in merged_edges[v]:
            edge_value = merged_edges[v][u][0] - merged_edges[v][u][1]
            edge_value = merged_edges[u][v][0] - merged_edges[u][v][1]
            if group[v] == group[u]:
                objective_value[group[v]] += edge_value
            else:
                objective_value['middle'] += edge_value
            if group[v] == 0 and group[u] == 1:
                move_score += 2 * edge_value
            elif group[v] == 1 and group[u] == 1:
                move_score -= 2 * edge_value
        moves_score[v] = move_score
        sorted_list.add((move_score, v))

    objective_value[0] /= 2
    objective_value[1] /= 2
    objective_value['middle'] /= 2
    print('highest move score(score, vertex): {s}'.format(s=sorted_list[-1]))
    print('group size: {s}'.format(s=group_size))
    print('objective_val: {o}, sum: {s}'.format(o=objective_value, s=-objective_value[0] + objective_value[1] - objective_value['middle']))
    return objective_value, moves_score, sorted_list


def get_truth_objective_value(truth_pos, vcf_dict, merged_edges):
    group_size = [0, 0]
    group = {}
    for v in merged_edges.keys():
        group[v] = vcf_dict[v] in truth_pos
        group_size[group[v]] += 1
    get_objective_value(group, merged_edges, group_size)
    return group, group_size


def accuracy(vcf_file, vcf_dict, incorrect_variants_pos, merged_edges, homo_variants, edges, output_true_variants_set):
    print('accuracy:')
    all_variants_pos = set([x for x in vcf_dict.values()])
    output_true_variants_pos = set([vcf_dict[x] for x in output_true_variants_set])
    all_variants_pos_to_idx = {pos: idx for idx, pos in vcf_dict.items()}

    truth_pos = set()
    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue
            if el[0] != 'chr20':
                continue
            if el[-1][2] == el[-1][0]:
                continue
            genomic_pos = int(el[1])-1
            truth_pos.add(genomic_pos)
    truth_pos = set([x for x in truth_pos if x <= max(max(incorrect_variants_pos), max(output_true_variants_pos))])
    get_truth_objective_value(truth_pos, vcf_dict, merged_edges)
    all_variants_intersect = all_variants_pos.intersection(truth_pos)
    statistics = [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]
    homo_counts = [[0, 0], [0, 0]]
    edges_by_vertex_type =  [[[], []], [[], []]]
    statistics_edges_by_vertex_type = [[[0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0]]]
    for pos_v in all_variants_pos:
        v = all_variants_pos_to_idx[pos_v]
        is_v_truth = pos_v in truth_pos
        if v in homo_variants[0]:
            homo_counts[is_v_truth][0] += 1
            continue
        elif v in homo_variants[1]:
            homo_counts[is_v_truth][1] += 1
            continue
        # for u in merged_edges[v]:
        #     pos_u = vcf_dict[u]
        #     is_u_truth = pos_u in truth_pos
        #     statistics[is_v_truth][is_u_truth][0] += merged_edges[v][u][0]
        #     statistics[is_v_truth][is_u_truth][1] += merged_edges[v][u][1]
        #     statistics_edges_by_vertex_type[is_v_truth][is_u_truth][0] += edges[v][u][0]
        #     statistics_edges_by_vertex_type[is_v_truth][is_u_truth][1] += edges[v][u][1]
        #     statistics_edges_by_vertex_type[is_v_truth][is_u_truth][2] += edges[v][u][2]
        #     statistics_edges_by_vertex_type[is_v_truth][is_u_truth][3] += edges[v][u][3]
        #     edges_by_vertex_type[is_v_truth][is_u_truth].append(edges[v][u])
        # edges_by_vertex_type[is_v_truth][0].append([-1, -1, -1, -1])
        # edges_by_vertex_type[is_v_truth][1].append([-1, -1, -1, -1])


    # for i in range(2):
    #     for j in range(2):
    #         with open('edges_by_vertex_type%d%d.txt' % (i,j), 'w') as file:
    #             file.write(str(edges_by_vertex_type[i][j]))
    # print(statistics)
    # print(np.array(statistics_edges_by_vertex_type)/[[[len(e) for _ in range(4)] for e in v] for v in edges_by_vertex_type])
    # print(homo_counts)

    print('''\t
          output_false in truth: {a} \t
          output_truth in truth: {b} \t
          output_homo in truth: {c} \t
          truth size: {d} \t
          output_false size: {e} \t
          output_truth size: {f} \t
          output_homo size: {g} \t'''.format(
        a=len(incorrect_variants_pos.intersection(truth_pos)),
        b=len(output_true_variants_pos.intersection(truth_pos)),
        c=homo_counts[1][1],
        d=len(truth_pos),
        e=len(incorrect_variants_pos),
        f=len(output_true_variants_pos),
        g=len(homo_variants[1]))
    )
