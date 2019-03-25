import numpy as np
from sortedcontainers import SortedList
import seaborn as sns
import matplotlib.pyplot as plt
from stats import get_barnard, get_fisher


def read_vcf_file(filename, chr_list=None, return_homo=True, return_hetero=True, from_pos=-1, to_pos=-1):
    snp_ix = 0
    vcf_dict = dict()
    with open(filename, 'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue
            if (chr_list is not None) and (el[0] not in chr_list):
                continue
            if (not return_homo) and el[-1][2] == el[-1][0]:
                continue
            if (not return_hetero) and el[-1][2] != el[-1][0]:
                continue
            genomic_pos = int(el[1]) - 1
            if to_pos != -1 and genomic_pos > to_pos:
                break
            if from_pos != -1 and genomic_pos < from_pos:
                continue
            vcf_dict[snp_ix] = genomic_pos
            snp_ix += 1
    return vcf_dict


def dict_to_set(d):
    return set(d.values())


def read_fragment_matrix(frag_matrix, vcf_file, max_num_reads, from_pos, to_pos):
    vcf_dict = read_vcf_file(vcf_file, from_pos=from_pos, to_pos=to_pos)
    min_var_idx = min(vcf_dict.keys())
    max_var_idx = max(vcf_dict.keys())
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
            alist = list(filter(lambda e: max_var_idx >= e[0] >= min_var_idx, alist))
            flist.append(alist)
    # print(cnt_quality)
    # print(flist[-1])
    return vcf_dict, flist


def keep_bed_intersect(vcf_dict, flist, intersect_file, from_pos, to_pos):
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


def keep_output_intersect(vcf_dict, flist, out_file, from_pos, to_pos):
    intersect_pos = dict_to_set(read_vcf_file(out_file, from_pos=from_pos, to_pos=to_pos))

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


def vertex_filter(vcf_dict, variant_count, min_hetero_reads):
    filtered_v = set()
    homo_variants = [set(), set()]
    for snp_idx, snp_count in variant_count.items():
        if snp_count[0] > min_hetero_reads and snp_count[1] > min_hetero_reads:
            filtered_v.add(snp_idx)
        else:
            is_variant = snp_count[1] > snp_count[0]
            homo_variants[is_variant].add(vcf_dict[snp_idx])
    return homo_variants, filtered_v


def merge_edge(edges, method, min_sum_edges, significance_threshold):
    print("Building edges for  {} vertexes".format(len(edges)))
    merged_edges = {}
    if method == 'barnard':
        f = get_barnard
    else:
        f = get_fisher
    for i, v in enumerate(edges):
        interval = int(len(edges) / 10)
        if i % interval == 0:
            frac = int(i / interval) * 10
            print("{}%...".format(frac))

        merged_edges[v] = {}
        for u in edges[v]:
            t = edges[v][u]
            # t.sort()
            # merged_edges[v][u] = (t[0]+t[1], t[2]+t[3])
            # merged_edges[v][u] = (t[0]*t[3] + t[1]*t[2], (t[0]*t[1] + t[2]*t[3] + t[0]*t[2] + t[1]*t[3])*2)
            # print(t)
            if sum(t) >= min_sum_edges:
                merged_edges[v][u] = -np.log(f(*t) + 0.000001) + np.log(significance_threshold)
                # merged_edges[v][u] = -np.log(barnard_test(*t)[0]+0.000001) + np.log(.01)
                # if merged_edges[v][u] < 0:
                #     print(v, u, t, -np.log(scstat.fisher_exact([[t[0], t[1]], [t[2], t[3]]])[-1]), + np.log(.10), -np.log(barnard_test(*t)[0]))
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


# return variants pos
def iterative_local_search(vcf_dict, merged_edges, threshold, iteration_num, random_change_factor, random_swap_factor, init_false_group_set):

    group, group_size, group_set = init_local_search(merged_edges, init_false_group_set)
    objective_value, moves_score, sorted_list = get_objective_value(group, merged_edges, group_size)
    # local_search
    for i in range(iteration_num):
        print('\n\n*******Iteration #{i}*********'.format(i=i))
        random_group_change(merged_edges, group, group_set, group_size, random_change_factor)
        random_group_swap(group, group_set, group_size, random_swap_factor)
        objective_value, moves_score, sorted_list = get_objective_value(group, merged_edges, group_size)
        local_search_step(objective_value, sorted_list, moves_score, group_size, threshold, group, merged_edges, group_set)
        get_objective_value(group, merged_edges, group_size)

    get_objective_value(group, merged_edges, group_size)

    variants = [set(), set()]
    for idx, is_true in group.items():
        variants[is_true].add(vcf_dict[idx])
    return variants[0], variants[1]


def init_local_search(merged_edges, init_false_group_set):
    group = {}  # [idx] = 0 -> error, 1 -> True
    group_size = [0, 0]
    group_set = (SortedList(), SortedList())
    # initialized
    for v in merged_edges:
        edge_score = 0
        for u in merged_edges[v]:
            edge_score += merged_edges[v][u]
        if init_false_group_set is not None:
            group[v] = v not in init_false_group_set
        else:
            group[v] = int(edge_score > 0)
            # group[v] = np.random.randint(2)
        group_set[group[v]].add(v)
        group_size[group[v]] += 1
    return group, group_size, group_set


# will mutate arguments
def random_group_change(merged_edges, group, group_set, group_size, random_change_factor):
    for v in merged_edges:
        if np.random.random_sample() < random_change_factor:
            change_vertex_group(v, group, group_set, group_size)
    return group, group_size


def random_group_swap(group, group_set, group_size, random_swap_factor):
    if min(group_size) == 0:
        return
    swap_size = int(min(group_size) * random_swap_factor) + 1
    swap_idx = np.random.choice(group_size[0], swap_size), np.random.choice(group_size[1], swap_size)
    for i in range(swap_size):
        v, u = group_set[0][swap_idx[0][i]], group_set[1][swap_idx[1][i]]
        change_vertex_group(v, group, group_set, group_size)
        change_vertex_group(u, group, group_set, group_size)


def change_vertex_group(v, group, group_set, group_size):
    group_size[group[v]] -= 1
    group_set[group[v]].remove(v)
    group[v] = 1 - group[v]
    group_size[group[v]] += 1
    group_set[group[v]].add(v)


# will mutate arguments
def local_search_step(objective_value, sorted_list, moves_score, group_size, threshold, group, merged_edges, group_set):
    while True:
        move_score, v = sorted_list.pop(-1)
        if move_score < threshold:
            break
        change_vertex_group(v, group, group_set, group_size)
        moves_score[v] = -move_score
        sorted_list.add((moves_score[v], v))
        for u in merged_edges[v]:
            edge_value = merged_edges[v][u]
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
            edge_value = merged_edges[v][u]
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


def draw_histogram(group, merged_edges):
    edges_by_group = {0: [], 1: [], 'middle': []}
    for v in merged_edges:
        for u in merged_edges[v]:
            edge_value = merged_edges[v][u]
            if group[v] == group[u]:
                edges_by_group[group[v]].append(edge_value)
            else:
                edges_by_group['middle'].append(edge_value)
    sns.distplot(edges_by_group[0])
    plt.show()
    sns.distplot(edges_by_group[1])
    plt.show()
    sns.distplot(edges_by_group['middle'])
    plt.show()


def get_truth_objective_value(truth_pos, vcf_dict, merged_edges):
    group_size = [0, 0]
    group = {}
    for v in merged_edges.keys():
        group[v] = vcf_dict[v] in truth_pos
        group_size[group[v]] += 1
    draw_histogram(group, merged_edges)
    get_objective_value(group, merged_edges, group_size)
    return group, group_size


def accuracy(merged_edges, output_variants_pos, homo_variants_pos, known_pos, hetero_truth_pos, homo_truth_pos, vcf_dict):
    print('accuracy:')

    get_truth_objective_value(hetero_truth_pos, vcf_dict, merged_edges)
    print('''\t
          output_hetero_false in real hetero truth: {} \t output_hetero_false(known) in real hetero truth: {} \t
          output_hetero_truth in real hetero truth: {} \t output_hetero_truth(known) in real hetero truth: {} \t
          output_homo_truth in real homo truth: {} \t
          real homo truth size: {} \t 
          real hetero truth size: {} \t
          output_false size: {} \t
          output_truth size: {} \t
          output_homo size: 0:{} 1:{}\t'''.format(
        len(output_variants_pos[0].intersection(hetero_truth_pos)),
        len(output_variants_pos[0].intersection(hetero_truth_pos).intersection(known_pos)),
        len(output_variants_pos[1].intersection(hetero_truth_pos)),
        len(output_variants_pos[1].intersection(hetero_truth_pos).intersection(known_pos)),
        len(homo_variants_pos[0].union(homo_variants_pos[1]).intersection(homo_truth_pos)),
        len(homo_truth_pos),
        len(hetero_truth_pos),
        len(output_variants_pos[0]),
        len(output_variants_pos[1]),
        len(homo_variants_pos[0]), len(homo_variants_pos[1])
        )
    )
