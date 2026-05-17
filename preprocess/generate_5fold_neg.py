import os
import numpy as np
import argparse

# import sys
# sys.path.append("../")

# from chinn import pair_features
from pair_generation import sample_from_neg_pairs
#min_dist=5000, max_dist=2000000
def load_pairs_as_dict(files, min_length=1000, min_dist=5000, max_dist=2000000, max_length=None):
    scores = {}
    t_dists = {}
    for inter_file in files:
        with open(inter_file, 'r') as f:
            for r in f:
                tokens = r.strip().split()
                if len(tokens) < 7:
                    tokens.append(0)
                for i in [1, 2, 4, 5, 6]:
                    try:
                        tokens[i] = int(tokens[i])
                    except:
                        tokens[i] = int(float(tokens[i]))
                if tokens[0] == tokens[3] and tokens[1] > tokens[4]:
                    temp = tokens[1]
                    tokens[1] = tokens[4]
                    tokens[4] = temp
                    temp = tokens[2]
                    tokens[2] = tokens[5]
                    tokens[5] = temp
                if max_length is not None and (tokens[2]-tokens[1] > max_length or tokens[5] - tokens[4] > max_length):
                    continue
                if min_length > 0:
                    for i,j in zip([1,4],[2,5]):
                        if tokens[j] - tokens[i] < min_length:
                            diff = min_length - (tokens[j]-tokens[i])
                            half_diff = diff // 2
                            tokens[i] -= half_diff
                            tokens[j] += diff - half_diff
                curr_dist = 0.5 * (tokens[4] + tokens[5] - tokens[1] - tokens[2])
                if min_dist <= curr_dist <= max_dist:
                    scores[tuple(tokens[:6])] = tokens[6]
                    if tokens[0] not in t_dists:
                        t_dists[tokens[0]] = []
                    t_dists[tokens[0]].append(curr_dist)

    return scores, t_dists


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Sampling 5x negative samples.")
    parser.add_argument("name", help="The prefix of the dataset. For example gm12878_ctcf.")
    parser.add_argument("datadir", help="The directory where the input and output are in.")
    args = parser.parse_args()

    num_bins = 50
    dist_range = (np.log10(5000), np.log10(2000000))

    pos_pairs, pos_dists = load_pairs_as_dict(
        [os.path.join(args.datadir, "{}.clustered_interactions.both_gene.bedpe".format(args.name))])

    neg_pairs, neg_dists = load_pairs_as_dict(
        [os.path.join(args.datadir, '{}.no_intra_all.negative_pairs.bedpe'.format(args.name)),
         os.path.join(args.datadir, '{}.random_tf_peak_pairs.filtered.bedpe'.format(args.name))],
        min_length=1000)
    other_neg_pairs, other_neg_dists = load_pairs_as_dict(
        [os.path.join(args.datadir, '{}.shuffled_neg_anchor.neg_pairs.filtered.tf_filtered.bedpe'.format(args.name))],
        min_length=1000)

    selected_neg_pairs = sample_from_neg_pairs(pos_dists, neg_pairs, 5, other_neg_pairs, num_bins, dist_range)
    selected_neg_dists = [0.5*(p[5]+p[4]-p[2]-p[1]) for p in selected_neg_pairs]
    with open(os.path.join(args.datadir,
                           "{}.neg_pairs_5x.from_singleton_inter_tf_random.bedpe".format(args.name)),'w') as out:
        for p in selected_neg_pairs:
            out.write("\t".join(map(str, p)) + "\n")