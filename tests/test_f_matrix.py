import unittest
from ete3 import Tree
import torch
from python.f_matrix import (
    make_all_isochronous_f_matrices,
    make_all_heterochronous_f_matrices,
    f_mat_to_d_mat,
    d_mat_to_f_mat,
    d_mat_to_e_mat,
    e_mat_to_d_mat,
    tree_to_e_mat,
    e_mat_to_tree,
    tree_to_f_mat,
    render_ranked_tree_shape,
    sample_heterochronous_f_matrix,
    heterochronous_ranked_tree_shape_from_time_tree,
)
from collections import defaultdict


class TestFMatrixMethods(unittest.TestCase):

    def test_isochronous_counts(self):
        """
        Verify that make_all_isochronous_f_matrices produces the correct number of
        F-matrices for small numbers of tips.
        """
        n_tips = range(3, 13)
        true_counts = [1, 2, 5, 16, 61, 272, 1385, 7936, 50521, 353792]
        mat_counts = list(map(len, map(make_all_isochronous_f_matrices, n_tips)))
        self.assertEqual(true_counts, mat_counts)

    def test_hetrochronous_counts(self):
        """
        Verify that make_all_heterochronous_f_matrices produces the correct number of
        F-matrices for small numbers of tips.
        """
        n_tips = range(3, 8)
        true_counts = [4, 34, 496, 11056, 349504]
        mat_counts = list(map(len, map(make_all_heterochronous_f_matrices, n_tips)))
        self.assertEqual(true_counts, mat_counts)

    def test_isochronous_f_mats(self):
        """
        Verify that make_all_isochronous_f_matrices produces the correct F-matrixes for
        very small numbers of tips.
        """
        n_tips = range(3, 6)
        three_tip_f_mats = {((2, 0), (1, 3))}
        four_tips_f_mats = {
            ((2, 0, 0), (1, 3, 0), (0, 2, 4)),
            ((2, 0, 0), (1, 3, 0), (1, 2, 4)),
        }
        five_tips_f_mats = {
            ((2, 0, 0, 0), (1, 3, 0, 0), (0, 2, 4, 0), (0, 1, 3, 5)),
            ((2, 0, 0, 0), (1, 3, 0, 0), (0, 2, 4, 0), (0, 2, 3, 5)),
            ((2, 0, 0, 0), (1, 3, 0, 0), (1, 2, 4, 0), (0, 1, 3, 5)),
            ((2, 0, 0, 0), (1, 3, 0, 0), (1, 2, 4, 0), (1, 1, 3, 5)),
            ((2, 0, 0, 0), (1, 3, 0, 0), (1, 2, 4, 0), (1, 2, 3, 5)),
        }
        true_f_mats = [three_tip_f_mats, four_tips_f_mats, five_tips_f_mats]
        f_mats = [
            {
                tuple(map(lambda x: tuple(map(int, x)), mat.numpy().tolist()))
                for mat in mats
            }
            for mats in map(make_all_isochronous_f_matrices, n_tips)
        ]
        self.assertEqual(true_f_mats, f_mats)

    def test_heterochronous_f_mats(self):
        """
        Verify that make_all_heterochronous_f_matrices produces the correct F-matrixes for
        very small numbers of tips.
        """
        n_tips = range(2, 4)
        two_tip_f_mats = {((2, 0), (1, 1))}
        three_tip_f_mats = {
            ((2, 0, 0, 0), (1, 1, 0, 0), (0, 0, 2, 0), (0, 0, 1, 1)),
            ((2, 0, 0, 0), (1, 3, 0, 0), (0, 2, 2, 0), (0, 1, 1, 1)),
            ((2, 0, 0, 0), (1, 3, 0, 0), (1, 2, 2, 0), (0, 1, 1, 1)),
            ((2, 0, 0, 0), (1, 3, 0, 0), (1, 2, 2, 0), (1, 1, 1, 1)),
        }
        true_f_mats = [two_tip_f_mats, three_tip_f_mats]
        f_mats = [
            {
                tuple(map(lambda x: tuple(map(int, x)), mat.numpy().tolist()))
                for mat in mats
            }
            for mats in map(make_all_heterochronous_f_matrices, n_tips)
        ]
        self.assertEqual(true_f_mats, f_mats)

    def test_matrix_bijections(self):
        """
        Verify that d_to_f_mat∘f_mat_to_d_mat, e_to_d_mat∘d_mat_to_e_mat, and
        tree_to_e_mat∘e_mat_to_tree are all the identity function on matrices for small
        numbers of tips.
        """
        methods_and_domains = [
            (True, make_all_isochronous_f_matrices, range(3, 10)),
            (False, make_all_heterochronous_f_matrices, range(3, 7)),
        ]
        for isochronous, method, domain in methods_and_domains:
            e_to_tree = lambda e: e_mat_to_tree(e, isochronous)
            tree_to_e = lambda tree: tree_to_e_mat(tree, isochronous)
            for n_tips in domain:
                f_mats = method(n_tips)
                d_mats = list(map(f_mat_to_d_mat, f_mats))
                e_mats = list(map(d_mat_to_e_mat, d_mats))
                trees = list(map(e_to_tree, e_mats))
                e_from_trees = list(map(tree_to_e, trees))
                d_from_e_mats = list(map(e_mat_to_d_mat, e_mats))
                f_from_d_mats = list(map(d_mat_to_f_mat, d_mats))

                c0 = all((m0 == m1).all() for m0, m1 in zip(f_mats, f_from_d_mats))
                c1 = all((m0 == m1).all() for m0, m1 in zip(d_mats, d_from_e_mats))
                c2 = all((m0 == m1).all() for m0, m1 in zip(e_mats, e_from_trees))
                for c in (c0, c1, c2):
                    self.assertTrue(c)

    def test_heterochronous_tree_and_e_bijections(self):
        """
        Verify that the heterochronous ranked tree shape:
              |
             _0_
           _1_  |
          |   | 2
         _3_  |
        |   4 |
        |     5
        6,
        produces the correct E matrix and vice versa.
        """
        nodes = [Tree() for _ in range(7)]
        for rank, node in enumerate(nodes):
            node.add_feature("rank_label", rank)
        n0, n1, n2, n3, n4, n5, n6 = nodes
        n0.add_child(n1, dist=1)
        n0.add_child(n2, dist=2)
        n1.add_child(n3, dist=2)
        n1.add_child(n5, dist=4)
        n3.add_child(n4, dist=1)
        n3.add_child(n6, dist=3)

        e_mat = [
            [1, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
        ]
        e_mat = torch.tensor(e_mat)

        e_from_tree = tree_to_e_mat(n0, isochronous=False)
        tree_from_e = e_mat_to_tree(e_mat, isochronous=False)
        mat0, _ = weighted_adjacency_matrix(n0)
        mat1, _ = weighted_adjacency_matrix(tree_from_e)
        self.assertTrue((e_from_tree == e_mat).all().item())
        self.assertEqual(mat0, mat1)

        # Uncomment the next line to visualize the tree.
        # render_ranked_tree_shape(n0)

    def test_isochronous_tree_and_e_bijections(self):
        """
        Verify that the isochronous ranked tree shape:
                   |
             ______0____
           __1__        |
          |     |     __2__
          |     |   _3_    |
         _4_    |  |   |   |
        *   *   *  *   *   *
        produces the correct E matrix and vice versa.
        """
        nodes = [Tree() for _ in range(5)]
        for rank, node in enumerate(nodes):
            node.add_feature("rank_label", rank)
        n0, n1, n2, n3, n4 = nodes
        n0.add_child(n1, dist=1)
        n0.add_child(n2, dist=2)
        n1.add_child(n4, dist=3)
        n2.add_child(n3, dist=1)

        n1.add_child(dist=4)
        n2.add_child(dist=3)
        n3.add_child(dist=2)
        n3.add_child(dist=2)
        n4.add_child(dist=1)
        n4.add_child(dist=1)

        e_mat = [
            [1, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 1, 1, 2, 2],
        ]
        e_mat = torch.tensor(e_mat)

        e_from_tree = tree_to_e_mat(n0, isochronous=True)
        tree_from_e = e_mat_to_tree(e_mat, isochronous=True)
        mat0, dists0 = weighted_adjacency_matrix(n0)
        mat1, dists1 = weighted_adjacency_matrix(tree_from_e)
        self.assertTrue((e_from_tree == e_mat).all().item())
        self.assertEqual(mat0, mat1)
        self.assertEqual(dists0, dists1)

        # Uncomment the next line to visualize the tree.
        # render_ranked_tree_shape(n0)

    def test_heterochronous_sampling(self):
        """
        Verify that the heterochronous sampling methods produces the correct matrices,
        the probabilities sum to 1, and the normalized frquencey counts are
        approximately the probabilities for small numbers of tips.
        """
        for n in range(3, 6):
            r = 10 ** (n + 1)
            samples = {}
            for _ in range(r):
                f_mat, prob = sample_heterochronous_f_matrix(n)
                f_mat = tuple(map(tuple, f_mat.tolist()))
                if f_mat not in samples:
                    samples[f_mat] = [prob, 1]
                else:
                    samples[f_mat][1] += 1

            all_mats = {
                tuple(map(tuple, mat.tolist()))
                for mat in make_all_heterochronous_f_matrices(n)
            }

            prob = sum(x[0] for x in samples.values())
            max_diff = max(x[0] - x[1] / r for x in samples.values())
            self.assertEqual(samples.keys(), all_mats)
            self.assertAlmostEqual(1, prob)
            self.assertAlmostEqual(0, max_diff, delta=0.01)

    def test_heterochronous_from_three_tip_time_trees(self):
        """
        Verify that the method heterochronous_ranked_tree_shape_from_time_tree breaks
        ties (nodes at the same distance from the root) and produces the correct
        heterochronous ranked tree shapes on various three tip trees and one four tip
        tree. Additionally, verify that when there are no ties among parent and child
        nodes, heterochronous_ranked_tree_shape_from_time_tree yields uniform sampling
        of the possible heterochronous ranked tree shapes.
        """
        # Map a newick string with branch lengths to bool for definitely uniform or
        # maybe not, the heterochronous ranked tree shapes, and the number of
        # replicates.
        nwk = "((a:{},b:{}):{},c:{});"
        F0 = ((2, 0, 0, 0), (1, 3, 0, 0), (1, 2, 2, 0), (0, 1, 1, 1))
        F1 = ((2, 0, 0, 0), (1, 3, 0, 0), (1, 2, 2, 0), (1, 1, 1, 1))
        F2 = ((2, 0, 0, 0), (1, 1, 0, 0), (0, 0, 2, 0), (0, 0, 1, 1))
        F3 = ((2, 0, 0, 0), (1, 3, 0, 0), (0, 2, 2, 0), (0, 1, 1, 1))
        r = 50000
        trees_to_shapes = {
            nwk.format(1, 1, 1, 3): [True, (F1,), r],
            nwk.format(1, 1, 0, 2): [True, (F1,), r],
            nwk.format(2, 2, 1, 2): [True, (F3,), r],
            nwk.format(2, 2, 0, 1): [True, (F3,), r],
            nwk.format(1, 1, 1, 2): [True, (F0, F1, F3), r],
            nwk.format(1, 1, 0, 1): [True, (F0, F1, F3), r],
            nwk.format(1, 2, 1, 2): [True, (F0, F3), r],
            nwk.format(2, 1, 1, 2): [True, (F0, F3), r],
            nwk.format(1, 2, 0, 1): [True, (F0, F3), r],
            nwk.format(2, 1, 0, 1): [True, (F0, F3), r],
            nwk.format(2, 1, 1, 3): [True, (F0, F1), r],
            nwk.format(1, 2, 1, 3): [True, (F0, F1), r],
            nwk.format(2, 1, 0, 2): [True, (F0, F1), r],
            nwk.format(1, 2, 0, 2): [True, (F0, F1), r],
            nwk.format(1, 1, 1, 1): [True, (F2, F3), r],
            nwk.format(1, 1, 0, 0): [True, (F2, F3), r],
            nwk.format(1, 2, 1, 1): [True, (F2, F3), r],
            nwk.format(2, 1, 1, 1): [True, (F2, F3), r],
            nwk.format(1, 2, 0, 0): [True, (F2, F3), r],
            nwk.format(2, 1, 0, 0): [True, (F2, F3), r],
            nwk.format(1, 1, 2, 1): [True, (F2,), r],
            nwk.format(1, 1, 1, 0): [True, (F2,), r],
            nwk.format(0, 1, 1, 3): [False, (F1,), r],
            nwk.format(1, 0, 1, 3): [False, (F1,), r],
            nwk.format(0, 1, 0, 2): [False, (F1,), r],
            nwk.format(1, 0, 0, 2): [False, (F1,), r],
            nwk.format(0, 1, 1, 2): [False, (F0, F1), r],
            nwk.format(1, 0, 1, 2): [False, (F0, F1), r],
            nwk.format(0, 1, 0, 1): [False, (F0, F1), r],
            nwk.format(1, 0, 0, 1): [False, (F0, F1), r],
            nwk.format(0, 2, 1, 2): [False, (F0,), r],
            nwk.format(2, 0, 1, 2): [False, (F0,), r],
            nwk.format(0, 2, 0, 1): [False, (F0,), r],
            nwk.format(2, 0, 0, 1): [False, (F0,), r],
            nwk.format(0, 0, 1, 2): [False, (F1,), r],
            nwk.format(0, 0, 0, 1): [False, (F1,), r],
            nwk.format(1, 2, 0, 0): [False, (F2, F3), r],
            nwk.format(2, 1, 0, 0): [False, (F2, F3), r],
            nwk.format(1, 2, 1, 1): [False, (F2, F3), r],
            nwk.format(2, 1, 0, 0): [False, (F2, F3), r],
            nwk.format(0, 0, 0, 0): [False, (F0, F1, F2, F3), r],
            nwk.format(0, 0, 1, 1): [False, (F0, F1, F2, F3), r],
            nwk.format(0, 1, 2, 1): [False, (F2,), r],
            nwk.format(0, 1, 1, 0): [False, (F2,), r],
            nwk.format(0, 0, 2, 1): [False, (F2,), r],
            nwk.format(0, 0, 1, 0): [False, (F2,), r],
            nwk.format(1, 2, 1, 1): [False, (F2, F3), r],
            nwk.format(2, 1, 1, 1): [False, (F2, F3), r],
            nwk.format(1, 2, 0, 0): [False, (F2, F3), r],
            nwk.format(2, 1, 0, 0): [False, (F2, F3), r],
        }
        nwk = "((a:1,b:1):1,(c:1,d:1):1);"
        f_mats = (
            (
                (2, 0, 0, 0, 0, 0),
                (1, 3, 0, 0, 0, 0),
                (0, 2, 4, 0, 0, 0),
                (0, 1, 3, 3, 0, 0),
                (0, 0, 2, 2, 2, 0),
                (0, 0, 1, 1, 1, 1),
            ),
            (
                (2, 0, 0, 0, 0, 0),
                (1, 3, 0, 0, 0, 0),
                (0, 2, 4, 0, 0, 0),
                (0, 2, 3, 3, 0, 0),
                (0, 1, 2, 2, 2, 0),
                (0, 0, 1, 1, 1, 1),
            ),
            (
                (2, 0, 0, 0, 0, 0),
                (1, 3, 0, 0, 0, 0),
                (0, 2, 4, 0, 0, 0),
                (0, 1, 3, 3, 0, 0),
                (0, 1, 2, 2, 2, 0),
                (0, 0, 1, 1, 1, 1),
            ),
            (
                (2, 0, 0, 0, 0, 0),
                (1, 3, 0, 0, 0, 0),
                (0, 2, 4, 0, 0, 0),
                (0, 2, 3, 3, 0, 0),
                (0, 2, 2, 2, 2, 0),
                (0, 1, 1, 1, 1, 1),
            ),
            (
                (2, 0, 0, 0, 0, 0),
                (1, 3, 0, 0, 0, 0),
                (0, 2, 4, 0, 0, 0),
                (0, 2, 3, 3, 0, 0),
                (0, 1, 2, 2, 2, 0),
                (0, 1, 1, 1, 1, 1),
            ),
            (
                (2, 0, 0, 0, 0, 0),
                (1, 3, 0, 0, 0, 0),
                (0, 2, 4, 0, 0, 0),
                (0, 1, 3, 3, 0, 0),
                (0, 1, 2, 2, 2, 0),
                (0, 1, 1, 1, 1, 1),
            ),
        )
        trees_to_shapes[nwk] = [True, f_mats, 100000]

        for nwk, (uniform, tree_shapes, r) in trees_to_shapes.items():
            tree = Tree(nwk)
            sample_counts = defaultdict(int)
            for _ in range(r):
                shape = heterochronous_ranked_tree_shape_from_time_tree(tree, False)
                f_mat = tree_to_f_mat(shape, False)
                f_mat = tuple(map(lambda row: tuple(map(int, row)), f_mat))
                sample_counts[f_mat] += 1
            # Verify we get the correct ranked tree shapes
            self.assertCountEqual(tree_shapes, sample_counts.keys())
            # Verify that if we should definitely have uniform sampling, then we have
            # uniform sampling.
            if uniform and len(tree_shapes) > 1:
                probs = [v / r for v in sample_counts.values()]
                max_diff = max(abs(probs[0] - p) for p in probs[1:])
                self.assertAlmostEqual(max_diff, 0, delta=0.01)


def weighted_adjacency_matrix(tree):
    """
    Given a ranked tree shape with N nodes labelled with ranks 0 to N-1 (so the tree has
    N nodes in the heterochronous case and 2N+1 nodes in the isochronous case), return
    a weighted adjacency matrix for the ranked nodes and a list of distances to unranked
    nodes.

    The matrix is a list of lists, where entry (i,j) is the length of the directed edge
    from the node with rank i to the node with rank j (the edge length is 0 when the
    rank j node is not a child of the rank i node). Note this follows the opposite
    convention of E-matrices for rows and columns.

    Entry i of the list of distances is a list containing the distances from the node
    with rank i to its child leaf nodes. Such an inner list is length 0, 1, or 2 and
    all distances in an inner list are equal. In the heterochronous case, the inner
    lists are all empty, as all nodes are ranked.
    """
    ranked_nodes = [n for n in tree.traverse() if hasattr(n, "rank_label")]
    n_nodes = len(ranked_nodes)
    adj_matrix = [[0] * n_nodes for _ in range(n_nodes)]
    leaf_distances = [[] for _ in range(n_nodes)]
    for parent in ranked_nodes:
        for child in parent.children:
            if hasattr(child, "rank_label"):
                adj_matrix[parent.rank_label][child.rank_label] = child.dist
            else:
                leaf_distances[parent.rank_label].append(child.dist)

    return adj_matrix, leaf_distances


if __name__ == "__main__":
    unittest.main()
