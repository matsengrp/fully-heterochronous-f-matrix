from collections import defaultdict
import numpy as np
from numpy.random import choice, random
import torch
from ete3 import Tree, TreeStyle, NodeStyle, AttrFace


def make_all_f_matrices(n_tips, isochronous=False):
    """
    Return a list of all F-matrices for isochronous or heterochronous ranked tree
    shapes.
    """
    if isochronous:
        return make_all_isochronous_f_matrices(n_tips)
    else:
        return make_all_heterochronous_f_matrices(n_tips)


def make_all_isochronous_f_matrices(n_tips):
    """
    Return a list of all F-matrices for trees with n_tips (all at identical sampling
    times). Memory and run time explodes for n_tips past 12.
    """
    # Matrix side length.
    m = n_tips - 1

    # Fill constant entries.
    starting_matrix = torch.zeros((m, m))
    for i in range(m - 1):
        starting_matrix[i, i] = i + 2
        starting_matrix[i + 1, i] = i + 1
    starting_matrix[m - 1, m - 1] = m + 1

    # Maintain a list of partially filled matrices.
    current_matrices = [starting_matrix]

    # Fill entries off the diagonal and subdiagonal.
    for i in range(2, m):
        for j in range(0, i - 1):
            more_matrices = []
            for mat in current_matrices:
                left, left_up = (0, 0) if j == 0 else (mat[i, j - 1], mat[i - 1, j - 1])
                up = mat[i - 1, j]

                lower, upper = valid_non_diagonal_values(
                    i, j, left, left_up, up, isochronous=True
                )
                if lower != upper:
                    # Make a new matrix for the second option.
                    mat0 = mat.detach().clone()
                    mat0[i, j] = upper
                    more_matrices.append(mat0)
                mat[i, j] = lower
            current_matrices.extend(more_matrices)
    return current_matrices


def make_all_heterochronous_f_matrices(n_tips):
    """
    Return a list of all F-matrices for trees with n_tips (all at different sampling
    times). Memory and run time explodes for n_tips past 7, but appears to scale with
    the number of F-matrices. Eight tips takes about two hours.
    """
    # Matrix side length.
    m = 2 * n_tips - 2

    # Fill constant entries.
    starting_matrix = torch.zeros((m, m))
    starting_matrix[0, 0] = 2
    starting_matrix[1, 0] = 1
    starting_matrix[m - 1, m - 1] = 1

    # Maintain a list of partially filled matrices.
    current_matrices = [starting_matrix]

    # Fill the diagonal, which forces values on the subdiagonal.
    for i in range(1, m - 1):
        more_matrices = []
        for mat in current_matrices:
            previous = mat[i - 1, i - 1]
            lower, upper = valid_heterochronous_diagonals(i, n_tips, previous)
            if lower != upper:
                # Make a new matrix for the second option.
                mat0 = mat.detach().clone()
                mat0[i, i] = upper
                mat0[i + 1, i] = mat0[i, i] - 1
                more_matrices.append(mat0)
            mat[i, i] = lower
            mat[i + 1, i] = mat[i, i] - 1
        current_matrices.extend(more_matrices)

    # Fill entries off the diagonal and subdiagonal.
    for i in range(2, m):
        for j in range(0, i - 1):
            more_matrices = []
            for mat in current_matrices:
                left, left_up = (0, 0) if j == 0 else (mat[i, j - 1], mat[i - 1, j - 1])
                up = mat[i - 1, j]
                prev_diag = mat[i - 1, i - 1]

                lower, upper = valid_non_diagonal_values(
                    i, j, left, left_up, up, prev_diag, isochronous=False
                )
                if lower != upper:
                    # Make a new matrix for the second option.
                    mat0 = mat.detach().clone()
                    mat0[i, j] = upper
                    more_matrices.append(mat0)
                mat[i, j] = lower
            current_matrices.extend(more_matrices)
    return current_matrices


def valid_values(i, j, dependencies, n=None, isochronous=False):
    """
    Return the tuple of possible values for entry (i,j), based on previous entries.
    There are exactly one or two possible values. When there is only one possible value,
    it is listed twice in the tuple.

    Parameters:
        i (int): The row index.
        j (int): The column index.
        dependencies (list): The dependency values for entry (i,j). In the isochronous
            case, dependencies are formatted like [value to the left, value to the
            above-left, value above]; in the heterochronous case the format is [value to
            the left, value to the above-left, value above, value of the previous
            diagonal]. For combinations  of i and j that do not require some of the
            dependency values (i.e., i=j of the heterochronous case and j=0 in either
            case), one can use any filler value.
        n (int): The number of tips; the matrix is (2n-2)x(2n-2). This parameter is
            only required for heterochronous diagonal entries.
        isochronous (bool): Isochronous or fully heterochronous F-matrices.
    """
    if i < j:
        return 0, 0
    elif i == j:
        if isochronous:
            return i + 2, i + 2
        else:
            return valid_heterochronous_diagonals(i, n, dependencies[-1])
    elif i == j + 1:
        return dependencies[2] - 1
    else:
        return valid_non_diagonal_values(
            i, j, *dependencies.tolist(), isochronous=isochronous
        )


def valid_heterochronous_diagonals(i, n, previous):
    """
    Return the tuple of possible values for diagonal entry i, based on previous
    diagonals. There are exactly one or two possible values. When there is only one
    possible value, it is listed twice in the tuple.

    Parameters:
        i (int): The diagonal index.
        n (int): number of tips
        previous (int): previous diagonal entry
    """
    if i == 0 or previous == 1:
        values = (2, 2)
    elif previous == 2 * n - i - 1:
        values = (previous - 1, previous - 1)
    else:
        values = (previous - 1, previous + 1)
    return values


def valid_non_diagonal_values(
    i, j, left, left_up, up, prev_diag=None, isochronous=False
):
    """
    Return the tuple of possible values for entry (i,j), based on previous entries.
    There are exactly one or two possible values. When there is only one possible value,
    it is listed twice in the tuple.

    Parameters:
        i (int): The row index.
        j (int): The column index.
        left (int): Neighboring entry at (i, j-1).
        left_up (int): Neighboring entry at (i-1, j-1).
        up (int): Neighboring entry at (i-1, j).
        prev_diag (int): Diagonal entry at (i-1,i-1). This is required only for the
            heterochronous case
        isochronous (bool): Isochronous or heterochronous.
    """
    if j == 0:
        # Inequalities for first column.
        lower = max(0, up - 1)
        upper = up
    else:
        # Generic inequalities.
        lower = max(left, up - 1, left + up - left_up - 1)
        upper = min(up, left + up - left_up)
    if not isochronous and up == prev_diag:
        forced_value = prev_diag - 1
        if lower <= forced_value <= upper:
            lower, upper = forced_value, forced_value
        else:
            raise ValueError(
                f"Problem with entry ({i},{j}) with forced value {forced_value}; bounds"
                + f"{lower} and {upper}; neighboring entries {left}, {left_up}, {up}; "
                + f"and previous diagonal value {prev_diag} appearing."
            )
    if 0 <= lower <= upper:
        return (lower, upper)
    else:
        raise ValueError(
            f"Problem with entry ({i},{j}) with bounds {lower} and {upper}, and "
            + f"neighboring entries {left}, {left_up}, {up}."
        )


def f_mat_to_d_mat(f_mat):
    """
    Return the D-matrix for an F-matrix. The D-matrix is a lower triangular matrix where
    column j of the D-matrix is column j minus column (j-1) of the F-Matrix (the first
    column of the D-matrix is the first column of the F-matrix). This type of matrix is
    discussed in the proof of Theorem 1 and labelled with a D.

    The nice thing: the (i,j) entry is the number of descendants of the node labelled j
    that are extant over the i-th time interval.
    """
    d_mat = f_mat.detach().clone()
    for j in range(1, d_mat.shape[1]):
        d_mat[:, j] -= f_mat[:, j - 1]
    d_mat += torch.diag(torch.diagonal(f_mat[:-1, :-1]), diagonal=1)
    return d_mat


def f_mat_to_e_mat(f_mat):
    """Return the E-matrix for the given F-Matrix."""
    return d_mat_to_e_mat(f_mat_to_d_mat(f_mat))


def d_mat_to_f_mat(d_mat):
    """Return the F-matrix for the given D-Matrix."""
    f_mat = d_mat.detach().clone()
    for j in range(1, d_mat.shape[1]):
        f_mat[j:, j] += f_mat[j:, j - 1]
    return f_mat


def d_mat_to_e_mat(d_mat):
    """
    Return the E-matrix for a D-matrix. The E-matrix is a lower triangular where row i
    of the E-matrix is row i minus row (i+1) of the D-matrix (the last row of E-matrix
    is the last row of the D-matrix). The row differences are discussed in the proof of
    Theorem 1, but aren't given a name, so for now we call them E-matrices.

    The nice thing: In the heterochronous case, each row i contains exactly one 1, at
    some column j. The node labelled j has a child that splits or is sampled at the
    (i+1)-th time event. In the isochronous case, this is true for all rows except the
    last. Instead, column j of the last row records the number of leaves that are
    children of the node labelled j.

    This is almost the adjacency matrix for the ranked tree shape. In particular, let A
    be the adjacency matrix for the directed graph of the ranked nodes of a ranked tree
    shape with n-tips (so A is size (2n-1)x(2n-1) in the heterochronous case and is size
    (n-1)x(n-1) in the isochronous case), with the convention that A[i,j]=1 if and only
    if there is an edge from j to i. In there heterochronous case, E[i,j]=A[i+1,j]
    for all i and j. In the isochronous case, E[i,j]=A[i+1,j] for all i and j with
    i < (n-2).
    """
    e_mat = d_mat.detach().clone()
    for i in range(d_mat.shape[0] - 1):
        e_mat[i, :] -= d_mat[i + 1, :]
    e_mat += torch.diag(torch.diagonal(d_mat[1:, 1:]), diagonal=1)
    return e_mat


def e_mat_to_d_mat(e_mat):
    """Return the D-matrix for the given E-matrix."""
    d_mat = e_mat.detach().clone()
    for i in range(d_mat.shape[0] - 1)[::-1]:
        d_mat[i, : i + 1] += d_mat[i + 1, : i + 1]
    return d_mat


def e_mat_to_f_mat(e_mat):
    """Return the F-matrix for the given E-Matrix."""
    return d_mat_to_f_mat(e_mat_to_d_mat(e_mat))


def f_mat_to_tree(f_mat, isochronous=False):
    """
    Return an ete3 tree from a valid F-matrix. The nodes of the tree have an integer
    valued feature "rank_label" and a boolean valued feature "is_sample". The branch
    lengths are set so that there is unit time between events of consecutively ranked
    nodes.
    """
    return d_mat_to_tree(f_mat_to_d_mat(f_mat), isochronous=isochronous)


def d_mat_to_tree(d_mat, isochronous=False):
    """
    Return an ete3 tree from a valid D-matrix. The nodes of the tree have an integer
    valued feature "rank_label". In the heterochronous case, nodes additionally have a
    boolean valued feature "is_sample". The branch lengths are set so that there is unit
    time between events of consecutively ranked nodes.
    """
    return e_mat_to_tree(d_mat_to_e_mat(d_mat), isochronous=isochronous)


def e_mat_to_tree(e_mat, isochronous=False):
    """
    Return the ranked tree shape from the given E-matrix.
    """
    root = Tree()
    root.add_feature("rank_label", 0)
    if not isochronous:
        root.add_feature("is_sample", False)
    root.add_child(dist=0), root.add_child(dist=0)
    labelled_nodes = [root]
    # The node at index i of labelled_nodes is labelled i
    unlabelled = lambda c: not hasattr(c, "rank_label")

    for i in range(e_mat.shape[0] - 1):
        longer_tips(root)
        j = e_mat[i, :].argmax()
        # An unlabelled child of the node labelled j now splits or is sampled.
        parent = labelled_nodes[j]
        # child = next(c for c in parent.children if c.is_leaf())
        child = next(c for c in parent.children if unlabelled(c))
        child.add_feature("rank_label", i + 1)
        if not isochronous:
            # Look at the column for the node labelled i+1. If the column is all zeros,
            # then this node is a sampling event. Otherwise, there is a one somewhere
            #  along the column, which indicates a child of this node will split later
            # (and so this node is a splitting event).
            if (e_mat[:, i + 1] == 0).all():
                child.add_feature("is_sample", True)
            else:
                child.add_feature("is_sample", False)
                child.add_child(dist=0), child.add_child(dist=0)
        else:
            child.add_child(dist=0), child.add_child(dist=0)
        labelled_nodes.append(child)
    longer_tips(root)

    if not isochronous:
        # Label last tip.
        j = e_mat[-1, :].argmax()
        parent = labelled_nodes[j]
        child = next(c for c in parent.children if unlabelled(c))
        child.add_feature("rank_label", e_mat.shape[0])
        child.add_feature("is_sample", True)

    return root


def longer_tips(tree, length=1):
    """Increase the branch length of all non-sampling event tips by length."""
    for leaf in tree.get_leaves():
        if not getattr(leaf, "is_sample", False):
            leaf.dist += length


def tree_to_f_mat(tree, isochronous=False):
    """
    Return the F-matrix for the given ranked tree shape.
    """
    return d_mat_to_f_mat(tree_to_d_mat(tree, isochronous=isochronous))


def tree_to_d_mat(tree, isochronous=False):
    """
    Return the D-matrix for the given ranked tree shape.
    """
    return e_mat_to_d_mat(tree_to_e_mat(tree, isochronous=isochronous))


def tree_to_e_mat(tree, isochronous=False):
    """
    Return the E-matrix for the given ranked tree shape.
    """
    nodes = [n for n in tree.get_descendants() if hasattr(n, "rank_label")]
    nodes.sort(key=lambda node: node.rank_label)
    m = len(nodes)
    if isochronous:
        m += 1
    e_mat = torch.zeros((m, m))
    for row, node in enumerate(nodes):
        column = node.up.rank_label
        e_mat[row, column] = 1
    if isochronous:
        for leaf in tree.get_leaves():
            column = leaf.up.rank_label
            e_mat[m - 1, column] += 1
    return e_mat


def render_ranked_tree_shape(tree, file_path=None, display=True, space=10):
    """
    Return a copy of tree with styles set for pretty printing. Optionally display the
    tree or write to file.

    When calling from a Jupyter notebook, do something like this to display as expected:
        tree, style = render_ranked_tree_shape(tree, display=False)
        display(tree.render("%%inline", tree_style=style))

    Parameters:
        tree (ete3.Tree): The tree.
        file_path (str): Optional file path of where to write the tree. Ete3 guesses the
            format from the file extension.
        display (bool): When true, display the tree on screen. This isn't what you want
            when calling from a Jupyter notebook.
        space (int): The minimum number of pixels between adjacent branches of the tree.
    """
    tree = tree.copy()
    ts = TreeStyle()
    ts.show_branch_length = False
    ts.rotation = 90
    ts.min_leaf_separation = space
    ts.show_scale = False
    ns = NodeStyle()
    ns["size"] = 0
    attr = AttrFace("rank_label")
    attr.rotation = -90
    for node in tree.traverse():
        node.set_style(ns)
        if hasattr(node, "rank_label"):
            node.add_face(attr, column=0, position="branch-top")
    if display:
        tree.show(tree_style=ts)
    if file_path is not None:
        tree.render(file_path, tree_style=ts)
    return tree, ts


def sample_heterochronous_f_matrix(n):
    """
    Randomly sample an F-matrix for a fully heterochronous ranked tree shape with n tips
    using a coalescent process on isochronous ranked tree shapes with 2n tips all
    forming cherries. Return a pair consisting of the F-matrix and the probability of
    the matrix under the coalescent model.

    Parameters:
        n (int): Number of tips in the heterochronous ranked tree shape.
    """
    # Initialize the F-matrix and fill the last column. In terms of ranked tree shapes,
    # this is merge two leaves for the isochronous tree and mark the merge with rank
    # 2n-2; this is sample the node with rank 2n-2 in the heterochronous tree.
    f_mat = torch.zeros((2 * n - 2, 2 * n - 2))
    f_mat[2 * n - 3, 2 * n - 3] = 1
    n_tips_to_merge = 2 * n - 2
    n_ranked_nodes_to_merge = 1
    ranked_nodes_to_merge = [2 * n - 2]

    prob = 1.0
    for rank in range(1, 2 * n - 2)[::-1]:
        # Decide if the node with the given rank is merging of leaves (heterochronous
        # sample event) or merging of internal nodes (heterochronous bifurcation event).
        # Based on this choice, fill in column (rank-1). The offsets follow from the
        # definition of an F-matrix. The first iteration is necessarilly a leaf merge.
        f_mat[rank:, rank - 1] = f_mat[rank:, rank]

        tip_weight = n_tips_to_merge * (n_tips_to_merge - 1) / 2
        internal_weight = n_ranked_nodes_to_merge * (n_ranked_nodes_to_merge - 1) / 2
        denom = tip_weight + internal_weight
        leaf_merge_prob = tip_weight / denom

        if random() < leaf_merge_prob:
            prob *= leaf_merge_prob
            n_tips_to_merge -= 2
            n_ranked_nodes_to_merge += 1

            f_mat[rank - 1, rank - 1] = f_mat[rank, rank] + 1
        else:
            prob *= 1 / denom
            n_ranked_nodes_to_merge -= 1
            left, right = sorted(choice(ranked_nodes_to_merge, size=2, replace=False))
            ranked_nodes_to_merge.remove(left)
            ranked_nodes_to_merge.remove(right)

            f_mat[rank - 1, rank - 1] = f_mat[rank, rank] - 1
            f_mat[rank:left, rank - 1] -= 2
            f_mat[left:right, rank - 1] -= 1
        ranked_nodes_to_merge.append(rank)

    return f_mat, prob


def heterochronous_f_matrix_coalescent_probability(f_mat):
    n = f_mat.shape[0] // 2 + 1
    n_tips_to_merge = 2 * n - 2
    n_ranked_nodes_to_merge = 1
    prob = 1.0
    for rank in range(1, 2 * n - 2)[::-1]:
        tip_weight = n_tips_to_merge * (n_tips_to_merge - 1) / 2
        internal_weight = n_ranked_nodes_to_merge * (n_ranked_nodes_to_merge - 1) / 2
        denom = tip_weight + internal_weight
        if f_mat[rank - 1, rank - 1] - f_mat[rank, rank] == 1:
            n_ranked_nodes_to_merge += 1
            n_tips_to_merge -= 2
            prob *= tip_weight / denom
        else:
            n_ranked_nodes_to_merge -= 1
            prob *= 1 / denom

    return prob


def heterochronous_ranked_tree_shape_from_time_tree(tree, inplace=True):
    """
    Sample a heterochronous ranked tree shape from a bifurcating tree with branch
    lengths. The branch lengths of nodes at a common height are modified by a small
    amount so that no nodes appear at the same height, which corresponds to a unique
    heterochronous ranked tree shape. When there are no branches of length zero, the
    sampling is uniform, but otherwise may not be uniform. See comments in the method
    body for details.

    Parameters:
        inplace (bool): Modify the tree in place or return a copy with modified branch
            lengths and nodes with the "rank_label" attribute.
    """
    # By setting all branch lengths to zero, we get a random sample from all
    # heterochronous ranked tree shapes on a common rooted topology. Doing so uniformly
    # is difficult.
    #
    # Uniform sampling when there are branches of length zero corresponds to uniform
    # sampling of linear extensions of a certain poset (the poset of the nodes to order,
    # with the partial ordering from the tree). For a general poset on n nodes, this can
    # be done with a MCMC setup where the distribution is mixed after O(n^3log(n))
    # steps. Improvements for restricted classes of posets is difficult. For posets of
    # height 2 with maximum degree \delta, the distribution is mixed after
    # O(n\delta^2log(n)) steps. Our posets are of arbitrary height. Counting the number
    # of linear extensions is also difficult. Igor Pak and his students work on stuff
    # like this.
    #
    # Redrawing permutations until we find a valid one would force uniform sampling,
    # but we don't know how long that will take on average unless we have a lower bound
    # on the number of linear extensions.
    #
    # When there are no branches of length zero, the posets are graphs with no edges,
    # which is easy.
    if not inplace:
        tree = tree.copy()
    tree.add_feature("rank_label", 0)
    if any(c.dist == 0 for c in tree.children):
        for c in tree.children:
            c.dist += 1

    common_lengths = defaultdict(set, {0: {tree}})
    for node in tree.get_descendants():
        common_lengths[tree.get_distance(node)].add(node)
    dists = sorted(common_lengths.keys(), reverse=True)
    if len(dists) == 1:
        epsilon = dists[0]
    else:
        epsilon = min(dists[i - 1] - dists[i] for i in range(1, len(dists)))
    if epsilon == 0:
        epsilon = 1

    rank = len(tree.get_descendants())
    for dist in dists:
        tied_nodes = common_lengths[dist]
        count = len(tied_nodes)
        epsilon0 = epsilon / count

        nodes_with_children = {
            n for n in tied_nodes if not tied_nodes.isdisjoint(n.children)
        }

        first_level = list(tied_nodes.difference(nodes_with_children))
        perm = np.random.permutation(len(first_level))
        first_level = [first_level[i] for i in perm]
        nodes_by_subtree_height = [first_level]
        # We partition tied nodes with nodes_by_subtree_height, where the set at index i
        # consists of the nodes whose rooted subtree at the node intersected with tied
        # nodes is height i+1.
        while len(nodes_with_children) != 0:
            this_level = [
                n
                for n in nodes_with_children
                if nodes_with_children.isdisjoint(n.children)
            ]
            # Nodes at this_level are uniformly randomly permuted.
            perm = np.random.permutation(len(this_level))
            this_level = [this_level[i] for i in perm]
            nodes_by_subtree_height.append(this_level)
            nodes_with_children.difference_update(this_level)

        # Create a random ordered of the tied nodes while ensuring children (larger
        # rank) appear after parents (smaller rank).
        ordered_nodes = nodes_by_subtree_height[0]
        for this_level in nodes_by_subtree_height[1:]:
            for node in this_level:
                last_index = min(
                    ordered_nodes.index(c) for c in node.children if c in ordered_nodes
                )
                index = np.random.randint(0, high=last_index + 1)
                ordered_nodes.insert(index, node)

        # Assign ranks and adjust branch lengths.
        for i, node in enumerate(ordered_nodes[::-1], 1):
            node.add_feature("rank_label", rank)
            rank -= 1
            node.dist += (count - i) * epsilon0
            for c in node.children:
                c.dist -= (count - i) * epsilon0

    return None if inplace else tree
