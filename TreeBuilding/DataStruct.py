

class Tree:
    def __init__(self, name):
        self.dist = 0  # distance to parent node
        self.name = name
        self.children = []

    def add_child(self, treenode):
        self.children.append(treenode)

    def is_leaf(self):
        """
        Return True if current node is a leaf.
        """
        return len(self.children) == 0

    def write(self, features=None, outfile=None, format=0, is_leaf_fn=None,
              format_root_node=False, dist_formatter=None, support_formatter=None,
              name_formatter=None):
        """
        Returns the newick representation of current node. Several
        arguments control the way in which extra data is shown for
        every node:
        :argument features: a list of feature names to be exported
          using the Extended Newick Format (i.e. features=["name",
          "dist"]). Use an empty list to export all available features
          in each node (features=[])
        :argument outfile: writes the output to a given file
        :argument format: defines the newick standard used to encode the
          tree. See tutorial for details.
        :argument False format_root_node: If True, it allows features
          and branch information from root node to be exported as a
          part of the newick text string. For newick compatibility
          reasons, this is False by default.
        :argument is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.

        **Example:**
        ::
             t.get_newick(features=["species","name"], format=1)
        """

        nw = write_newick(self, features=features,
                          format=format,
                          is_leaf_fn=is_leaf_fn,
                          format_root_node=format_root_node,
                          dist_formatter=dist_formatter,
                          support_formatter=support_formatter,
                          name_formatter=name_formatter)

        if outfile is not None:
            open(outfile, "w").write(nw)
        else:
            return nw

    def iter_prepostorder(self, is_leaf_fn=None):
        """
        Iterate over all nodes in a tree yielding every node in both
        pre and post order. Each iteration returns a postorder flag
        (True if node is being visited in postorder) and a node
        instance.
        """
        to_visit = [self]
        if is_leaf_fn is not None:
            _leaf = is_leaf_fn
        else:
            _leaf = self.__class__.is_leaf

        while to_visit:
            node = to_visit.pop(-1)
            try:
                node = node[1]
            except TypeError:
                # PREORDER ACTIONS
                yield (False, node)
                if not _leaf(node):
                    # ADD CHILDREN
                    to_visit.extend(reversed(node.children + [[1, node]]))
                else:
                    #POSTORDER ACTIONS
                    yield (True, node)



def write_newick(rootnode, features=None, format=1, format_root_node=True,
                 is_leaf_fn=None, dist_formatter=None, support_formatter=None,
                 name_formatter=None):
    """ Iteratively export a tree structure and returns its NHX
    representation. """
    newick = []
    leaf = is_leaf_fn if is_leaf_fn else lambda n: not bool(n.children)
    for postorder, node in rootnode.iter_prepostorder(is_leaf_fn=is_leaf_fn):
        if postorder:
            newick.append(")")
            if node.up is not None or format_root_node:
                newick.append(format_node(node, "internal", format,
                                          dist_formatter=dist_formatter,
                                          support_formatter=support_formatter,
                                          name_formatter=name_formatter))
                newick.append(_get_features_string(node, features))
        else:
            if node is not rootnode and node != node.up.children[0]:
                newick.append(",")

            if leaf(node):
                safe_name = re.sub("["+_ILEGAL_NEWICK_CHARS+"]", "_", \
                               str(getattr(node, "name")))
                newick.append(format_node(node, "leaf", format,
                              dist_formatter=dist_formatter,
                              support_formatter=support_formatter,
                              name_formatter=name_formatter))
                newick.append(_get_features_string(node, features))
            else:
                newick.append("(")

    newick.append(";")
    return ''.join(newick)
