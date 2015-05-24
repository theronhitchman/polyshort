r"""
A polygon class to support modeling of different polygon shortening flows.

<Paragraph description>

AUTHORS:

- Theron J Hitchman (2015-05-31): initial version

EXAMPLES::

<Lots and lots of examples>
"""

#*****************************************************************************
#       Copyright (C) 2015 Theron J Hitchman
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


class PlanarPolygon(SageObject):
    def __init__(self, vertices):
        """Construct a planar polygon."""
        assert isinstance(vertices, list), \
        "input must be a list of vertices, given as coordinate tuples"
        for vert in vertices:
            assert isinstance(vert, tuple) and len(vert) == 2, \
            "each vertex should be a coordinate tuple of length 2"
        self.vertex_list = vertices
        self.number_sides = len(vertices)

    def _repr_(self):
        return "A %d sided polygon with vertices at " % self.number_sides \
        + ", ".join(str(vert) for vert in self.vertex_list)

    def side_as_vector(self,k):
        """
        Return side joining vertex k-1 to vertex k as a vector.

        Indexing is cyclic, so k can be any integer.
        """
        assert isinstance(k, Integer), "Argument must be an integer"
        m = k % self.number_sides
        x1, y1 = self.vertex_list[m][0], self.vertex_list[m][1]
        x0, y0 = self.vertex_list[m-1][0], self.vertex_list[m-1][1]
        xcoord, ycoord = x1 - x0, y1 - y0
        return vector([xcoord, ycoord])

    def side_length(self,k):
        """
        Return length of side joining vertex k-1 to vertex k.

        Indexing is cyclic, so k can be any integer.
        """
        assert isinstance(k, Integer), "Argument must be an integer"
        return self.side_as_vector(k).norm()