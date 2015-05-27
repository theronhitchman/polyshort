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

    def vertex_as_vector(self,k):
        """
        Return vertex with given index as a vector.

        Indexing is cyclic, so input can be any integer.
        """
        m = k % self.number_sides
        return vector(self.vertex_list[m])

    def side_as_vector(self,k):
        """
        Return side joining vertex k-1 to vertex k as a vector.

        Indexing is cyclic, so k can be any integer.
        """
        return self.vertex_as_vector(k) - self.vertex_as_vector(k-1)

    def side_length(self,k):
        """
        Return length of side joining vertex k-1 to vertex k.

        Indexing is cyclic, so k can be any integer.
        """
        #assert isinstance(k, Integer), "Argument must be an integer"
        return self.side_as_vector(k).norm()

    def perimeter(self):
        """
        Return the sum of lengths of sides of the polygon.
        """
        return sum(self.side_length(k) for k in range(self.number_sides))

    def tangent_vector(self,k):
        """
        Return the unit vector in the direction of side with given index.
        """
        return self.side_as_vector(k).normalized()

    def exterior_angle(self,k):
        """
        Return the radian measure of the exterior angle at the given vertex.
        """
        return arccos(self.tangent_vector(k+1).inner_product(self.tangent_vector(k)))

    def exterior_angle_degrees(self,k,num_digits=10):
        """
        Return the degree measure of the exterior angle at the given vertex.
        """
        return ((180*self.exterior_angle(k))/pi).n(digits=num_digits)

    def curvature_vector(self,k):
        """
        Return the curvature vector at the given vertex.
        """
        return 2*(self.tangent_vector(k+1) - self.tangent_vector(k))/(self.side_length(k+1) + self.side_length(k))

    def curvature(self,k):
        """
        Return the curvature at the given vertex.
        """
        return norm(self.curvature_vector(k))

    def visual(self):
        """
        Return graphics object which contains a plot of the polygon.
        """
        return line(self.vertex_list + [self.vertex_list[0]])

    def visual_with_curvatures(self):
        """
        Return graphics object which contains a plot of the polygon and its
        curvature vectors, located at the appropriate vertices.
        """
        G = Graphics()
        G += self.visual()
        for k in range(self.number_sides):
            G += arrow(self.vertex_list[k],
                       vector(self.vertex_list[k]) + self.curvature_vector(k),
                       color='red', linestyle=":")
        return G

    def side_normal(self,k):
        """
        Return vector normal to given side.
        """
        w = self.side_as_vector(k)
        return vector([-w[1],w[0]])


    def side_functional(self,k):
        """
        Return linear functional that vanishes along side k.
        """
        return lambda x,y : self.side_normal(k).dot_product(vector([x,y])-vector(self.vertex_list[k-1]))

    def side_is_support(self,k):
        """
        Return a Boolean indicating if the given side supports the polygon.
        """
        f = self.side_functional(k)
        for l in range(self.number_sides):
            x = self.vertex_list[l]
            for m in range(l,self.number_sides):
                y = self.vertex_list[m]
                if f(x[0],x[1]) * f(y[0],y[1]) < 0:
                    return False
        return True

    def is_convex(self):
        """
        Returns a Boolean indicating if the polygon is convex.
        """
        for k in range(self.number_sides):
            if not self.side_is_support(k):
                return False
        return True

    def vertices_to_list(self):
        """
        Returns a flat list of coordinates for the vertices of the polygon.
        """
        flat_list = []
        for vert in self.vertex_list:
            flat_list += [vert[0],vert[1]]
        return flat_list

    def bundled_curvature_list(self):
        curv_list = []
        for m in range(self.number_sides):
            curv_list += [self.curvature_vector(m)[0], self.curvature_vector(m)[1]]
        return curv_list

def list_to_two_tuples(lst):
    """
    Process a flat list to a list of 2-tuples.
    """
    assert(len(lst) % 2 == 0)
    verts = []
    if len(lst) == 2:
        return [(lst[0],lst[1])]
    else:
        return [(lst[0],lst[1])] + list_to_two_tuples(lst[2:])

curvature_coordinate0 = \
'-2*((y[{0}] - y[{2}])/sqrt((y[{0}] - y[{2}])^2 + (y[{1}] - y[{3}])^2) + (y[{0}] - y[{4}])/sqrt((y[{0}] - y[{4}])^2 + (y[{1}] - y[{5}])^2))/(sqrt((y[{0}] - y[{2}])^2 + (y[{1}] - y[{3}])^2) + sqrt((y[{0}] - y[{4}])^2 + (y[{1}] - y[{5}])^2))'

curvature_coordiate1 = \
'-2*((y[{1}] - y[{3}])/sqrt((y[{0}] - y[{2}])^2 + (y[{1}] - y[{3}])^2) + (y[{1}] - y[{5}])/sqrt((y[{0}] - y[{4}])^2 + (y[{1}] - y[{5}])^2))/(sqrt((y[{0}] - y[{2}])^2 + (y[{1}] - y[{3}])^2) + sqrt((y[{0}] - y[{4}])^2 + (y[{1}] - y[{5}])^2))'
