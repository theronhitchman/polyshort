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

    def visual(self, aspect_ratio=1, **kwds):
        """
        Return graphics object which contains a plot of the polygon.
        """
        return line(self.vertex_list + [self.vertex_list[0]], aspect_ratio=aspect_ratio, **kwds)

    def visual_with_curvatures(self, **kwds):
        """
        Return graphics object which contains a plot of the polygon and its
        curvature vectors, located at the appropriate vertices.
        """
        G = Graphics()
        G += self.visual(**kwds)
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

    def compute_flow(self, t_span=[0,1],num_points=1000,**kwds):
        """
        Return a a list of the form [(t0,[y0,...,yn]),...] containing solution
        data from running the incenter curvature flow on the polygon.
        """
        ODE = ode_solver()
        ODE.y_0 = self.vertices_to_list()
        ODE.function = sageobj('lambda t,y: '+incenter_odes_str(self.number_sides))
        ODE.ode_solve(t_span=t_span, num_points=num_points,**kwds)
        return ODE.solution

    def animate_flow(self, t_span=[0,1],num_points=1000,**kwds):
        """
        Returns an animation of the incenter flow.
        """
        raw_sols = self.compute_flow(t_span=t_span, num_points=num_points)
        xs = [x for (x,y) in self.vertex_list]
        ys = [y for (x,y) in self.vertex_list]
        xM, xm = max(x for x in xs), min(x for x in xs)
        yM, ym = max(y for y in ys), min(y for y in ys)
        pics = [PlanarPolygon(list_to_two_tuples(raw_sols[k][1])).visual_with_curvatures(xmin=xm,xmax=xM,ymin=ym,ymax=yM) for k in range(num_points)]
        the_show = animate(pics)
        return the_show

# Helper Functions used to set up the solution of differential equations.

def list_to_two_tuples(lst):
    """
    Process a flat list to a list of 2-tuples.
    """
    assert(len(lst) % 2 == 0), "list must have an even number of elements"
    verts = []
    if len(lst) == 2:
        return [(lst[0],lst[1])]
    else:
        return [(lst[0],lst[1])] + list_to_two_tuples(lst[2:])

curvature_coordinate0 = \
'2*((y[{4}] - y[{2}])/sqrt((y[{4}] - y[{2}])^2 + (y[{5}] - y[{3}])^2) - (y[{2}] - y[{0}])/sqrt((y[{2}] - y[{0}])^2 + (y[{3}] - y[{1}])^2))/(sqrt((y[{4}] - y[{2}])^2 + (y[{5}] - y[{3}])^2) + sqrt((y[{2}] - y[{0}])^2 + (y[{3}] - y[{1}])^2))'

curvature_coordinate1 = \
'2*((y[{5}] - y[{3}])/sqrt((y[{4}] - y[{2}])^2 + (y[{5}] - y[{3}])^2) - (y[{3}] - y[{1}])/sqrt((y[{2}] - y[{0}])^2 + (y[{3}] - y[{1}])^2))/(sqrt((y[{4}] - y[{2}])^2 + (y[{5}] - y[{3}])^2) + sqrt((y[{2}] - y[{0}])^2 + (y[{3}] - y[{1}])^2))'

curvature_vector_pair = curvature_coordinate0 + ', ' + curvature_coordinate1

def incenter_odes_str(num_verts):
    """
    Returns a string used to construct the ODE's defining the incenter flow for
    a planar polygon with `num` sides.
    """
    m = 2*num_verts
    ode_list = '[' + curvature_vector_pair.format(str(m-2),str(m-1),str(0),str(1),str(2),str(3))
    for ind in range(1,num_verts-1):
        ode_list += ', ' + curvature_vector_pair.format(str(2*ind-2),str(2*ind-1),str(2*ind),str(2*ind+1),str(2*ind+2),str(2*ind+3))
    ode_list += ', ' + curvature_vector_pair.format(str(m-4),str(m-3),str(m-2),str(m-1),str(0),str(1)) + ']'
    return ode_list
