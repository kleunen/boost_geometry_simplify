# boost_geometry_simplify

The following project provides an implementation of a self-intersection aware implementation of the simplify operation in boost geometry:
https://www.boost.org/doc/libs/1_76_0/libs/geometry/doc/html/geometry/reference/algorithms/simplify/simplify_3.html

The current implementation in boost geometry may generate polygons and multipolygons with self-intersection. This implementation makes sure, that when a polygon or multi-polygon is simplified that: 

* The simplified rings do not contain self-intersection (inners and outers)
* The outer does not overlap the inners
* In case one inner overlaps another inner, it is combined into a single inner using union
* In case one polygon in a multipolygon overlaps another polygon, it is combined into a single polygon using a union

For simplification of the ring, this library uses the approach as described in: https://www.jasondavies.com/simplify/
The end result is a polygon/multipolygon which is simplified and is valid. 
