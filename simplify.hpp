#ifndef __SIMPLIFY_H__
#define __SIMPLIFY_H__

#include <boost/geometry.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>

namespace impl {

// Combine overlapping elements by performing a union
template<typename C, typename T>
static inline void simplify_combine(C &result, T &&new_element)
{
    result.push_back(new_element);

   	for(std::size_t i = 0; i < result.size() - 1; ) {
        if(!boost::geometry::intersects(result[i], result.back())) {
            ++i;
            continue;
        }

        std::vector<T> union_result;
        boost::geometry::union_(result[i], result.back(), union_result);

        if(union_result.size() != 1) {
			++i;
			continue;
		}

       	result.back() = std::move(union_result[0]);
       	result.erase(result.begin() + i);
    } 
}

template<
	typename point_t = boost::geometry::model::d2::point_xy<double>, 
	typename ring_t = boost::geometry::model::ring<point_t>,
	typename segment_t = boost::geometry::model::segment<point_t>,
	typename rtree_t = boost::geometry::index::rtree<segment_t, boost::geometry::index::quadratic<16>>
	>
static inline ring_t simplify_ring(ring_t const &input, double distance, bool is_closed = true, rtree_t const &outer_rtree = rtree_t())
{        
	ring_t output;

    std::deque<std::size_t> nodes(input.size());
    for(std::size_t i = 0; i < input.size(); ++i) 
        nodes[i] = i;

	rtree_t rtree(
		boost::irange<std::size_t>(0, input.size() - 1)
		| boost::adaptors::transformed([&input](std::size_t i) {
			return segment_t(input[i], input[i+1]);
		}));

	for(std::size_t pq = input.size() - 2; pq--; ) {
        auto entry = pq;
        
        auto start = nodes[entry];
        auto middle = nodes[entry + 1];
        auto end = nodes[entry + 2];

        segment_t line(input[start], input[end]);

        double max_comp_distance = 0.0;
		std::size_t max_comp_i = start + 1;
	
        for(auto i = start + 1; i < end; ++i) { 
			auto comp_distance = boost::geometry::comparable_distance(line, input[i]);
            if(comp_distance > max_comp_distance) {
				max_comp_distance = comp_distance;
				max_comp_i = i;
			}
		}
 
        if(boost::geometry::distance(line, input[max_comp_i]) < distance) {
			std::size_t query_count = 0;
            for(auto const &result: rtree | boost::geometry::index::adaptors::queried(boost::geometry::index::intersects(line)))
				++query_count;
            for(auto const &result: outer_rtree | boost::geometry::index::adaptors::queried(boost::geometry::index::intersects(line)))
				++query_count;
			
			std::size_t expected_count = std::min<std::size_t>(4, nodes.size() - 1);

			if(!is_closed && ((start == 0 || end == input.size() - 1))) 
				expected_count = std::min<std::size_t>(expected_count, 2);

            if(query_count == expected_count) {
                nodes.erase(nodes.begin() + entry + 1);
                rtree.remove(segment_t(input[start], input[middle]));
                rtree.remove(segment_t(input[middle], input[end]));
                rtree.insert(line);
            }
        }
    }
 
	output.resize(nodes.size());
	for(std::size_t i = 0; i < nodes.size(); ++i)
		output[i] = input[nodes[i]];   
	return output;
}

}

template<
	typename point_t = boost::geometry::model::d2::point_xy<double>, 
	typename segment_t = boost::geometry::model::segment<point_t>,
	typename rtree_t = boost::geometry::index::rtree<segment_t, boost::geometry::index::quadratic<16>>
	>
static inline void simplify(boost::geometry::model::polygon<point_t> const &p, boost::geometry::model::polygon<point_t> &output, double max_distance) 
{
	boost::geometry::model::polygon<point_t> result;

	rtree_t outer_rtree(
		boost::irange<std::size_t>(0, p.outer().size() - 1)
		| boost::adaptors::transformed([&p](std::size_t i) { 
			return segment_t(p.outer()[i], p.outer()[i+1]); 
		}));

	for(auto const &inner: p.inners()) {
		auto new_inner = impl::simplify_ring(inner, max_distance, true, outer_rtree);

		std::reverse(new_inner.begin(), new_inner.end());
		if(new_inner.size() > 3 && boost::geometry::perimeter(new_inner) > 3 * max_distance) {
			impl::simplify_combine(result.inners(), std::move(new_inner));
		}
	}

	rtree_t inners_rtree;
	for(auto &inner: result.inners()) {
		std::reverse(inner.begin(), inner.end());

		inners_rtree.insert(
			boost::irange<std::size_t>(0, inner.size() - 1)
			| boost::adaptors::transformed([&inner](std::size_t i) {
				return segment_t(inner[i], inner[i+1]);
			}));
	} 

	result.outer() = impl::simplify_ring(p.outer(), max_distance, true, inners_rtree);
	if(result.outer().size() > 3 && boost::geometry::perimeter(result.outer()) > 3 * max_distance) {
		output = std::move(result);
	}
}

template<
	typename point_t = boost::geometry::model::d2::point_xy<double>, 
	typename polygon_t = boost::geometry::model::polygon<point_t>
	>
static inline void simplify(boost::geometry::model::multi_polygon<polygon_t> const &mp, boost::geometry::model::multi_polygon<polygon_t> &result, double max_distance) 
{
	for(auto const &p: mp) {
		polygon_t new_p;
		simplify(p, new_p, max_distance);
    	if(!new_p.outer().empty()) {
			impl::simplify_combine(result, std::move(new_p));
		}
	}
}

template<
	typename point_t = boost::geometry::model::d2::point_xy<double>
	>
static inline void simplify(boost::geometry::model::linestring<point_t> const &ls, boost::geometry::model::linestring<point_t> &result, double max_distance) 
{
	if(ls.empty()) true;
	bool is_closed = boost::geometry::equals(ls.front(), ls.back());
	result = impl::simplify_ring(ls, max_distance, is_closed);
}

template<
	typename point_t = boost::geometry::model::d2::point_xy<double>
	>
static inline void simplify(boost::geometry::model::multi_linestring<point_t> const &mls, boost::geometry::model::multi_linestring<point_t> &result, double max_distance) 
{
	for(auto const &ls: mls) {
		if(ls.empty()) continue;

		bool is_closed = boost::geometry::equals(ls.front(), ls.back());
		auto new_ls = impl::simplify_ring(ls, max_distance, is_closed);

    	if(!new_ls.empty()) {
			if(boost::geometry::intersects(new_ls, result)) {
				boost::geometry::model::multi_linestring<point_t> output;
				boost::geometry::union_(new_ls, result, output);
				result = std::move(output);
			} else {
				result.push_back(std::move(new_ls));
			}
		}
	}
}
#endif
