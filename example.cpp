
#include <boost/geometry.hpp>

typedef boost::geometry::model::d2::point_xy<double> Point; 
typedef boost::geometry::model::polygon<Point> Polygon;
typedef boost::geometry::ring_type<Polygon>::type Ring;
typedef boost::geometry::model::multi_polygon<Polygon> MultiPolygon;

typedef boost::geometry::model::segment<Point> simplify_segment;
typedef boost::geometry::index::rtree<simplify_segment, boost::geometry::index::quadratic<16>> simplify_rtree;

// Combine overlapping elements by performing a union
template<typename C, typename T>
void simplify_combine(C &result, T &&new_element)
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

struct simplify_rtree_counter
{
	using value_type = simplify_segment;
	std::size_t count = 0;
	void push_back(value_type const &) { ++count; }
	std::size_t size() const { return count; }
};

template<typename GeometryType>
void simplify(GeometryType const &input, GeometryType &output, double max_distance, simplify_rtree const &outer_rtree = simplify_rtree())
{        
	simplify_rtree rtree;

    std::deque<std::size_t> nodes(input.size());
    for(std::size_t i = 0; i < input.size(); ++i) 
        nodes[i] = i;
    for(std::size_t i = 0; i < input.size() - 1; ++i)
        rtree.insert({ input[i], input[i + 1] });    

    std::priority_queue<std::size_t, std::vector<size_t>> pq;
    for(std::size_t i = 0; i < input.size() - 2; ++i) 
        pq.push(i);      
        
    while(!pq.empty()) {
        auto entry = pq.top();
        pq.pop();
        
        auto start = nodes[entry];
        auto middle = nodes[entry + 1];
        auto end = nodes[entry + 2];

        simplify_segment line(input[start], input[end]);
        double distance = 0.0;
        for(auto i = start + 1; i < end; ++i) 
            distance = std::max(distance, boost::geometry::distance(line, input[i]));          
    
        if(distance < max_distance) {
            simplify_rtree_counter result;
            boost::geometry::index::query(rtree, boost::geometry::index::intersects(line), std::back_inserter(result));
            boost::geometry::index::query(outer_rtree, boost::geometry::index::intersects(line), std::back_inserter(result));

            std::size_t query_expected = ((start == 0 || end == input.size() - 1) ? 2 : 4);
            if(result.size() == query_expected) {
                nodes.erase(nodes.begin() + entry + 1);
                rtree.remove(simplify_segment(input[start], input[middle]));
                rtree.remove(simplify_segment(input[middle], input[end]));
                rtree.insert(line);
        
                if(entry + 2 < nodes.size()) {
                    pq.push(start);             
                }
            }
        }
    }
    
    for(auto i: nodes)
        boost::geometry::append(output, input[i]);
}

void simplify(Polygon const &p, Polygon &result, double max_distance) 
{
	simplify_rtree outer_rtree;
	for(std::size_t j = 0; j < p.outer().size() - 1; ++j) 
		outer_rtree.insert({ p.outer()[j], p.outer()[j + 1] });    

	std::vector<Ring> new_inners;
	for(auto const &inner: p.inners()) {
		Ring new_inner;
		simplify(inner, new_inner, max_distance, outer_rtree);

		std::reverse(new_inner.begin(), new_inner.end());
		if(boost::geometry::area(new_inner) > max_distance * max_distance) {
			simplify_combine(new_inners, std::move(new_inner));
		}
	}

	simplify_rtree inners_rtree;
	for(auto const &inner: new_inners) {
		for(std::size_t z = 0; z < inner.size() - 1; ++z) 
			inners_rtree.insert({ inner[z], inner[z + 1] });    
	} 

	simplify(p.outer(), result.outer(), max_distance, inners_rtree);
	if(boost::geometry::area(result.outer()) < max_distance * max_distance) {
		return;
	}

	for(auto&& r: new_inners) {
		std::reverse(r.begin(), r.end());
		result.inners().push_back(r);
	} 
}

void simplify(MultiPolygon const &mp, MultiPolygon &result, double max_distance) 
{
	for(auto const &p: mp) {
		Polygon new_p;
		simplify(p, new_p, max_distance);
    	if(!new_p.outer().empty()) {
			simplify_combine(result, std::move(new_p));
		}
	}
}

#include <iostream>

int main()
{
	MultiPolygon mp;
	boost::geometry::read_wkt("MULTIPOLYGON(((0 10,10 10,10 0,0 0,0 10),(1.1 1.1,3.1 1.9,4.9 1.1,3.1 3.1,2.5 2.1,1.1 1.1)))", mp);

    // Simplify it, using distance of 0.5 units
    MultiPolygon simplified;
    boost::geometry::simplify(mp, simplified, 0.5);

    std::cout
        << "  original: " << boost::geometry::wkt(mp) << ", valid: " << std::boolalpha << boost::geometry::is_valid(mp) << std::endl
        << "simplified: " << boost::geometry::wkt(simplified) << ", valid: " << std::boolalpha << boost::geometry::is_valid(simplified) << std::endl; 
    
	return 0;
}
