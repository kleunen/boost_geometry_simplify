#include "simplify.hpp"
#include <iostream>

typedef boost::geometry::model::d2::point_xy<double> Point;
typedef boost::geometry::model::polygon<Point> Polygon;
typedef boost::geometry::model::multi_polygon<Polygon> MultiPolygon;

int main()
{
	MultiPolygon mp;
	boost::geometry::read_wkt("MULTIPOLYGON(((-3.69338090482 55.715692269,-3.6926079 55.7159519,-3.6924337 55.7151986,-3.6930661 55.715117,-3.6930157 55.7147726,-3.69338091231 55.7147132692,-3.69338090482 55.715692269)),((-3.69338090688 55.7142211068,-3.6933768 55.7141871,-3.6929634 55.7141725,-3.6929954 55.7145025,-3.6923885 55.7145554,-3.6923617 55.7147792,-3.6914998 55.7150295,-3.6912315 55.7134982,-3.6893145 55.7140656,-3.6885759 55.7143373,-3.6878468 55.7144903,-3.6885649 55.7161566,-3.6879954 55.7164484,-3.6873303 55.7153805,-3.686997 55.7156142,-3.6864469 55.716,-3.686103 55.716241,-3.6862139 55.7163757,-3.6858363 55.7166538,-3.6853115 55.7167809,-3.6849156 55.7164513,-3.6845733 55.7165856,-3.6839978 55.7159007,-3.6839527 55.7158282,-3.6835296 55.7160357,-3.6808072 55.7170818,-3.6807411 55.7168675,-3.6805504 55.7150641,-3.682455 55.7150385,-3.6834742 55.7143695,-3.6845373 55.7140717,-3.6845125 55.7139033,-3.684369 55.7127279,-3.6806985 55.7123158,-3.6810755 55.7133089,-3.6813647 55.7139107,-3.681652 55.7137697,-3.6820175 55.7147227,-3.6805049 55.7147351,-3.6795311 55.7129405,-3.6778782 55.7129041,-3.6777658 55.714827,-3.6760112 55.7159293,-3.6752767 55.7124876,-3.6715262 55.7125441,-3.6712776 55.7122602,-3.6708886 55.711493,-3.6708159 55.7113861,-3.6706013 55.7110893,-3.670957 55.7111165,-3.6731425 55.7111427,-3.6742991 55.709201,-3.6753334 55.7095268,-3.6758595 55.710859,-3.6788138 55.7112687,-3.676502 55.7081632,-3.6780904 55.7067399,-3.6788205 55.7067314,-3.6792276 55.707739,-3.6793662 55.7086275,-3.6806176 55.7088675,-3.6804964 55.7104918,-3.6803571 55.7111529,-3.6809193 55.7114274,-3.6806673 55.7120776,-3.6843319 55.712513,-3.6853143 55.7126239,-3.6900574 55.7130953,-3.6914546 55.71325,-3.6915457 55.7128876,-3.6919092 55.712801,-3.6913866 55.7112063,-3.6933809 55.7105007,-3.69338090688 55.7142211068)),((-3.6933809 55.7105007,-3.6933589 55.7104093,-3.6925489 55.710479,-3.6924738 55.7097065,-3.6933482 55.7096368,-3.69338090149 55.7095931641,-3.6933809 55.7105007)),((-3.69338091082 55.7090016493,-3.692978 55.7091164,-3.692286 55.7091593,-3.6917755 55.7087839,-3.6906731 55.7089272,-3.6902922 55.7073596,-3.6908162 55.7071423,-3.6910434 55.7069466,-3.6912775 55.7066595,-3.6908162 55.7064289,-3.690583 55.7055368,-3.6920467 55.7047551,-3.6918633 55.7046247,-3.6926388 55.7037173,-3.69338090692 55.7037081763,-3.69338091082 55.7090016493)))", mp);

    double distance = 0.0006;

    MultiPolygon simplified;
    simplify(mp, simplified, distance);

    MultiPolygon boost_simplified;
    boost::geometry::simplify(mp, boost_simplified, distance);

	std::string message;
	if(!boost::geometry::is_valid(boost_simplified, message)) {
		std::cout << "Boost simplify generated invalid geometry: " << message << std::endl;
	}

    std::cout << std::setprecision(12) << std::boolalpha 
        << "simplified: valid: " << boost::geometry::is_valid(simplified) << ", wkt: " << boost::geometry::wkt(simplified) << ", area: " << boost::geometry::area(simplified) << std::endl
        << "boost     : valid: " << boost::geometry::is_valid(boost_simplified) << ", wkt: " << boost::geometry::wkt(boost_simplified) << ", area: " << boost::geometry::area(boost_simplified) << std::endl;

	return 0;
}
