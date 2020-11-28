#ifndef CGAL_CLUSTERGRABBER_H
#define CGAL_CLUSTERGRABBER_H

#include <boost/function_output_iterator.hpp>

#include <CGAL/Polygon_2.h>

// cluster grabbing class:
// point p -> cluster {p}
// segment pq -> cluster {p,q}
// polygon -> cluster with the same vertices as polygon

namespace CGAL{

namespace internal{

template <class Kernel, class output_iterator>
class Cluster_grabber{
  output_iterator out;

private:
  typedef CGAL::Polygon_2<Kernel> Cluster_2;

public:
  Cluster_grabber(output_iterator it):out(it){}

  void operator()(const typename Kernel::Point_2& p){
    // a cluster with one point
    Cluster_2 cluster;
    cluster.push_back(p);
    *out++ = cluster;
  }

  void operator()(const typename Kernel::Segment_2& s){
    // a cluster with two points
    Cluster_2 cluster;
    cluster.push_back(s[0]);
    cluster.push_back(s[1]);
    *out++ = cluster;
  }

  template<class Container>
  void operator()(const CGAL::Polygon_2<Kernel,Container>& p){
    // a cluster with more than two points
    Cluster_2 cluster;
    for(typename CGAL::Polygon_2<Kernel,Container>::Vertex_iterator it=
        p.vertices_begin();it!=p.vertices_end();++it) {
      cluster.push_back(*it);
    }
    *out++ = cluster;
  }
};

template<class Kernel,class output_iterator>
boost::function_output_iterator<Cluster_grabber<Kernel,output_iterator> >
cluster_grabber(output_iterator it){
  return boost::make_function_output_iterator(
      Cluster_grabber<Kernel,output_iterator>(it));
}


}//internal
}//CGAL
#endif //CGAL_GRABBER_H
