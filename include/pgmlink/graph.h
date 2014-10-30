/**
   @file
   @ingroup util
   @brief lemon graph extensions
*/

#ifndef GRAPH_H
#define GRAPH_H

#include <map>
#include <stdexcept>
#include <string>
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace pgmlink {
  ////
  //// Graph Properties
  ////
  /** 
   Property map traits class
  */
  template<typename PropertyTag, typename Graph>
    struct property_map {
    typedef void type;
    static const std::string name;
  };

  
  
  ////
  //// Graph
  ////
  /**
     Graph with associated property maps
  */
  template <typename Graph>
    class PropertyGraph : public Graph {
  public:
    typedef Graph base_graph;

    std::map<std::string, boost::any> getProperties() const {
		return properties_;
	}

    template <typename PropertyTag>
      typename property_map<PropertyTag, Graph>::type &
      get(PropertyTag) const;
    
    template<typename PropertyTag>
      PropertyGraph& 
      add(PropertyTag);

    template<typename PropertyTag>
      bool 
      has_property(PropertyTag) const;

    /**
     * Attach existing property_map to graph.
     *
     * PropertyGraph takes ownership of the pointer.
     */
    template<typename PropertyTag>
      PropertyGraph&
      insert(PropertyTag, typename property_map<PropertyTag, Graph>::type*);

    static void copy(PropertyGraph<Graph> &src, PropertyGraph<Graph> &dest);
    
  private:
    typedef std::map<std::string, boost::any> properties_map;  
    properties_map properties_;
  };



  /******************/
  /* Implementation */
  /******************/
  
  template <typename Graph>
    template <typename PropertyTag>
    typename property_map<PropertyTag, Graph>::type &
    PropertyGraph<Graph>::get(PropertyTag) const {
  
    std::string name = property_map<PropertyTag, Graph>::name;
    
    if(properties_.count(name) > 0) {
      // internally, boost::any stores a (stable) pointer to its content
      // therefore, the following cast is safe, even though the std::map guarantees only
      // stability of iterators to elements and not actual memory addresses
      return *boost::any_cast<boost::shared_ptr<typename property_map<PropertyTag, Graph>::type> >(properties_.find(name)->second);
    } else {
      throw std::runtime_error("PropertyGraph::get(): property " + name + " not found");
    }
  }
  
  template <typename Graph>
    template <typename PropertyTag>
    PropertyGraph<Graph> & 
    PropertyGraph<Graph>::add(PropertyTag) {

    std::string name = property_map<PropertyTag, Graph>::name;
    // we have to use (shared) ptrs here, because some property maps have private copy constructors
    // and boost::any needs a copy of value during construction 
    boost::any property = boost::make_shared<typename property_map<PropertyTag, Graph>::type >(*this);
    properties_.insert(properties_map::value_type(name, property));
    return *this;
  }

  template <typename Graph>
    template <typename PropertyTag>
    bool 
    PropertyGraph<Graph>::has_property(PropertyTag) const {
    std::string name = property_map<PropertyTag, Graph>::name;
    return properties_.count(name) ? true : false;
  }
  
  template <typename Graph>
    template<typename PropertyTag>
    PropertyGraph<Graph> &
    PropertyGraph<Graph>::insert(PropertyTag, typename property_map<PropertyTag, Graph>::type* m) {

    std::string name = property_map<PropertyTag, Graph>::name;
    boost::any property = boost::shared_ptr<typename property_map<PropertyTag, Graph>::type>(m);
    properties_.insert(properties_map::value_type(name, property));
    return *this;	
    }

} /* namespace pgmlink */
#endif /* GRAPH_H */
