#ifndef IMAGE_INPUT_H
#define IMAGE_INPUT_H

// stl
#include <string>
#include <vector>

// vigra
#include <vigra/multi_array.hxx> // vigra::MultiArray
#include <vigra/impex.hxx> // 

// boost
#include <boost/shared_ptr.hpp>


namespace pgmlink {
class ImageSelectorBase;

typedef boost::shared_ptr<ImageSelectorBase> ImageSelectorPtr;

  ////
  //// ImageSelectorBase
  ////
  class ImageSelectorBase {
  public:
    virtual boost::shared_ptr<std::vector<std::string> > select(std::string directory) = 0;
    virtual ~ImageSelectorBase();
  };


  ////
  //// ImageSelectorAll
  ////
  class ImageSelectorAll {
  public:
    virtual boost::shared_ptr<std::vector<std::string> > select(std::string directory);
    virtual ~ImageSelectorAll();
  };


  ////
  //// ImageSelectorType
  ////
  class ImageSelectorType {
  public:
    virtual boost::shared_ptr<std::vector<std::string> > select(std::string directory);
    ImageSelectorType(std::string type);
    virtual ~ImageSelectorType;
  private:
    std::string type_;
  };
}

#endif /* IMAGE_INPUT_H */
