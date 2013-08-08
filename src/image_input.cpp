// stl
#include <string>

// boost
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

// vigra

// pgmlink
#include <pgmlink/image_input.h>


namespace pgmlink {
  namespace fs = boost::filesystem;

  ////
  //// ImageSelectorAll
  ////
  ImageSelectorAll::~ImageSelectorAll() {

  }

  
  FilenameListPtr ImageSelectorAll::select(const std::string& directory) {
    FilenameListPtr filenames(new FilenameList);
    fs::path directory_path(directory);
    if (!fs::exists(directory_path)) {
      return filenames;
    }
    fs::directory_iterator end;
    for (fs::directory_iterator it(directory_path);
         it != end;
         ++it) {
      filenames->push_back(it->path().string());
    }
    return filenames;
  }



  ////
  //// ImageSelectorType
  ////
  ImageSelectorType::ImageSelectorType(const std::string& type) :
    type_(type) {
    
  }


  ImageSelectorType::~ImageSelectorType() {

  }


  FilenameListPtr ImageSelectorType::select(const std::string& directory) {
    FilenameListPtr filenames(new FilenameList);
    fs::path directory_path(directory);
    fs::path valid_extension(type_);
    if (!fs::exists(directory_path)) {
      return filenames;
    }
    fs::directory_iterator end;
    for (fs::directory_iterator it(directory_path);
         it != end;
         ++it) {
      fs::path extension = fs::extension(it->path());
      if (extension == valid_extension) {
        filenames->push_back(it->path().string());
      }
    }
    return filenames;
  }

}
