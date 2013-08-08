#ifndef IMAGE_INPUT_H
#define IMAGE_INPUT_H

// stl
#include <string>
#include <vector>

// vigra
#include <vigra/multi_array.hxx> // vigra::MultiArray
#include <vigra/impex.hxx> // image import and export
#include <vigra/imageinfo.hxx> // image import and export info

// boost
#include <boost/shared_ptr.hpp>


namespace pgmlink {
  class ImageSelectorBase;

  typedef boost::shared_ptr<ImageSelectorBase> ImageSelectorPtr;

  typedef std::vector<std::string> FilenameList;

  typedef boost::shared_ptr<FilenameList > FilenameListPtr;
  
  template <typename T, int N>
  class ImageRetrieverBase;

  template <typename T, int N>
  struct ImageRetrieverPtr {
    typedef boost::shared_ptr<ImageRetrieverBase<T, N> > type;
  };

  class ImageWriterBase;

  typedef boost::shared_ptr<ImageWriterBase> ImageWriterPtr;
  


  ////
  //// ImageSelectorBase
  ////
  class ImageSelectorBase {
  public:
    virtual FilenameListPtr select(const std::string& directory) = 0;
    virtual ~ImageSelectorBase() {};
  };


  ////
  //// ImageSelectorAll
  ////
  class ImageSelectorAll : public ImageSelectorBase{
  public:
    virtual FilenameListPtr select(const std::string& directory);
    virtual ~ImageSelectorAll();
  };


  ////
  //// ImageSelectorType
  ////
  class ImageSelectorType : public ImageSelectorBase{
  public:
    virtual FilenameListPtr select(const std::string& directory);
    ImageSelectorType(const std::string& type);
    virtual ~ImageSelectorType();
  private:
    std::string type_;
  };


  ////
  //// ImageRetrieverBase
  ////
  template <typename T, int N>
  class ImageRetrieverBase {
  public:
    virtual vigra::MultiArrayView<N, T> retrieve() = 0;
    virtual void set_current_image(const std::string& image) = 0;
    virtual ~ImageRetrieverBase() {}
  };


  ////
  //// VigraReader
  ////
  template <typename T, int N>
  class VigraReader : public ImageRetrieverBase<T, N> {
  public:
    virtual vigra::MultiArrayView<N, T> retrieve();
    virtual void set_current_image(const std::string& image);
    VigraReader(const std::string& filename);
    virtual ~VigraReader();
  private:
    std::string filename_;
  };


  ////
  //// HDF5Reader
  ////
  template <typename T, int N>
  class HDF5Reader : public ImageRetrieverBase<T, N> {
  public:
    virtual vigra::MultiArrayView<N, T> retrieve();
    virtual void set_current_image(const std::string& image);
    HDF5Reader(const std::string& filename,
               const std::string& path);
    virtual ~HDF5Reader();
  private:
    std::string filename_;
    std::string path_;
  };


  ////
  //// ImageWriterBase
  ////
  class ImageWriterBase {
  public:
    virtual void write() = 0;
    virtual void set_current_image(const std::string& filename) = 0;
    ~ImageWriterBase() {}
  };


  ////
  //// VigraWriter
  ////
  template <typename T, int N>
  class VigraWriter : public ImageWriterBase {
  public:
    virtual void write();
    virtual void set_current_image(const std::string& filename);
    VigraWriter(vigra::MultiArrayView<N, T> image,
                const std::string& filename);
    virtual ~VigraWriter();
  private:
    vigra::MultiArrayView<N, T> image_;
    std::string filename_;
  };


  ////
  //// HDF5Writer
  ////
  template <typename T, int N>
  class HDF5Writer : public ImageWriterBase {
  public:
    virtual void write();
    virtual void set_current_image(const std::string& filename) = 0;
    HDF5Writer(vigra::MultiArrayView<N, T> image,
               const std::string& filename,
               const std::string& path);
    virtual ~HDF5Writer();
  private:
    vigra::MultiArrayView<N, T> image_;
    std::string filename_;
    std::string path_;
  };
    


  /* IMPLEMENTATIONS */

  ////
  //// VigraReader
  ////
  template <typename T, int N>
  vigra::MultiArrayView<N, T> VigraReader<T, N>::retrieve() {
    vigra::ImageImportInfo info(filename_.c_str());
    vigra::MultiArray<N, T> image(info.shape());
    vigra::importImage(info, destImage(image));
    return vigra::MultiArrayView<N, T>(image);
  }


  template <typename T, int N>
  void VigraReader<T, N>::set_current_image(const std::string& image) {
    filename_ = image;
  }


  template <typename T, int N>
  VigraReader<T, N>::VigraReader(const std::string& filename) :
    filename_(filename) {
    // nothing to be done here
  }


  template <typename T, int N>
  VigraReader<T, N>::~VigraReader() {
    // nothing to be done here
  }


  ////
  //// VigraWriter
  ////
  template <typename T, int N>
  void VigraWriter<T, N>::write() {
    vigra::ImageExportInfo info(filename_.c_str());
    vigra::exportImage(image_, info);
  }


  template <typename T, int N>
  void VigraWriter<T, N>::set_current_image(const std::string& filename) {
    filename_ = filename;
  }


  template <typename T, int N>
  VigraWriter<T, N>::VigraWriter(vigra::MultiArrayView<N, T> image,
                                 const std::string& filename) :
    image_(image),
    filename_(filename) {
    // nothing to be done here
  }


  template <typename T, int N>
  VigraWriter<T, N>::~VigraWriter() {
    // nothing to be done here
  }
    

}

#endif /* IMAGE_INPUT_H */
