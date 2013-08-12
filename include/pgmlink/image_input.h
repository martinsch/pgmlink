#ifndef IMAGE_INPUT_H
#define IMAGE_INPUT_H

// stl
#include <string>
#include <vector>
#include <map>

// vigra
#include <vigra/multi_array.hxx> // vigra::MultiArray
#include <vigra/impex.hxx> // image import and export
#include <vigra/imageinfo.hxx> // image import and export info

// boost
#include <boost/shared_ptr.hpp>

// pgmlink
#include <pgmlink/log.h>


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

  template <typename T, int N>
  class ImageWriterBase;

  template <typename T, int N>
  struct ImageWriterPtr {
    typedef boost::shared_ptr<ImageWriterBase<T, N> > type;
  };
  




  ////
  //// ImageOptions
  ////
  class ImageOptions {
  public:
    inline void set(const std::string& key,
                    const std::string& value) {
      options[key] = value;
    }
    inline std::string get(const std::string& key) {
      return options[key];
    }
  private:
    std::map<std::string, std::string> options;
  };



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
    virtual vigra::MultiArray<N, T> retrieve() = 0;
    ImageRetrieverBase() {}
    ImageRetrieverBase(const ImageOptions& options) :
      options_(options) {}
    virtual ~ImageRetrieverBase() {}
    ImageOptions options_;
  };


  ////
  //// VigraReader
  ////
  template <typename T, int N>
  class VigraReader : public ImageRetrieverBase<T, N> {
  public:
    virtual vigra::MultiArray<N, T> retrieve();
    VigraReader(const ImageOptions& options);
    virtual ~VigraReader();
  };


  ////
  //// HDF5Reader
  ////
  template <typename T, int N>
  class HDF5Reader : public ImageRetrieverBase<T, N> {
  public:
    virtual vigra::MultiArray<N, T> retrieve();
    HDF5Reader(const ImageOptions& options);
    virtual ~HDF5Reader();
  };


  ////
  //// ImageWriterBase
  ////
  template <typename T, int N>
  class ImageWriterBase {
  public:
    virtual void write(const vigra::MultiArrayView<N, T> image) = 0;
    ImageWriterBase() {}
    ImageWriterBase(const ImageOptions& options) :
      options_(options) {}
    ~ImageWriterBase() {}
    ImageOptions options_;
  };


  ////
  //// VigraWriter
  ////
  template <typename T, int N>
  class VigraWriter : public ImageWriterBase<T, N> {
  public:
    virtual void write(const vigra::MultiArrayView<N, T> image);
    VigraWriter(const ImageOptions& options);
    virtual ~VigraWriter();
  };


  ////
  //// HDF5Writer
  ////
  template <typename T, int N>
  class HDF5Writer : public ImageWriterBase<T, N> {
  public:
    virtual void write(const vigra::MultiArrayView<N, T> image);
    HDF5Writer(const ImageOptions& options);
    virtual ~HDF5Writer();
  };
    


  /* IMPLEMENTATIONS */

  ////
  //// VigraReader
  ////
  template <typename T, int N>
  vigra::MultiArray<N, T> VigraReader<T, N>::retrieve() {
    LOG(logDEBUG1) << "VigraReader<T, " << N << ">::retrieve() -- reading "
                   << ImageRetrieverBase<T, N>::options_.get("filename");
    vigra::ImageImportInfo info(ImageRetrieverBase<T, N>::options_.get("filename").c_str());
    vigra::MultiArray<N, T> image(info.shape());
    vigra::importImage(info, destImage(image));
    return image;
    // return vigra::MultiArrayView<N, T>(image);
  }


  template <typename T, int N>
  VigraReader<T, N>::VigraReader(const ImageOptions& options) :
    ImageRetrieverBase<T, N>(options) {
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
  void VigraWriter<T, N>::write(vigra::MultiArrayView<N, T> image) {
    vigra::ImageExportInfo info(ImageWriterBase<T, N>::options_.get("filename").c_str());
    vigra::exportImage(vigra::srcImageRange(image), info);
  }


  template <typename T, int N>
  VigraWriter<T, N>::VigraWriter(const ImageOptions& options) :
    ImageWriterBase<T, N>(options) {
    // nothing to be done here
  }


  template <typename T, int N>
  VigraWriter<T, N>::~VigraWriter() {
    // nothing to be done here
  }
    

}

#endif /* IMAGE_INPUT_H */
