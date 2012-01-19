#ifndef __WeightedVotingLabelFusionImageFilter_h_
#define __WeightedVotingLabelFusionImageFilter_h_

#include "itkImageToImageFilter.h"
#include "itkConstNeighborhoodIterator.h"

template <class TInputImage, class TOutputImage>
class WeightedVotingLabelFusionImageFilter : public itk::ImageToImageFilter <TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef WeightedVotingLabelFusionImageFilter Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage>  Superclass;
  typedef itk::SmartPointer<Self>         Pointer;
  typedef itk::SmartPointer<const Self>   ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToImageFilter,ImageSource);

  itkNewMacro(Self);

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename Superclass::OutputImagePixelType  OutputImagePixelType;

  /** Some convenient typedefs. */
  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::PixelType      InputImagePixelType;

  typedef typename InputImageType::RegionType     RegionType;
  typedef typename InputImageType::SizeType       SizeType;
  typedef typename InputImageType::IndexType      IndexType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
                                                                               
  /** Set target image */
  void SetTargetImage(InputImageType *image)
    { this->SetInput(0, image); }

  /** Set atlas images */
  void SetNumberOfAtlases(int nAtlases)
    {
    this->SetNumberOfInputs(nAtlases * 2 + 1);
    }

  void SetAtlas(int i, InputImageType *grey, InputImageType *seg)
    {
    SetInput(1 + 2 * i, grey);
    SetInput(2 + 2 * i, seg);
    }

  /** Set the parameters */
  itkSetMacro(SearchRadius, SizeType);
  itkGetMacro(SearchRadius, SizeType);

  itkSetMacro(PatchRadius, SizeType);
  itkGetMacro(PatchRadius, SizeType);

  itkSetMacro(Alpha, double);
  itkGetMacro(Alpha, double);

  itkSetMacro(Beta, double);
  itkGetMacro(Beta, double);

  /** Set the requested region */
  void GenerateInputRequestedRegion();


  void GenerateData();
 
protected:

  WeightedVotingLabelFusionImageFilter() { m_Alpha=0.01; m_Beta=2; }
  ~WeightedVotingLabelFusionImageFilter() {}

private:

  typedef itk::Neighborhood<InputImagePixelType, InputImageDimension> HoodType;
  typedef itk::ConstNeighborhoodIterator<InputImageType> NIter;

  double PatchSimilarity(
    const InputImagePixelType *psearch, const InputImagePixelType *pnormtrg, 
    size_t n, int *offsets, InputImagePixelType &psearchSum, InputImagePixelType &psearchSSQ);

  void ComputeOffsetTable(
    const InputImageType *image, const SizeType &radius, 
    int **offset, size_t &nPatch, int **manhattan = NULL);

  void PatchStats(const InputImagePixelType *p, size_t n, int *offsets, InputImagePixelType &mean, InputImagePixelType &sd);

  double JointErrorEstimate(const InputImagePixelType *t, const InputImagePixelType *a1, const InputImagePixelType *a2, size_t n, int *offsets);

  SizeType m_SearchRadius, m_PatchRadius;

  double m_Alpha, m_Beta;

};


#endif
