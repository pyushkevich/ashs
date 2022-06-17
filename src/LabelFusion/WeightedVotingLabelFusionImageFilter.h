/*===================================================================

  Program:   ASHS (Automatic Segmentation of Hippocampal Subfields)
  Module:    $Id: WeightedVotingLabelFusionImageFilter.h 114 2016-02-25 17:11:29Z yushkevich $
  Language:  C++ program
  Copyright (c) 2012 Paul A. Yushkevich, University of Pennsylvania
  
  This file is part of ASHS

  ASHS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details. 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  =================================================================== 
  
  CITATION:

    This program implements the method described in the paper

    H. Wang, J. W. Suh, J. Pluta, M. Altinay, and P. Yushkevich, 
       (2011) "Computing Optimal Weights for Label Fusion based 
       Multi-Atlas Segmentation," in Proc. Information Processing 
       in Medical Imaging (IPMI), 2011


  =================================================================== */
  
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
    { m_Target = image; UpdateInputs(); }

  /** Add an atlas */
  void AddAtlas(InputImageType *grey, InputImageType *seg)
    {
    m_Atlases.push_back(grey);
    m_AtlasSegs.push_back(seg);
    UpdateInputs();
    }

  /** Add an atlas without labels */
  void AddAtlas(InputImageType *grey)
    {
    m_Atlases.push_back(grey);
    UpdateInputs();
    }

  void AddExclusionMap(InputImagePixelType label, InputImageType *excl)
    {
    m_Exclusions[label] = excl;
    UpdateInputs();
    }

  /** Set the mask image. A mask image explicitly specifies where voting is performed */
  void SetMaskImage(InputImageType *mask)
    {
    m_MaskImage = mask;
    UpdateInputs();
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

  /** 
   * Whether the posterior maps should be retained. This can have a negative effect
   * on memory use, so it should only be done if one wishes to save the posterior
   * maps. The posterior maps given the probability of each voxel in the target image
   * belonging to each label.
   */
  itkSetMacro(RetainPosteriorMaps, bool)
  itkGetMacro(RetainPosteriorMaps, bool)

  /**
   * Whether per-atlas weight maps should be generated. This really only makes sense
   * when the search radius is zero, because the weight corresponds to the atlas 
   * patch that was found to be the best match to the target patch, and this information
   * is not currently exported in any way
   */
  itkSetMacro(GenerateWeightMaps, bool)
  itkGetMacro(GenerateWeightMaps, bool)

  typedef itk::Image<float, InputImageDimension> PosteriorImage;
  typedef typename PosteriorImage::Pointer PosteriorImagePtr;
  typedef typename std::map<InputImagePixelType, PosteriorImagePtr> PosteriorMap;

  typedef itk::Image<float, InputImageDimension> WeightMapImage;
  typedef typename WeightMapImage::Pointer WeightMapImagePtr;
  typedef typename std::vector<WeightMapImagePtr> WeightMapArray;
                                                                    
  /**
   * Get the posterior maps (if they have been retained)
   */
  const PosteriorMap &GetPosteriorMaps()
    { return m_PosteriorMap; }

  /**
   * Get the weight image for each atlas
   */
  WeightMapImage* GetWeightMap(int iAtlas) const
    { return m_WeightMapArray[iAtlas]; }



  void BeforeThreadedGenerateData();
  void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId);
  void AfterThreadedGenerateData();

 
protected:

  WeightedVotingLabelFusionImageFilter() 
    { 
    m_Alpha=0.01; 
    m_Beta=2; 
    m_RetainPosteriorMaps = false;
    m_GenerateWeightMaps = false;
    }
  ~WeightedVotingLabelFusionImageFilter() {}

private:

  typedef itk::Neighborhood<InputImagePixelType, InputImageDimension> HoodType;
  typedef itk::ConstNeighborhoodIterator<InputImageType> NIter;

  double PatchSimilarity(
    const InputImagePixelType *psearch, const InputImagePixelType *pnormtrg, 
    size_t n, const int *offsets, InputImagePixelType &psearchSum, InputImagePixelType &psearchSSQ);

  void ComputeOffsetTable(
    const InputImageType *image, const SizeType &radius, 
    int **offset, size_t &nPatch, int **manhattan = NULL);

  void UpdateInputs();

  void PatchStats(const InputImagePixelType *p, size_t n, const int *offsets, 
                  InputImagePixelType &mean, InputImagePixelType &sd);

  double JointErrorEstimate(const InputImagePixelType *t, const InputImagePixelType *a1, const InputImagePixelType *a2, size_t n, int *offsets);

  SizeType m_SearchRadius, m_PatchRadius;

  double m_Alpha, m_Beta;

  typedef std::vector<InputImagePointer> InputImageList;
  typedef std::map<InputImagePixelType, InputImagePointer> ExclusionMap;

  // Posterior maps
  PosteriorMap m_PosteriorMap;

  // Whether they are retained
  bool m_RetainPosteriorMaps;

  // Whether weight maps are computed
  bool m_GenerateWeightMaps;

  // Optional weight map array
  WeightMapArray m_WeightMapArray;

  // Array of weight map data pointers - for faster access
  float **m_WeightMapArrayBuffer;

  // Organized lists of inputs
  InputImagePointer m_Target, m_MaskImage;
  InputImageList m_AtlasSegs, m_Atlases;
  ExclusionMap m_Exclusions;

  // Stuff used internally
  int *m_OffPatchTarget, **m_OffPatchAtlas, **m_OffPatchSeg, **m_OffSearchAtlas, **m_OffSearchSeg;
  int *m_Manhattan;

  // Mask - may be maskimage or may be internal
  InputImagePointer m_Mask;

  // Neighborhood sizes
  size_t m_NPatch, m_NSearch;

  // Set of labels
  std::set<InputImagePixelType> m_LabelSet;

  // Counter map
  PosteriorImagePtr m_CounterMap;

  // Thread-specific data
  struct ThreadData
    {
    std::vector<int> m_SearchHisto;
    };

  std::vector<ThreadData> m_ThreadData;

};


#endif
