/*===================================================================

  Program:   ASHS (Automatic Segmentation of Hippocampal Subfields)
  Module:    $Id: WeightedVotingLabelFusionImageFilter.txx 115 2016-02-25 17:34:42Z yushkevich $
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
  
#include <itkNeighborhoodIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkBinaryFunctorImageFilter.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_cholesky.h>

#include <set>
#include <map>
#include <vector>

template <class TInput1, class TInput2, class TOutput>
class NormalizeFunctor
{
public:
  bool operator !=(const NormalizeFunctor &other) { return false; }

  TOutput operator()(const TInput1 &val, const TInput2 &scale)
    {
    if(scale < 0.1)
      return val;
    else
      return val / scale;
    }
};

template <class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::UpdateInputs()
{
  char buffer[64];

  // Set the target as the primary input
  this->itk::ProcessObject::SetInput("Primary", m_Target);

  // Set the atlases and their segmentations as secondary inputs
  for(size_t i = 0; i < m_Atlases.size(); i++)
    {
    sprintf(buffer, "atlas_%04d", (int) i);
    this->itk::ProcessObject::SetInput(buffer, m_Atlases[i]);

    if(m_AtlasSegs.size())
      {
      sprintf(buffer, "atseg_%04d", (int) i);
      this->itk::ProcessObject::SetInput(buffer, m_AtlasSegs[i]);
      }
    }

  for(typename ExclusionMap::iterator it = m_Exclusions.begin(); it != m_Exclusions.end(); ++it)
    {
    sprintf(buffer, "excl_%04f", it->first);
    this->itk::ProcessObject::SetInput(buffer, it->second);
    }

  // If the mask is defined, add it as input
  if(m_MaskImage)
    this->itk::ProcessObject::SetInput("mask", m_MaskImage);
}




template <class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // Get the output requested region
  RegionType outRegion = this->GetOutput()->GetRequestedRegion();

  // Pad this region by the search window and patch size
  outRegion.PadByRadius(m_SearchRadius);
  outRegion.PadByRadius(m_PatchRadius);

  // Iterate over all the inputs to this filter
  itk::ProcessObject::DataObjectPointerArray inputs = this->GetInputs();
  for(size_t i = 0; i < inputs.size(); i++)
    {
    // Get i-th input
    InputImageType *input = dynamic_cast<InputImageType *>(inputs[i].GetPointer());
    RegionType region = outRegion;
    region.Crop(input->GetLargestPossibleRegion());
    input->SetRequestedRegion(region);
    }
}

template<class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::ComputeOffsetTable(
  const InputImageType *image, 
  const SizeType &radius, 
  int **offset, 
  size_t &nPatch,
  int **manhattan)
{
  // Use iterators to construct offset tables
  RegionType r = image->GetBufferedRegion();
  NIter itTempPatch(radius, image, r);

  // Position the iterator in the middle to avoid problems with boundary conditions
  IndexType iCenter;
  for(size_t i = 0; i < InputImageDimension; i++)
    iCenter[i] = r.GetIndex(i) + r.GetSize(i)/2;
  itTempPatch.SetLocation(iCenter);

  // Compute the offsets 
  nPatch = itTempPatch.Size();
  (*offset) = new int[nPatch];
  if(manhattan)
    (*manhattan) = new int[nPatch];
  for(size_t i = 0; i < nPatch; i++)
  {
    (*offset)[i] = itTempPatch[i] - itTempPatch.GetCenterPointer();
    if(manhattan)
      {
      typename NIter::OffsetType off = itTempPatch.GetOffset(i);
      (*manhattan)[i] = 0;
      for(int d = 0; d < InputImageDimension; d++)
        (*manhattan)[i] += abs((int) off[d]);
      }
  }
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  // Allocate the output
  this->GetOutput()->SetBufferedRegion( this->GetOutput()->GetRequestedRegion() );
  this->GetOutput()->Allocate();

  // Get the target image
  InputImageType *target = m_Target;

  // Create a neighborhood iterator for the target image
  NIter itTarget(m_PatchRadius, target, this->GetOutput()->GetRequestedRegion());

  // Get the number of atlases
  int n = m_Atlases.size();

  // Do we have segmentation inputs
  bool have_segs = m_AtlasSegs.size() == n;

  // Construct offset tables for all the images (these can be different because they
  // depend on the buffered region)
  m_OffPatchAtlas = new int *[n];
  m_OffPatchSeg = new int *[n];
  m_OffSearchAtlas = new int *[n];
  m_OffSearchSeg = new int *[n];

  // Compute the offset table for the target image
  ComputeOffsetTable(target, m_PatchRadius, &m_OffPatchTarget, m_NPatch);

  // Find all unique labels in the requested region
  m_LabelSet.clear();

  for(int i = 0; i < n; i++)
    {
    // Compute the offset table for that atlas
    ComputeOffsetTable(m_Atlases[i], m_PatchRadius, m_OffPatchAtlas+i, m_NPatch);
    ComputeOffsetTable(m_Atlases[i], m_SearchRadius, m_OffSearchAtlas+i, m_NSearch, &m_Manhattan);

    // If there are segmentation inputs, process them
    if(have_segs)
      {
      // Find all the labels. This is fast enough to not require threading
      const InputImageType *seg = m_AtlasSegs[i];
      itk::ImageRegionConstIteratorWithIndex<InputImageType> it(seg, seg->GetRequestedRegion());
      bool have_last_label = false;
      InputImagePixelType last_label;
      for(; !it.IsAtEnd(); ++it)
        {
        InputImagePixelType label = it.Get();
        if(!have_last_label || last_label != label)
          {
          m_LabelSet.insert(label);
          last_label = label;
          have_last_label = true;
          }
        }

      ComputeOffsetTable(m_AtlasSegs[i], m_PatchRadius, m_OffPatchSeg+i, m_NPatch);
      ComputeOffsetTable(m_AtlasSegs[i], m_SearchRadius, m_OffSearchSeg+i, m_NSearch, &m_Manhattan);
      }
    }

  // Initialize the posterior maps
  m_PosteriorMap.clear();

  // Allocate posterior images for the different labels
  if(have_segs)
    {
    for(typename std::set<InputImagePixelType>::iterator sit = m_LabelSet.begin();
      sit != m_LabelSet.end(); ++sit)
      {
      m_PosteriorMap[*sit] = PosteriorImage::New();
      m_PosteriorMap[*sit]->SetLargestPossibleRegion(target->GetLargestPossibleRegion());
      m_PosteriorMap[*sit]->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
      m_PosteriorMap[*sit]->SetBufferedRegion(this->GetOutput()->GetRequestedRegion());
      m_PosteriorMap[*sit]->Allocate();
      m_PosteriorMap[*sit]->FillBuffer(0.0f);
      }
    }

  // Generate the optional weight maps
<<<<<<< HEAD
  if(m_GenerateWeightMaps)
    {
    m_WeightMapArray.resize(n);
    m_WeightMapArrayBuffer = new float *[n];
=======
  m_WeightMapArray.resize(n);
  if(m_GenerateWeightMaps)
    {
>>>>>>> 515ff7c2f50928adabc4e64bded9a7e76fc750b1
    for(int i = 0; i < n; i++)
      {
      m_WeightMapArray[i] = WeightMapImage::New();
      m_WeightMapArray[i]->SetLargestPossibleRegion(target->GetLargestPossibleRegion());
      m_WeightMapArray[i]->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
      m_WeightMapArray[i]->SetBufferedRegion(this->GetOutput()->GetRequestedRegion());
      m_WeightMapArray[i]->Allocate();
      m_WeightMapArray[i]->FillBuffer(1.0f / n);
<<<<<<< HEAD
      m_WeightMapArrayBuffer[i] = m_WeightMapArray[i]->GetBufferPointer();
      }
    }
=======
      }
    }

  int iter = 0;

  // We need an array of absolute patch differences between target image and atlases
  // (apd - atlas patch difference)
  InputImagePixelType **apd = new InputImagePixelType*[n];
  for(int i = 0; i < n; i++)
    apd[i] = new InputImagePixelType[nPatch];

  // Also an array of pointers to the segmentations of different atlases
  const InputImagePixelType **patchSeg = new const InputImagePixelType*[n]; 
>>>>>>> 515ff7c2f50928adabc4e64bded9a7e76fc750b1

  // Create a counter map -- needed if we have weights or posteriors - so always
  m_CounterMap = PosteriorImage::New();
  m_CounterMap->SetLargestPossibleRegion(target->GetLargestPossibleRegion());
  m_CounterMap->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
  m_CounterMap->SetBufferedRegion(this->GetOutput()->GetRequestedRegion());
  m_CounterMap->Allocate();
  m_CounterMap->FillBuffer(0.0f);

  // In earlier code, we iterated over all the voxels in the target image. But this is often
  // not necessary because much of the image is just background. Now we allow the user to 
  // provide a flag to automatically mask the iterated region by the dilated union of all
  // the segmentations. This should have no effect on the output segmentation, but will affect
  // the posterior maps
  m_Mask = NULL;
  if(m_MaskImage.IsNull() && have_segs)
    {
    std::cout << "Computing mask based on input segmentations" << std::endl;

    // Create a mask from all the segmentations
    m_Mask = InputImageType::New();
    m_Mask->CopyInformation(this->GetOutput());
    m_Mask->SetRegions(this->GetOutput()->GetBufferedRegion());
    m_Mask->Allocate();
    m_Mask->FillBuffer(0);

    unsigned int nmasked = 0, nskipped = 0;

    // The mask image will be constructed as follows. For each voxel, we search for all the 
    // different labels in each of the input segmentation images, using the search radius 
    // specified. As soon as we find more than one unique voxel, we mask the voxel as being
    // useful for labeling. Otherwise, we just assign the common label to the output and 
    // ignore this voxel in the label fusion.
    typedef itk::ImageRegionIteratorWithIndex<InputImageType> MaskIter;
    for(MaskIter it(m_Mask, m_Mask->GetBufferedRegion()); !it.IsAtEnd(); ++it)
      {
      // Start by recording the label of the first atlas at this location
      InputImagePixelType uniq = m_AtlasSegs[0]->GetPixel(it.GetIndex());
      bool is_uniqie = true;

      // Compare this value to all the other segmentation values in the search radius
      for(int i = 0; i < n; i++)
        {
        // Get the i-th segmentation and its offset table for search
        const InputImageType *seg = m_AtlasSegs[i];
        int *offSearch = m_OffSearchSeg[i];

        // Find the current voxel in the atlas seg
        const InputImagePixelType *pSeg = seg->GetBufferPointer() + seg->ComputeOffset(it.GetIndex());
        for(unsigned int k = 0; k < m_NSearch; k++)
          {
          InputImagePixelType intensity = pSeg[offSearch[k]];
          if(intensity != uniq)
            {
            is_uniqie = false;
            break;
            }
          }

        // No need to search if intensity is not unique
        if(!is_uniqie)
          break;
        }

      // If the pixel has unique possible label, assign that label to the output
      if(is_uniqie)
        {
        this->GetOutput()->SetPixel(it.GetIndex(), uniq);
        it.Set(0);
        nskipped++;
        }
      else
        {
        it.Set(1);
        nmasked++;
        }
      }

      std::cout << "  Skipping " << nskipped << " out of " << nskipped+nmasked << " voxels." << std::endl;
    }
  else if(m_MaskImage.IsNotNull())
    {
    m_Mask = m_MaskImage;
    std::cout << "  Using mask provided by user" << std::endl;
    }
  else
    {
    m_Mask = NULL;
    std::cout << "  No mask supplied, using whole image" << std::endl;
    }

  // Initialize thread data
  m_ThreadData.resize(this->GetNumberOfThreads());
}

template <class T>
T* allocate_aligned(int elements)
{
  void* pointer;
  posix_memalign(&pointer, 16, elements * sizeof(T));
  return static_cast<T *>(pointer)  ;
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
  // Get the target image
  InputImageType *target = m_Target;

  // Create a neighborhood iterator for the target image
  NIter itTarget(m_PatchRadius, target, outputRegionForThread);

  // Get the number of atlases
  int n = m_Atlases.size();
  bool have_segs = m_AtlasSegs.size() == n;

  // Allocate Mx
  typedef vnl_matrix<double> MatrixType;
  MatrixType Mx(n, n);

  // Define a vector of all ones
  vnl_vector<double> ones(n, 1.0);

  // Solve for the weights
  vnl_vector<double> W(n, 0.0);

  // Collect search statistics
  m_ThreadData[threadId].m_SearchHisto.resize(100, 0);

  // Keep track of iterations
  int iter = 0;

  // For faster code, we will 16-byte align the arrays and also allocate them as factor of 4
  int n_PatchRnd = (m_NPatch % 4 == 0) ? m_NPatch : ((m_NPatch >> 2) + 1) << 2;

  // We need an array of absolute patch differences between target image and atlases
  // (apd - atlas patch difference)
  InputImagePixelType **apd = new InputImagePixelType*[n];
  for(int i = 0; i < n; i++)
    {
    apd[i] = allocate_aligned<float>(n_PatchRnd);
    for(int j = 0; j < n_PatchRnd; j++)
      apd[i][j] = 0.0f;
    }

  // Also an array of pointers to the segmentations of different atlases
  const InputImagePixelType **patchSeg = new const InputImagePixelType*[n]; 

  // Create an array for storing the normalized target patch to save more time
  InputImagePixelType *xNormTargetPatch = new InputImagePixelType[m_NPatch];

  // Iterate over voxels in the output region
  typedef itk::ImageRegionIteratorWithIndex<TOutputImage> OutIter;
  for(OutIter it(this->GetOutput(), outputRegionForThread); !it.IsAtEnd(); ++it)
    {
    // If this point is outside of the mask, skip it for posterior computation
    if(m_Mask && m_Mask->GetPixel(it.GetIndex()) == 0)
      continue;

    // Point the target iterator to the output location
    itTarget.SetLocation(it.GetIndex());
    InputImagePixelType *pTargetCurrent = target->GetBufferPointer() + target->ComputeOffset(it.GetIndex());

    // Compute stats for the target patch
    InputImagePixelType mu, sigma;
    PatchStats(pTargetCurrent, m_NPatch, m_OffPatchTarget, mu, sigma);
    for(unsigned int i = 0; i < m_NPatch; i++)
      xNormTargetPatch[i] = (*(pTargetCurrent + m_OffPatchTarget[i]) - mu) / sigma;

    // In each atlas, search for a patch that matches our patch
    for(int i = 0; i < n; i++)
      {
      const InputImageType *atlas = m_Atlases[i];
      int *offPatch = m_OffPatchAtlas[i], *offSearch = m_OffSearchAtlas[i];

      // Search over neighborhood
      const InputImagePixelType *pAtlasCurrent = atlas->GetBufferPointer() + atlas->ComputeOffset(it.GetIndex());
      double bestMatch = 1e100;
      const InputImagePixelType *bestMatchPtr = NULL;
      InputImagePixelType bestMatchSum = 0, bestMatchSSQ = 0;
      int bestK = 0;
      for(unsigned int k = 0; k < m_NSearch; k++)
        {
        // Pointer to the voxel at the center of the search
        const InputImagePixelType *pSearchCenter = pAtlasCurrent + offSearch[k];
        InputImagePixelType matchSum = 0, matchSSQ = 0;
        double match = this->PatchSimilarity(pSearchCenter, xNormTargetPatch, m_NPatch, offPatch,
                                             matchSum, matchSSQ);
        if(k == 0 || match < bestMatch)
          {
          bestMatch = match;
          bestMatchPtr = pSearchCenter;
          bestMatchSum = matchSum;
          bestMatchSSQ = matchSSQ;
          bestK = k;
          }
        }

      // Update the manhattan distance histogram
      m_ThreadData[threadId].m_SearchHisto[m_Manhattan[bestK]]++;

      // Once the patch has been found, compute the absolute difference with target image
      InputImagePixelType bestMatchMean = bestMatchSum / m_NPatch;
      InputImagePixelType bestMatchVar = 
        (bestMatchSSQ - m_NPatch * bestMatchMean * bestMatchMean) / (m_NPatch - 1);
      if(bestMatchVar < 1.0e-12)
        bestMatchVar = 1.0e-12;
      InputImagePixelType bestMatchSD = sqrt(bestMatchVar);

      for(unsigned int m = 0; m < m_NPatch; m++)
        {
        InputImagePixelType x = *(bestMatchPtr + offPatch[m]);
        apd[i][m] = fabs(xNormTargetPatch[m] - (x - bestMatchMean) / bestMatchSD);
        }

      // Store the best found neighborhood
      if(have_segs)
        {
        const InputImageType *seg = m_AtlasSegs[i];
        patchSeg[i] = (bestMatchPtr - atlas->GetBufferPointer()) + seg->GetBufferPointer();
        }
      }

    // Now we can compute Mx
    for(int i = 0; i < n; i++) 
      {
      float *apdi = apd[i];
      for(int k = 0; k <= i; k++) 
        {
        float *apdk = apd[k];

        // Multiply through the apd arrays - this is slow C code
 #ifdef _NO_SSE_

        InputImagePixelType mxval = 0.0;
        for(unsigned int m = 0; m < n_PatchRnd; m+=4)
          {
          mxval += apdi[m] * apdk[m];
          mxval += apdi[m+1] * apdk[m+1];
          mxval += apdi[m+2] * apdk[m+2];
          mxval += apdi[m+3] * apdk[m+3];
          }

#else
        // Fast multiplication
        __m128 acc, x, y; 

        // Zero out the accumulator
        acc = _mm_set_ps(0.0f, 0.0f, 0.0f, 0.0f);

        for(unsigned int m = 0; m < n_PatchRnd; m+=4)
          {
          x = _mm_load_ps(apdi + m);
          y = _mm_load_ps(apdk + m);
          acc = _mm_add_ps(acc, _mm_mul_ps(x, y));
          }

        InputImagePixelType mxval;
        __m128 shuf   = _mm_shuffle_ps(acc, acc, _MM_SHUFFLE(2, 3, 0, 1));  // [ C D | B A ]
        __m128 sums   = _mm_add_ps(acc, shuf);      // sums = [ D+C C+D | B+A A+B ]
        shuf          = _mm_movehl_ps(shuf, sums);      //  [   C   D | D+C C+D ]  // let the compiler avoid a mov by reusing shuf
        sums          = _mm_add_ss(sums, shuf);
        mxval         = _mm_cvtss_f32(sums);

#endif

        mxval /= (m_NPatch - 1);
        
        if(m_Beta == 2)
          mxval *= mxval;
        else
          mxval = pow(mxval, m_Beta);

        Mx(i,k) = Mx(k,i) = mxval;
        }

      // Add alpha
      Mx(i,i) += m_Alpha;
      }

    // Now we can compute the weights by solving for the inverse of Mx
    vnl_cholesky cholesky(Mx, vnl_cholesky::estimate_condition);
    if(cholesky.rcond() > vnl_math::sqrteps)
      {
      // Well-conditioned matrix
      cholesky.solve(ones, &W);
      }
    else
      {
      // Matrix badly conditioned
      W = vnl_svd<double>(Mx).solve(ones);
      }

    // Normalize the weights
    W *= 1.0 / dot_product(W, ones);

    // Compute the sum of the weights (shouldn't this always be one?)
    float Wsum = 0.0;
    for(int i = 0; i < n; i++)
      Wsum += W[i];

    /*
    # Debugging placeholder - for verifying weights
    if(it.GetIndex()[0] == 193 && it.GetIndex()[1] == 78 && it.GetIndex()[2] == 17)
      {
      std::cout << "Mx:" << std::endl;
      std::cout << Mx << std::endl;
      std::cout << "W:" << std::endl;
      std::cout << W << std::endl;
      }
    */
    
    // Reduce the number of std::map lookups for speed
    bool have_last = false;
    InputImagePixelType last_label;
    typename PosteriorImage::PixelType *last_posterior_buffer = NULL;

    // Counter map buffer - direct access
    typename PosteriorImage::PixelType *countermap_buffer = m_CounterMap->GetBufferPointer();

    // Output the weight maps if desired
    if(m_GenerateWeightMaps)
      {
      IndexType idx = it.GetIndex();
      for(int q = 0; q < n; q++)
        m_WeightMapArray[q]->SetPixel(idx, W[q]);
      }

    // Perform voting using Hongzhi's averaging scheme. Iterate over all segmentation patches
    for(unsigned int ni = 0; ni < m_NPatch; ni++)
      {
      // The index of the patch voxel. This index may fall outside of the thread's output
      // region. In this case, we must use a mutex to ensure that two threads are not writing
      // to the same location at the same time. Hopefully this will not create a bottleneck!
      IndexType idx = itTarget.GetIndex(ni);

      // Outside of the overall region - ignore
      if(!this->GetOutput()->GetRequestedRegion().IsInside(idx))
        continue;

      // Outside of the threaded region - need to have exclusivity. However, the chances 
      // of two threads trying to write to the same location at once are next to nil, so
      // for now we will just sweep it under the rug!
        
      // To save some time, we can convert this index into an offset since all the images
      // below use the same regions
      typename InputImageType::OffsetValueType idx_offset = this->GetOutput()->ComputeOffset(idx);

      for(int i = 0; i < n; i++)
        {
        // Update the posteriors - if they exist!
        if(have_segs)
          {
          // The segmentation at the corresponding patch location in atlas i
          InputImagePixelType label = *(patchSeg[i] + m_OffPatchSeg[i][ni]);

          // Update the posterior - reduce number of map lookups
          if(!have_last || label != last_label)
            {
            last_label = label;
            last_posterior_buffer = m_PosteriorMap[label]->GetBufferPointer();
            have_last = true;
            }

          // Add that weight the posterior map for voxel at idx
          last_posterior_buffer[idx_offset] += W[i];
          }

        // Add the weight to the weight map too
        if(m_GenerateWeightMaps)
          {
          m_WeightMapArrayBuffer[i][idx_offset] += W[i];
          }
        }

      // Add the weight to the counter
      countermap_buffer[idx_offset] += Wsum;
      }

      if(++iter % 1000 == 0)
        {
        static double t = clock();
        std::cout << "." << std::flush;
        }
    }
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::AfterThreadedGenerateData()
{
  /* Who cares?
  std::cout << std::endl << "Search Manhattan Distance Histogram " << std::endl;
  for(size_t i = 0; i < searchHisto.size() && searchHisto[i] > 0; i++)
    std::cout << "    " << i << "\t" << searchHisto[i] << std::endl;
  */

  std::cout << std::endl << "VOTING " << std::endl;

  // Filter type for normalizing by the counter
  typedef NormalizeFunctor<float, float, float> FloatNormalizeFunctor;
  typedef itk::BinaryFunctorImageFilter<PosteriorImage, PosteriorImage, PosteriorImage, FloatNormalizeFunctor> NormFilter;

  // Perform voting at each voxel
  if(m_AtlasSegs.size() == m_Atlases.size())
    {
    typedef itk::ImageRegionIteratorWithIndex<TOutputImage> OutIter;
    for(OutIter it(this->GetOutput(), this->GetOutput()->GetBufferedRegion()); !it.IsAtEnd(); ++it)
      {
      // If this point is outside of the mask, skip it for posterior computation
      if(m_Mask && m_Mask->GetPixel(it.GetIndex()) == 0)
        continue;

      double wmax = 0;
      InputImagePixelType winner = 0;

      for(typename std::set<InputImagePixelType>::iterator sit = m_LabelSet.begin();
        sit != m_LabelSet.end(); ++sit)
        {
        double posterior = m_PosteriorMap[*sit]->GetPixel(it.GetIndex());

        // check if the label is excluded
        typename ExclusionMap::iterator xit = m_Exclusions.find(*sit);
        bool excluded = (xit != m_Exclusions.end() && xit->second->GetPixel(it.GetIndex()) != 0);

        // Vote!
        if (wmax < posterior && !excluded)
          {
          wmax = posterior;
          winner = *sit;
          }
        }

      it.Set(winner);
      }

    // Clear posterior maps
    if(!m_RetainPosteriorMaps)
      {
      m_PosteriorMap.clear();
      }
    else
      {
      for(typename PosteriorMap::const_iterator itp = m_PosteriorMap.begin();
        itp != m_PosteriorMap.end(); ++itp)
        {
        typename NormFilter::Pointer norm = NormFilter::New();
        norm->SetInput1(itp->second);
        norm->SetInput2(m_CounterMap);
        norm->GraftOutput(itp->second);
        norm->Update();
        }
      }
    }

  // Normalize weight maps
  if(m_GenerateWeightMaps)
    {
    for(int i = 0; i < m_WeightMapArray.size(); i++)
      {
      typename NormFilter::Pointer norm = NormFilter::New();
      norm->SetInput1(m_WeightMapArray[i]);
      norm->SetInput2(m_CounterMap);
      norm->GraftOutput(m_WeightMapArray[i]);
      norm->Update();
      }
    }
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::PatchStats(const InputImagePixelType *p, size_t n, const int *offsets,
             InputImagePixelType &mean, InputImagePixelType &std)
{
  InputImagePixelType sum = 0, ssq = 0;
  for(unsigned int i = 0; i < n; i++)
    {
    InputImagePixelType v = *(p + offsets[i]);
    sum += v;
    ssq += v * v;
    }

  mean = sum / n;
  std = sqrt((ssq - n * mean * mean) / (n - 1));

  // Check for very small values or NaN
  if(std < 1e-6 || std != std) 
    std = 1e-6;
}

/**
 * This function computes similarity between a normalized patch (normtrg) and a patch
 * that has not been normalized (psearch). It can be shown that the sum of squared 
 * differences between a normalized patch u and a unnormalized patch v is equal to
 *
 * 2 [ (n-1) - (\Sum u_i v_i ) / \sigma_v ]
 *
 * Since we are only interested in finding the patch with the smallest SSD, we can simplify
 * this expression further to minimizing - [ (\Sum u_i v_i ) / \sigma_v ] ^ 2. To further
 * simpify computation, we return
 *
 *        - (\Sum u_i v_i)^2 / z,   where z = sigma_v^2 * (n-1)
 */
template <class TInputImage, class TOutputImage>
double
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::PatchSimilarity(
  const InputImagePixelType *psearch, 
  const InputImagePixelType *normtrg, 
  size_t n, 
  const int *offsets,
  InputImagePixelType &sum_psearch,
  InputImagePixelType &ssq_psearch)
{
  // Here the patch normtrg should already be normalized.
  // We simultaneously compute the patch stats and solve the problem
  InputImagePixelType sum_uv = 0;
  for(unsigned int i = 0; i < n; i++)
    {
    InputImagePixelType u = *(psearch + offsets[i]);
    InputImagePixelType v = normtrg[i];
    sum_psearch += u;
    ssq_psearch += u * u;
    sum_uv += u * v;
    }

  InputImagePixelType var_u_unnorm = ssq_psearch - sum_psearch * sum_psearch / n;
  if(var_u_unnorm < 1.0e-6)
    var_u_unnorm = 1.0e-6;

  if(sum_uv > 0)
    return - (sum_uv * sum_uv) / var_u_unnorm;
  else
    return (sum_uv * sum_uv) / var_u_unnorm;

  // InputImagePixelType sd_u = std::max(1e-6, sqrt((ssq_u - sum_u * sum_u / n) / (n - 1)));
  // return 2 * ((n - 1) - sum_uv / sd_u);
}

template <class TInputImage, class TOutputImage>
double
WeightedVotingLabelFusionImageFilter<TInputImage, TOutputImage>
::JointErrorEstimate(const InputImagePixelType *t, const InputImagePixelType *a1, const InputImagePixelType *a2, size_t n, int *offsets)
{
  InputImagePixelType mu_t, sigma_t, mu1, sigma1, mu2, sigma2;
  PatchStats(t, n, offsets, mu_t, sigma_t);
  PatchStats(a1, n, offsets, mu1, sigma1);
  PatchStats(a2, n, offsets, mu2, sigma2);

  // What should we return when patches have zero variance?
  if(sigma1 == 0 || sigma2 == 0 || sigma_t == 0)
    return 0;

  double Mxval = 0.0;
  for(unsigned int i = 0; i < n; i++)
    {
    int off = offsets[i];
    InputImagePixelType ft = (*(t + off) - mu_t) / sigma_t;
    InputImagePixelType f1 = (*(a1 + off) - mu1) / sigma1;
    InputImagePixelType f2 = (*(a2 + off) - mu2) / sigma2;

    Mxval += fabs(ft-f1) * fabs(ft-f2);
    }

  return pow(Mxval, -m_Beta);
}

