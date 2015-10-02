/*===================================================================

  Program:   ASHS (Automatic Segmentation of Hippocampal Subfields)
  Module:    $Id$
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
  

#include "WeightedVotingLabelFusionImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <iostream>

#include "itkMirrorPadImageFilter.h"
#include "itkCropImageFilter.h"
#include "WeightedVotingLabelFusionImageFilter.txx"

using namespace std;

int usage()
{
  cout << "label_fusion: " << endl;
  cout << "usage: " << endl;
  cout << "  label_fusion [dim] [options] target_image output_image" << endl;
  cout << endl; 
  cout << "required options:" << endl;
  cout << "  dim                             Image dimension (2 or 3)" << endl;
  cout << "  -g atlas1.nii ... atlasN.nii    Atlas intensity images" << endl;
  cout << "  -l label1.nii ... labelN.nii    Atlas segmentation images" << endl;
  cout << "  -m <method> [parameters]        Select voting method." << endl;
  cout << "                                  Options: Gauss (Gaussian Weighting), " << endl;
  cout << "                                           Joint (Joint Label Fusion) " << endl;
  cout << "                                  May be followed by optional parameters" << endl;
  cout << "                                  in brackets, e.g., -m Gauss[0.5] or -m Joint[0.01,2]" << endl;
  cout << "                                  See below for parameters" << endl;
  cout << "other options: " << endl;
  cout << "  -rp radius                      Patch radius for similarity measures " << endl;
  cout << "                                  scalar or vector (AxBxC) " << endl;
  cout << "                                  Default: 3x3x3" << endl;
  cout << "  -rs radius                      Search radius for correcting registration." << endl;
  cout << "                                  Default: 3x3x3" << endl;
  cout << "  -pd radius                      Additional boundary padding for the images (use only if " << endl;
  cout << "                                  the segmentation extends all the way to image boundaries." << endl;
  cout << "  -x label image.nii              Specify an exclusion region for the given label. " << endl;
  cout << "                                  If a voxel has non-zero value in an exclusion image," << endl;
  cout << "                                  the corresponding label is not allowed at that voxel." << endl;
  cout << "  -p filenamePattern              Save the posterior voting maps (probability that each " << endl;
  cout << "                                  voxel belongs to each label) as images. The number of " << endl;
  cout << "                                  images saved equals the number of labels. The filename" << endl;
  cout << "                                  pattern must be in C printf format, e.g. posterior%04d.nii.gz" << endl; 
  cout << "  -w filenamePattern              Save weight maps corresponding to the atlases." << endl;
  cout << "                                  The pattern should be like weight%04d.nii.gz" << endl;
  cout << "                                  This really only makes sense for -rs 0x0x0" << endl;
  cout << "Parameters for -m Gauss option:" << endl;
  cout << "  sigma                           Standard deviation of Gaussian" << endl;
  cout << "                                  Default: X.XX" << endl;
  cout << "Parameters for -m Joint option:" << endl;
  cout << "  alpha                           Regularization term added to matrix Mx for inverse" << endl;
  cout << "                                  Default: X.XX" << endl;
  cout << "  beta                            Exponent for mapping intensity difference to joint error" << endl;
  cout << "                                  Default: X.XX" << endl;

  return -1;
}

enum LFMethod
{
  JOINT, GAUSSIAN, INVERSE 
};

template<unsigned int VDim> 
struct LFParam
{
  vector<string> fnAtlas;
  vector<string> fnLabel;
  string fnTarget;
  string fnOutput;
  string fnPosterior;
  string fnWeight;
  LFMethod method;

  map<int, string> fnExclusion;

  double alpha, beta, sigma;
  itk::Size<VDim> r_patch, r_search;
  
  bool padding;
  itk::Size<VDim> paddingSize;

  LFParam()
    {
    alpha = 0.01;
    beta = 2;
    sigma = 0.5;
    r_patch.Fill(3);
    r_search.Fill(3);
    method = JOINT;
    padding = false;
    }

  void Print(std::ostream &oss)
    {
    oss << "Target image: " << fnTarget << endl;
    oss << "Output image: " << fnOutput << endl;
    oss << "Atlas images: " << endl;
    for(size_t i = 0; i < fnAtlas.size(); i++)
      {
      oss << "    " << i << "\t" << fnAtlas[i] << " | " << fnLabel[i] << endl;
      }
    if(method == GAUSSIAN)
      {
      oss << "Method:     Gaussian" << endl;
      oss << "    Sigma:  " << sigma << endl;
      }
    else if(method == JOINT)
      {
      oss << "Method:     Joint" << endl;
      oss << "    Alpha:  " << alpha << endl;
      oss << "    Beta:  " << beta << endl;
      }
    else if(method == INVERSE)
      {
      oss << "Method:     Inverse" << endl;
      oss << "    Beta:  " << beta << endl;
      }
    oss << "Search Radius: " << r_search << endl;
    oss << "Patch Radius: " << r_patch << endl;
    if(fnPosterior.size())
      oss << "Posterior Filename Pattern: " << fnPosterior << endl;
    if(fnWeight.size())
      oss << "Weight Map Filename Pattern: " << fnWeight << endl;

    oss << "Padding Enabled: " << padding << endl;
    if(padding)
      oss << "Padding Radius: " << paddingSize << endl;
    }
};


template <unsigned int VDim>
void ExpandRegion(itk::ImageRegion<VDim> &r, bool &isinit, const itk::Index<VDim> &idx)
{
  if(!isinit)
    {
    for(size_t d = 0; d < VDim; d++)
      {
      r.SetIndex(d, idx[d]);
      r.SetSize(d,1);
      }
    isinit = true;
    }
  else
    {
    for(size_t d = 0; d < VDim; d++)
      {
      int x = r.GetIndex(d), s = r.GetSize(d);
      if(idx[d] < x)
        {
        r.SetSize(d, s + x - idx[d]);
        r.SetIndex(d, idx[d]);
        }
      else if(idx[d] >= x + s)
        {
        r.SetSize(d, 1 + idx[d] - x);
        } 
      }
    }
}


template<unsigned int VDim>
bool
parse_vector(char *text, itk::Size<VDim> &s)
{
  char *t = strtok(text,"x");
  size_t i = 0;
  while(t && i < VDim)
    {
    s[i++] = atoi(t);
    t = strtok(NULL, "x");
    }

  if(i == VDim)
    return true;
  else if(i == 1)
    {
    s.Fill(s[0]);
    return true;
    }
  return false;
}

template <unsigned int VDim>
int lfapp(int argc, char *argv[])
{
  // Parameter vector
  LFParam<VDim> p;

  // Read the parameters from command line
  p.fnOutput = argv[argc-1];
  p.fnTarget = argv[argc-2];
  int argend = argc-2;

  for(int j = 2; j < argc-2; j++)
    {
    string arg = argv[j];
    
    if(arg == "-g")
      {
      // Read the following options as images
      while(argv[j+1][0] != '-' && j < argend-1)
        p.fnAtlas.push_back(argv[++j]);
      }

    else if(arg == "-l")
      {
      // Read the following options as images
      while(argv[j+1][0] != '-' && j < argend-1)
        p.fnLabel.push_back(argv[++j]);
      }

    else if(arg == "-p")
      {
      p.fnPosterior = argv[++j];
      }

    else if(arg == "-w")
      {
      p.fnWeight = argv[++j];
      }

    else if(arg == "-x")
      {
      int label = atoi(argv[++j]);
      string image = argv[++j];
      p.fnExclusion[label] = image;
      }

    else if(arg == "-m" && j < argend-1)
      {
      char *parm = argv[++j];
      if(!strncmp(parm, "Joint", 5))
        {
        float alpha, beta;
        switch(sscanf(parm, "Joint[%f,%f]", &alpha, &beta))
          {
        case 2:
          p.alpha = alpha; p.beta = beta; break;
        case 1:
          p.alpha = alpha; break;
          }
        p.method = JOINT;
        }
      else if(!strncmp(parm, "Gauss", 5))
        {
        float sigma;
        if(sscanf(parm, "Gauss[%f]", &sigma) == 1)
          p.sigma = sigma;
        p.method = GAUSSIAN;
        }
      else
        {
        cerr << "Unknown method specification " << parm << endl;
        return -1;
        }
      }

    else if(arg == "-rs" && j < argend-1)
      {
      if(!parse_vector<VDim>(argv[++j], p.r_search))
        {
        cerr << "Bad vector spec " << argv[j] << endl;
        return -1;
        }
      }
    
    else if(arg == "-rp" && j < argend-1)
      {
      if(!parse_vector<VDim>(argv[++j], p.r_patch))
        {
        cerr << "Bad vector spec " << argv[j] << endl;
        return -1;
        }
      }

    else if(arg == "-pd" && j < argend-1)
      {
      if(!parse_vector<VDim>(argv[++j], p.paddingSize))
        {
        cerr << "Bad vector spec " << argv[j] << endl;
        return -1;
        }
      p.padding = false;
      for(int d = 0; d < VDim; d++)
        if(p.paddingSize[d] > 0)
          p.padding = true;
      }
    else
      {
      cerr << "Unknown option " << arg << endl;
      return -1;
      }
    }

  // We have the parameters now. Check for validity
  if(p.fnAtlas.size() != p.fnLabel.size())
    {
    cerr << "Number of atlases and segmentations does not match" << endl;
    return -1;
    }

  if(p.fnAtlas.size() < 2)
    {
    cerr << "Too few atlases" << endl;
    return -1;
    }

  // Check the posterior filename pattern
  if(p.fnPosterior.size())
    {
    char buffer[4096];
    sprintf(buffer, p.fnPosterior.c_str(), 100);
    if(strcmp(buffer, p.fnPosterior.c_str()) == 0)
      {
      cerr << "Invalid filename pattern " << p.fnPosterior << endl;
      return -1;
      }
    }

  if(p.fnWeight.size())
    {
    char buffer[4096];
    sprintf(buffer, p.fnWeight.c_str(), 100);
    if(strcmp(buffer, p.fnWeight.c_str()) == 0)
      {
      cerr << "Invalid filename pattern " << p.fnWeight << endl;
      return -1;
      }

    itk::Size<VDim> zeroSize; zeroSize.Fill(0);
    if(p.r_search != zeroSize)
      {
      cerr << "Weight maps can only be used with zero search radius! Use -rs 0x0x0" << endl;
      return -1;
      }
    }

  // Print parametes
  cout << "LABEL FUSION PARAMETERS:" << endl;
  p.Print(cout);

  // Configure the filter
  typedef itk::Image<float, VDim> ImageType;
  typedef WeightedVotingLabelFusionImageFilter<ImageType, ImageType> VoterType;
  typename VoterType::Pointer voter = VoterType::New();

  // Set inputs
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  typename ReaderType::Pointer rTarget = ReaderType::New();
  rTarget->SetFileName(p.fnTarget.c_str());
  rTarget->Update();
  typename ImageType::Pointer target = rTarget->GetOutput();

  if(p.padding)
    {
    typedef itk::MirrorPadImageFilter<ImageType, ImageType> PadFilter;
    typename PadFilter::Pointer padTarget = PadFilter::New();
    padTarget->SetInput(target);
    padTarget->SetPadLowerBound(p.paddingSize);
    padTarget->SetPadUpperBound(p.paddingSize);
    padTarget->Update();
    target = padTarget->GetOutput();
    }

  voter->SetTargetImage(target);

  // Compute the output region by merging all segmentations
  itk::ImageRegion<VDim> rMask;
  bool isMaskInit = false;

  std::vector<typename ReaderType::Pointer> rAtlas, rLabel;
  for(size_t i = 0; i < p.fnAtlas.size(); i++)
    {
    typename ReaderType::Pointer ra, rl;
    ra = ReaderType::New();
    ra->SetFileName(p.fnAtlas[i].c_str());
    ra->Update();
    rAtlas.push_back(ra);
    rl = ReaderType::New();
    rl->SetFileName(p.fnLabel[i].c_str());
    rl->Update();
    rLabel.push_back(rl);

    // Apply padding if requested
    typename ImageType::Pointer imgAtlas = ra->GetOutput(), imgLabel = rl->GetOutput();

    if(p.padding)
      {
      typedef itk::MirrorPadImageFilter<ImageType, ImageType> PadFilter;
      typename PadFilter::Pointer padAtlas = PadFilter::New();
      padAtlas->SetInput(imgAtlas);
      padAtlas->SetPadLowerBound(p.paddingSize);
      padAtlas->SetPadUpperBound(p.paddingSize);
      padAtlas->Update();
      imgAtlas = padAtlas->GetOutput();

      typename PadFilter::Pointer padLabel = PadFilter::New();
      padLabel->SetInput(imgLabel);
      padLabel->SetPadLowerBound(p.paddingSize);
      padLabel->SetPadUpperBound(p.paddingSize);
      padLabel->Update();
      imgLabel = padLabel->GetOutput();
      }

    voter->AddAtlas(imgAtlas, imgLabel);

    // Update the mask region
    for(itk::ImageRegionIteratorWithIndex<ImageType> 
      it(imgLabel, imgLabel->GetBufferedRegion()); 
      !it.IsAtEnd(); ++it)
      {
      if(it.Get())
        {
        ExpandRegion<VDim>(rMask, isMaskInit, it.GetIndex());
        }
      }

    // Unload the image if needed (so that memory is freed up)
    // PY: I DON'T UNDERSTAND THIS - IT'S CAUSING A BUG!
    // rl->GetOutput()->ReleaseData();
    }

  // Make sure the region is inside bounds
  itk::ImageRegion<VDim> rOut = target->GetLargestPossibleRegion();
  for(int d = 0; d < 3; d++)
    {
    rOut.SetIndex(d, p.r_patch[d] + p.r_search[d] + rOut.GetIndex(d));
    rOut.SetSize(d, rOut.GetSize(d) - 2 * (p.r_search[d] + p.r_patch[d]));
    }
  rMask.Crop(rOut);

  voter->SetPatchRadius(p.r_patch);
  voter->SetSearchRadius(p.r_search);
  voter->SetAlpha(p.alpha);
  voter->SetBeta(p.beta);

  // The posterior maps
  if(p.fnPosterior.size())
    voter->SetRetainPosteriorMaps(true);

  if(p.fnWeight.size())
    voter->SetGenerateWeightMaps(true);

  std::cout << "Output Requested Region: " << rMask.GetIndex() << ", " << rMask.GetSize() << " ("
    << rMask.GetNumberOfPixels() << " pixels)" << std::endl;
  /*
  rr.SetIndex(0,275); rr.SetSize(0,10);
  rr.SetIndex(1,210); rr.SetSize(1,10);
  rr.SetIndex(2,4); rr.SetSize(2,10);
  */

  // Set the exclusions in the atlas
  for(typename map<int,string>::iterator xit = p.fnExclusion.begin(); xit != p.fnExclusion.end(); ++xit)
    {
    typename ReaderType::Pointer rx;
    rx = ReaderType::New();
    rx->SetFileName(xit->second.c_str());
    rx->Update();
    typename ImageType::Pointer xmap = rx->GetOutput();

    if(p.padding)
      {
      typedef itk::MirrorPadImageFilter<ImageType, ImageType> PadFilter;
      typename PadFilter::Pointer padXMap = PadFilter::New();
      padXMap->SetInput(xmap);
      padXMap->SetPadLowerBound(p.paddingSize);
      padXMap->SetPadUpperBound(p.paddingSize);
      padXMap->Update();
      xmap = padXMap->GetOutput();
      }

    voter->AddExclusionMap(xit->first, xmap);
    }

  voter->GetOutput()->SetRequestedRegion(rMask);
  voter->Update();

  // Convert to an output image
  target->FillBuffer(0.0);
  for(itk::ImageRegionIteratorWithIndex<ImageType> it(voter->GetOutput(), rMask);
    !it.IsAtEnd(); ++it)
    {
    target->SetPixel(it.GetIndex(), it.Get());
    }

  // If the image was padded, crop it
  typename ImageType::Pointer outSeg = target;
  if(p.padding)
    {
    typedef itk::CropImageFilter<ImageType, ImageType> CropFilter;
    typename CropFilter::Pointer fltCrop = CropFilter::New();
    fltCrop->SetInput(target);
    fltCrop->SetUpperBoundaryCropSize(p.paddingSize);
    fltCrop->SetLowerBoundaryCropSize(p.paddingSize);
    fltCrop->Update();
    outSeg = fltCrop->GetOutput();
    }


  // Create writer
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(outSeg);
  writer->SetFileName(p.fnOutput.c_str());
  writer->Update();

  // Store the weight maps
  if(p.fnWeight.size())
    {
    for(int i = 0; i < p.fnAtlas.size(); i++)
      {
      // Create a full-size image
      typename VoterType::WeightMapImagePtr wout = VoterType::WeightMapImage::New();
      wout->SetRegions(target->GetBufferedRegion());
      wout->CopyInformation(target);
      wout->Allocate();
      wout->FillBuffer(1.0f / p.fnAtlas.size());

      // Replace the portion affected by the filter
      for(itk::ImageRegionIteratorWithIndex<typename VoterType::WeightMapImage> qt(voter->GetWeightMap(i), rMask);
        !qt.IsAtEnd(); ++qt)
        {
        wout->SetPixel(qt.GetIndex(), qt.Get());
        }

      // Deal with padding
      typename VoterType::WeightMapImagePtr imgWeightOut = wout;

      if(p.padding)
        {
        typedef itk::CropImageFilter<typename VoterType::WeightMapImage, typename VoterType::WeightMapImage> CropFilter;
        typename CropFilter::Pointer fltCrop = CropFilter::New();
        fltCrop->SetInput(wout);
        fltCrop->SetUpperBoundaryCropSize(p.paddingSize);
        fltCrop->SetLowerBoundaryCropSize(p.paddingSize);
        fltCrop->Update();
        imgWeightOut = fltCrop->GetOutput();
        }

      // Get the filename
      char buffer[4096];
      sprintf(buffer, p.fnWeight.c_str(), i);

      // Create writer
      typedef itk::ImageFileWriter<typename VoterType::WeightMapImage> WeightWriter;
      typename WeightWriter::Pointer writer = WeightWriter::New();
      writer->SetInput(imgWeightOut);
      writer->SetFileName(buffer);
      writer->Update();
      }
    }

  // Finally, store the posterior maps
  if(p.fnPosterior.size())
    {
    // Get the posterior maps
    const typename VoterType::PosteriorMap &pm = voter->GetPosteriorMaps();

    // Create a full-size image
    typename VoterType::PosteriorImagePtr pout = VoterType::PosteriorImage::New();
    pout->SetRegions(target->GetBufferedRegion());
    pout->CopyInformation(target);
    pout->Allocate();

    // Iterate over the labels
    typename VoterType::PosteriorMap::const_iterator it;
    for(it = pm.begin(); it != pm.end(); it++)
      {
      // Get the filename
      char buffer[4096];
      sprintf(buffer, p.fnPosterior.c_str(), (int) it->first);

      // Initialize to zeros
      pout->FillBuffer(0.0);

      // Replace the portion affected by the filter
      for(itk::ImageRegionIteratorWithIndex<typename VoterType::PosteriorImage> qt(it->second, rMask);
        !qt.IsAtEnd(); ++qt)
        {
        pout->SetPixel(qt.GetIndex(), qt.Get());
        }

      typename ImageType::Pointer imgPostOut = pout;

      if(p.padding)
        {
        typedef itk::CropImageFilter<ImageType, ImageType> CropFilter;
        typename CropFilter::Pointer fltCrop = CropFilter::New();
        fltCrop->SetInput(pout);
        fltCrop->SetUpperBoundaryCropSize(p.paddingSize);
        fltCrop->SetLowerBoundaryCropSize(p.paddingSize);
        fltCrop->Update();
        imgPostOut = fltCrop->GetOutput();
        }

      // Create writer
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput(imgPostOut);
      writer->SetFileName(buffer);
      writer->Update();
      }
    }

  return 0;
}


int main(int argc, char *argv[])
{
  // Parse user input
  if(argc < 5) return usage();

  // Get the first option
  int dim = atoi(argv[1]);
  
  // Call the templated method
  if(dim == 2)
    return lfapp<2>(argc, argv);
  else if(dim == 3)
    return lfapp<3>(argc, argv);
  else
    {
    cerr << "Dimension " << argv[1] << " is not supported" << endl;
    return -1;
    }
}
