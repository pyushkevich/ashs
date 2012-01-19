#include "WeightedVotingLabelFusionImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <iostream>

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
  cout << "  -x label image.nii              Specify an exclusion region for the given label. " << endl;
  cout << "                                  If a voxel has non-zero value in an exclusion image," << endl;
  cout << "                                  the corresponding label is not allowed at that voxel." << endl;
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
  LFMethod method;

  map<int, string> fnExclusion;

  double alpha, beta, sigma;
  itk::Size<VDim> r_patch, r_search;
  
  LFParam()
    {
    alpha = 0.01;
    beta = 2;
    sigma = 0.5;
    r_patch.Fill(3);
    r_search.Fill(3);
    method = JOINT;
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

  // Check if the exclusion pattern makes sense
  char test[4096];
  if(sprintf(test, fnExclusionPattern.c_str(), 1) < 0)
    {
    cerr << "The exclusion pattern %s does not have proper C printf format\n";
    return -1;
    }

  // Print parametes
  cout << "LABEL FUSION PARAMETERS:" << endl;
  p.Print(cout);

  // Configure the filter
  typedef itk::Image<float, VDim> ImageType;
  typedef WeightedVotingLabelFusionImageFilter<ImageType, ImageType> VoterType;
  typename VoterType::Pointer voter = VoterType::New();
  voter->SetNumberOfAtlases(p.fnAtlas.size());

  // Set inputs
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  typename ReaderType::Pointer rTarget = ReaderType::New();
  rTarget->SetFileName(p.fnTarget.c_str());
  rTarget->Update();
  typename ImageType::Pointer target = rTarget->GetOutput();
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
    rAtlas.push_back(ra);
    rl = ReaderType::New();
    rl->SetFileName(p.fnLabel[i].c_str());
    rLabel.push_back(rl);
    voter->SetAtlas(i, ra->GetOutput(), rl->GetOutput());

    // Update the mask region
    rl->Update();
    for(itk::ImageRegionIteratorWithIndex<ImageType> 
      it(rl->GetOutput(), rl->GetOutput()->GetBufferedRegion()); 
      !it.IsAtEnd(); ++it)
      {
      if(it.Get())
        {
        ExpandRegion<VDim>(rMask, isMaskInit, it.GetIndex());
        }
      }

    // Unload the image if needed (so that memory is freed up)
    rl->GetOutput()->ReleaseData();
    }

  // Make sure the region is inside bounds
  itk::ImageRegion<VDim> rOut = rTarget->GetOutput()->GetLargestPossibleRegion();
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

  std::cout << "Output Requested Region: " << rMask.GetIndex() << ", " << rMask.GetSize() << " ("
    << rMask.GetNumberOfPixels() << " pixels)" << std::endl;
  /*
  rr.SetIndex(0,275); rr.SetSize(0,10);
  rr.SetIndex(1,210); rr.SetSize(1,10);
  rr.SetIndex(2,4); rr.SetSize(2,10);
  */
  voter->GetOutput()->SetRequestedRegion(rMask);
  voter->Update();

  // Convert to an output image
  target->FillBuffer(0.0);
  for(itk::ImageRegionIteratorWithIndex<ImageType> it(voter->GetOutput(), rMask);
    !it.IsAtEnd(); ++it)
    {
    target->SetPixel(it.GetIndex(), it.Get());
    }

  // Create writer
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(target);
  writer->SetFileName(p.fnOutput.c_str());
  writer->Update();

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
