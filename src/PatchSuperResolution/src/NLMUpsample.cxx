//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%
//%   Jose V. Manjon - jmanjon@fis.upv.es
//%   Universidad Politecinca de Valencia, Spain
//%   Pierrick Coupe - pierrick.coupe@gmail.com
//%   Brain Imaging Center, Montreal Neurological Institute.
//%   Mc Gill University
//%
//%   Copyright (C) 2010 Jose V. Manjon and Pierrick Coupe
//%
//%    Usage of NLMUpsample:
//%
//%    fima=NLMUpsample(ima,f)
//%
//%    ima: LR volume
//%    f: Magnification factor in each dimension (for example [2,2,2])
//%    fima: HR upsampled volume
//%
//%**************************************************************************


#include <iostream>
#include <string>
#include "math.h"
#include <cstdlib>
#include "CommandLineHelper.h"
#include "NLMUpsample.h"
#include "util/itkOrientedRASImage.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkNoiseImageFilter.h"

/* Multithreading stuff*/
#ifdef _WIN32
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#endif

#define pi 3.1415926535

using namespace std;

int usage()
{
  cout << "NLMUpsample: none local upsample for medical images" << endl;
  cout << "Usage:" << endl;
  cout << "  NLMUpsample [options]" << endl;
  cout << "Required Options:" << endl;
  cout << "  -i input.nii.gz            : input image" << endl;
  cout << "  -o output.nii.gz           : output image" << endl;
  cout << "Additional Options" << endl;
  cout << "  -v search radius           : search radius (3)" << endl;
  cout << "  -f patch radius            : patch radius (1)" << endl;
  cout << "  -lf upsample factor        : upsample factor in x, y and z (2 2 2)" << endl;
  return -1;
}

void check(bool condition, char *format,...)
{
  if(!condition)
    {
      char buffer[256];
      va_list args;
      va_start (args, format);
      vsprintf (buffer,format, args);
      va_end (args);

      cerr << buffer << endl;
      exit(-1);
    }
}

struct UpsampleParameters
{
  string fnImage;
  string fnOutput;
  int v;
  int f;
  int lf_x;
  int lf_y;
  int lf_z;

  UpsampleParameters() :
    v(3), f(1), lf_x(2), lf_y(2), lf_z(2){}
};


typedef struct{
  int rows;
  int cols;
  int slices;
  double * in_image;
  double * out_image;
  double * mean_image;
  double * pesos;
  int ini;
  int fin;
  int radio;
  int f;
  /*int th;*/
  double * sigma;
}myargument;


template <class TFloat, unsigned int VDim>
class NLMUpsampleProblem
{
public:
  typedef itk::Image<TFloat, VDim> ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  typedef typename ImageReaderType::Pointer  ImageReaderPointer;
  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  typedef typename ImageWriterType::Pointer  ImageWriterPointer;

  static double distancia(double* ima,int x,int y,int z,int nx,int ny,int nz,int f,int sx,int sy,int sz);

  static void* ThreadFunc( void* pArguments );

  static void PatchCorrection(double* ima, double* fima, int* dims, double* h, UpsampleParameters param);

  static ImagePointer Resample(ImagePointer ITKima, UpsampleParameters param);

  static void MeanCorrection(double * bima, double * ref, double * fima, int* dims, UpsampleParameters param);

  static ImagePointer ComputeLevel(ImagePointer input, int radius);

  static void main_loop(UpsampleParameters param);

private:
};


template <class TFloat, unsigned int VDim>
double
NLMUpsampleProblem<TFloat, VDim>
::distancia(double* ima,int x,int y,int z,int nx,int ny,int nz,int f,int sx,int sy,int sz)
{
  double d,acu,distancetotal,inc;
  int i,j,k,ni1,nj1,ni2,nj2,nk1,nk2;

  distancetotal=0;

  for(k=-f;k<=f;k++)
    {
      nk1=z+k;
      nk2=nz+k;
      if(nk1<0) nk1=-nk1;
      if(nk2<0) nk2=-nk2;
      if(nk1>=sz) nk1=2*sz-nk1-1;
      if(nk2>=sz) nk2=2*sz-nk2-1;

      for(j=-f;j<=f;j++)
        {
          nj1=y+j;
          nj2=ny+j;
          if(nj1<0) nj1=-nj1;
          if(nj2<0) nj2=-nj2;
          if(nj1>=sy) nj1=2*sy-nj1-1;
          if(nj2>=sy) nj2=2*sy-nj2-1;

          for(i=-f;i<=f;i++)
            {
              ni1=x+i;
              ni2=nx+i;
              if(ni1<0) ni1=-ni1;
              if(ni2<0) ni2=-ni2;
              if(ni1>=sx) ni1=2*sx-ni1-1;
              if(ni2>=sx) ni2=2*sx-ni2-1;

              distancetotal = distancetotal + ((ima[nk1*(sx*sy)+(nj1*sx)+ni1]-ima[nk2*(sx*sy)+(nj2*sx)+ni2])*(ima[nk1*(sx*sy)+(nj1*sx)+ni1]-ima[nk2*(sx*sy)+(nj2*sx)+ni2]));
            }
        }
    }

  acu=(2*f+1)*(2*f+1)*(2*f+1);
  d=distancetotal/acu;

  return d;

}

template <class TFloat, unsigned int VDim>
void*
NLMUpsampleProblem<TFloat, VDim>
::ThreadFunc( void* pArguments )
{
  double *ima,*fima,*medias,*pesos,w,d,*h,t1;
  int ii,jj,kk,ni,nj,nk,i,j,k,ini,fin,rows,cols,slices,v,p,p1,f,rc;

  myargument arg;
  arg=*(myargument *) pArguments;

  rows=arg.rows;
  cols=arg.cols;
  slices=arg.slices;
  ini=arg.ini;
  fin=arg.fin;
  ima=arg.in_image;
  fima=arg.out_image;
  medias=arg.mean_image;
  pesos=arg.pesos;
  v=arg.radio;
  f=arg.f;
  /*th=arg.th;*/
  h=arg.sigma;
  rc=rows*cols;

  /* filter*/
  for(k=ini;k<fin;k++)
    {
      for(j=0;j<rows;j++)
        {
          for(i=0;i<cols;i++)
            {
              p=k*rc+j*cols+i;
              if(h[p]<1) continue;

              for(kk=0;kk<=v;kk++)
                {
                  nk=k+kk;
                  for(ii=-v;ii<=v;ii++)
                    {
                      ni=i+ii;
                      for(jj=-v;jj<=v;jj++)
                        {
                          nj=j+jj;

                          if(kk==0 && jj<0) continue;
                          if(kk==0 && jj==0 && ii<=0) continue;

                          if(ni>=0 && nj>=0 && nk>=0 && ni<cols && nj<rows && nk<slices)
                            {
                              p1=nk*rc+nj*cols+ni;

                              t1 = abs(medias[p]-medias[p1]);

                              if(t1>0.6*h[p]) continue;

                              d=distancia(ima,i,j,k,ni,nj,nk,f,cols,rows,slices);

                              if(h[p]>0) d=d/(2*h[p]*h[p])-1;
                              if(d<0) d=0;

                              w = exp(-d);

                              fima[p] = fima[p] + w*ima[p1];
                              pesos[p] = pesos[p] + w;

                              fima[p1] = fima[p1] + w*ima[p];
                              pesos[p1] = pesos[p1] + w;
                            }
                        }
                    }
                }
            }
        }
    }

#ifdef _WIN32
  _endthreadex(0);
#else
  pthread_exit(0);
#endif

  return 0;
}

template <class TFloat, unsigned int VDim>
void
NLMUpsampleProblem<TFloat, VDim>
::PatchCorrection(double* ima, double* fima, int* dims, double* h, UpsampleParameters param)
{
  // read inputs
  int v = param.v;
  int f = param.f;
  int fac[3];
  fac[0] = param.lf_x;
  fac[1] = param.lf_y;
  fac[2] = param.lf_z;
  int nelement = dims[2] * dims[1] * dims[0];

  // multi thread
  myargument *ThreadArgs;
#ifdef _WIN32
  HANDLE *ThreadList; /* Handles to the worker threads*/
#else
  pthread_t * ThreadList;
#endif

  // allocate memory for intermediate matrices
  double *tmp = new double [nelement];
  double *medias = new double [nelement];
  double *pesos = new double [nelement];

  /*Declarations*/
  double off,media;
  int ini,fin,i,j,k,ii,jj,kk,ni,nj,nk,indice,Nthreads,rc,ft;
  bool salir;


  //if(nrhs<5)
  //  {
  //    printf("Wrong number of arguments!!!\r");
  //   exit;
  //}

  /*Copy input pointer x*/
  //xData = prhs[0];

  /*Get matrix x*/
  //ima = mxGetPr(xData);

  //ndim = mxGetNumberOfDimensions(prhs[0]);
  //dims= mxGetDimensions(prhs[0]);

  //pv = prhs[1];
  //v = (int)(mxGetScalar(pv));
  //pv = prhs[2];
  //f = (int)(mxGetScalar(pv));
  //xData = prhs[3];
  //h = (double*) mxGetPr(xData);
  //pv = prhs[4];
  //lf = (double*)(mxGetPr(pv));
  //for(i=0;i<3;i++) fac[i]=(int)lf[i];

  /*Allocate memory and assign output pointer*/

  //plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);
  //Mxmedias = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);
  //xtmp = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);
  //tmp = mxGetPr(xtmp);

  /*Get a pointer to the data space in our newly allocated memory*/
  //fima = mxGetPr(plhs[0]);
  //medias = mxGetPr(Mxmedias);

  //Mxpesos = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);
  //pesos = mxGetPr(Mxpesos);



  // start computation
  rc=dims[0]*dims[1];
  for(k=0;k<dims[2];k++)
    {
      for(j=0;j<dims[1];j++)
        {
          for(i=0;i<dims[0];i++)
            {
              media=0;
              for(ii=-1;ii<=1;ii++)
                {
                  ni=i+ii;
                  if(ni<0) ni=-ni;
                  if(ni>=dims[0]) ni=2*dims[0]-ni-1;
                  for(jj=-1;jj<=1;jj++)
                    {
                      nj=j+jj;
                      if(nj<0) nj=-nj;
                      if(nj>=dims[1]) nj=2*dims[1]-nj-1;
                      for(kk=-1;kk<=1;kk++)
                        {
                          nk=k+kk;
                          if(nk<0) nk=-nk;
                          if(nk>=dims[2]) nk=2*dims[2]-nk-1;

                          media = media + ima[nk*rc+nj*dims[0]+ni];
                        }
                    }
                }
              medias[k*rc+j*dims[0]+i]=media/27;
            }
        }
    }

  ft=fac[2]*fac[1]*fac[0];
  for(k=0;k<dims[2]/fac[2];k++)
    for(j=0;j<dims[1]/fac[1];j++)
      for(i=0;i<dims[0]/fac[0];i++)
        {
          media=0;
          for (kk=0;kk<fac[2];kk++)
            for (jj=0;jj<fac[1];jj++)
              for (ii=0;ii<fac[0];ii++)
                media+=ima[(k*fac[2]+kk)*rc + (j*fac[1]+jj)*dims[0] + i*fac[0]+ii];
          tmp[k*rc+(j*dims[0])+i]=media/ft;
        }

  for(k=0;k<dims[2]*dims[1]*dims[0];k++)
    {
      fima[k]=ima[k];
      pesos[k]=1.0;
    }
  /*th=0.6*h;  //3/sqrt(27) */

  // start multi-thread computation
  //Nthreads=dims[2]<8?dims[2]:8;
  int MaxThreadN = (int) (itk::MultiThreader::GetGlobalDefaultNumberOfThreads() * 0.8);
  if( MaxThreadN < 1 )
  {
    MaxThreadN = 1;
  }
  //int MaxThreadN = 8;
  Nthreads = dims[2]<MaxThreadN?dims[2]:MaxThreadN;
  //Nthreads = (int) (itk::MultiThreader::GetGlobalMaximumNumberOfThreads() * 0.5);

#ifdef _WIN32

  // Reserve room for handles of threads in ThreadList
  ThreadList = (HANDLE*)malloc(Nthreads* sizeof( HANDLE ));
  ThreadArgs = (myargument*) malloc( Nthreads*sizeof(myargument));

  for (i=0; i<Nthreads; i++)
    {
      // Make Thread Structure
      ini=(i*dims[2])/Nthreads;
      fin=((i+1)*dims[2])/Nthreads;

      ThreadArgs[i].cols=dims[0];
      ThreadArgs[i].rows=dims[1];
      ThreadArgs[i].slices=dims[2];
      ThreadArgs[i].in_image=ima;
      ThreadArgs[i].out_image=fima;
      ThreadArgs[i].mean_image=medias;
      ThreadArgs[i].pesos=pesos;
      ThreadArgs[i].ini=ini;
      ThreadArgs[i].fin=fin;
      ThreadArgs[i].radio=v;
      ThreadArgs[i].f=f;
      ThreadArgs[i].sigma=h;

      ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &ThreadFunc, &ThreadArgs[i] , 0, NULL );
    }

  for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
  for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }

#else

  /* Reserve room for handles of threads in ThreadList*/
  ThreadList = (pthread_t *) calloc(Nthreads,sizeof(pthread_t));
  ThreadArgs = (myargument*) calloc( Nthreads,sizeof(myargument));

  for (i=0; i<Nthreads; i++)
    {
      /* Make Thread Structure */
      ini=(i*dims[2])/Nthreads;
      fin=((i+1)*dims[2])/Nthreads;

      ThreadArgs[i].cols=dims[0];
      ThreadArgs[i].rows=dims[1];
      ThreadArgs[i].slices=dims[2];
      ThreadArgs[i].in_image=ima;
      ThreadArgs[i].out_image=fima;
      ThreadArgs[i].mean_image=medias;
      ThreadArgs[i].pesos=pesos;
      ThreadArgs[i].ini=ini;
      ThreadArgs[i].fin=fin;
      ThreadArgs[i].radio=v;
      ThreadArgs[i].f=f;
      ThreadArgs[i].sigma=h;
    }
  for (i=0; i<Nthreads; i++)
    {
      if(pthread_create(&ThreadList[i], NULL, ThreadFunc,&ThreadArgs[i]))
        {
          printf("Threads cannot be created\n");
          exit(1);
        }
    }

  for (i=0; i<Nthreads; i++)
    {
      pthread_join(ThreadList[i],NULL);
    }

#endif
  free(ThreadArgs);
  free(ThreadList);

  for(i=0;i<dims[0]*dims[1]*dims[2];i++) fima[i]/=pesos[i];


  /* apply mean preservation constraint*/
  for(k=0;k<dims[2];k=k+fac[2])
    for(j=0;j<dims[1];j=j+fac[1])
      for(i=0;i<dims[0];i=i+fac[0])
        {
          salir=false;

          media=0;
          for (kk=0;kk<fac[2];kk++)
            for (jj=0;jj<fac[1];jj++)
              for (ii=0;ii<fac[0];ii++)
                media+=fima[(k+kk)*rc + (j+jj)*dims[0] + i+ii];
          media=media/ft;

          off=tmp[(k/fac[2])*rc+(j/fac[1])*dims[0]+(i/fac[0])]-media;

          for (kk=0;kk<fac[2];kk++)
            for (jj=0;jj<fac[1];jj++)
              for (ii=0;ii<fac[0];ii++)
                {
                  fima[(k+kk)*rc + (j+jj)*dims[0] + (i+ii)]+=off;
                }
        }

  // free memory
  delete tmp;
  delete medias;
  delete pesos;

  //return;
}

template <class TFloat, unsigned int VDim>
typename NLMUpsampleProblem<TFloat, VDim>::ImagePointer
NLMUpsampleProblem<TFloat, VDim>
::Resample(ImagePointer ITKima, UpsampleParameters param)
{
  // code from Paul's c3d

  // Get the size of the image in voxels
  typename ImageType::SizeType sz;
  sz[0] = (unsigned long)(ITKima->GetBufferedRegion().GetSize(0) * param.lf_x + 0.5);
  sz[1] = (unsigned long)(ITKima->GetBufferedRegion().GetSize(1) * param.lf_y + 0.5);
  sz[2] = (unsigned long)(ITKima->GetBufferedRegion().GetSize(2) * param.lf_z + 0.5);

  // Build the resampling filter
  typedef itk::ResampleImageFilter<ImageType,ImageType> ResampleFilterType;
  typename ResampleFilterType::Pointer fltSample = ResampleFilterType::New();

  // Specify interpolator
  typedef itk::BSplineInterpolateImageFunction<ImageType,TFloat> CubicInterpolatorType;

  // Initialize the resampling filter with an identity transform
  fltSample->SetInput(ITKima);
  fltSample->SetTransform(itk::IdentityTransform<TFloat,VDim>::New());
  fltSample->SetInterpolator(CubicInterpolatorType::New());

  // Compute the spacing of the new image
  typename ImageType::SpacingType spc_pre = ITKima->GetSpacing();
  typename ImageType::SpacingType spc_post = spc_pre;
  for(size_t i = 0; i < 3; i++)
    spc_post[i] *= ITKima->GetBufferedRegion().GetSize()[i] * 1.0 / sz[i];

  // Get the bounding box of the input image
  typename ImageType::PointType origin_pre = ITKima->GetOrigin();

  // Recalculate the origin. The origin describes the center of voxel 0,0,0
  // so that as the voxel size changes, the origin will change as well.
  typename ImageType::SpacingType off_pre = (ITKima->GetDirection() * spc_pre) * 0.5;
  typename ImageType::SpacingType off_post = (ITKima->GetDirection() * spc_post) * 0.5;
  typename ImageType::PointType origin_post = origin_pre - off_pre + off_post;

  // Set the image sizes and spacing.
  fltSample->SetSize(sz);
  fltSample->SetOutputSpacing(spc_post);
  fltSample->SetOutputOrigin(origin_post);
  fltSample->SetOutputDirection(ITKima->GetDirection());

  // Perform resampling
  fltSample->UpdateLargestPossibleRegion();
  return fltSample->GetOutput();
}

template <class TFloat, unsigned int VDim>
void
NLMUpsampleProblem<TFloat, VDim>
::MeanCorrection(double * bima, double * ref, double * fima, int* dims, UpsampleParameters param)
{
  /* Declarations */
  double media, off;
  int i,j,k,ii,jj,kk;

  // compute size
  int lf[3];
  lf[0] = param.lf_x;
  lf[1] = param.lf_y;
  lf[2] = param.lf_z;
  int sx=dims[0]/lf[0];
  int sy=dims[1]/lf[1];
  int sz=dims[2]/lf[2];

  /*  mean correction */
  for(k=0;k<dims[2];k=k+lf[2])
    for(j=0;j<dims[1];j=j+lf[1])
      for(i=0;i<dims[0];i=i+lf[0])
        {
          media=0;
          for (kk=0;kk<lf[2];kk++)
            for (jj=0;jj<lf[1];jj++)
              for (ii=0;ii<lf[0];ii++)
                media=media+bima[(k+kk)*dims[0]*dims[1] + (j+jj)*dims[0] + i+ii];
          media=media/(lf[0]*lf[1]*lf[2]);

          off=ref[(int)((k/lf[2])*sx*sy+(j/lf[1])*sx+(i/lf[0]))]-media;

          for (kk=0;kk<lf[2];kk++)
            for (jj=0;jj<lf[1];jj++)
              for (ii=0;ii<lf[0];ii++)
                fima[(k+kk)*dims[0]*dims[1] + (j+jj)*dims[0] + i+ii]=bima[(k+kk)*dims[0]*dims[1] + (j+jj)*dims[0] + i+ii]+off;
        }
}

template <class TFloat, unsigned int VDim>
typename NLMUpsampleProblem<TFloat, VDim>::ImagePointer
NLMUpsampleProblem<TFloat, VDim>
::ComputeLevel(ImagePointer input, int radius)
{
  // compute local std map
  typedef itk::NoiseImageFilter< ImageType, ImageType > NoiseImageFilterType;
  typename NoiseImageFilterType::Pointer noiseImageFilter = NoiseImageFilterType::New();
  noiseImageFilter->SetInput(input);
  noiseImageFilter->SetRadius(radius);
  noiseImageFilter->Update();

  // perform mean filter
  typedef itk::MeanImageFilter< ImageType, ImageType > MeanImageFilterType;
  typename MeanImageFilterType::Pointer meanImageFilter = MeanImageFilterType::New();
  meanImageFilter->SetInput(noiseImageFilter->GetOutput());
  meanImageFilter->SetRadius(radius);
  meanImageFilter->Update();

  return meanImageFilter->GetOutput();
}


template <class TFloat, unsigned int VDim>
void
NLMUpsampleProblem<TFloat, VDim>
::main_loop(UpsampleParameters param)
{
  // dimension
  const unsigned int ndim = 3;

  // read image
  ImageReaderPointer reader = ImageReaderType::New();
  reader->SetFileName(param.fnImage.c_str());
  try { reader->Update(); }
  catch(itk::ExceptionObject &exc)
  {
    cerr << "Error reading image " << param.fnImage.c_str() << endl;
    cerr << "ITK Exception: " << exc << endl;
    throw -1;
  }
  ImagePointer ITKima = reader->GetOutput();
  double *ima = ITKima->GetBufferPointer();

  // get image size
  itk::Size<ndim> ITKdims = ITKima->GetBufferedRegion().GetSize();
  int* dims = new int[3];
  dims[0] = ITKdims[0];
  dims[1] = ITKdims[1];
  dims[2] = ITKdims[2];
  int n = dims[2] * dims[1] * dims[0];

  // normalize signal to [0 256]
  double iMax = ima[0];
  for(size_t i = 0; i < n; i++)
    iMax = (iMax > ima[i]) ? iMax : ima[i];
  if(iMax == 0)
    {
      cerr << "Maximum of image equal to 0." << endl;
      throw -1;
    }
  for(size_t i = 0; i < n; i++)
    ima[i] = 256 * (ima[i]/iMax);

  // perform initial interpolation
  ImagePointer ITKbima = Resample(ITKima, param);
  double* bima = ITKbima->GetBufferPointer();

  // update dimension info
  ITKdims = ITKbima->GetBufferedRegion().GetSize();
  dims[0] = ITKdims[0];
  dims[1] = ITKdims[1];
  dims[2] = ITKdims[2];
  n = dims[2] * dims[1] * dims[0];

  // multi-thread
  int MaxThreadN = (int) (itk::MultiThreader::GetGlobalDefaultNumberOfThreads() * 0.8);
  if( MaxThreadN < 1 )
  {
    MaxThreadN = 1;
  }
  int Nthreads = dims[2]<MaxThreadN?dims[2]:MaxThreadN;
  //int Nthreads = (int) (itk::MultiThreader::GetGlobalMaximumNumberOfThreads() * 0.5);
  cout << "Using " << Nthreads << " threads." << endl;

  // allocate memory for output image
  ImagePointer ITKfima = ImageType::New();
  ITKfima->SetRegions(ITKbima->GetBufferedRegion());
  ITKfima->SetSpacing(ITKbima->GetSpacing());
  ITKfima->SetOrigin(ITKbima->GetOrigin());
  ITKfima->SetDirection(ITKbima->GetDirection());
  ITKfima->SetMetaDataDictionary(ITKbima->GetMetaDataDictionary());
  ITKfima->Allocate();
  double *fima = ITKfima->GetBufferPointer();

  // perform mean correction
  MeanCorrection(bima, ima, fima, dims, param);

  // compute level (sigma
  int radius = 1;
  ImagePointer ITKlevel = ComputeLevel(ITKfima, radius);
  double* level = ITKlevel->GetBufferPointer();

  // perform iteration
  int maxiter = 1000;
  double mean[maxiter];
  double mean_down[maxiter];
  double tol = 1.2;
  int down = 0;
  int ii = 1, iii = 1;
  double* tmpima = new double [n];
  for(size_t i = 0; i < n; i++)
    {
      bima[i] = fima[i];
      tmpima[i] = fima[i];
    }

  while(1)
    {
      // perform upsample
      PatchCorrection(bima, fima, dims, level, param);

      // compute average differences between ima and fima
      mean[ii] = 0;
      for(size_t i = 0; i < n; i++)
        mean[ii] += abs(bima[i] - fima[i]);
      mean[ii] /= n;
      cout << "Iter " << ii << ": down fac = " << iii
           << "; abs mean diff = " << mean[ii] << "." << endl;

      // determine continue or not
      if(ii>1)
        {
          if(mean[ii-1]/mean[ii] < tol && down == 0)
            {
              down = 1;
              for(size_t i = 0; i < n; i++)
                {
                  level[i] = level[i]/2;
                  mean_down[iii] += abs(tmpima[i] - fima[i]);
                  tmpima[i] = fima[i];
                }
              if(iii > 1)
                if(mean_down[iii-1]/mean_down[iii] < tol)
                  break;
              iii += 1;
            }
          else
            down = 0;

          if(mean[ii] <= 0.001)
            break;

          if(ii >= maxiter)
            break;
        }

      // prepare for the next iteration
      for(size_t i = 0; i < n; i++)
        bima[i] = fima[i];
      ii += 1;
    }

  // rescale it back to original range
  for(size_t i = 0; i < n; i++)
    fima[i] = iMax * (fima[i]/256);

  // write output image
  ImageWriterPointer writer = ImageWriterType::New();
  writer->SetInput(ITKfima);
  writer->SetFileName(param.fnOutput.c_str());
  try { writer->Update(); }
  catch (itk::ExceptionObject &exc) {
    cerr << "Error writing image to " << param.fnOutput.c_str() << endl;
    cerr << "ITK Exception: " << exc << endl;
    throw -1;
  }
}



int main(int argc, char *argv[])
{
  if(argc < 2)
    return usage();

  UpsampleParameters param;

  // Process parameters
  CommandLineHelper cl(argc, argv);

  while(!cl.is_at_end())
    {
      // Read the next command
      std::string arg = cl.read_command();

      if(arg == "-i")
        {
          param.fnImage = cl.read_existing_filename();
        }
      else if(arg == "-o")
        {
          param.fnOutput = cl.read_output_filename();
        }
      else if(arg == "-v")
        {
          param.v = cl.read_integer();
        }
      else if(arg == "-f")
        {
          param.f = cl.read_integer();
        }
      else if(arg == "-lf")
        {
          param.lf_x = cl.read_integer();
          param.lf_y = cl.read_integer();
          param.lf_z = cl.read_integer();
        }
      else
        {
          cerr << "Unknown option " << arg << endl;
          return -1;
        }
    }

  // Check parameters
  check(param.v >= 1, "Search radius (v) must a positive integer.");
  check(param.f >= 1, "Patch radius (f) must a positive integer.");
  check(param.lf_x >= 1, "Upsample factor in x (lf_x) must a positive integer.");
  check(param.lf_y >= 1, "Upsample factor in y (lf_y) must a positive integer.");
  check(param.lf_z >= 1, "Upsample factor in z (lf_z) must a positive integer.");

  NLMUpsampleProblem<double, 3>::main_loop(param);
  //main_loop(param);

  cout << "Done NLM upsampling." << endl;

}
