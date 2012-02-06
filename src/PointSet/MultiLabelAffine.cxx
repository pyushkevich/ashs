#include <itkImageFileReader.h>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkOrientedRASImage.h>
#include "itkVTKImageExport.h"
#include "vtkImageImport.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkMarchingCubes.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkImplicitFunction.h"
#include "vtkQuadricClustering.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolyDataWriter.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkCellArray.h"
#include "vtkVertex.h"
#include <vtkMatrix4x4.h>
#include "vnl/vnl_matrix_fixed.h"
#include "itk_to_nifti_xform.h"
#include <string>
#include <vector>
#include <set>

using namespace std;
using namespace itk;

// Global datatypes
typedef itk::OrientedRASImage<short,3> LabelImageType;
typedef itk::ImageFileReader<LabelImageType> LabelImageReader;

struct Parameters
{
  // Filenames
  string fnMoving, fnTarget, fnOutput;
};

int usage()
{
  cout << "ml_affine: multi-label image affine registration" << endl
       << "usage:" << endl
       << "  ml_affine [options] target.nii moving.nii output.txt" << endl;
  return -1;
}

template<class TImage>
void ConnectITKToVTK(itk::VTKImageExport<TImage> *fltExport,vtkImageImport *fltImport)
{
  fltImport->SetUpdateInformationCallback( fltExport->GetUpdateInformationCallback());
  fltImport->SetPipelineModifiedCallback( fltExport->GetPipelineModifiedCallback());
  fltImport->SetWholeExtentCallback( fltExport->GetWholeExtentCallback());
  fltImport->SetSpacingCallback( fltExport->GetSpacingCallback());
  fltImport->SetOriginCallback( fltExport->GetOriginCallback());
  fltImport->SetScalarTypeCallback( fltExport->GetScalarTypeCallback());
  fltImport->SetNumberOfComponentsCallback( fltExport->GetNumberOfComponentsCallback());
  fltImport->SetPropagateUpdateExtentCallback( fltExport->GetPropagateUpdateExtentCallback());
  fltImport->SetUpdateDataCallback( fltExport->GetUpdateDataCallback());
  fltImport->SetDataExtentCallback( fltExport->GetDataExtentCallback());
  fltImport->SetBufferPointerCallback( fltExport->GetBufferPointerCallback());
  fltImport->SetCallbackUserData( fltExport->GetCallbackUserData());
}

void GetLabeledBoundaryMesh(LabelImageType *image, set<int> labels,
                            vnl_matrix<double> &P, vnl_vector<double> &Pl)
{
  // Read the input image
  typedef itk::OrientedRASImage<float, 3> FloatImageType;

  // For each label, threshold the image
  FloatImageType::Pointer ifloat = FloatImageType::New();
  ifloat->SetRegions(image->GetBufferedRegion());
  ifloat->CopyInformation(image);
  ifloat->Allocate();

  // Create a point array to hold all the points
  vtkPoints *allPts = vtkPoints::New();
  allPts->Initialize();

  // Create a data array
  vtkFloatArray *allLabels = vtkFloatArray::New();
  allLabels->SetName("Label");
  allLabels->Initialize();

  vtkCellArray *allCells = vtkCellArray::New();
  allCells->Initialize();

  // Iterate over the labels
  for(set<int>::iterator ilab = labels.begin(); ilab != labels.end(); ilab++)
    {
    // Set up the floating point image
    ImageRegionIterator<LabelImageType> itsrc(image, image->GetBufferedRegion());
    ImageRegionIterator<FloatImageType> ittrg(ifloat, ifloat->GetBufferedRegion());
    for(; !itsrc.IsAtEnd(); ++itsrc, ++ittrg)
      {
      ittrg.Set(itsrc.Get() == *ilab ? 1.0 : -1.0);
      }

    // Create an importer and an exporter in VTK
    typedef itk::VTKImageExport<FloatImageType> ExporterType;
    ExporterType::Pointer fltExport = ExporterType::New();
    fltExport->SetInput(ifloat);
    vtkImageImport *fltImport = vtkImageImport::New();
    ConnectITKToVTK(fltExport.GetPointer(), fltImport);

    // Run marching cubes on the input image
    vtkMarchingCubes *fltMarching = vtkMarchingCubes::New();
    fltMarching->SetInput(fltImport->GetOutput());
    fltMarching->ComputeScalarsOff();
    fltMarching->ComputeGradientsOff();
    fltMarching->ComputeNormalsOff();
    fltMarching->SetNumberOfContours(1);
    fltMarching->SetValue(0,0.0);
    fltMarching->Update();

    // Create the transform filter
    vtkTransformPolyDataFilter *fltTransform = vtkTransformPolyDataFilter::New();
    fltTransform->SetInput(fltMarching->GetOutput());

    // Compute the transform from VTK coordinates to NIFTI/RAS coordinates
    vnl_matrix_fixed<double, 4, 4> vtk2nii =
      ConstructVTKtoNiftiTransform(
        image->GetDirection().GetVnlMatrix(),
        image->GetOrigin().GetVnlVector(),
        image->GetSpacing().GetVnlVector());

    // Update the VTK transform to match
    vtkTransform *transform = vtkTransform::New();
    transform->SetMatrix(vtk2nii.data_block());
    fltTransform->SetTransform(transform);
    fltTransform->Update();

    // Get final output
    vtkPolyData *mesh = fltTransform->GetOutput();

    int poff = allPts->GetNumberOfPoints();

    for(int i = 0; i < mesh->GetNumberOfPoints(); i++)
      {
      double *p = mesh->GetPoint(i);
      allPts->InsertNextPoint(p[0],p[1],p[2]);
      allLabels->InsertNextTuple1(*ilab);
      }

    for(int i = 0; i < mesh->GetNumberOfCells(); i++)
      {
      vtkCell *cell = mesh->GetCell(i);
      for(int j = 0; j < cell->GetNumberOfPoints(); j++)
        cell->GetPointIds()->SetId(j, cell->GetPointId(j) + poff);
      allCells->InsertNextCell(cell);
      }
    }

  // Create a polydata
  vtkPolyData *outgrid = vtkPolyData::New();
  outgrid->SetPoints(allPts);
  outgrid->SetPolys(allCells);
  outgrid->GetPointData()->SetScalars(allLabels);

  // Cluster points
  vtkQuadricClustering *cluster = vtkQuadricClustering::New();
  cluster->SetNumberOfDivisions(16,16,4);
  cluster->SetInput(outgrid);
  cluster->SetUseInputPoints(1);
  cluster->Update();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName("dump.vtk");
  writer->SetInput(cluster->GetOutput());
  writer->Update();

  // Construct the output data
  vtkPolyData *result = cluster->GetOutput();
  P.set_size(result->GetNumberOfPoints(),3);
  Pl.set_size(result->GetNumberOfPoints());
  for(int i = 0; i < result->GetNumberOfPoints(); i++)
    {
    P(i,0) = result->GetPoint(i)[0];
    P(i,1) = result->GetPoint(i)[1];
    P(i,2) = result->GetPoint(i)[2];
    Pl[i] = result->GetPointData()->GetScalars()->GetTuple1(i);
    }


}

void normalize_pointset(
    vnl_matrix<double> &X,
    vnl_vector_fixed<double,3> &center)
{
  // Just compute the mean of the points
  center.fill(0.0);
  for(unsigned int i = 0; i < X.rows(); i++)
    {
    center += X.get_row(i);
    }
  center = center / (1.0 * X.rows());

  // Now subtract the mean
  for(unsigned int i = 0; i < X.rows(); i++)
    for(unsigned int j = 0; j < 3; j++)
      X[i][j] -= center[j];
}

// This is the iterative soft-assign algorithm
void softassign(
    vnl_matrix<double> X,
    vnl_matrix<double> Y,
    vnl_vector<double> Xlab,
    vnl_vector<double> Ylab,
    vnl_matrix_fixed<double,4,4> &xform)
{
  typedef vnl_vector_fixed<double, 3> Vec3;
  typedef vnl_matrix_fixed<double, 3, 3> Mat3;

  // Normalize both datasets
  Vec3 ctrX, ctrY;

  normalize_pointset(X,ctrX);
  normalize_pointset(Y,ctrY);

  // Allocate the matrix M
  int nx = X.rows(), ny = Y.rows();
  vnl_matrix<double> M(nx, ny);
  vnl_vector<double> C(nx), R(ny), w(nx);


  // Initialize the annealing parameter
  double temp = 4.0; // BE SMARTER!
  double temptarget = 0.1;
  double tempscale = 0.93;

  // Apply affine transform
  vnl_matrix<double> GX = X, V(nx,3);
  vnl_matrix_fixed<double, 4, 4> G;

  // Main annealing loop
  int iter = 0;
  while (temp > temptarget)
    {
    // Scaling factors
    double scale_exp = -1.0 / (2 * temp * temp);
    double scale_outer = 1.0 / (temp * sqrt(2 * vnl_math::pi));

    // Row and column min square distance
    vnl_vector<double> md_row(ny), md_col(nx);
    md_row.fill(1e100); md_col.fill(1e100);

    // Compute the pairwise distances between points
    for(int i = 0; i < nx; i++)
      {
      double x0 = GX(i,0), x1 = GX(i,1), x2 = GX(i,2);
      for(int j = 0; j < ny; j++)
        {
        if(Xlab(i) == Ylab(j))
          {
          double y0 = Y(j,0), y1 = Y(j,1), y2 = Y(j,2);
          double dst_sq = (x0-y0)*(x0-y0) + (x1-y1)*(x1-y1) + (x2-y2)*(x2-y2);

          double z = scale_exp * dst_sq;
          double m = (z > -16.0) ? exp(z) : 0.0;
          M(i,j) = m * scale_outer;

          if(md_col(i) > dst_sq)
            md_col(i) = dst_sq;
          if(md_row(j) > dst_sq)
            md_row(j) = dst_sq;
          }
        else
          {
          M(i,j) = 0;
          }
        }
      }

    // Compute RMS distance
    double rmsd = sqrt(0.5 * md_col.mean() + 0.5 * md_row.mean());

    // Initialize the outlier vectors C and R
    C.fill(0.01); R.fill(0.01);

    // Iteratively normalize the matrix
    for(int p = 0; p < 200; p++)
      {
      double maxerr = 0.0;
      for(int i = 0; i < nx; i++)
        {
        double row_sum = C(i);
        for(int j = 0; j < ny; j++)
          row_sum += M(i,j);
        if(row_sum < 0.0001)
          row_sum = 0.0001;
        if(row_sum > 0)
          {
          double scale = 1.0 / row_sum;
          for(int j = 0; j < ny; j++)
            M(i,j) *= scale;
          C(i) *= scale;
          }
        double del = fabs(row_sum - 1.0);
        maxerr = (del > maxerr) ? del : maxerr;
        }
      for(int j = 0; j < ny; j++)
        {
        double col_sum = R(j);
        for(int i = 0; i < nx; i++)
          col_sum += M(i,j);
        if(col_sum < 0.0001)
          col_sum = 0.0001;
        if(col_sum > 0)
          {
          double scale = 1.0 / col_sum;
          for(int i = 0; i < nx; i++)
            M(i,j) *= scale;
          R(j) *= scale;
          }
        double del = fabs(col_sum - 1.0);
        maxerr = (del > maxerr) ? del : maxerr;
        }

      if(maxerr < 1e-5)
        break;
      }

    // Compute V and w
    V = M * Y;
    for(int i = 0; i < nx; i++)
      {
      double scale = C(i) < 0.9999 ? 1.0 / (1 - C(i)) : 0.0;
      for(int k = 0; k < 3; k++)
        V(i,k) *= scale;
      }

    // The weights
    w = 1.0 - C;

    // Given the correspondence, compute the optimal transform
    vnl_matrix_fixed<double,4,4> Gxx,Gvx;
    Gxx.fill(0.0); Gvx.fill(0.0);
    double Xh[4], Vh[4];
    for(int i = 0; i < nx; i++)
      {
      Xh[0] = X(i,0); Xh[1] = X(i,1); Xh[2] = X(i,2); Xh[3] = 1.0;
      Vh[0] = V(i,0); Vh[1] = V(i,1); Vh[2] = V(i,2); Vh[3] = 1.0;
      for(int k = 0; k < 4; k++)
        {
        for(int l = 0; l < 4; l++)
          {
          Gxx(k,l) += Xh[k] * Xh[l] * w(i);
          Gvx(k,l) += Vh[k] * Xh[l] * w(i);
          }
        }
      }

    // Solve for the matrix
    G = vnl_svd<double>(Gxx.transpose()).solve(Gvx.transpose()).transpose();
    // std::cout << G << std::endl;

    // Update GX and compute the difference between GX and V
    double delta = 0.0;
    for(int i = 0; i < nx; i++)
      {
      double dsq = 0.0;
      for(int k = 0; k < 3; k++)
        {
        GX(i,k) = G(k,0) * X(i,0) + G(k,1) * X(i,1) + G(k,2) * X(i,2) + G(k,3);
        dsq += (GX(i,k) - V(i,k)) * (GX(i,k) - V(i,k));
        }
      delta += w(i) * dsq;
      }

    // Compute the fraction of outliers
    double foutl = 0.0;
    for(int i = 0; i < nx; i++)
      foutl += w(i) < 0.1 ? 1.0 : 0.0;
    foutl /= nx;

    // Compute the difference
    printf("Iter %03d \t Temp %8.2f \t Delta %8.2f \t OutFr %8.4f \t RMSD %10.4f\n",
           ++iter, temp, delta, foutl, rmsd);


    // Quit if the outlier fraction is too large. This means annealing is
    // becoming counterproductive. In reality this is a problem with the XP
    // approach, but seems like this works in practice
    if(foutl > 0.2)
      {
      cout << "Terminating due to outlier fraction too large" << endl;
      break;
      }

    temp = temp * tempscale;
    }

  // Compute the transform (account for centers)
  xform = G;
  vnl_matrix<double> A = G.extract(3,3);
  vnl_vector<double> b = G.get_column(3).extract(3);
  vnl_vector<double> bn = b + ctrY - A * ctrX;
  xform(0,3) = bn(0); xform(1,3) = bn(1); xform(2,3) = bn(2);
}

int main(int argc, char *argv[])
{
  // Scan the parameters
  if(argc < 3)
    return usage();

  Parameters p;
  p.fnTarget = argv[argc-3];
  p.fnMoving = argv[argc-2];
  p.fnOutput = argv[argc-1];

  // Read the input datasets
  typedef ImageFileReader<LabelImageType> ReaderType;
  ReaderType::Pointer readerTrg = ReaderType::New();
  readerTrg->SetFileName(p.fnTarget);
  readerTrg->Update();
  LabelImageType::Pointer target = readerTrg->GetOutput();

  ReaderType::Pointer readerMov = ReaderType::New();
  readerMov->SetFileName(p.fnMoving);
  readerMov->Update();
  LabelImageType::Pointer moving = readerMov->GetOutput();

  // Get the list of labels in each dataset
  set<int> lt, lm, lcommon;

  // Generate labeled boundary vertices for each dataset
  typedef ImageRegionIterator<LabelImageType> IterType;
  for(IterType it(target,target->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    lt.insert(it.Get());
  for(IterType im(moving,moving->GetBufferedRegion()); !im.IsAtEnd(); ++im)
    lm.insert(im.Get());

  // Get the common list of labels
  for(set<int>::iterator ic = lt.begin(); ic!=lt.end(); ic++)
    if(lm.find(*ic) != lm.end() && *ic > 0)
      {
      printf("Common label: %d\n", *ic);
      lcommon.insert(*ic);
      }

  // Generate point data for the labels
  vnl_matrix<double> X, Y;
  vnl_vector<double> lX, lY;
  GetLabeledBoundaryMesh(target, lcommon, X, lX);
  GetLabeledBoundaryMesh(moving, lcommon, Y, lY);

  // Now execute the softassign algorithm
  vnl_matrix_fixed<double,4,4> G;
  softassign(X,Y,lX,lY,G);

  // Store the matrix in the text file
  ofstream fout(p.fnOutput.c_str());
  fout << G << endl;
  fout.close();
}
