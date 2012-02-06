#ifndef __itkOrientedRASImage_h_
#define __itkOrientedRASImage_h_

#include "itkOrientedImage.h"

namespace itk {

/** 
 * Oriented image with RAS physical coordinates (as opposed to LPS)
 */
template <class TPixel, unsigned int VImageDimension>
class ITK_EXPORT OrientedRASImage : public OrientedImage<TPixel, VImageDimension>
{
public:
  /** Standard class typedefs */
  typedef OrientedRASImage               Self;
  typedef OrientedImage<TPixel, VImageDimension>  Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  typedef WeakPointer<const Self>  ConstWeakPointer;
  typedef Matrix<double, VImageDimension+1, VImageDimension+1> TransformMatrixType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(OrientedRASImage, OrientedImage);

  /** Index typedef support. An index is used to access pixel values. */
  typedef typename Superclass::IndexType  IndexType;

  /** Direction typedef support. The direction cosines of the image. */
  typedef typename Superclass::DirectionType  DirectionType;

  /** Spacing typedef support.  Spacing holds the size of a pixel.  The
   * spacing is the geometric distance between image samples. */
  typedef typename Superclass::SpacingType SpacingType;

  typedef typename Superclass::AccessorType        AccessorType;
  typedef typename Superclass::AccessorFunctorType AccessorFunctorType;
  typedef typename Superclass::IOPixelType         IOPixelType;

  /** Tyepdef for the functor used to access a neighborhood of pixel pointers.*/
  typedef NeighborhoodAccessorFunctor< Self > 
                                            NeighborhoodAccessorFunctorType;

  /** Return the NeighborhoodAccessor functor. This method is called by the 
   * neighborhood iterators. */
  NeighborhoodAccessorFunctorType GetNeighborhoodAccessor() 
    { return NeighborhoodAccessorFunctorType(); }
  
  /** Return the NeighborhoodAccessor functor. This method is called by the 
   * neighborhood iterators. */
  const NeighborhoodAccessorFunctorType GetNeighborhoodAccessor() const
    { return NeighborhoodAccessorFunctorType(); }
  

  /** \brief Get the continuous index from a physical point
   *
   * Returns true if the resulting index is within the image, false otherwise.
   * \sa Transform */
  template<class TCoordRep>
  bool TransformRASPhysicalPointToContinuousIndex(
              const Point<TCoordRep, VImageDimension>& point,
              ContinuousIndex<TCoordRep, VImageDimension>& index   ) const
    {
    Point<TCoordRep, VImageDimension> p_lps = point;
    p_lps[0] = -point[0]; p_lps[1] = -point[1];
    return Superclass::TransformPhysicalPointToContinuousIndex(p_lps, index);
    }

  /** Get the index (discrete) from a physical point.
   * Floating point index results are truncated to integers.
   * Returns true if the resulting index is within the image, false otherwise
   * \sa Transform */
  template<class TCoordRep>
  bool TransformRASPhysicalPointToIndex(
            const Point<TCoordRep, VImageDimension>& point,
            IndexType & index                                ) const
    {
    Point<TCoordRep, VImageDimension> p_lps = point;
    p_lps[0] = -point[0]; p_lps[1] = -point[1];
    return Superclass::TransformPhysicalPointToIndex(p_lps, index);
    }

  /** Get a physical point (in the space which
   * the origin and spacing infomation comes from)
   * from a continuous index (in the index space)
   * \sa Transform */
  template<class TCoordRep>
  void TransformContinuousIndexToRASPhysicalPoint(
            const ContinuousIndex<TCoordRep, VImageDimension>& index,
            Point<TCoordRep, VImageDimension>& point        ) const
    {
    Superclass::TransformContinuousIndexToPhysicalPoint(index, point);
    point[0] = -point[0];
    point[1] = -point[1];
    }

  /** Get a physical point (in the space which
   * the origin and spacing infomation comes from)
   * from a discrete index (in the index space)
   *
   * \sa Transform */
  template<class TCoordRep>
  void TransformIndexToRASPhysicalPoint(
                      const IndexType & index,
                      Point<TCoordRep, VImageDimension>& point ) const
    {
    Superclass::TransformIndexToPhysicalPoint(index, point);
    point[0] = -point[0];
    point[1] = -point[1];
    }

  /** Take a vector or covariant vector that has been computed in the
   * coordinate system parallel to the image grid and rotate it by the
   * direction cosines in order to get it in terms of the coordinate system of
   * the image acquisition device.  This implementation in the OrientedImage
   * multiply the array (vector or covariant vector) by the matrix of Direction
   * Cosines. The arguments of the method are of type FixedArray to make
   * possible to use this method with both Vector and CovariantVector.
   * The Method is implemented differently in the itk::Image.
   *
   * \sa Image
   */ 
  template<class TCoordRep>
  void TransformLocalVectorToRASPhysicalVector(
    const FixedArray<TCoordRep, VImageDimension> & inputGradient,
          FixedArray<TCoordRep, VImageDimension> & outputGradient ) const
    {
    Superclass::TransformLocalVectorToPhysicalVector(inputGradient, outputGradient);
    outputGradient[0] = -outputGradient[0];
    outputGradient[1] = -outputGradient[1];
    }

  /** 
   * Get a matrix that maps points voxel coordinates to RAS coordinates
   */
  TransformMatrixType GetVoxelSpaceToRASPhysicalSpaceMatrix()
    {
    TransformMatrixType mat;
    mat.SetIdentity();

    for(size_t i = 0; i < VImageDimension; i++)
      {
      double ras_flip = (i < 2) ? -1 : 1;
      for(size_t j = 0; j < VImageDimension; j++)
        {
        mat[i][j] = ras_flip * this->GetDirection()(i,j) * this->GetSpacing()[i];
        }
      mat[i][VImageDimension] = ras_flip * this->GetOrigin()[i];
      }

    return mat;
    };

  /** 
   * Get a matrix that maps points in the x * spacing + origin space to
   * the RAS space
   */
  TransformMatrixType GetSpacingOriginPhysicalSpaceToRASPhysicalSpaceMatrix()
    {
    TransformMatrixType mat;
    mat.SetIdentity();

    for(size_t i = 0; i < VImageDimension; i++)
      {
      double ras_flip = (i < 2) ? -1 : 1;
      mat[i][VImageDimension] = ras_flip * this->GetOrigin()[i];
      for(size_t j = 0; j < VImageDimension; j++)
        {
        mat[i][j] = ras_flip * this->GetDirection()(i,j);
        mat[i][VImageDimension] -= ras_flip * this->GetDirection()(i,j) * this->GetOrigin()[i];
        }      
      }

    return mat;
    }


protected:
  OrientedRASImage() {};
  virtual ~OrientedRASImage() {};

private:
  OrientedRASImage(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} //namespace itk

#endif

