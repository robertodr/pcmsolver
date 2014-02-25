#ifndef MATHUTILS_HPP
#define MATHUTILS_HPP

#include <iostream>
//#include <iosfwd>
#include <stdexcept>
#include <vector>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

/*! \file MathUtils.hpp
 *  \brief Some math utilities.
 *  \author Roberto Di Remigio
 *  \date 2014
 */

/*! \brief Returns an Eigen matrix of type T, with dimensions _rows*_columns.
 *  \param _rows the number of rows.
 *  \param _columns the number of columns.
 *  \param _inData the raw data buffer.
 *  \tparam T the data type of the matrix to be returned.
 *
 *  Warning! This function assumes that the raw buffer is in column-major order
 *  as in Fortran.
 */

template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getFromRawBuffer(size_t _rows, size_t _columns, void * _inData)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _outData;
	_outData.resize(_rows, _columns);
	T * data = reinterpret_cast<T*>(_inData);
	Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > mapped_data(data, _rows, _columns);
 	_outData = mapped_data;
	return _outData;
};

/*! \brief Packs a block-diagonal matrix.
 *  
 * 
 *  We currently assume that all the blocks are square and have the same dimensionality.
 *  This data type is to be used in conjuction with symmetry handling in the solver.
 */

template <typename T, int nrBlocks_, int dimBlock_>
class BlockDiagonalMatrix
{
	private:
		typedef Eigen::Matrix<T, dimBlock_, dimBlock_, Eigen::ColMajor> EigenBlock;
		typedef Eigen::Matrix<T, (nrBlocks_*dimBlock_), (nrBlocks_*dimBlock_), Eigen::ColMajor> EigenFull;
		typedef std::vector<EigenBlock> blockVector;
	private:
		blockVector blockedMatrix_;
		EigenFull fullMatrix_;
		int fullDim_;
		std::ostream & printObject(std::ostream & os)
		{
                	os << "Matrix has " << nrBlocks_ << " square blocks of dimension " << dimBlock_ << std::endl;
                	for (int i = 0; i < nrBlocks_; ++i)
                	{
                		os << "Block number " << i << std::endl;
                		os << blockedMatrix_[i] << std::endl;
                	}
                	return os;
                }
	public:
		BlockDiagonalMatrix() {}
		/*!
		 *  Packs a square block diagonal matrix given the full matrix.
		 */
		BlockDiagonalMatrix(const EigenFull & matrix) : fullMatrix_(matrix), fullDim_(nrBlocks_*dimBlock_)
		{	
			// The dimension of the full matrix is (nrBlocks * dimBlock)
			// The full matrix is assumed to be square, with square blocks on the diagonal
			// all with the same dimension.
			int j = 0;
			for (int i = 0; i < nrBlocks_; ++i)
			{
				blockedMatrix_.push_back(fullMatrix_.block(j, j, dimBlock_, dimBlock_));
				j += dimBlock_;
			}
		}
		/// Copy constructor
		BlockDiagonalMatrix(const BlockDiagonalMatrix & other) : blockedMatrix_(other.blockedMatrix_), fullMatrix_(other.fullMatrix_), fullDim_(other.fullDim_) {}
		friend inline void swap(BlockDiagonalMatrix & left, BlockDiagonalMatrix & right);
		inline void swap(BlockDiagonalMatrix & other);
		/// Addition-assignment operator
		BlockDiagonalMatrix & operator+=(const BlockDiagonalMatrix & other);
		/// Subtraction-assignment operator
		BlockDiagonalMatrix & operator-=(const BlockDiagonalMatrix & other);
		~BlockDiagonalMatrix() {}
         	friend std::ostream & operator<<(std::ostream & os, BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & bd)
	 	{
			return bd.printObject(os);
         	}

		int nrBlocks() const { return nrBlocks_; }
		int dimBlock() const { return dimBlock_; }
		int fullDim() const { return fullDim_; }

		const EigenFull & fullMatrix() const { return fullMatrix_; }
		void fullMatrix(const EigenFull & matrix) const { matrix = fullMatrix_; }

		const EigenBlock & block(int i) const { return blockedMatrix_[i]; }
		void block(const EigenBlock & matrix, int i) const { matrix = blockedMatrix_[i]; }
};

template <typename T, int nrBlocks_, int dimBlock_>
BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & BlockDiagonalMatrix<T, nrBlocks_, dimBlock_>::operator+=(const BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & other)
{
	for (int i = 0; i < nrBlocks_; ++i)
	{
		blockedMatrix_[i] += other.blockedMatrix_[i];
	}
	return *this;
}

template <typename T, int nrBlocks_, int dimBlock_>
BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & BlockDiagonalMatrix<T, nrBlocks_, dimBlock_>::operator-=(const BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & other)
{
	for (int i = 0; i < nrBlocks_; ++i)
	{
		blockedMatrix_[i] -= other.blockedMatrix_[i];
	}
	return *this;
}

/*!
 * \fn inline BlockDiagonalMatrix operator+(BlockDiagonalMatrix left, const BlockDiagonalMatrix & right)
 * \brief Addition operator
 * \param left the left hand side of the addition
 * \param right the right hand side of the addition
 */
template <typename T, int nrBlocks_, int dimBlock_>
inline BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> operator+(BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> left, const BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & right)
{
	left += right;
	return left;
};

/*!
 * \fn inline BlockDiagonalMatrix operator-(BlockDiagonalMatrix left, const BlockDiagonalMatrix & right)
 * \brief Subtraction operator
 * \param left the left hand side of the subtraction
 * \param right the right hand side of the subtraction
 */
template <typename T, int nrBlocks_, int dimBlock_>
inline BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> operator-(BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> left, const BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & right)
{
	left -= right;
	return left;
};

/*!
 * \fn inline BlockDiagonalMatrix operator*(BlockDiagonalMatrix left, const BlockDiagonalMatrix & right)
 * \brief Multiplication operator: uniform scaling of BlockDiagonalMatrix version
 * \param scaling the scaling factor
 * \param object the surface function to be scaled
 */
template <typename T, int nrBlocks_, int dimBlock_>
inline BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> operator*(double scaling, BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & object)
{
	object *= scaling;
	return object;
};

#endif // MATHUTILS_HPP
