// This is a slight modification of the Umfpack solver in the Eigen library. Especially we take care of the long int type

#ifndef EIGEN_UMFPACKSUPPORT_NEW_H
#define EIGEN_UMFPACKSUPPORT_NEW_H

namespace Eigen {

/* TODO extract L, extract U, compute det, etc... */

// generic double/complex<double> wrapper functions:


inline void umfpack_defaults_new(double control[UMFPACK_CONTROL], double) { umfpack_dl_defaults(control); }
inline void umfpack_defaults_new(double control[UMFPACK_CONTROL], std::complex<double>){ umfpack_zl_defaults(control); }

inline void umfpack_report_info_new(double control[UMFPACK_CONTROL], double info[UMFPACK_INFO], double){ umfpack_dl_report_info(control, info);}
inline void umfpack_report_info_new(double control[UMFPACK_CONTROL], double info[UMFPACK_INFO], std::complex<double>){ umfpack_zl_report_info(control, info);}

inline void umfpack_report_status_new(double control[UMFPACK_CONTROL], int status, double){ umfpack_dl_report_status(control, status);}
inline void umfpack_report_status_new(double control[UMFPACK_CONTROL], int status, std::complex<double>){ umfpack_zl_report_status(control, status);}

inline void umfpack_report_control_new(double control[UMFPACK_CONTROL], double){ umfpack_dl_report_control(control);}
inline void umfpack_report_control_new(double control[UMFPACK_CONTROL], std::complex<double>){ umfpack_zl_report_control(control);}

inline void umfpack_free_numeric_new(void **Numeric, double){ umfpack_dl_free_numeric(Numeric); *Numeric = 0; }
inline void umfpack_free_numeric_new(void **Numeric, std::complex<double>){ umfpack_zl_free_numeric(Numeric); *Numeric = 0; }

inline void umfpack_free_symbolic_new(void **Symbolic, double){ umfpack_dl_free_symbolic(Symbolic); *Symbolic = 0; }
inline void umfpack_free_symbolic_new(void **Symbolic, std::complex<double>){ umfpack_zl_free_symbolic(Symbolic); *Symbolic = 0; }

inline long int umfpack_symbolic_new(long int n_row,long int n_col, const long int Ap[], const long int Ai[], const double Ax[], void **Symbolic, const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]){
  return umfpack_dl_symbolic(n_row,n_col,Ap,Ai,Ax,Symbolic,Control,Info);
}

inline long int umfpack_symbolic_new(long int n_row,long int n_col, const long int Ap[], const long int Ai[], const std::complex<double> Ax[], void **Symbolic, const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]){
  return umfpack_zl_symbolic(n_row,n_col,Ap,Ai,&numext::real_ref(Ax[0]),0,Symbolic,Control,Info);
}

inline long int umfpack_numeric_new( const long int Ap[], const long int Ai[], const double Ax[],
                            void *Symbolic, void **Numeric,
                            const double Control[UMFPACK_CONTROL],double Info [UMFPACK_INFO])
{
  return umfpack_dl_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info);
}

inline long int umfpack_numeric_new( const long int Ap[], const long int Ai[], const std::complex<double> Ax[],
                            void *Symbolic, void **Numeric,
                            const double Control[UMFPACK_CONTROL],double Info [UMFPACK_INFO])
{
  return umfpack_zl_numeric(Ap,Ai,&numext::real_ref(Ax[0]),0,Symbolic,Numeric,Control,Info);
}

inline long int umfpack_solve_new( long int sys, const long int Ap[], const long int Ai[], const double Ax[],
                          double X[], const double B[], void *Numeric,
                          const double Control[UMFPACK_CONTROL], double Info[UMFPACK_INFO])
{
  return umfpack_dl_solve(sys,Ap,Ai,Ax,X,B,Numeric,Control,Info);
}

inline long int umfpack_solve_new( long int sys, const long int Ap[], const long int Ai[], const std::complex<double> Ax[],
                          std::complex<double> X[], const std::complex<double> B[], void *Numeric,
                          const double Control[UMFPACK_CONTROL], double Info[UMFPACK_INFO])
{
  return umfpack_zl_solve(sys,Ap,Ai,&numext::real_ref(Ax[0]),0,&numext::real_ref(X[0]),0,&numext::real_ref(B[0]),0,Numeric,Control,Info);
}

inline long int umfpack_get_lunz(long int *lnz, long int *unz, long int *n_row, long int *n_col, long int *nz_udiag, void *Numeric, double)
{
  return umfpack_dl_get_lunz(lnz,unz,n_row,n_col,nz_udiag,Numeric);
}

inline long int umfpack_get_lunz(long int *lnz, long int *unz, long int *n_row, long int *n_col, long int *nz_udiag, void *Numeric, std::complex<double>)
{
  return umfpack_zl_get_lunz(lnz,unz,n_row,n_col,nz_udiag,Numeric);
}

inline long int umfpack_get_numeric(long int Lp[], long int Lj[], double Lx[], long int Up[], long int Ui[], double Ux[],
                               long int P[], long int Q[], double Dx[], long int *do_recip, double Rs[], void *Numeric)
{
  return umfpack_dl_get_numeric(Lp,Lj,Lx,Up,Ui,Ux,P,Q,Dx,do_recip,Rs,Numeric);
}

inline long int umfpack_get_numeric(long int Lp[], long int Lj[], std::complex<double> Lx[], long int Up[], long int Ui[], std::complex<double> Ux[],
                               long int P[], long int Q[], std::complex<double> Dx[], long int *do_recip, double Rs[], void *Numeric)
{
  double& lx0_real = numext::real_ref(Lx[0]);
  double& ux0_real = numext::real_ref(Ux[0]);
  double& dx0_real = numext::real_ref(Dx[0]);
  return umfpack_zl_get_numeric(Lp,Lj,Lx?&lx0_real:0,0,Up,Ui,Ux?&ux0_real:0,0,P,Q,
                                Dx?&dx0_real:0,0,do_recip,Rs,Numeric);
}

inline long int umfpack_get_determinant_new(double *Mx, double *Ex, void *NumericHandle, double User_Info [UMFPACK_INFO])
{
  return umfpack_dl_get_determinant(Mx,Ex,NumericHandle,User_Info);
}

inline long int umfpack_get_determinant_new(std::complex<double> *Mx, double *Ex, void *NumericHandle, double User_Info [UMFPACK_INFO])
{
  double& mx_real = numext::real_ref(*Mx);
  return umfpack_zl_get_determinant(&mx_real,0,Ex,NumericHandle,User_Info);
}


/** \ingroup UmfPackSupport_Module
  * \brief A sparse LU factorization and solver based on UmfPack
  *
  * This class allows to solve for A.X = B sparse linear problems via a LU factorization
  * using the UmfPack library. The sparse matrix A must be squared and full rank.
  * The vectors or matrices X and B can be either dense or sparse.
  *
  * \warning The input matrix A should be in a \b compressed and \b column-major form.
  * Otherwise an expensive copy will be made. You can call the inexpensive makeCompressed() to get a compressed matrix.
  * \tparam _MatrixType the type of the sparse matrix A, it must be a SparseMatrix<>
  *
  * \implsparsesolverconcept
  *
  * \sa \ref TutorialSparseSolverConcept, class SparseLU
  */
template<typename _MatrixType>
class UmfPackLU_NEW : public SparseSolverBase<UmfPackLU_NEW<_MatrixType> >
{
  protected:
    typedef SparseSolverBase<UmfPackLU_NEW<_MatrixType> > Base;
    using Base::m_isInitialized;
  public:
    using Base::_solve_impl;
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename MatrixType::StorageIndex StorageIndex;
    typedef Matrix<Scalar,Dynamic,1> Vector;
    typedef Matrix<long int, 1, MatrixType::ColsAtCompileTime> IntRowVectorType;
    typedef Matrix<long int, MatrixType::RowsAtCompileTime, 1> IntColVectorType;
    typedef SparseMatrix<Scalar> LUMatrixType;
    typedef SparseMatrix<Scalar,ColMajor,long int> UmfpackMatrixType;
    typedef Ref<const UmfpackMatrixType, StandardCompressedFormat> UmfpackMatrixRef;
    enum {
      ColsAtCompileTime = MatrixType::ColsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };

  public:

    typedef Array<double, UMFPACK_CONTROL, 1> UmfpackControl;
    typedef Array<double, UMFPACK_INFO, 1> UmfpackInfo;

    UmfPackLU_NEW()
      : m_dummy(0,0), mp_matrix(m_dummy)
    {
      init();
    }

    template<typename InputMatrixType>
    explicit UmfPackLU_NEW(const InputMatrixType& matrix)
      : mp_matrix(matrix)
    {
      init();
      compute(matrix);
    }

    ~UmfPackLU_NEW()
    {
      if(m_symbolic) umfpack_free_symbolic_new(&m_symbolic,Scalar());
      if(m_numeric)  umfpack_free_numeric_new(&m_numeric,Scalar());
    }

    inline Index rows() const { return mp_matrix.rows(); }
    inline Index cols() const { return mp_matrix.cols(); }

    /** \brief Reports whether previous computation was successful.
      *
      * \returns \c Success if computation was succesful,
      *          \c NumericalIssue if the matrix.appears to be negative.
      */
    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return m_info;
    }

    inline const LUMatrixType& matrixL() const
    {
      if (m_extractedDataAreDirty) extractData();
      return m_l;
    }

    inline const LUMatrixType& matrixU() const
    {
      if (m_extractedDataAreDirty) extractData();
      return m_u;
    }

    inline const IntColVectorType& permutationP() const
    {
      if (m_extractedDataAreDirty) extractData();
      return m_p;
    }

    inline const IntRowVectorType& permutationQ() const
    {
      if (m_extractedDataAreDirty) extractData();
      return m_q;
    }

    /** Computes the sparse Cholesky decomposition of \a matrix
     *  Note that the matrix should be column-major, and in compressed format for best performance.
     *  \sa SparseMatrix::makeCompressed().
     */
    template<typename InputMatrixType>
    void compute(const InputMatrixType& matrix)
    {
      if(m_symbolic) umfpack_free_symbolic_new(&m_symbolic,Scalar());
      if(m_numeric)  umfpack_free_numeric_new(&m_numeric,Scalar());
      grab(matrix.derived());
      analyzePattern_impl();
      factorize_impl();
    }

    /** Performs a symbolic decomposition on the sparcity of \a matrix.
      *
      * This function is particularly useful when solving for several problems having the same structure.
      *
      * \sa factorize(), compute()
      */
    template<typename InputMatrixType>
    void analyzePattern(const InputMatrixType& matrix)
    {
      if(m_symbolic) umfpack_free_symbolic_new(&m_symbolic,Scalar());
      if(m_numeric)  umfpack_free_numeric_new(&m_numeric,Scalar());

      grab(matrix.derived());

      analyzePattern_impl();
    }

    /** Provides the return status code returned by UmfPack during the numeric
      * factorization.
      *
      * \sa factorize(), compute()
      */
    inline int umfpackFactorizeReturncode() const
    {
      eigen_assert(m_numeric && "UmfPackLU_NEW: you must first call factorize()");
      return m_fact_errorCode;
    }

    /** Provides access to the control settings array used by UmfPack.
      *
      * If this array contains NaN's, the default values are used.
      *
      * See UMFPACK documentation for details.
      */
    inline const UmfpackControl& umfpackControl() const
    {
      return m_control;
    }

    /** Provides access to the control settings array used by UmfPack.
      *
      * If this array contains NaN's, the default values are used.
      *
      * See UMFPACK documentation for details.
      */
    inline UmfpackControl& umfpackControl()
    {
      return m_control;
    }

    /** Performs a numeric decomposition of \a matrix
      *
      * The given matrix must has the same sparcity than the matrix on which the pattern anylysis has been performed.
      *
      * \sa analyzePattern(), compute()
      */
    template<typename InputMatrixType>
    void factorize(const InputMatrixType& matrix)
    {
      eigen_assert(m_analysisIsOk && "UmfPackLU_NEW: you must first call analyzePattern()");
      if(m_numeric)
        umfpack_free_numeric_new(&m_numeric,Scalar());

      grab(matrix.derived());

      factorize_impl();
    }

    /** Prints the current UmfPack control settings.
      *
      * \sa umfpackControl()
      */
    void umfpackReportControl()
    {
      umfpack_report_control_new(m_control.data(), Scalar());
    }

    /** Prints statistics collected by UmfPack.
      *
      * \sa analyzePattern(), compute()
      */
    void umfpackReportInfo()
    {
      eigen_assert(m_analysisIsOk && "UmfPackLU_NEW: you must first call analyzePattern()");
      umfpack_report_info_new(m_control.data(), m_umfpackInfo.data(), Scalar());
    }

    /** Prints the status of the previous factorization operation performed by UmfPack (symbolic or numerical factorization).
      *
      * \sa analyzePattern(), compute()
      */
    void umfpackReportStatus() {
      eigen_assert(m_analysisIsOk && "UmfPackLU_NEW: you must first call analyzePattern()");
      umfpack_report_status_new(m_control.data(), m_fact_errorCode, Scalar());
    }

    /** \internal */
    template<typename BDerived,typename XDerived>
    bool _solve_impl(const MatrixBase<BDerived> &b, MatrixBase<XDerived> &x) const;

    Scalar determinant() const;

    void extractData() const;

  protected:

    void init()
    {
      m_info                  = InvalidInput;
      m_isInitialized         = false;
      m_numeric               = 0;
      m_symbolic              = 0;
      m_extractedDataAreDirty = true;

      umfpack_defaults_new(m_control.data(), Scalar());
    }

    void analyzePattern_impl()
    {
      m_fact_errorCode = umfpack_symbolic_new(internal::convert_index<long int>(mp_matrix.rows()),
                                          internal::convert_index<long int>(mp_matrix.cols()),
                                          mp_matrix.outerIndexPtr(), mp_matrix.innerIndexPtr(), mp_matrix.valuePtr(),
                                          &m_symbolic, m_control.data(), m_umfpackInfo.data());

      m_isInitialized = true;
      m_info = m_fact_errorCode ? InvalidInput : Success;
      m_analysisIsOk = true;
      m_factorizationIsOk = false;
      m_extractedDataAreDirty = true;
    }

    void factorize_impl()
    {

      m_fact_errorCode = umfpack_numeric_new(mp_matrix.outerIndexPtr(), mp_matrix.innerIndexPtr(), mp_matrix.valuePtr(),
                                         m_symbolic, &m_numeric, m_control.data(), m_umfpackInfo.data());

      m_info = m_fact_errorCode == UMFPACK_OK ? Success : NumericalIssue;
      m_factorizationIsOk = true;
      m_extractedDataAreDirty = true;
    }

    template<typename MatrixDerived>
    void grab(const EigenBase<MatrixDerived> &A)
    {
      mp_matrix.~UmfpackMatrixRef();
      ::new (&mp_matrix) UmfpackMatrixRef(A.derived());
    }

    void grab(const UmfpackMatrixRef &A)
    {
      if(&(A.derived()) != &mp_matrix)
      {
        mp_matrix.~UmfpackMatrixRef();
        ::new (&mp_matrix) UmfpackMatrixRef(A);
      }
    }

    // cached data to reduce reallocation, etc.
    mutable LUMatrixType m_l;
    int m_fact_errorCode;
    UmfpackControl m_control;
    mutable UmfpackInfo m_umfpackInfo;

    mutable LUMatrixType m_u;
    mutable IntColVectorType m_p;
    mutable IntRowVectorType m_q;

    UmfpackMatrixType m_dummy;
    UmfpackMatrixRef mp_matrix;

    void* m_numeric;
    void* m_symbolic;

    mutable ComputationInfo m_info;
    int m_factorizationIsOk;
    int m_analysisIsOk;
    mutable bool m_extractedDataAreDirty;

  private:
    UmfPackLU_NEW(const UmfPackLU_NEW& ) { }
};


template<typename MatrixType>
void UmfPackLU_NEW<MatrixType>::extractData() const
{
  if (m_extractedDataAreDirty)
  {
    // get size of the data
    long int lnz, unz, rows, cols, nz_udiag;
    umfpack_get_lunz(&lnz, &unz, &rows, &cols, &nz_udiag, m_numeric, Scalar());

    // allocate data
    m_l.resize(rows,(std::min)(rows,cols));
    m_l.resizeNonZeros(lnz);

    m_u.resize((std::min)(rows,cols),cols);
    m_u.resizeNonZeros(unz);

    m_p.resize(rows);
    m_q.resize(cols);

    // extract
    umfpack_get_numeric(m_l.outerIndexPtr(), m_l.innerIndexPtr(), m_l.valuePtr(),
                        m_u.outerIndexPtr(), m_u.innerIndexPtr(), m_u.valuePtr(),
                        m_p.data(), m_q.data(), 0, 0, 0, m_numeric);

    m_extractedDataAreDirty = false;
  }
}

template<typename MatrixType>
typename UmfPackLU_NEW<MatrixType>::Scalar UmfPackLU_NEW<MatrixType>::determinant() const
{
  Scalar det;
  umfpack_get_determinant_new(&det, 0, m_numeric, 0);
  return det;
}

template<typename MatrixType>
template<typename BDerived,typename XDerived>
bool UmfPackLU_NEW<MatrixType>::_solve_impl(const MatrixBase<BDerived> &b, MatrixBase<XDerived> &x) const
{
  Index rhsCols = b.cols();
  eigen_assert((BDerived::Flags&RowMajorBit)==0 && "UmfPackLU_NEW backend does not support non col-major rhs yet");
  eigen_assert((XDerived::Flags&RowMajorBit)==0 && "UmfPackLU_NEW backend does not support non col-major result yet");
  eigen_assert(b.derived().data() != x.derived().data() && " Umfpack does not support inplace solve");

  int errorCode;
  Scalar* x_ptr = 0;
  Matrix<Scalar,Dynamic,1> x_tmp;
  if(x.innerStride()!=1)
  {
    x_tmp.resize(x.rows());
    x_ptr = x_tmp.data();
  }
  for (long int j=0; j<rhsCols; ++j)
  {
    if(x.innerStride()==1) x_ptr = &x.col(j).coeffRef(0);
    errorCode = umfpack_solve_new(UMFPACK_A, mp_matrix.outerIndexPtr(), mp_matrix.innerIndexPtr(), mp_matrix.valuePtr(), x_ptr, &b.const_cast_derived().col(j).coeffRef(0), m_numeric, m_control.data(), m_umfpackInfo.data());
    if(x.innerStride()!=1) x.col(j) = x_tmp;
    if (errorCode!=0) return false;
  }

  return true;
}

} // end namespace Eigen

#endif // EIGEN_UMFPACKSUPPORT_NEW_H
