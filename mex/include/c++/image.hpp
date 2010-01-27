#ifndef IMAGE_HPP
# define IMAGE_HPP

# include <cassert>
# include <cstring>

# include <vector.hpp>

// namespace detail
// {
//   template <typename T>
//   void allocate_data(T* & data,
// 		     T** & idx_data,
// 		     int width,
// 		     int height,
// 		     int margin)
//   {
//     int width_eff = width + (margin << 1);
//     int height_eff = height + (margin << 1);

//     data = new T[width_eff * height_eff];

//     idx_data = new T*[width_eff];

//     T* ptr = data + margin;

//     for (int x = 0; x < width_eff; ++x)
//       {
// 	idx_data[x] = ptr;
// 	ptr += height_eff;
//       }

//     idx_data += margin;
//   }

//   template <typename T>
//   void desallocate_data(T* & data, T** & idx_data, int margin)
//   {
//     idx_data -= margin;
//     delete[] idx_data;
//     idx_data = 0;
//     delete[] data;
//     data = 0;
//   }
// }

// // Important note: image<T> class decribes an array of 2-dimensional
// // data which is column-wise stored in memory. Consquently, although
// // in Matlab, accessing to the point at the ith row and jth column is
// // given by I(i, j), elements in image class are accessible via
// // I(j, i) or more traditionnaly I(x, y). We choose this since matlab
// // array are also column-wise in memory and the copy between Matlab
// // arrays and images class are therefore faster to perform.

// template <typename T>
// class image
// {
// public:
//   image() :
//     width_(0), height_(0), margin_(0), data_(0), idx_data_(0)
//   {
//   }

//   image(int width, int height, int margin = 1) :
//     width_(width), height_(height), margin_(margin), data_(0), idx_data_(0)
//   {
//     detail::allocate_data(data_, idx_data_, width_, height_, margin_);
//   }

//   // src is a pointer to an array of width * height elements
//   // src needs to be column-wise stored image.
//   image(T* src, int width, int height, int margin = 1) :
//     width_(width), height_(height), margin_(margin), data_(0), idx_data_(0)
//   {
//     detail::allocate_data(data_, idx_data_, width_, height_, margin_);

//     T* ptr = src;

//     for (int x = 0; x < width_; ++x)
//       {
// 	std::memcpy(idx_data_[x], ptr, sizeof(T) * height_);
// 	ptr += height_;
//       }
//   }

//   image(const image & rhs) :
//     width_(rhs.width()), height_(rhs.height()), margin_(rhs.margin()),
//     data_(0), idx_data_(0)
//   {
//     detail::allocate_data(data_, idx_data_, width_, height_, margin_);
    
//     int width_eff = width_ + (margin_ << 1);
//     int height_eff = height_ + (margin_ << 1);

//     std::memcpy(data_, rhs.data(), sizeof(T) * width_eff * height_eff);
//   }

//   ~image()
//   {
//     detail::desallocate_data(data_, idx_data_, margin_);
//   }

//   image& operator=(const image& rhs)
//   {
//     if (&rhs == this)
//       return *this;
    
//     if (rhs.width() != width_ ||
// 	rhs.height() != height_ ||
// 	rhs.margin() != margin_)
//       {
// 	detail::desallocate_data(data_, idx_data_, margin_);
// 	width_ = rhs.width();
// 	height_ = rhs.height();
// 	margin_ = rhs.margin();
// 	detail::allocate_data(data_, idx_data_, width_, height_, margin_);
//       }

//     int width_eff = width_ + (margin_ << 1);
//     int height_eff = height_ + (margin_ << 1);

//     std::memcpy(data_, rhs.data(), sizeof(T) * width_eff * height_eff);

//     return *this;
//   }

// public:
//   int width() const { return width_; }
//   int height() const { return height_; }
//   int margin() const { return margin_; }

//   const T* data() const { return data_; }

//   template <typename U>
//   const T& operator[](const vector<2,U> & p) const
//   {
//     return (*this)(p[0], p[1]);
//   }

//   template <typename U>
//   T& operator[](const vector<2,U> & p)
//   {
//     return (*this)(p[0], p[1]);
//   }

//   const T& operator()(int x, int y) const
//   {
//     assert(contains_large(x, y));

//     return idx_data_[x][y];
//   }

//   T& operator()(int x, int y)
//   {
//     assert(contains_large(x, y));

//     return idx_data_[x][y];
//   }

//   T operator()(double x, double y) const
//   {
//     int i = (int) floor(x);
//     int j = (int) floor(y);
    
//     assert(contains_large(i, j) && contains_large(i + 1, j + 1));

//     double dx = x - i;
//     double dy = y - j;
    
//     double v00 = idx_data_[i][j];
//     double v10 = idx_data_[i + 1][j];
//     double v01 = idx_data_[i][j + 1];
//     double v11 = idx_data_[i + 1][j + 1];
    
//     return (dx * (v11 * dy - v10 * (dy - 1.0)) -
// 	    (dx - 1.0) * (v01 * dy - v00 * (dy - 1.0)));
//   }

// public:
//   bool contains(int x, int y) const
//   {
//     return (x >= 0 && x < width_ && y >= 0 && y < height_);
//   }
  
//   bool contains_large(int x, int y) const
//   {
//     return (x >= -margin_ && x < width_ + margin_ &&
// 	    y >= -margin_ && y < height_ + margin_);
//   }

// public:
//   bool operator==(const image & rhs) const
//   {
//     if (this == &rhs)
//       return true;

//     if (rhs.width() != width_ ||
// 	rhs.height() != height_ ||
// 	rhs.margin() != margin_)
//       return false;
	
//     bool same = true;

//     int width_eff = width_ + (margin_ << 1);
//     int height_eff = height_ + (margin_ << 1);

//     const T* rhs_data = rhs.data();

//     for (int i = 0; i < width_eff * height_eff; ++i)
//       same &= (data_[i] == rhs_data[i]);

//     return same;
//   }

//   bool operator!=(const image & rhs) const
//   {
//     return !(this->operator==(rhs));
//   }

// public:
//   void fill(T value)
//   {
//     int width_eff = width_ + (margin_ << 1);
//     int height_eff = height_ + (margin_ << 1);

//     std::fill_n(data_, width_eff * height_eff, value);
//   }

// public:
//   void border_replicate(int new_margin) const
//   {
//     if (new_margin != margin_)
//       const_cast<image<T>*>(this)->border_reallocate_and_copy_(new_margin);

//     const int imax = width_ - 1;
//     const int jmax = height_ - 1;

//     for (int i = - margin_; i; ++i)
//       for (int j = 0; j <= jmax; ++j)
// 	{
// 	  idx_data_[i][j] = idx_data_[0][j];
// 	  idx_data_[imax - i][j] = idx_data_[imax][j];
// 	}

//     for (int i = - margin_; i <= imax + margin_; ++i)
//       for (int j = - margin_; j; ++j)
// 	{
// 	  idx_data_[i][j] = idx_data_[i][0];
// 	  idx_data_[i][jmax - j] = idx_data_[i][jmax];
// 	}
//   }
      
//   void border_mirror(int new_margin) const
//   {
//     if (new_margin != margin_)
//       const_cast<image<T>*>(this)->border_reallocate_and_copy_(new_margin);

//     const int imax = width_ - 1;
//     const int jmax = height_ - 1;

//     for (int i = - margin_; i; ++i)
//       for (int j = 0; j <= jmax; ++j)
// 	{
// 	  idx_data_[i][j] = idx_data_[-i][j];
// 	  idx_data_[imax - i][j] = idx_data_[imax + i][j];
// 	}

//     for (int i = - margin_; i <= imax + margin_; ++i)
//       for (int j = - margin_; j; ++j)
// 	{
// 	  idx_data_[i][j] = idx_data_[i][-j];
// 	  idx_data_[i][jmax - j] = idx_data_[i][jmax + j];
// 	}
//   }

//   void border_assign(int new_margin, T val) const
//   {
//     if (new_margin != margin_)
//       const_cast<image<T>*>(this)->border_reallocate_and_copy_(new_margin);

//     const int imax = width_ - 1;
//     const int jmax = height_ - 1;

//     for (int i = - margin_; i; ++i)
//       for (int j = 0; j <= jmax; ++j)
// 	{
// 	  idx_data_[i][j] = val;
// 	  idx_data_[imax - i][j] = val;
// 	}

//     for (int i = - margin_; i <= imax + margin_; ++i)
//       for (int j = - margin_; j; ++j)
// 	{
// 	  idx_data_[i][j] = val;
// 	  idx_data_[i][jmax - j] = val;
// 	}
//   }

// public:
//   void image2mxArray(mxArray*& dst) const
//   {
//     if (dst)
//       {
// 	const mwSize* size = mxGetDimensions(dst);

// 	if (mxGetNumberOfDimensions(dst) != 2 ||
// 	    size[0] != height_ ||
// 	    size[1] != width_)
// 	  {
// 	    mxDestroyArray(dst);
// 	    dst = mxCreateDoubleMatrix(height_, width_, mxREAL);
// 	  }
//       }
//     else
//       dst = mxCreateDoubleMatrix(height_, width_, mxREAL);
    
//     double* ptr = mxGetPr(dst);
    
//     for (int x = 0; x < width_; ++x)
//       {
// 	double* offset = ptr + x * height_;
	
// 	for (int y = 0; y < height_; ++y)
// 	  offset[y] = idx_data_[x][y];
//       }
//   }

// private:
//   void border_reallocate_and_copy_(int new_margin)
//   {
//     T* data = 0;

//     T** idx_data = 0;

//     detail::allocate_data(data, idx_data, width_, height_, new_margin);

//     for (int x = 0; x < width_; ++x)
//       std::memcpy(idx_data[x], & idx_data_[x][0], sizeof(T) * height_);
      
//     detail::desallocate_data(data_, idx_data_, margin_);

//     margin_ = new_margin;
//     data_ = data;
//     idx_data_ = idx_data;
//   }

// private:
//   int width_, height_, margin_;

//   T* data_;
//   T** idx_data_;
// };


// New implementation

namespace detail
{
  template <int n, typename T>
  class image_;

  template <typename T>
  class image_<2, T>
  {
  public:
    int width() const { return size_[0]; }
    int height() const { return size_[1]; }
    int margin() const { return margin_; }

    template <typename U>
    const T& operator[](const vector<2,U> & p) const
    {
      return (*this)(p[0], p[1]);
    }

    template <typename U>
    T& operator[](const vector<2,U> & p)
    {
      return (*this)(p[0], p[1]);
    }

    const T& operator()(int x, int y) const
    {
      assert(contains_large(x, y));

      return idx_data_[x][y];
    }

    T& operator()(int x, int y)
    {
      assert(contains_large(x, y));

      return idx_data_[x][y];
    }

    double operator()(double x, double y) const
    {
      int i = (int) floor(x);
      int j = (int) floor(y);
    
      assert(contains_large(i, j) && contains_large(i + 1, j + 1));

      double dx = x - i;
      double dy = y - j;
    
      double v00 = idx_data_[i][j];
      double v10 = idx_data_[i + 1][j];
      double v01 = idx_data_[i][j + 1];
      double v11 = idx_data_[i + 1][j + 1];
    
      return (dx * (v11 * dy - v10 * (dy - 1.0)) -
	      (dx - 1.0) * (v01 * dy - v00 * (dy - 1.0)));
    }

    bool contains(int x, int y) const
    {
      return (x >= 0 && x < size_[0] && y >= 0 && y < size_[1]);
    }
  
    bool contains_large(int x, int y) const
    {
      return (x >= -margin_ && x < size_[0] + margin_ &&
	      y >= -margin_ && y < size_[1] + margin_);
    }
    
    // Import export functions:
    // src is a pointer to an array of width * height elements
    // src needs to be column-wise stored image.
    void fill(const T * src)
    {
      const T* ptr = src;

      for (int x = 0; x < size_[0]; ++x)
	{
	  std::memcpy(idx_data_[x], ptr, sizeof(T) * size_[1]);
	  ptr += size_[1];
	}
    }

    template <typename U>
    void raw_data(U* dst) const
    {
      U* ptr = dst;
      
      for (int x = 0; x < size_[0]; ++x)
	{
	  const T* ptr2 = idx_data_[x];

	  for (int y = 0; y < size_[1]; ++y)
	    ptr[y] = ptr2[y];
	  ptr += size_[1];
	}
    }

    void border_replicate(int new_margin) const
    {
      if (new_margin != margin_)
	const_cast< image_* >(this)->border_reallocate_and_copy_(new_margin);

      const int imax = size_[0] - 1;
      const int jmax = size_[1] - 1;

      for (int i = - margin_; i; ++i)
	for (int j = 0; j <= jmax; ++j)
	  {
	    idx_data_[i][j] = idx_data_[0][j];
	    idx_data_[imax - i][j] = idx_data_[imax][j];
	  }

      for (int i = - margin_; i <= imax + margin_; ++i)
	for (int j = - margin_; j; ++j)
	  {
	    idx_data_[i][j] = idx_data_[i][0];
	    idx_data_[i][jmax - j] = idx_data_[i][jmax];
	  }
    }
    
    void border_mirror(int new_margin) const
    {
      if (new_margin != margin_)
	const_cast< image_* >(this)->border_reallocate_and_copy_(new_margin);
      
      const int imax = size_[0] - 1;
      const int jmax = size_[1] - 1;
      
      for (int i = - margin_; i; ++i)
	for (int j = 0; j <= jmax; ++j)
	  {
	    idx_data_[i][j] = idx_data_[-i][j];
	    idx_data_[imax - i][j] = idx_data_[imax + i][j];
	  }

      for (int i = - margin_; i <= imax + margin_; ++i)
	for (int j = - margin_; j; ++j)
	  {
	    idx_data_[i][j] = idx_data_[i][-j];
	    idx_data_[i][jmax - j] = idx_data_[i][jmax + j];
	  }
    }

    void border_assign(int new_margin, T val) const
    {
      if (new_margin != margin_)
	const_cast< image_* >(this)->border_reallocate_and_copy_(new_margin);
      
      const int imax = size_[0] - 1;
      const int jmax = size_[1] - 1;

      for (int i = - margin_; i; ++i)
	for (int j = 0; j <= jmax; ++j)
	  {
	    idx_data_[i][j] = val;
	    idx_data_[imax - i][j] = val;
	  }

      for (int i = - margin_; i <= imax + margin_; ++i)
	for (int j = - margin_; j; ++j)
	  {
	    idx_data_[i][j] = val;
	    idx_data_[i][jmax - j] = val;
	  }
    }

  private:
    void border_reallocate_and_copy_(int new_margin)
    {
      T* new_data = 0;
      T** new_idx_data = 0;

      allocate_(new_data, new_idx_data, size_, new_margin);
      
      for (int x = 0; x < size_[0]; ++x)
	std::memcpy(new_idx_data[x], & idx_data_[x][0], sizeof(T) * size_[1]);
      
      desallocate_(data_, idx_data_, size_, margin_);

      margin_ = new_margin;
      data_ = new_data;
      idx_data_ = new_idx_data;
    }

  protected:
    image_() : margin_(0), data_(0), idx_data_(0) {}

    static void allocate_(T* & data,
			  T** & idx_data,
			  int size[2],
			  int margin)
    {
      int width_eff = size[0] + (margin << 1);
      int height_eff = size[1] + (margin << 1);
      
      data = new T[width_eff * height_eff];
      
      idx_data = new T*[width_eff];
      
      T* ptr = data + margin;
      
      for (int x = 0; x < width_eff; ++x)
	{
	  idx_data[x] = ptr;
	  ptr += height_eff;
	}
      
      idx_data += margin;
    }

    static void desallocate_(T* & data,
			     T** & idx_data,
			     int size[2],
			     int margin)
    {
      delete[] data;
      data = 0;

      idx_data -= margin;
      delete[] idx_data;
      idx_data = 0;
    }

  protected:
    int size_[2];
    int margin_;

    T* data_;
    T** idx_data_;
  };

  template <typename T>
  class image_<3, T>
  {
  public:
    int width() const { return size_[0]; }
    int height() const { return size_[1]; }
    int depth() const { return size_[2]; }
    int margin() const { return margin_; }

    template <typename U>
    const T& operator[](const vector<3,U> & p) const
    {
      return (*this)(p[0], p[1], p[2]);
    }

    template <typename U>
    T& operator[](const vector<3,U> & p)
    {
      return (*this)(p[0], p[1], p[2]);
    }

    const T& operator()(int x, int y, int z) const
    {
      assert(contains_large(x, y, z));

      return idx_data_[z][x][y];
    }

    T& operator()(int x, int y, int z)
    {
      assert(contains_large(x, y, z));

      return idx_data_[z][x][y];
    }

    double operator()(double x, double y, double z) const
    {
      int i = (int) floor(x);
      int j = (int) floor(y);
      int k = (int) floor(z);

      assert(contains_large(i, j, k) &&
	     contains_large(i + 1, j + 1, k + 1));

      double dx = x - i;
      double dy = y - j;
      double dz = z - k;

      double i1 = idx_data_[k][i][j] * (1 - dz) +
	idx_data_[k + 1][i][j] * dz;

      double i2 = idx_data_[k][i][j + 1] * (1 - dz) +
	idx_data_[k + 1][i][j + 1] * dz;

      double j1 = idx_data_[k][i + 1][j] * (1 - dz) +
	idx_data_[k + 1][i + 1][j] * dz;

      double j2 = idx_data_[k][i + 1][j + 1] * (1 - dz) +
	idx_data_[k + 1][i + 1][j + 1] * dz;
      
      double w1 = i1 * (1 - dy) + i2 * dy;
      double w2 = j1 * (1 - dy) + j2 * dy;

      return w1 * (1 - dx) + w2 * dx;
    }

    bool contains(int x, int y, int z) const
    {
      return (x >= 0 && x < size_[0] &&
	      y >= 0 && y < size_[1] &&
	      z >= 0 && z < size_[2]);
    }
  
    bool contains_large(int x, int y, int z) const
    {
      return (x >= -margin_ && x < size_[0] + margin_ &&
	      y >= -margin_ && y < size_[1] + margin_ &&
	      z >= -margin_ && z < size_[2] + margin_);
    }

    // Import export functions: 
    // src is a pointer to an array of width * height * depth
    // elements. src needs to be column-wise sorted image.
    void fill(const T * src)
    {
      const T* ptr = src;

      for (int z = 0; z < size_[2]; ++z)
	{
	  T** ptr2 = idx_data_[z];

	  for (int x = 0; x < size_[0]; ++x)
	    {
	      std::memcpy(ptr2[x], ptr, sizeof(T) * size_[1]);
	      ptr += size_[1];
	    }
	}
    }
      
    template <typename U>
    void raw_data(U* dst) const
    {
      U* ptr = dst;

      for (int z = 0; z < size_[2]; ++z)
	{
	  T** ptr2 = idx_data_[z];

	  for (int x = 0; x < size_[0]; ++x)
	    {
	      const T* ptr3 = ptr2[x];

	      for (int y = 0; y < size_[1]; ++y)
		ptr[y] = ptr3[y];
	      ptr += size_[1];
	    }
	}
    }

    static void allocate_(T* & data,
			  T*** & idx_data,
			  int size[3],
			  int margin)
    {
      int width_eff = size[0] + (margin << 1);
      int height_eff = size[1] + (margin << 1);
      int depth_eff = size[2] + (margin << 1);

      data = new T[width_eff * height_eff * depth_eff];

      idx_data = new T**[depth_eff];

      T* ptr = data + margin;

      for (int z = 0; z < depth_eff; ++z)
	{
	  T** tmp = new T*[width_eff];
	  idx_data[z] = tmp;
	  for (int x = 0; x < width_eff; ++x)
	    {
	      idx_data[z][x] = ptr;
	      ptr += height_eff;
	    }
	  idx_data[z] += margin;
	}
      idx_data += margin;
    }

    static void desallocate_(T* & data,
			     T*** & idx_data,
			     int size[3],
			     int margin)
    {
      delete[] data;
      data = 0;

      idx_data -= margin;

      int depth_eff = size[2] + (margin << 1);

      for (int z = 0; z < depth_eff; ++z)
	{
	  idx_data[z] -= margin;
	  delete[] idx_data[z];
	}

      idx_data = 0;
    }

  protected:
    int size_[3];
    int margin_;

    T* data_;
    T*** idx_data_; // order: z,x,y
  };
}

template <int n, typename T>
class image : public detail::image_<n, T>
{
public:
  typedef typename detail::image_<n, T> super_type;

protected:  
  using super_type::size_;
  using super_type::margin_;
  using super_type::data_;
  using super_type::idx_data_;

public:
  image()
  {
    for (int i = 0; i < n; ++i)
      size_[i] = 0;
  }

  image(const int size[n], int margin = 1)
  {
    margin_ = margin;

    for (int i = 0; i < n; ++i)
      size_[i] = size[i];

    super_type::allocate_(data_, idx_data_, size_, margin_);
  }

  image(const image & rhs)
  {
    margin_ = rhs.margin();

    const int* rhs_size = rhs.size();

    for (int i = 0; i < n; ++i)
      size_[i] = rhs_size[i];

    super_type::allocate_(data_, idx_data_, size_, margin_);
 
    unsigned long int s = 1;

    for (int i = 0; i < n; ++i)
      s *= (size_[i] + (margin_ << 1));

    std::memcpy(data_, rhs.data(), sizeof(T) * s);
  }

  ~image()
  {
    super_type::desallocate_(data_, idx_data_, size_, margin_);
  }

  image& operator=(const image& rhs)
  {
    if (&rhs == this)
      return *this;
   
    bool same = true;

    const int* rhs_size = rhs.size();

    for (unsigned i = 0; i < n; ++i)
      same &= size_[i] == rhs_size[i];
 
    if (!same || rhs.margin() != margin_)
      {
	super_type::desallocate_(data_, idx_data_, size_, margin_);
	for (int i = 0; i < n; ++i)
	  size_[i] = rhs.size(i);
	margin_ = rhs.margin();
	super_type::allocate_(data_, idx_data_, size_, margin_);
      }

    unsigned long int s = 1;

    for (int i = 0; i < n; ++i)
      s *= (size_[i] + (margin_ << 1));

    std::memcpy(data_, rhs.data(), sizeof(T) * s);

    return *this;
  }

  const T* data() const { return data_; }
  const int*  size() const { return size_; }

  bool operator==(const image & rhs) const
  {
    if (this == &rhs)
      return true;

    bool same = true;

    for (unsigned i = 0; i < n; ++i)
      same &= size_[i] == rhs.size(i);
 
    if (!same || rhs.margin() != margin_)
      return false;
	
    unsigned long int s = 1;

    for (int i = 0; i < n; ++i)
      s *= (size_[i] + (margin_ << 1));

    const T* rhs_data = rhs.data();

    same = true;

    for (unsigned long int i = 0; i < s; ++i)
      same &= (data_[i] == rhs_data[i]);

    return same;
  }

  bool operator!=(const image & rhs) const
  {
    return !(this->operator==(rhs));
  }

  void fill(T value)
  {
    unsigned long int s = 1;

    for (int i = 0; i < n; ++i)
      s *= (size_[i] + (margin_ << 1));
    
    std::fill_n(data_, s, value);
  }

  void fill(const T* src)
  {
    static_cast<super_type*>(this)->fill(src);
  }
};

#endif /* IMAGE_HPP */
