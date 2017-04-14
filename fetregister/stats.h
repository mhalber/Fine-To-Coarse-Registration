/* 
  A simple set of functions for common statistical operations related to
  arrays. Implemented on top of std::vector.
 */

#ifndef STATS_H
#define STATS_H

template<typename T> T stats_mean( const std::vector<T> &vals );
template<typename T> T stats_sdev( const std::vector<T> & vals );
template<typename T> T stats_sdev( const std::vector<T> & vals, T mean );
template<typename T> std::vector<T> stats_pdf( std::vector<T> density );
template<typename T> std::vector<T> stats_cdf( const std::vector<T> &pdf );
template<typename T> std::vector<T> stats_icdf( const std::vector<T> & cdf );
template<typename T> void stats_print_vector( const std::vector<T> &v );

#endif /* STATS_H */

#ifdef STATS_IMPLEMENTATION

template<typename T>
T stats_mean( const std::vector<T> &vals )
{
  return std::accumulate( vals.begin(), vals.end(), 0.0 ) / (double)vals.size();
}

template<typename T>
T stats_sdev( const std::vector<T> & vals )
{
  T mean = stats_mean( vals );
  T sq_sum = std::inner_product( vals.begin(), vals.end(), vals.begin(), 0.0);
  return std::sqrt( sq_sum / (double)vals.size() - mean * mean );
}

template<typename T>
T stats_sdev( const std::vector<T> & vals, T mean )
{
  T sq_sum = std::inner_product( vals.begin(), vals.end(), vals.begin(), 0.0);
  return std::sqrt( sq_sum / (double)vals.size() - mean * mean );
}

template<typename T> 
std::vector<T> stats_pdf( std::vector<T> density )
{
  T sum = std::accumulate( density.begin(), density.end(), 0.0 );
  std::for_each( density.begin(), density.end(), [sum](T &n){ n /= sum; } );
  return density;
}

template<typename T> 
std::vector<T> stats_cdf( const std::vector<T> & pdf )
{
  std::vector<T> cdf( pdf.size(), 0 );
  std::partial_sum( pdf.begin(), pdf.end(), cdf.begin() );
  return cdf;
}

/*note: given random number n in (0,1) you should convert it 
        to the range in your distribution (0,n) */
template<typename T> 
std::vector<T> stats_icdf( const std::vector<T> & cdf )
{
  size_t size = cdf.size();
  std::vector<T> icdf( size, -1 );
  
  /* flip the existing values / probabilities */
  for ( size_t i = 0 ; i < size ; ++i )
  {
    T prob = cdf[i];
    size_t idx = round( prob * (size - 1) );
    icdf[ idx ] = i;
  }

  /* interpolate missing values */
  int first = 0, last = size - 1;
  for ( size_t i = 0 ; i < size ; ++i )
  {
    int prev = -1, next = -1;
    int prev_idx = i, next_idx = i;
    if ( icdf[i] == -1 )
    {
      while ( prev_idx >= 0 && prev == -1 )
      {
        prev = icdf[--prev_idx];
      }
      if ( prev == -1 ) { prev = first; prev_idx = 0; }

      while ( next_idx < size && next == -1 )
      {
        next = icdf[++next_idx];
      }
      if ( next == -1 ) { next = last; next_idx = 0; }

      float distance = next_idx - prev_idx;
      while( i != next_idx )
      {
        float val = ((next_idx - i)/distance) * prev + 
                    ((i - prev_idx)/distance) * next;
        icdf[i++] = round(val);
      }
    }
  }
  return icdf;
}

template<typename T> void stats_print_vector( const std::vector<T> &v )
{
  std::cout << "[";
  for( int i = 0 ; i < v.size() ; ++i )
  {
    std::cout << std::setfill(' ') << std::setw(3) << v[i] << " ";
  } 
  std::cout << std::setw(0) << "]" << std::endl;
  
}

#endif /* STATS_IMPLEMENTATION */


/* For a c library, at some point 
float stats_accum( float *start, float *end, float init );
float stats_inner_product( float *start1, float *end1, 
                           float *start2, float init);
float stats_mean( float *vals, size_t vals );
float stats_sdev( float *vals, size_t vals, float mean  );

float stats_cumsum(  );
*/