

#ifndef SIMDCompressionAndIntersection_INTERSECTION_H_
#define SIMDCompressionAndIntersection_INTERSECTION_H_

#include <common.h>

#ifndef CONST_IF
#  ifdef __cpp_if_constexpr
#    define CONST_IF(...) if constexpr(#__VA_ARGS__)
#  else
#    define CONST_IF(...) if(#__VA_ARGS__)
#  endif
#endif


namespace SIMDCompressionLib {

/**
 * This is often called galloping or exponential search.
 *
 * Used by frogintersectioncardinality below
 *
 * Based on binary search...
 * Find the smallest integer larger than pos such
 * that array[pos]>= min.
 * If none can be found, return array.length.
 * From code by O. Kaser.
 */
template<typename T>
size_t __frogadvanceUntil(const T *array, const size_t pos,
                          const size_t length, const size_t min) {
  size_t lower = pos + 1;

  // special handling for a possibly common sequential case
  if ((lower >= length) or (array[lower] >= min)) {
    return lower;
  }

  size_t spansize = 1; // could set larger
  // bootstrap an upper limit

  while ((lower + spansize < length) and (array[lower + spansize] < min))
    spansize *= 2;
  size_t upper = (lower + spansize < length) ? lower + spansize : length - 1;

  if (array[upper] < min) { // means array has no item >= min
    return length;
  }

  // we know that the next-smallest span was too small
  lower += (spansize / 2);

  // else begin binary search
  size_t mid = 0;
  while (lower + 1 != upper) {
    mid = (lower + upper) / 2;
    if (array[mid] == min) {
      return mid;
    } else if (array[mid] < min)
      lower = mid;
    else
      upper = mid;
  }
  return upper;
}

template<typename T>
size_t onesidedgallopingintersection(const T *smallset,
                                     const size_t smalllength,
                                     const T *largeset,
                                     const size_t largelength, T *out) {
  if (largelength < smalllength)
    return onesidedgallopingintersection(largeset, largelength, smallset,
                                         smalllength, out);
  if (0 == smalllength)
    return 0;
  const T *const initout(out);
  size_t k1 = 0, k2 = 0;
  while (true) {
    if (largeset[k1] < smallset[k2]) {
      k1 = __frogadvanceUntil(largeset, k1, largelength, smallset[k2]);
      if (k1 == largelength)
        break;
    }
  midpoint:
    if (smallset[k2] < largeset[k1]) {
      ++k2;
      if (k2 == smalllength)
        break;
    } else {
      *out++ = smallset[k2];
      ++k2;
      if (k2 == smalllength)
        break;
      k1 = __frogadvanceUntil(largeset, k1, largelength, smallset[k2]);
      if (k1 == largelength)
        break;
      goto midpoint;
    }
  }
  return out - initout;
}

/**
 * Fast scalar scheme designed by N. Kurz.
 */
template<typename T, bool emit_output=true>
size_t scalar(const T *A, const size_t lenA, const T *B,
              const size_t lenB, T *out=nullptr) {
  const T *const initout(out);
  if (lenA == 0 || lenB == 0)
    return 0;
  size_t n = 0;

  const T *endA = A + lenA;
  const T *endB = B + lenB;

  for(;;) {
    while (*A < *B) {
    SKIP_FIRST_COMPARE:
      if (++A == endA)
        return (out - initout);
    }
    while (*A > *B) {
      if (++B == endB)
        return (out - initout);
    }
    if (*A == *B) {
      CONST_IF(emit_output) *out++ = *A;
      else                          ++n;
      if (++A == endA || ++B == endB) {
        CONST_IF(emit_output) return out - initout;
        else                              return n;
      }
    } else {
      goto SKIP_FIRST_COMPARE;
    }
  }

#ifdef __GNUC__ /* clang defines __GNUC__ */
  __builtin_unreachable();
#endif
  return (out - initout); // NOTREACHED
}

template<typename T>
size_t match_scalar(const T *A, const size_t lenA, const T *B,
                    const size_t lenB, T *out) {

  const T *initout = out;
  if (lenA == 0 || lenB == 0)
    return 0;

  const T *endA = A + lenA;
  const T *endB = B + lenB;

  while (1) {
    while (*A < *B) {
    SKIP_FIRST_COMPARE:
      if (++A == endA)
        goto FINISH;
    }
    while (*A > *B) {
      if (++B == endB)
        goto FINISH;
    }
    if (*A == *B) {
      *out++ = *A;
      if (++A == endA || ++B == endB)
        goto FINISH;
    } else {
      goto SKIP_FIRST_COMPARE;
    }
  }

FINISH:
  return (out - initout);
}
using namespace std;
/*
 * Given two arrays, this writes the intersection to out. Returns the
 * cardinality of the intersection.
 */
typedef size_t (*intersectionfunction)(const uint32_t *set1,
                                       const size_t length1,
                                       const uint32_t *set2,
                                       const size_t length2, uint32_t *out);

/*
 * Given two arrays, this writes the intersection to out. Returns the
 * cardinality of the intersection.
 *
 * This is a mix of very fast vectorized intersection algorithms, several
 * designed by N. Kurz, with adaptations by D. Lemire.
 */
size_t SIMDintersection(const uint32_t *set1, const size_t length1,
                        const uint32_t *set2, const size_t length2,
                        uint32_t *out);

#ifdef __AVX2__
#include <immintrin.h>

/*
 * Straight port of SIMDintersection to AVX2.
 */
size_t SIMDintersection_avx2(const uint32_t *set1, const size_t length1,
                        const uint32_t *set2, const size_t length2,
                        uint32_t *out);

#endif
/*
 * Given two arrays, this writes the intersection to out. Returns the
 * cardinality of the intersection.
 *
 * This is a well-written, but otherwise unsophisticated function.
 * Written by N. Kurz.
size_t nate_scalar(const uint32_t *set1, const size_t length1,
                   const uint32_t *set2, const size_t length2, uint32_t *out);
 */

/*
 * Given two arrays, this writes the intersection to out. Returns the
 * cardinality of the intersection.
 *
 * This applies a state-of-the-art algorithm. First coded by O. Kaser, adapted
 * by D. Lemire.
size_t onesidedgallopingintersection(const uint32_t *smallset,
                                     const size_t smalllength,
                                     const uint32_t *largeset,
                                     const size_t largelength, uint32_t *out);
 */

class IntersectionFactory {
public:
  static std::map<std::string, intersectionfunction> intersection_schemes;

  static vector<string> allNames() {
    vector<string> ans;
    for (auto i = intersection_schemes.begin(); i != intersection_schemes.end();
         ++i) {
      ans.push_back(i->first);
    }
    return ans;
  }

  static string getName(intersectionfunction v) {
    for (auto i = intersection_schemes.begin(); i != intersection_schemes.end();
         ++i) {
      if (i->second == v)
        return i->first;
    }
    return "UNKNOWN";
  }

  static bool valid(string name) {
    return (intersection_schemes.find(name) != intersection_schemes.end());
  }

  static intersectionfunction getFromName(string name) {
    if (intersection_schemes.find(name) == intersection_schemes.end()) {
      cerr << "name " << name << " does not refer to an intersection procedure."
           << endl;
      cerr << "possible choices:" << endl;
      for (auto i = intersection_schemes.begin();
           i != intersection_schemes.end(); ++i) {
        cerr << static_cast<string>(i->first)
             << endl; // useless cast, but just to be clear
      }
      return NULL;
    }
    return intersection_schemes[name];
  }
};

} // namespace SIMDCompressionLib

#endif /* SIMDCompressionAndIntersection_INTERSECTION_H_ */
