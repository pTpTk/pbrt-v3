#ifndef PBRT_UTILS_H
#define PBRT_UTILS_H

#include <utility>

namespace pbrt{
namespace gpu{
namespace utils{

template<typename T>
__both__
int get_buffer_size(T const * const ptr) {
  int count = 0;
  for (auto it = ptr; it != nullptr; ++it)
    ++count;
  return count;
}

template<typename T>
__both__
void swap(T& x, T& y) {
  T temp = std::move(x);
  x = std::move(y);
  y = std::move(temp);
}

}  // namespace utils
}  // namespace gpu
}  // namespace pbrt
#endif
