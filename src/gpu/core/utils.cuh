#ifndef PBRT_UTILS_H
#define PBRT_UTILS_H

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

}  // namespace utils
}  // namespace gpu
}  // namespace pbrt
#endif
