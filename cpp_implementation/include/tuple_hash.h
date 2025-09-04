#ifndef TUPLE_HASH_H
#define TUPLE_HASH_H

#include <tuple>
#include <functional>

namespace std {
  template<>
  struct hash<std::tuple<int, float, float>> {
    std::size_t operator()(const std::tuple<int, float, float>& t) const {
      auto h1 = std::hash<int>{}(std::get<0>(t));
      auto h2 = std::hash<float>{}(std::get<1>(t));
      auto h3 = std::hash<float>{}(std::get<2>(t));
      // FNV-1a inspired hash combining with golden ratio for better distribution
      std::size_t result = h1;
      result ^= h2 + 0x9e3779b9 + (result << 6) + (result >> 2);
      result ^= h3 + 0x9e3779b9 + (result << 6) + (result >> 2);
      return result;
    }
  };
  
  template<>
  struct hash<std::tuple<int, float>> {
    std::size_t operator()(const std::tuple<int, float>& t) const {
      auto h1 = std::hash<int>{}(std::get<0>(t));
      auto h2 = std::hash<float>{}(std::get<1>(t));
      // FNV-1a inspired hash combining
      std::size_t result = h1;
      result ^= h2 + 0x9e3779b9 + (result << 6) + (result >> 2);
      return result;
    }
  };
  
  template<>
  struct hash<std::tuple<int, int, int>> {
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
      auto h1 = std::hash<int>{}(std::get<0>(t));
      auto h2 = std::hash<int>{}(std::get<1>(t));
      auto h3 = std::hash<int>{}(std::get<2>(t));
      // FNV-1a inspired hash combining
      std::size_t result = h1;
      result ^= h2 + 0x9e3779b9 + (result << 6) + (result >> 2);
      result ^= h3 + 0x9e3779b9 + (result << 6) + (result >> 2);
      return result;
    }
  };
  
  template<>
  struct hash<std::pair<int, int>> {
    std::size_t operator()(const std::pair<int, int>& p) const {
      auto h1 = std::hash<int>{}(p.first);
      auto h2 = std::hash<int>{}(p.second);
      // FNV-1a inspired hash combining
      std::size_t result = h1;
      result ^= h2 + 0x9e3779b9 + (result << 6) + (result >> 2);
      return result;
    }
  };

  // Hash for nested pair of tuples (used in dynamic version)
  template<>
  struct hash<std::pair<std::tuple<int, float, float>, std::tuple<int, float, float>>> {
    std::size_t operator()(const std::pair<std::tuple<int, float, float>, std::tuple<int, float, float>>& p) const {
      auto h1 = std::hash<std::tuple<int, float, float>>{}(p.first);
      auto h2 = std::hash<std::tuple<int, float, float>>{}(p.second);
      // FNV-1a inspired hash combining
      std::size_t result = h1;
      result ^= h2 + 0x9e3779b9 + (result << 6) + (result >> 2);
      return result;
    }
  };
}

#endif // TUPLE_HASH_H
